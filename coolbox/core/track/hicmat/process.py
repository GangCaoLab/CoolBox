import abc
import re
from typing import Tuple

import numpy as np
from scipy.linalg import toeplitz
from scipy.ndimage import gaussian_filter
from scipy.signal import convolve2d

from coolbox.utilities.logtools import get_logger

log = get_logger(__name__)


class ProcessHicMat(object):
    SMALL_VALUE = 1e-12

    def process_matrix(self, arr: np.ndarray) -> np.ndarray:
        properties = self.properties
        transform = properties['transform']
        normalize = properties['normalize']
        gaussian_sigma = properties['gaussian_sigma']
        process_func = properties['process_func']
        # process the matrix
        if transform != 'no':
            arr = self.transform_matrix(arr, transform)
        if normalize != 'no':
            arr = self.normalize_matrix(arr, normalize)
        if gaussian_sigma != 'no':
            arr = self.gaussian_matrix(arr, gaussian_sigma)
        if process_func != 'no':
            # user-defined process function
            func = process_func
            try:
                if callable(func):
                    arr = func(arr)
                elif isinstance(func, str):
                    func = eval(func)
                    arr = func(arr)
                else:
                    raise ValueError("process_func")
            except Exception as e:
                log.error(str(e))
                raise ValueError(
                    "process_func should a one argument function "
                    "receive a matrix return a processed matrix.")
        return arr

    def fill_zero_nan(self, arr: np.ndarray) -> np.ndarray:
        # fill zero and nan with small value
        small = self.SMALL_VALUE
        zero_indices = arr == 0
        nan_indices = np.isnan(arr)
        self.zero_indices = zero_indices
        self.nan_indices = nan_indices
        arr[zero_indices] = small
        arr[nan_indices] = small
        return arr

    @staticmethod
    def transform_matrix(arr: np.ndarray, method='log') -> np.ndarray:
        if method == 'log10':
            arr = np.log10(arr)
        elif method == 'log2':
            arr = np.log2(arr)
        elif method == 'log':
            arr = np.log(arr)
        return arr

    @staticmethod
    def gaussian_matrix(arr: np.ndarray, sigma=1.0) -> np.ndarray:
        arr = gaussian_filter(arr, sigma)
        return arr

    @classmethod
    def normalize_matrix(cls, mat: np.ndarray, cis=True, method='zscore') -> np.ndarray:
        if cis:
            res = cls.normalize_cis(mat, method)
        else:
            res = cls.normalize_trans(mat, method)
        return res

    @classmethod
    def normalize_cis(cls, mat: np.ndarray, method: str) -> np.ndarray:
        res = mat
        if method == 'total':
            total = np.sum(mat)
            if total != 0:
                res = mat / total
        elif method == 'expect':
            means = cls.diagonal_mean(mat)
            expect = toeplitz(means)
            res = mat / expect
        elif method == 'zscore':
            means, stds = cls.diagonal_mean_std(mat)
            mat_mean = toeplitz(means)
            mat_std = toeplitz(stds)
            res = (mat - mat_mean) / mat_std
        elif re.match("hiccups:.:.", method):
            p, w = method.strip("hiccups:").split(":")
            p, w = int(w), int(p)
            kernel = cls.__donut_kernel(p, w)

            def apply_donut(m):
                m_ext = np.zeros((m.shape[0] + 2 * (w - 1), m.shape[0] + 2 * (w - 1)))
                idx_center = slice(w, w + m.shape[0]), slice(w, w + m.shape[1])
                m_ext[idx_center] = m
                m_f = convolve2d(m_ext, kernel, mode='same')
                m_f = m_f[idx_center]
                return m_f

            means = cls.diagonal_mean(mat)
            exp_decay = toeplitz(means)
            m_donut = apply_donut(mat)
            exp_donut = apply_donut(exp_decay)
            exp = (m_donut / exp_donut) * exp_decay
            res = mat / exp
        else:
            log.warning(f"Cis-matrix does not support the {method} normalize.")
        return res

    @classmethod
    def normalize_trans(cls, mat: np.ndarray, method: str) -> np.ndarray:
        res = mat
        if method == 'total':
            total = np.sum(mat)
            if total != 0:
                res = mat / total
        elif method == 'zscore':
            mean = np.mean(mat)
            std = np.std(mat)
            res = (mat - mean) / std
        else:
            log.warning(f"Trans-matrix doest not support the {method} normalize.")
        return res

    @staticmethod
    def diagonal_mean(mat: np.ndarray) -> np.ndarray:
        return np.array([np.diagonal(mat, i).mean() for i in range(mat.shape[0])])

    @staticmethod
    def diagonal_mean_std(mat: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        means = []
        stds = []
        for i in range(mat.shape[0]):
            diagonal = np.diagonal(mat, i)
            means.append(diagonal.mean())
            stds.append(diagonal.std())
        stds = np.array(stds)
        try:
            stds[stds == 0] = stds[stds > 0].min()
        except Exception:
            log.error("Error when executing 'stds[stds > 0].min()'")
            pass
        return np.array(means), stds


    @staticmethod
    def __donut_kernel(p: int, w: int) -> np.ndarray:
        k1 = np.ones((2 * w + 1, 2 * w + 1))
        k2 = np.zeros((2 * w + 1, 2 * w + 1))
        k2[w - p:w + p + 1, w - p:w + p + 1] = 1
        k3 = np.zeros((2 * w + 1, 2 * w + 1))
        k3[:w - p, w] = 1
        k3[w, :w - p] = 1
        k3[w, w + 2 * p - 1:2 * w + 1] = 1
        k3[w + 2 * p - 1:2 * w + 1, w] = 1
        k = k1 - k2 - k3
        return k
