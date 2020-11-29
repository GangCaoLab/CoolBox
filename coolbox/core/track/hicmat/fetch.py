import abc

import numpy as np
from scipy.linalg import toeplitz
from scipy.ndimage import gaussian_filter

from coolbox.utilities.logtools import get_logger

log = get_logger(__name__)


class FetchHiC(abc.ABC):
    SMALL_VALUE = 1e-12

    def fetch_data(self, genome_range1, genome_range2=None, resolution='auto'):
        """
        Parameters
        ----------
        genome_range1 : {str, GenomeRange}

        genome_range2 : {str, GenomeRange}, optional.

        resolution : {int, 'auto'}, optional

        Return
        ------
        matrix : np.array
            Hi-C contact matrix.
        """
        return self.fetch_matrix(genome_range1, genome_range2, resolution=resolution)

    @abc.abstractmethod
    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        pass

    def fetch_matrix(self, genome_range, genome_range2=None, resolution='auto'):
        """
        Fetch the matrix for plot.

        Parameters
        ----------
        genome_range : coolbox.utilities.GenomeRange
            The genome range to fetch.

        genome_range2 : coolbox.utilities.GenomeRange, optional
            Second genome range to fetch.

        resolution : {'auto', int}
            The matrix resolution, for multi-resolution(.hic or multi-cool) file.
            Use 'auto' to infer the resolution automatically.
            default 'auto'
        """
        from coolbox.utilities.hic.wrap import StrawWrap, CoolerWrap

        path = self.properties['file']

        from coolbox.utilities.hic.tools import file_type
        if file_type(self.properties['file']) == '.hic':
            wrap = StrawWrap(path, normalization=self.balance, binsize=resolution)
        else:
            wrap = CoolerWrap(path, balance=self.balance, binsize=resolution)

        arr = wrap.fetch(genome_range, genome_range2)

        self.fetched_binsize = wrap.fetched_binsize  # expose fetched binsize

        # fill zero and nan with small value
        small = self.SMALL_VALUE
        arr[arr == 0] = small
        arr[np.isnan(arr)] = small

        # process the matrix
        if 'transform' in self.properties and self.properties['transform'] != 'no':
            arr = self.__transform_matrix(arr)
        if 'normalize' in self.properties and self.properties['normalize'] != 'no':
            arr = self.__normalize_matrix(arr)
        if 'gaussian_sigma' in self.properties and self.properties['gaussian_sigma'] != 'no':
            arr = self.__gaussian_matrix(arr)
        if 'process_func' in self.properties and self.properties['process_func'] != 'no':
            # user-defined process function
            func = self.properties['process_func']
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

    def __normalize_matrix(self, mat):
        norm_mth = self.properties['normalize']
        res = mat
        if norm_mth == 'total':
            total = np.sum(mat)
            if total != 0:
                res = mat / total
        elif norm_mth == 'expect':
            means = [np.diagonal(mat, i).mean() for i in range(mat.shape[0])]
            expect = toeplitz(means)
            res = mat / expect
        elif norm_mth == 'zscore':
            means = []
            stds = []
            for i in range(mat.shape[0]):
                diagonal = np.diagonal(mat, i)
                means.append(diagonal.mean())
                stds.append(diagonal.std())
            stds = np.array(stds)
            stds[stds == 0] = stds[stds > 0].min()
            mat_mean = toeplitz(means)
            mat_std = toeplitz(stds)
            res = (mat - mat_mean) / mat_std
        return res

    def __transform_matrix(self, arr):
        if self.properties['transform'] == 'log10':
            arr = np.log10(arr)
        elif self.properties['transform'] == 'log2':
            arr = np.log2(arr)
        elif self.properties['transform'] == 'log':
            arr = np.log(arr)
        return arr

    def __gaussian_matrix(self, arr):
        sigma = self.properties['gaussian_sigma']
        arr = gaussian_filter(arr, sigma)
        return arr
