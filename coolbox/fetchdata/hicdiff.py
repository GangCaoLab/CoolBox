import numpy as np
from scipy.linalg import toeplitz

from .base import FetchTrackData


class FetchHiCDiff(FetchTrackData):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def fetch_related_tracks(self, genome_range, resolution=None):
        if resolution:
            reso = resolution
        else:
            reso = self.properties['resolution']
        hic1 = self.properties['hic1']
        hic2 = self.properties['hic2']
        mat1 = hic1.fetch_matrix(genome_range, reso)
        mat2 = hic2.fetch_matrix(genome_range, reso)
        return mat1, mat2

    def __normalize_data(self, mat):
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

    def __diff_data(self, mat1, mat2):
        diff_mth = self.properties['diff_method']
        if diff_mth == 'log2fc':
            return np.log2((mat1 + 1)/(mat2 + 1))
        else:
            return mat1 - mat2

    def fetch_data(self, genome_range, resolution=None):
        mat1, mat2 = self.fetch_related_tracks(genome_range, resolution)
        mat1, mat2 = self.__normalize_data(mat1), self.__normalize_data(mat2)
        diff = self.__diff_data(mat1, mat2)
        return diff
