import abc

import numpy as np


class FetchHiC(abc.ABC):
    SMALL_VALUE = 1e-12

    def fetch_data(self, genome_range1, genome_range2=None):
        """
        Parameters
        ----------
        genome_range1 : {str, GenomeRange}

        genome_range2 : {str, GenomeRange}, optional.

        Return
        ------
        pixels : pandas.core.frame.DataFrame
            Hi-C pixels table.
            The pixel table contains the non-zero upper triangle entries of the contact map.
        """
        return self.fetch_array(genome_range1, genome_range2)

    @abc.abstractmethod
    def fetch_array(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        pass

    @abc.abstractmethod
    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        pass

    def fetch_matrix(self, genome_range, resolution='auto'):
        """
        Fetch the matrix for plot.

        Parameters
        ----------
        genome_range : coolbox.utilities.GenomeRange
            The genome range to fetch.

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

        arr = wrap.fetch(genome_range)

        self.fetched_binsize = wrap.fetched_binsize  # expose fetched binsize

        # fill zero and nan with small value
        small = self.SMALL_VALUE
        arr[arr == 0] = small
        arr[np.isnan(arr)] = small

        if 'transform' in self.properties and self.properties['transform'] != 'no':
            arr = self.__transform_matrix(arr)

        return arr

    def __transform_matrix(self, arr):
        if self.properties['transform'] == 'log10':
            arr = np.log10(arr)
        elif self.properties['transform'] == 'log2':
            arr = np.log2(arr)
        elif self.properties['transform'] == 'log':
            arr = np.log(arr)
        return arr

