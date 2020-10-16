from coolbox.fetchdata.base import FetchTrackData


class FetchCool(FetchTrackData):

    def fetch_data(self, genome_range1, genome_range2=None):
        """
        Parameters
        ----------
        genome_range1 : {str, GenomeRange}

        genome_range2 : {str, GenomeRange}, optional.

        Return
        ------
        arr : numpy.ndarray
        """
        return self.fetch_array(genome_range1, genome_range2)

    def fetch_array(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        genome_range2 : {str, GenomeRange}, optional.

        balance : bool, optional
            balance matrix or not,
            default `self.is_balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        Return
        ------
        arr : numpy.ndarray
        """
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        if balance is None:
            balance = self.is_balance
        wrap = CoolerWrap(path, balance=balance, binsize=resolution)

        arr = wrap.fetch(genome_range, genome_range2)
        return arr

    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, resolution='auto', join=True):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        genome_range2 : {str, GenomeRange}, optional.

        balance : bool, optional
            balance matrix or not,
            default `self.is_balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        join : bool
            whether to expand the bin ID columns
            into (chrom, start, end).
            default True

        Return
        ------
        pixels : pandas.core.frame.DataFrame
            Hi-C pixels table.
            The pixel table contains the non-zero upper triangle entries of the contact map.
        """
        from coolbox.utilities.hic.wrap import CoolerWrap

        path = self.properties['file']
        if balance is None:
            balance = self.is_balance
        wrap = CoolerWrap(path, balance=balance, binsize=resolution)

        pixels = wrap.fetch_pixels(genome_range, genome_range2, join=join)
        return pixels