from coolbox.fetchdata.base import FetchTrackData


class FetchDotHiC(FetchTrackData):

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

    def fetch_array(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
            matrix balance method,
            default `self.balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        Return
        ------
        arr : numpy.ndarray
        """
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        if balance is None:
            balance = self.balance
        wrap = StrawWrap(path, normalization=balance, binsize=resolution)

        arr = wrap.fetch(genome_range, genome_range2)
        return arr

    def fetch_pixels(self, genome_range, genome_range2=None, balance=None, resolution='auto'):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}
            Intervals within input chromosome range.

        balance : {bool, 'KR', 'VC', 'VC_SQRT'}, optional
            matrix balance method,
            default `self.balance`.

        resolution : {'auto', int}
            resolution of the data. for example 5000.
            'auto' for calculate resolution automatically.
            default 'auto'

        Return
        ------
        pixels : pandas.core.frame.DataFrame
        """
        from coolbox.utilities.hic.wrap import StrawWrap

        path = self.properties['file']
        if balance is None:
            balance = self.balance
        wrap = StrawWrap(path, normalization=balance, binsize=resolution)

        pixels = wrap.fetch_pixels(genome_range, genome_range2)
        return pixels