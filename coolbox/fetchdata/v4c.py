from coolbox.fetchdata.base import FetchTrackData


class FetchVirtual4C(FetchTrackData):

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        mean_arr : numpy.ndarray
        """
        return self.fetch_mean_arr(genome_range)