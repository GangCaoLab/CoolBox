import joblib
import pandas as pd
from scipy.sparse import coo_matrix
from peakachu.scoreUtils import Chromosome
from peakachu import peakacluster
from coolbox.api import *
from coolbox.core.track.hicmat import Cool
from coolbox.utilities import GenomeRange, correspond_track

model = joblib.load('../../tests/test_data/down100.ctcf.pkl')


class Peakachu(Cool):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def fetch_data(self, genome_range1, **kwargs):
        matrix = super().fetch_data(genome_range1, **kwargs) # reuse Cool.fetch_data to get the raw contacts
        ch = Chromosome(coo_matrix(matrix), model) # call peakachu
        pval_ma, raw_ma = ch.score()
        self.pval_ma = pval_ma.tocoo()
        # should return a np.ndarray representing matrix
        return (pval_ma + pval_ma.T).toarray()

@track_to_coverage
class PeakaChuLoops(ArcsBase):
    def __init__(self, peakachu, **kwargs):
        super().__init__(style="hicpeaks", **kwargs)
        peakachu = correspond_track(peakachu)
        self.peakachu = peakachu

    def fetch_data(self, gr: GenomeRange, **kwargs):
        ma = self.peakachu.pval_ma
        mask = ma.data > 0.87 # threshold to filter pvals
        pixels = dict(zip(zip(ma.row[mask], ma.col[mask]), 1 / ma.data[mask]))
        rows = [(x, x + 1, y, y + 1) # call peakachu
                for [x, y], *_ in peakacluster.local_clustering(pixels, res=10000)]
        # should return a pd.DataFrame with columns ['start1', 'end1', 'start2', 'end2']
        peaks = pd.DataFrame(rows, columns=['start1', 'end1', 'start2', 'end2']) * self.peakachu.fetched_binsize + gr.start
        return peaks
