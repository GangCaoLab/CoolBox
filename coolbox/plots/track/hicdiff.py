from .hicmatrix import PlotHiCMatrix


class PlotHiCDiff(PlotHiCMatrix):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.properties['transform'] = 'no'
        self.properties['norm'] = 'no'

    def fetch_matrix(self, genome_range, resolution='auto'):
        diff = self.fetch_data(genome_range, None)
        try:
            self.small_value = diff[diff > 0].min()
        except:
            self.small_value = 1e-12
        return diff

