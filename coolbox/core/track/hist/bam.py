from .base import HistBase, HistData, GenomeRange
from coolbox.utilities.bam import process_bam, coverage_by_samtools


class BAMCov(HistBase):
    """
    Alignment reads coverage track.

    Parameters
    ----------
    file: str
        File path of bam file.

    num_bins: int, optional
        Number of bins to plot hist fig.


    """

    DEFAULT_PROPERTIES = {
        "height": 3,
        "style": HistBase.STYLE_FILL,
        "color": "#6688ff",
        "num_bins": 200
    }

    def __init__(self, file, **kwargs):
        properties = BAMCov.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs,
        })
        super().__init__(**properties)
        self.indexed_bam = process_bam(file)

    def fetch_data(self, gr: GenomeRange, **kwargs) -> HistData:
        bins = self.properties.get("num_bins", 200)
        scores_per_bin = self.fetch_coverage(gr, bins)
        return scores_per_bin

    def fetch_coverage(self, genome_range: GenomeRange, bins=100):
        scores_per_bin = coverage_by_samtools(
            self.indexed_bam,
            str(genome_range),
            bins
        )
        return scores_per_bin
