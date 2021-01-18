import numpy as np
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord

from coolbox.utilities import (
    get_logger, GenomeRange, split_genome_range
)
from coolbox.utilities.bam import process_bam, query_bam
from .base import Track

log = get_logger(__name__)


class BAM(Track):
    """
    BAM alignment track for plotting reads.

    Parameters
    ----------
    file : str
        Path to .bam .sam file.

    length_ratio_thresh : float
        Length ratio threshold of show alignments. (Default 0.01)

    """

    DEFAULT_PROPERTIES = {
        "height": 3,
        "length_ratio_thresh": 0.05,
        "color": "#6688ff"
    }

    def __init__(self, file, **kwargs):
        properties = BAM.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs,
        })
        properties.update(kwargs)
        super().__init__(properties)
        self.indexed_bam = process_bam(file)

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            Sam interval table.
            The DataFrame table should has columns like:

            columns = ["qname", "flag", "rname", "pos", "mapq", "cigar",
                      "rnext", "pnext", "tlen", "seq", "qual", "options"]
        """
        return self.fetch_intervals(gr)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.plot_align(ax, gr)

    def plot_align(self, ax, gr: GenomeRange):
        assert isinstance(gr, GenomeRange), "The input gr should be type GenomeRange"
        df = self.fetch_plot_data(gr)
        df_ = df[np.bitwise_and(df['flag'], 0b100) == 0]
        len_thresh = self.properties["length_ratio_thresh"]
        df_ = df_[df_['seq'].str.len() > (gr.length * len_thresh)]
        if df_.shape[0] <= 0:
            return
        rev_flag = np.bitwise_and(df['flag'], 0b10000) != 0
        features = []
        for idx, row in df_.iterrows():
            start = row['pos'] - gr.start
            end = row['pos'] + len(row['seq']) - gr.start
            strand = -1 if rev_flag.iloc[idx] else 1
            gf = GraphicFeature(
                start=start,
                end=end,
                strand=strand,
                color=self.properties['color'],
            )
            features.append(gf)
        record = GraphicRecord(sequence_length=gr.length, features=features)
        record.plot(
            ax=ax,
            with_ruler=False,
            draw_line=False
        )

    def fetch_intervals(self, genome_range: GenomeRange):
        chrom, start, end = split_genome_range(genome_range)
        rows = []
        for row_items in query_bam(self.indexed_bam, chrom, start, end, split=True):
            rows.append(row_items)
        # https://samtools.github.io/hts-specs/SAMv1.pdf
        fields = ["qname", "flag", "rname", "pos", "mapq", "cigar",
                  "rnext", "pnext", "tlen", "seq", "qual", "options"]
        df = pd.DataFrame(rows, columns=fields)
        if df.shape[0] > 0:
            df['flag'] = df['flag'].astype(int)
            df['pos'] = df['pos'].astype(int)
            df['mapq'] = df['mapq'].astype(int)
        return df
