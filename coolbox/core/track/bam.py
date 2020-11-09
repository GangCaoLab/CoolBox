import pandas as pd
import numpy as np
from dna_features_viewer import GraphicFeature, GraphicRecord

from .base import Track
from .hist.plot import CoveragePlot
from coolbox.utilities import (
    get_logger, GenomeRange, split_genome_range
)
from coolbox.utilities.bam import *


log = get_logger(__name__)


class BAM(Track, CoveragePlot):

    """
    BAM alignment track.

    Parameters
    ----------
    file_ : str
        Path to .gtf(or .gtf.bgz) file.

    length_ratio_thresh : float
        Length ratio threshold of show alignments. (Default 0.01)

    height : float, optional
        The height of Spacer track. (Default: BAM.DEFAULT_HEIGHT)

    plot_type : str
        Plot type, 'alignment' or 'coverage'. (Default 'coverage')

    style : str, optional
        Track graph style(for 'coverage'), format {'fill', 'line:`size`', 'points:`size`'},
        example: 'line:2', 'points:0.5'. default: 'fill'

    color : {str}
        Plot color.

    alhpa : float.
        Plot alpha. (Default 1.0)

    bins : int
        Number of bins when plot coverage. (Default 200)

    data_range_style : str
        Data range show style(y-axis or text), when plot coverage. (Default 'y-axis')

    max_value : {float, 'auto'}, optional
        Max value of track. 'auto' for specify max value automatically, default 'auto'.

    min_value : {float, 'auto'}, optional
        Min value of track. 'auto' for specify max value automatically, default 'auto'.

    title : str, optional
        Label text, default ''.

    name : str, optional
        Track's name.

    """
    DEFAULT_HEIGHT = 3
    DEFAULT_COLOR = "#6688ff"

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "length_ratio_thresh": 0.005,
            "height": BAM.DEFAULT_HEIGHT,
            "plot_type": "coverage",
            "style": 'fill',
            "title": '',
            "color": BAM.DEFAULT_COLOR,
            "alpha": 1.0,
            "bins": 200,
            "data_range_style": 'y-axis',
            "max_value": "auto",
            "min_value": "auto",
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.indexed_bam = process_bam(file_)

    def plot(self, ax, chrom_region, start_region, end_region):
        gr = GenomeRange(chrom_region, start_region, end_region)
        ptype = self.properties.get("plot_type", "alignment")
        self.ax = ax
        if ptype == "alignment":
            self.plot_align(ax, gr)
        else:
            self.plot_coverage(ax, gr)

    def plot_align(self, ax, genome_range):
        gr = genome_range
        df = self.fetch_intervals(gr)
        df_ = df[np.bitwise_and(df['flag'], 0b100) == 0]
        len_thresh = self.properties.get("length_ratio_thresh", 0.005)
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

    def plot_coverage(self, ax, genome_range):
        gr = genome_range
        bins = self.properties.get("bins", 40)
        scores_per_bin = self.fetch_coverage(str(gr), bins)
        super().plot_coverage(ax, genome_range, scores_per_bin)
        self.plot_label()

    def fetch_data(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Sam interval table.
        """
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
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

    def fetch_coverage(self, genome_range, bins=100):
        scores_per_bin = coverage_by_samtools(
            self.indexed_bam,
            str(genome_range),
            bins
        )
        return scores_per_bin
