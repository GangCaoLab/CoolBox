import random
import re

import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord

from .base import Track
from coolbox.utilities.bed import (
    build_gtf_index, tabix_query
)
from coolbox.utilities import (
    split_genome_range, get_logger, GenomeRange
)


log = get_logger(__name__)


class GTF(Track):

    """
    GTF gene annotation track.

    Parameters
    ----------
    file_ : str
        Path to .gtf(or .gtf.bgz) file.

    row_filter : str, optional
        Filter rows, only keep the rows for draw. (Default 'feature == "gene"')

    length_ratio_thresh : float
        Length ratio threshold of features, (Default 0.01)

    height : float, optional
        The height of track. (Default: GTF.DEFAULT_HEIGHT)

    color : {str, List[str]}
        Annotation color. (Default: 'random')

    title : str, optional
        Label text, default ''.

    name : str, optional
        Track's name.

    """
    DEFAULT_HEIGHT = 4
    DEFAULT_COLOR = '#2855d8'
    RANDOM_COLORS = [
        "#ffcccc",
        "#ffd700",
        "#cffccc",
        "#ff9c9c",
        "#66ccff"
    ]

    def __init__(self, file_, **kwargs):
        properties_dict = {
            "file": file_,
            "row_filter": 'feature == "gene"',
            "length_ratio_thresh": 0.005,
            "height": GTF.DEFAULT_HEIGHT,
            "title": '',
            "color": 'random',
        }
        properties_dict.update(kwargs)
        super().__init__(properties_dict)
        self.bgz_file = build_gtf_index(file_)
        color = self.properties['color']
        if (type(color) is str) and (color.startswith('#')):
            self.colors = [color]
        elif type(color) is list:
            self.colors = [c for c in color if (type(c) is str) and c.startswith('#')]
            if len(self.colors) == 0:
                self.colors = GTF.RANDOM_COLORS
        else:
            self.colors = GTF.RANDOM_COLORS

    def fetch_data(self, genome_range):
        return self.fetch_intervals(genome_range)

    def fetch_intervals(self, genome_range):
        """
        Parameters
        ----------
        genome_range : {str, GenomeRange}

        Return
        ------
        intervals : pandas.core.frame.DataFrame
            Annotation interval table.
        """
        chrom, start, end = split_genome_range(genome_range)
        rows = []
        for row in tabix_query(self.bgz_file, chrom, start, end):
            rows.append(row)
        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        df = pd.DataFrame(rows, columns=columns)
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        df['gene_name'] = df['attribute'].str.extract(".*gene_name (.*?) ").iloc[:, 0].str.strip('\";')
        df['gene_name'][df['gene_name'].isna()] = ""
        return df

    def plot(self, ax, chrom_region, start_region, end_region):
        self.ax = ax
        genome_range = GenomeRange(chrom_region, start_region, end_region)
        itv_df = self.fetch_intervals(genome_range)
        df = itv_df
        if self.has_prop("row_filter"):
            filters = self.properties["row_filter"]
            for filter_ in filters.split(";"):
                try:
                    op_idx = list(re.finditer("[=><!]", filter_))[0].start()
                    l_ = filter_[:op_idx].strip()
                    r_ = filter_[op_idx:]
                    df = eval(f'df[df["{l_}"]{r_}]')
                except IndexError:
                    log.warning(f"row filter {filter_} is not valid.")
        region_length = end_region - start_region
        if self.has_prop("length_ratio_thresh"):
            len_ratio_th = self.properties["length_ratio_thresh"]
            df = df[(df["end"] - df["start"]) > region_length*len_ratio_th]
        features = []
        for _, row in df.iterrows():
            gf = GraphicFeature(
                start=row['start']-start_region,
                end=row['end']-start_region,
                strand=(1 if row['strand'] == '+' else -1),
                label=row['gene_name'],
                color=random.choice(self.colors),
            )
            features.append(gf)
        record = GraphicRecord(sequence_length=end_region-start_region, features=features)
        record.plot(
            ax=ax,
            with_ruler=False,
            draw_line=False
        )
        self.plot_label()

