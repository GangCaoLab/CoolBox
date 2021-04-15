import random
import re

import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord

from coolbox.utilities import (
    split_genome_range, get_logger, GenomeRange
)
from coolbox.utilities.bed import (
    build_gtf_index, tabix_query
)
from .base import Track

log = get_logger(__name__)


class GTF(Track):
    """
    GTF gene annotation track.

    Parameters
    ----------
    file : str
        Path to .gtf(or .gtf.bgz) file.

    row_filter : str, optional
        Row filter expression, only keep the rows for draw. (Default 'feature == "gene"')

    length_ratio_thresh : float
        Length ratio threshold of features, (Default 0.01)

    color : {str, 'random'}, optional
        When the color is random, color for each gene will be randomly selected.

    name_attribute : {'auto', 'gene_name', 'gene_id', str}, optional
        Use which attribute to show feature's name.
        Default use 'auto'(try 'gene_name' -> 'gene_id' -> 'position_string')
    """
    RANDOM_COLORS = [
        "#ffcccc",
        "#ffd700",
        "#cffccc",
        "#ff9c9c",
        "#66ccff"
    ]

    DEFAULT_PROPERTIES = {
        "height": 4,
        "color": "random",
        "row_filter": 'feature == "gene"',
        "length_ratio_thresh": 0.005,
        "name_attribute": "auto",
    }

    def __init__(self, file, **kwargs):
        properties = GTF.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(properties)
        self.bgz_file = build_gtf_index(file)
        color = self.properties['color']
        if (type(color) is str) and (color.startswith('#')):
            self.colors = [color]
        elif type(color) is list:
            self.colors = [c for c in color if (type(c) is str) and c.startswith('#')]
            if not self.colors:
                self.colors = GTF.RANDOM_COLORS
        else:
            self.colors = GTF.RANDOM_COLORS

    def fetch_data(self, gr: GenomeRange, **kwargs):
        """
        Returns
        -------
        df: pandas.DataFrame
            should be with the format like:

            columns = ['seqname', 'source', 'feature', 'start', 'end',
                        'score', 'strand', 'frame', 'attribute', 'feature_name']

        """
        return self.fetch_intervals(gr)

    def fetch_intervals(self, gr: GenomeRange):
        """

        Parameters
        ----------
        gr : {str, GenomeRange}

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            Annotation interval table.
        """
        rows = [row for row in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end)]
        if not rows:
            gr.change_chrom_names()
            for row in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end):
                rows.append(row)

        columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        df = pd.DataFrame(rows, columns=columns)
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        name_attr = self.properties.get("name_attr", "auto")
        if name_attr == "auto":
            gene_name = df['attribute'].str.extract(".*gene_name (.*?) ").iloc[:, 0].str.strip('\";')
            if gene_name.hasnans:
                gene_id = df['attribute'].str.extract(".*gene_id (.*?) ").iloc[:, 0].str.strip('\";')
                gene_name.fillna(gene_id, inplace=True)
                if gene_name.hasnans:
                    pos_str = df['seqname'].astype(str) + ":" +\
                              df['start'].astype(str) + "-" +\
                              df['end'].astype(str)
                    gene_name.fillna(pos_str, inplace=True)
            df['feature_name'] = gene_name
        else:
            df['feature_name'] = df['attribute'].str.extract(f".*{name_attr} (.*?) ").iloc[:, 0].str.strip('\";')
        return df

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        df = self.fetch_plot_data(gr)
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
        region_length = gr.end - gr.start
        len_ratio_th = self.properties["length_ratio_thresh"]
        df = df[(df["end"] - df["start"]) > region_length * len_ratio_th]
        features = []
        for _, row in df.iterrows():
            gf = GraphicFeature(
                start=row['start'],
                end=row['end'],
                strand=(1 if row['strand'] == '+' else -1),
                label=row['feature_name'],
                color=random.choice(self.colors),
            )
            features.append(gf)
        record = GraphicRecord(
            sequence_length=gr.end - gr.start,
            features=features,
            first_index=gr.start
        )
        record.plot(
            ax=ax,
            with_ruler=False,
            draw_line=False
        )
        self.plot_label()
