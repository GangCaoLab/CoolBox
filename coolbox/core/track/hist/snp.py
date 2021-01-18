import os.path as osp
import subprocess as subp

import numpy as np
import pandas as pd

from coolbox.utilities import GenomeRange
from coolbox.utilities.bed import tabix_query, build_snp_index
from .base import HistBase


class SNP(HistBase):
    """
    Track for show SNPs Manhattan plot.

    Input file is a tab-split file contains SNP's chrom, position, pvalue information.
    You should specify it's column indexes by `col_chrom`, `col_pos`, `col_pval` arguments.

    Parameters
    ----------
    file : str
        Path to input .snp/.vcf file.

    col_chrom : int
        Column index of seqname(chromosome).

    col_pos : int
        Column index of snp position.

    col_pval : int
        Column index of snp p-value.

    pval_transform : {'-log2', '-log10'}
        Transform the p value. Default '-log10'.
    """
    COL_CHROM = 0
    COL_POS = 2
    COL_PVAL = 9

    DEFAULT_PROPERTIES = {
        "style": HistBase.STYLE_SCATTER,
        "color": "grey",
        "threshold_color": "#ff9c9c",
        "threshold": 0.05,
        "alpha": 0.5,
        "size": 10,
        # TODO if min_value is set to 'auto', the min_value is not 0, different to the original codes
        "min_value": "0",
        "max_value": "auto",
        "height": 5.0,
        # processing
        "pval_transform": "-log10",
        "col_chrom": COL_CHROM,
        "col_pos": COL_POS,
        "col_pval": COL_PVAL,
    }

    def __init__(self, file, **kwargs):
        properties = SNP.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.bgz_file = build_snp_index(
            file,
            self.properties['col_chrom'],
            self.properties['col_pos']
        )
        # TODO what does this mean?
        self.properties['threshold'] = self.transform_fn()(self.properties['threshold'])

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        df = self.fetch_data(gr, **kwargs)
        df['score'] = self.transform_fn()(df['score'])
        return df

    def fetch_data(self, gr: GenomeRange, **kwargs):
        ix_chrom = self.properties['col_chrom']
        ix_pos = self.properties['col_pos']
        ix_pval = self.properties['col_pval']
        rows = self.load_range(gr)
        if len(rows) == 0:
            gr.change_chrom_names()
            rows = self.load_range(gr)
        df = pd.DataFrame(rows)
        if df.shape[0] > 0:
            columns = [f'col_{i}' for i in range(df.shape[1])]
            columns[ix_chrom] = "chrom"
            columns[ix_pos] = "pos"
            columns[ix_pval] = "score"
            df.columns = columns
        return df

    def load_range(self, gr):
        rows = []
        ix_pos = self.properties['col_pos']
        ix_pval = self.properties['col_pval']
        for items in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end):
            items[ix_pos] = int(items[ix_pos])
            items[ix_pval] = float(items[ix_pval])
            rows.append(items)
        return rows

    def transform_fn(self):
        method = self.properties.get('pval_transform', '-log2')
        if method == "-log2":
            return lambda x: -np.log2(x)
        elif method == "-log10":
            return lambda x: -np.log10(x)
        else:
            return lambda x: x
