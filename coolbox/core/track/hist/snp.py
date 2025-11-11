import numpy as np

from coolbox.utilities import GenomeRange
from coolbox.utilities.reader.tab import get_indexed_tab_reader
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
    FIELDS = [
        "chrom", "rsid", "pos", "a1", "a2", "n", "maf", "beta", "se", "pval"
    ]

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
        "fields": FIELDS,
    }

    def __init__(self, file, **kwargs):
        properties = SNP.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(**properties)
        self.reader = get_indexed_tab_reader(file, columns=properties['fields'])
        # TODO what does this mean?
        self.properties['threshold'] = self.transform_fn()(self.properties['threshold'])

    def fetch_plot_data(self, gr: GenomeRange, **kwargs):
        df = self.fetch_data(gr, **kwargs)
        df['pos'] = df['pos'].astype(int)
        df['score'] = self.transform_fn()(df['pval'].astype(float))
        return df

    def fetch_data(self, gr: GenomeRange, **kwargs):
        df = self.reader.query_var_chr(gr)
        return df

    def transform_fn(self):
        method = self.properties.get('pval_transform', '-log2')
        if method == "-log2":
            return lambda x: -np.log2(x)
        elif method == "-log10":
            return lambda x: -np.log10(x)
        else:
            return lambda x: x
