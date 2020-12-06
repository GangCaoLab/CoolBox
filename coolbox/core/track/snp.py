import os.path as osp
import subprocess as subp

import numpy as np
import pandas as pd

from coolbox.utilities import (
    to_gr, GenomeRange
)
from coolbox.utilities.bed import tabix_query
from .base import Track


def build_indexed_snp(file, col_chrom, col_pos):
    c = col_chrom + 1
    p = col_pos + 1
    if file.endswith(".bgz"):
        bgz_file = file
    elif osp.exists(file + '.bgz'):
        bgz_file = file + '.bgz'
    else:
        bgz_file = file + '.bgz'
        if file.endswith('.gz'):
            cmd = "zcat"
        else:
            cmd = "cat"
        cmd += f" {file} | sort -k{c},{c} -k{p},{p}n | bgzip > {bgz_file}"
        subp.check_call(cmd, shell=True)
    index_file = bgz_file + '.tbi'
    if not osp.exists(index_file):
        cmd = ['tabix', '-s', str(c), '-b', str(p), '-e', str(p), bgz_file]
        subp.check_call(cmd)
    return bgz_file


class SNP(Track):
    """
    Track for show SNPs Manhattan plot.
    Input file is a tab-split file,
    contain SNP's chrom, position, pvalue information.
    You should specify it's column indexes by
    `col_chrom`, `col_pos`, `col_pval` arguments.

    Parameters
    ----------
    file : str
        Path to input SNP file.

    col_chrom : int
        Column index of seqname(chromosome).

    col_pos : int
        Column index of snp position.

    col_pval : int
        Column index of snp p-value.

    color : str
        Point color. Default '#66ccff'

    alpha : float
        Point alpha. Default 0.5

    size : float
        Point size. Default 10.

    pval_thresh : float
        p-value threshold. Default 0.05

    pval_transform : {'-log2', '-log10'}
        Transform the p value. Default '-log10'.

    sig_color : str
        Significant points color. Default '#ff9c9c'

    show_line : bool
        Show threshold line or not. Default True

    line_color : str
        Threshold line color. Default '#000000'

    line_width : float
        Threshold line width. Default 2.0

    line_style : str
        Threshold line style. Default '--'

    min_value : {'auto', float}
        Track min value.

    max_value : {'auto', float}
        Track max value.

    height : float
        Track height.

    title : str
        Track title.
    """
    COL_CHROM = 0
    COL_POS = 2
    COL_PVAL = 9

    def __init__(self, file, **kwargs):
        properties_dict = {
            "file": file,
            "col_chrom": self.COL_CHROM,
            "col_pos": self.COL_POS,
            "col_pval": self.COL_PVAL,
            "color": "grey",
            "alpha": 0.5,
            "size": 10,
            "pval_thresh": 0.05,
            "pval_transform": "-log10",
            "sig_color": "#ff9c9c",
            "show_line": True,
            "line_color": "#000000",
            "line_width": 2.0,
            "line_style": "--",
            "max_value": "auto",
            "min_value": "auto",
            "height": 5.0,
            "title": "",
        }
        super().__init__(properties_dict)
        self.bgz_file = build_indexed_snp(
            file,
            properties_dict['col_chrom'],
            properties_dict['col_pos']
        )

    def plot(self, ax, region_chrom, region_start, region_end):
        self.ax = ax
        gr = GenomeRange(region_chrom, region_start, region_end)
        df = self.fetch_data(gr)
        if df.shape[0] == 0:
            return
        thresh = self.properties['pval_thresh']
        df, thresh = self.__transform_pval(df, thresh)
        df_sig = df[df['pval'] > thresh]
        df_nosig = df[df['pval'] <= thresh]
        if df_sig.shape[0] > 0:
            ax.scatter(
                df_sig['position'],
                df_sig['pval'],
                s=self.properties['size'],
                alpha=self.properties['alpha'],
                c=self.properties['sig_color']
            )
        if df_nosig.shape[0] > 0:
            ax.scatter(
                df_nosig['position'],
                df_nosig['pval'],
                s=self.properties['size'],
                alpha=self.properties['alpha'],
                c=self.properties['color']
            )
        if self.properties['show_line'] != "no":
            ax.hlines(
                y=thresh,
                xmin=gr.start,
                xmax=gr.end,
                colors=[self.properties['line_color']],
                linestyles=[self.properties['line_style']],
                linewidths=[self.properties['line_width']],
            )
        self.__adjust_plot(ax, gr, df)
        self.plot_label()
        if hasattr(self, 'y_ax'):
            self.plot_y_axis(ax, self.y_ax)

    def __adjust_plot(self, ax, gr, df):
        ax.set_xlim(gr.start, gr.end)
        max_ = self.properties['max_value']
        min_ = self.properties['min_value']
        if max_ == "auto":
            max_ = df['pval'].max()
        if min_ == "auto":
            min_ = df['pval'].min()
        ax.set_ylim(min_, max_)

    def __transform_pval(self, df, thresh):
        method = self.properties['pval_transform']
        if method == '-log2':
            df['pval'] = -np.log2(df['pval'])
            thresh = -np.log2(thresh)
        else:
            df['pval'] = -np.log10(df['pval'])
            thresh = -np.log10(thresh)
        return df, thresh

    def fetch_data(self, genome_range):
        gr = to_gr(genome_range)
        ix_chrom = self.properties['col_chrom']
        ix_pos = self.properties['col_pos']
        ix_pval = self.properties['col_pval']
        rows = self.__load(gr)
        if len(rows) == 0:
            gr.change_chrom_names()
            rows = self.__load(gr)
        df = pd.DataFrame(rows)
        if df.shape[0] > 0:
            columns = [f'col_{i}' for i in range(df.shape[1])]
            columns[ix_chrom] = "chrom"
            columns[ix_pos] = "position"
            columns[ix_pval] = "pval"
            df.columns = columns
        return df

    def __load(self, gr):
        rows = []
        ix_pos = self.properties['col_pos']
        ix_pval = self.properties['col_pval']
        for items in tabix_query(self.bgz_file, gr.chrom, gr.start, gr.end):
            items[ix_pos] = int(items[ix_pos])
            items[ix_pval] = float(items[ix_pval])
            rows.append(items)
        return rows
