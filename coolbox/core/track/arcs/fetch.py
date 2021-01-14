import pandas as pd
from coolbox.utilities.genome import GenomeRange
from coolbox.utilities.bed import pairix_query


class FetchParix(object):
    @staticmethod
    def fetch_intervals(bgz_file, gr: GenomeRange, gr2: GenomeRange = None) -> pd.DataFrame:
        rows = list(pairix_query(bgz_file, gr, second=gr2, split=True))
        if len(rows) == 0:
            gr.change_chrom_names()
            if gr2:
                gr2.change_chrom_names()
            rows = list(pairix_query(bgz_file, gr, second=gr2, split=True))

        return pd.DataFrame(rows)
