import pandas as pd
from coolbox.utilities.genome import GenomeRange
from coolbox.utilities.bed import pairix_query


class FetchParix(object):
    def fetch_intervals(self, bgz_file, gr: GenomeRange, gr2: GenomeRange = None) -> pd.DataFrame:
        open_region = self.properties.get("open_region") == "yes"
        rows = list(pairix_query(bgz_file, gr, second=gr2, open_region=open_region, split=True))
        if len(rows) == 0:
            gr.change_chrom_names()
            if gr2:
                gr2.change_chrom_names()
            rows = list(pairix_query(bgz_file, gr, second=gr2, open_region=open_region, split=True))

        return pd.DataFrame(rows)
