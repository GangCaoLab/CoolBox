import random
import re

from dna_features_viewer import GraphicFeature, GraphicRecord

from coolbox.plots.track.base import TrackPlot
from coolbox.utilities import (
    get_logger, GenomeRange
)


log = get_logger(__name__)


class PlotGTF(TrackPlot):

    DEFAULT_COLOR = '#2855d8'
    RANDOM_COLORS = [
        "#ffcccc",
        "#ffd700",
        "#cffccc",
        "#ff9c9c",
        "#66ccff"
    ]

    def __init__(self, *args, **kwargs):
        TrackPlot.__init__(self, *args, **kwargs)
        color = self.properties['color']
        if (type(color) is str) and (color.startswith('#')):
            self.colors = [color]
        elif type(color) is list:
            self.colors = [c for c in color if (type(c) is str) and c.startswith('#')]
            if len(self.colors) == 0:
                self.colors = PlotGTF.RANDOM_COLORS
        else:
            self.colors = PlotGTF.RANDOM_COLORS

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

