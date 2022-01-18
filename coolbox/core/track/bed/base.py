from typing import Union

import pandas as pd
import matplotlib

from coolbox.utilities import get_logger
from coolbox.utilities.bed import build_bed_index
from coolbox.utilities.genome import GenomeRange
from coolbox.core.track.base import Track

log = get_logger(__name__)


class BedBase(Track):
    """
    BED Base track.

    Parameters
    ----------
    file: str
        The file path of `.bed` file.

    color : str, optional
        Track color, 'bed_rgb' for auto specify color according to bed record.
        (Default: 'bed_rgb')

    border_color : str, optional
        Border_color of gene. (Default: 'black')

    max_value : float, optional
        Max score. (Default: inf)

    min_value : float, optional
        Min score. (Default: -inf)

    """

    COLOR = "#1f78b4"

    DEFAULT_PROPERTIES = {
        'color': "bed_rgb",
        'border_color': "#1f78b4",
        'min_score': '-inf',
        'max_score': 'inf',
        'bed_type': None,
    }

    def __init__(self, file, **kwargs):
        properties = BedBase.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs
        })
        super().__init__(properties)
        self.bgz_file = build_bed_index(file)

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            BED interval table. The table should be in format like:

            bed_fields = ['chromosome', 'start', 'end',
                          'name', 'score', 'strand',
                          'thick_start', 'thick_end',
                          'rgb', 'block_count',
                          'block_sizes', 'block_starts']
            The table can be in bed6/bed9/bed12 format and the trailing columns can be omited.

        """
        return self.fetch_intervals(self.bgz_file, gr)

    def init_colormap(self):
        self.colormap = None
        if not matplotlib.colors.is_color_like(self.properties['color']) and self.properties['color'] != 'bed_rgb':
            if self.properties['color'] not in matplotlib.cm.datad:
                log.debug("*WARNING* color: '{}' for Track {} is not valid. Color has "
                          "been set to {}".format(self.properties['color'], self.properties['name'],
                                                    self.COLOR))
                self.properties['color'] = self.COLOR
            else:
                self.colormap = self.properties['color']

    def set_colormap(self, df):
        """As min_score and max_score change every plot, we compute them for every plot"""
        props = self.properties
        min_score, max_score = props['min_score'], props['max_score']
        has_score_col = props['bed_type'] in ('bed6', 'bed9', 'bed12')
        if has_score_col and (df.shape[0] > 0):
            min_score = (min_score != 'inf') or df['score'].min()
            max_score = (max_score != '-inf') or df['score'].max()
        min_score, max_score = float(min_score), float(max_score)
        # set colormap
        if self.colormap is not None:
            norm = matplotlib.colors.Normalize(vmin=min_score, vmax=max_score)
            cmap = matplotlib.cm.get_cmap(props['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        if props['color'] == 'bed_rgb' and props['bed_type'] not in ['bed12', 'bed9']:
            log.debug("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                      "been set to {}".format(self.COLOR))
            self.properties['color'] = self.COLOR
            self.colormap = None

    def get_rgb_and_edge_color(self, bed):
        # TODO need simplification
        props = self.properties
        rgb = props['color']
        edgecolor = props['border_color']

        if self.colormap:
            # translate value field (in the example above is 0 or 0.2686...) into a color
            rgb = self.colormap.to_rgba(bed.score)

        # for tad coverage
        if props.get('border_only', 'no') == 'yes':
            rgb = 'none'
        elif props['color'] == 'bed_rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if props['bed_type'] in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                    if 'border_color' in props:
                        edgecolor = props['border_color']
                    else:
                        edgecolor = props['color']
                except IndexError:
                    rgb = self.COLOR
            else:
                rgb = self.COLOR
        return rgb, edgecolor

    @staticmethod
    def infer_bed_type(df: pd.DataFrame) -> Union[str, None]:
        #  bed_type of dataframe are store in dataframe's __dict__ in FetchBed.fetch_intervals
        if 'bed_type' in df.__dict__:
            bed_type = df.bed_type
        else:
            bed_types = {
                12: 'bed12',
                9: 'bed9',
                6: 'bed6',
                3: 'bed3'
            }
            num_col = len(df.columns)
            bed_type = bed_types[num_col] if num_col in bed_types else 'bed3'
            if bed_type == 'bed3' and num_col < 3:
                raise ValueError(f"Invalid dataframe for bed3 with columns: {df.columns}")
        return bed_type

