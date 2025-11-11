from typing import Union

import pandas as pd
import matplotlib

from coolbox.utilities import get_logger
from coolbox.utilities.reader.tab import get_indexed_tab_reader
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
        Track color, 'rgb' for auto specify color according to bed record.
        (Default: 'rgb')

    border_color : str, optional
        Border_color of gene. (Default: 'black')

    max_value : float, optional
        Max score. (Default: inf)

    min_value : float, optional
        Min score. (Default: -inf)

    """

    COLOR = "#1f78b4"

    DEFAULT_PROPERTIES = {
        'color': "rgb",
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
        self.reader = get_indexed_tab_reader(file)

    def fetch_data(self, gr: GenomeRange, **kwargs) -> pd.DataFrame:
        """

        Returns
        -------
        intervals : pandas.core.frame.DataFrame
            BED interval table. The table should be in format like::

                bed_fields = ['chromosome', 'start', 'end',
                              'name', 'score', 'strand',
                              'thick_start', 'thick_end',
                              'rgb', 'block_count',
                              'block_sizes', 'block_starts']

            The table can be in bed6/bed9/bed12 format and the trailing columns can be omited.

        """
        return self.fetch_intervals(gr)

    def init_colormap(self):
        self.colormap = None
        if not matplotlib.colors.is_color_like(self.properties['color']) and self.properties['color'] != 'rgb':
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
        if props['color'] == 'rgb' and props['bed_type'] not in ['bed12', 'bed9']:
            log.debug("*WARNING* Color set to 'rgb', but bed file does not have the rgb field. The color has "
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
        elif props['color'] == 'rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if props['bed_type'] in ['bed9', 'bed12']:
                try:
                    if isinstance(bed.rgb, str):
                        rgb = [int(c.strip()) for c in bed.rgb.split(',') if c]
                        assert len(rgb) == 3, "rgb must be a list or tuple of length 3"
                    else:
                        assert isinstance(bed.rgb, (list, tuple)), "rgb must be a list or tuple"
                        assert len(bed.rgb) == 3, "rgb must be a list or tuple of length 3"
                        assert all(0 <= x <= 255 for x in bed.rgb), "rgb must be a list or tuple of integers between 0 and 255"
                        rgb = bed.rgb
                    rgb = [float(x) / 255 for x in rgb]
                    if 'border_color' in props:
                        edgecolor = props['border_color']
                    else:
                        edgecolor = props['color']
                except (IndexError, AssertionError):
                    rgb = self.COLOR
            else:
                rgb = self.COLOR
        return rgb, edgecolor

    @staticmethod
    def infer_bed_type(df: pd.DataFrame) -> Union[str, None]:
        #  bed_type of dataframe are store in dataframe's __dict__ in FetchBed.fetch_intervals
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

    def fetch_intervals(self, gr: GenomeRange) -> pd.DataFrame:
        """
        Fetch intervals within input chromosome range.
        """
        df = self.reader.query_var_chr(gr)
        return df
