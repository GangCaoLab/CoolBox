import abc
from abc import ABC
from typing import Union

from coolbox.utilities import GenomeRange, get_logger
from ..base import Track
from .plot import CoveragePlot

log = get_logger(__name__)


class HistBase(Track, CoveragePlot, ABC):

    def __init__(self, **kwargs):
        properties_dict = {
            'max_value': 'auto',
            'min_value': 'auto',
            'show_data_range': True,
            'data_range_style': 'y-axis',
            'style': 'line',
            'title': '',
        }
        properties_dict.update(kwargs)
        properties_dict['type'] = properties_dict['style']  # change key word
        super().__init__(properties_dict)
        self.genome_range = None

    @abc.abstractmethod
    def fetch_data(self, genome_range: Union[str, GenomeRange]):
        pass

    @abc.abstractmethod
    def fetch_plot_data(self, genome_range: GenomeRange):
        pass

    def plot(self, ax, chrom_region, start_region, end_region):
        log.debug("plotting {}".format(self.properties['file']))

        self.ax = ax
        # if genome range would change?
        self.genome_range = GenomeRange(chrom_region, start_region, end_region)
        plot_data = self.fetch_plot_data(self.genome_range)
        if plot_data is not None:
            if isinstance(plot_data, tuple):
                scores_per_bin, x_values = plot_data
            else:
                scores_per_bin, x_values = plot_data, None
            self.plot_coverage(ax, self.genome_range, scores_per_bin, x_values)
        self.plot_label()
