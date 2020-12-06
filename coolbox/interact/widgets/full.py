from collections import OrderedDict

from ipywidgets import (Tab, VBox, HTML)

from .base import WidgetsBox
from .navigation import NavigationBar


def compose_track_config_panel(widgets_panel):
    return VBox([])


class FullWidgets(WidgetsBox):

    def __init__(self, browser, *args, **kwargs):
        chromosomes = browser.chrom_lengthes.keys()
        frame_widget = HTML()
        self.navigation_bar = NavigationBar(chromosomes, frame_widget)
        super().__init__(browser, frame_widget, *args, **kwargs)

    def get_widgets_dict(self):
        navigation_bar = self.navigation_bar.widgets
        track_config = OrderedDict([

        ])
        widgets_dict = OrderedDict([
            ("navigation_bar", navigation_bar),
            ("track_config", track_config),
        ])
        return widgets_dict

    def compose_panel(self, widgets_dict):
        # compose navigation_bar
        navigation_bar = self.navigation_bar.panel
        track_config = compose_track_config_panel(widgets_dict['track_config'])
        panel = Tab()
        panel.children = [navigation_bar, track_config]
        panel.set_title(0, "Navigation")
        panel.set_title(1, "Tracks")
        return panel

    def register_events_handler(self):
        self.navigation_bar.register_events_handler(self.browser)

    def refresh_widgets(self, who=None):
        self.navigation_bar.refresh_widgets(self.browser, who)
