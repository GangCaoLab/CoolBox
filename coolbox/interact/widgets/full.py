from .base import *

from ipywidgets import (Tab)


def compose_track_config_panel(widgets_panel):
    return []


class FullWidgets(WidgetsBox):

    def get_widgets_dict(self):
        frame = Image()
        navigation_bar = compose_navigation_bar_panel(frame)
        track_config = OrderedDict([

        ])
        widgets_dict = OrderedDict([
            ("navigation_bar",  navigation_bar),
            ("track_config",  track_config),
        ])
        return widgets_dict

    def compose_panel(self, widgets_dict):
        # compose navigation_bar
        navigation_bar = compose_navigation_bar_panel(widgets_dict['navigation_bar'])
        track_config = compose_track_config_panel(widgets_dict['track_config'])
        panel = Tab()
        panel.children = [navigation_bar, track_config]
        panel.set_title(0, "Navigation")
        panel.set_title(1, "Tracks")
        return panel
