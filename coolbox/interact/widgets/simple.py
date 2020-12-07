import logging

from ipywidgets import HTML

from .base import WidgetsBox
from .navigation import NavigationBar

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class SimpleWidgets(WidgetsBox):
    """
    Simplest widgets panel design.

    -----

    """

    def __init__(self, browser, *args, **kwargs):
        frame_widget = HTML()
        self.navigation_bar = NavigationBar(browser, frame_widget)
        super().__init__(browser, frame_widget, *args, **kwargs)

    def get_widgets_dict(self):
        widgets = self.navigation_bar.widgets
        return widgets

    def refresh_widgets(self, who=None):
        self.navigation_bar.refresh_widgets(self.browser, who)

    def compose_panel(self, widgets_dict):
        panel = self.navigation_bar.panel
        return panel

    def register_events_handler(self):
        self.navigation_bar.register_events_handler(self.browser)
