import abc


class WidgetsBox(abc.ABC):
    """
    Widgets panel base class.
    """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        browser : BrowserBase
            Browser object.

        frame_widget : ipywidgets.Image
            frame widget object.
        """
        self.browser = args[0]
        self.frame_widget = args[1]

        self.widgets_dict = self.get_widgets_dict()
        self.refresh_widgets()
        self.panel = self.compose_panel(self.widgets_dict)
        self.register_events_handler()

    def get_widgets_dict(self):
        pass

    def refresh_widgets(self):
        pass

    def compose_panel(self, widgets_dict):
        pass

    def register_events_handler(self):
        pass
