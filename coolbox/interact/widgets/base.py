from ipywidgets import (
    HBox, VBox, Label, Dropdown, Button, Label,
    Checkbox, FloatText,
    Text, Image, IntRangeSlider, Layout
)

from coolbox.utilities import (
    GenomeRange
)

from collections import OrderedDict


class WidgetsBox(object):
    """
    Widgets panel base class.
    """
    def __init__(self, *args, **kwargs):
        """
        Args:
            browser (:obj:`BrowserBase`): Browser object.
        """
        self.browser = args[0]

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


def navigation_bar_widgets(chromosomes, frame=None):
    if frame is None:
        frame = Image()
    widgets = OrderedDict([
        ("chromosomes_list",          Dropdown(options=chromosomes)),
        ("left_button",               Button(icon="arrow-left")),
        ("right_button",              Button(icon="arrow-right")),
        ("zoom_out_button",           Button(icon="search-minus")),
        ("zoom_in_button",            Button(icon="search-plus")),
        ("range_textbox",             Text(placeholder="genome range like: 'chr1:10000-20000'")),
        ("go_button",                 Button(description="Go")),

        ("range_slider",              IntRangeSlider(continuous_update=False, readout=False, layout=Layout(width='90%'))),
        ("range_min_label",           Label("", layout=Layout(width='2%'))),
        ("range_max_label",           Label("", layout=Layout(width='20%'))),

        ("auto_check_box",            Checkbox(value=True, description="Auto Range",
                                               layout=Layout(width='120px'),
                                               style={'description_width': 'initial'})),
        ("track_min_val_float_text",  FloatText(value=0,  description="track's min value:", step=0.5, disabled=True,
                                                layout=Layout(width='30%'),
                                                style={'description_width': 'initial'})),
        ("track_max_val_float_text",  FloatText(value=10, description="track's max value:", step=0.5, disabled=True,
                                                layout=Layout(width='30%'),
                                                style={'description_width': 'initial'})),

        ("frame",                     frame)
    ])
    return widgets


def compose_navigation_bar_panel(widgets_dict):
    panel = VBox([
        VBox([
            HBox(list(widgets_dict.values())[:7]),
            HBox([
                widgets_dict["range_min_label"],
                widgets_dict["range_slider"],
                widgets_dict["range_max_label"],
            ]),
            HBox([
                widgets_dict["auto_check_box"],
                widgets_dict["track_min_val_float_text"],
                widgets_dict["track_max_val_float_text"],
            ], layout=Layout(justify_content="flex-start")),
        ], layout=Layout(border='solid 2px')),
        widgets_dict["frame"],
    ])
    return panel
