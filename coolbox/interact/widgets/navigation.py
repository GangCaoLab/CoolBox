from collections import OrderedDict

from ipywidgets import (
    HBox, VBox, Dropdown, Button, Label,
    Checkbox, FloatText,
    Text, IntRangeSlider, Layout, HTML
)

from coolbox.utilities import (
    GenomeRange
)

ALL_BW_MARK = "ALL(BW/BG)"


class NavigationBar(object):
    """
    Layout:

    -----------------------------------------------------------------------------------------
    left_button | right_button | zoom_out_button | zoom_in_button | range_textbox | go_button
    -----------------------------------------------------------------------------------------
    range_min_label | range_slider | range_max_label
    --------------------------------------------------------------------
    auto_check_box | track_min_val_float_text | track_max_val_float_text
    --------------------------------------------------------------------
    frame
    """

    def __init__(self, browser, frame=None):
        chromsomes = list(browser.chrom_lengthes.keys())
        self.widgets = self.__get_widgets(chromsomes, browser, frame)
        self.panel = self.__compose_panel(self.widgets)
        self.selected_tracks = self.__get_tracks_name(browser)

    def __get_tracks_name(self, browser):
        return [t.name for t in browser.tracks.values() if 'min_value' in t.properties]

    def refresh_widgets(self, browser, who=None):
        range_ = browser.current_range
        max_ = browser.chrom_lengthes[range_.chrom]
        # update range_textbox
        if who != "go_button":
            self.widgets['range_textbox'].value = str(range_)
        # update range_slider
        if who != "range_slider":
            slider = self.widgets['range_slider']
            slider.min = 1
            slider.max = max_
            slider.value = (range_.start, range_.end)
        # update chromosome list
        if who != "chromosomes_list":
            self.widgets['chromosomes_list'].value = range_.chrom
        # update range max min label
        self.widgets['range_min_label'].value = str(1)
        self.widgets['range_max_label'].value = str(max_)

    def register_events_handler(self, browser):
        """
        Parameters
        ----------
        browser : `BrowserBase`
            Browser object.
        """

        # chromosomes_list value change
        def chrom_dropdown_val_change(change):
            new_chrom = change['new']
            current_range = browser.current_range
            # only change chromosome
            range_ = GenomeRange(new_chrom, current_range.start, current_range.end)
            range_ = browser.chrom_lengthes.bound_range(range_)
            browser.goto(range_, who='chromosomes_list')
            browser.refresh()

        self.widgets['chromosomes_list'].observe(chrom_dropdown_val_change, names="value")

        # left_button click
        def left_button_click(b):
            browser.go_left()
            browser.refresh()

        #            browser.preload_imgs('left')
        self.widgets['left_button'].on_click(left_button_click)

        # right_button click
        def right_button_click(b):
            browser.go_right()
            browser.refresh()

        #            browser.preload_imgs('right')
        self.widgets['right_button'].on_click(right_button_click)

        # zoom_in_button click
        def zoom_in_button_click(b):
            browser.zoom_in()
            browser.refresh()

        #            browser.preload_imgs('zoom-in')
        self.widgets['zoom_in_button'].on_click(zoom_in_button_click)

        # zoom_out_button click
        def zoom_out_button_click(b):
            browser.zoom_out()
            browser.refresh()

        #            browser.preload_imgs('zoom-out')
        self.widgets['zoom_out_button'].on_click(zoom_out_button_click)

        # go_button click
        def go_button_click(b):
            range_str = self.widgets['range_textbox'].value.strip("'")
            range_ = GenomeRange(range_str)
            browser.goto(range_, who='go_button')
            browser.refresh()

        self.widgets['go_button'].on_click(go_button_click)

        # range_slider value change
        def range_slider_val_change(change):
            start_old, end_old = change['old']
            length_old = end_old - start_old

            start, end = change['new']
            chrom = browser.current_range.chrom
            if end - start <= 0:
                end = start + length_old
            new_range = GenomeRange(chrom, start, end)
            new_range = browser.chrom_lengthes.bound_range(new_range)
            browser.goto(new_range, who='range_slider')
            browser.refresh()

        self.widgets['range_slider'].observe(range_slider_val_change, names="value")

        # auto_check_box value change
        def auto_check_box_val_change(change):
            if change['new'] == True:
                self.widgets['track_min_val_float_text'].disabled = True
                self.widgets['track_max_val_float_text'].disabled = True
                self.widgets['track_dropdown'].disabled = True
                browser.frame.set_tracks_min_max('auto', 'auto')
            else:
                self.widgets['track_min_val_float_text'].disabled = False
                self.widgets['track_max_val_float_text'].disabled = False
                self.widgets['track_dropdown'].disabled = False
                min_ = self.widgets['track_min_val_float_text'].value
                max_ = self.widgets['track_max_val_float_text'].value
                browser.frame.set_tracks_min_max(min_, max_)
            browser.clear_fig_cache()
            browser.refresh()

        self.widgets['auto_check_box'].observe(auto_check_box_val_change, names="value")

        # track_max_value_float_text value change and
        #   track_min_value_float_text value change
        def track_float_text_val_change(change):
            min_ = self.widgets['track_min_val_float_text'].value
            max_ = self.widgets['track_max_val_float_text'].value
            if self.widgets['track_dropdown'].value == ALL_BW_MARK:
                browser.frame.set_tracks_min_max(min_, max_)
            else:
                for tname in self.selected_tracks:
                    browser.frame.set_tracks_min_max(min_, max_, tname)
            browser.clear_fig_cache()
            browser.refresh()

        self.widgets['track_min_val_float_text'].observe(track_float_text_val_change, names="value")
        self.widgets['track_max_val_float_text'].observe(track_float_text_val_change, names="value")

        def track_dropdown_change(change):
            if change['new'] == ALL_BW_MARK:
                self.selected_tracks = self.__get_tracks_name(browser)
            else:
                self.selected_tracks = [change['new']]

        self.widgets['track_dropdown'].observe(track_dropdown_change, names="value")

    def __get_widgets(self, chromosomes, browser, frame=None):
        if frame is None:
            frame = HTML()
        tracks = self.__get_tracks_name(browser)
        widgets = OrderedDict([
            ("chromosomes_list", Dropdown(options=chromosomes)),
            ("left_button", Button(icon="arrow-left")),
            ("right_button", Button(icon="arrow-right")),
            ("zoom_out_button", Button(icon="search-minus")),
            ("zoom_in_button", Button(icon="search-plus")),
            ("range_textbox", Text(placeholder="genome range like: 'chr1:10000-20000'")),
            ("go_button", Button(description="Go")),

            ("range_slider", IntRangeSlider(continuous_update=False, readout=False, layout=Layout(width='90%'))),
            ("range_min_label", Label("", layout=Layout(width='2%'))),
            ("range_max_label", Label("", layout=Layout(width='20%'))),

            ("auto_check_box", Checkbox(value=True, description="Auto Range",
                                        layout=Layout(width='120px'),
                                        style={'description_width': 'initial'})),
            ("track_min_val_float_text",
             FloatText(value=0.0001, description="Track's min value:", step=0.5, disabled=True,
                       layout=Layout(width='30%'),
                       style={'description_width': 'initial'})),
            ("track_max_val_float_text", FloatText(value=10, description="Track's max value:", step=0.5, disabled=True,
                                                   layout=Layout(width='30%'),
                                                   style={'description_width': 'initial'})),
            ("track_dropdown", Dropdown(options=[ALL_BW_MARK] + tracks,
                                        value=ALL_BW_MARK,
                                        description="Select track",
                                        disabled=True,
                                        )),
            ("frame", frame)
        ])
        return widgets

    def __compose_panel(self, widgets_dict):
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
                    widgets_dict["track_dropdown"],
                ], layout=Layout(justify_content="flex-start")),
            ], layout=Layout(border='solid 2px')),
            widgets_dict["frame"],
        ])
        return panel
