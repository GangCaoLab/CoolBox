from IPython.display import display

from ipywidgets import (
    HBox, VBox, Label, Dropdown, Button, Label,
    Checkbox, FloatText,
    Text, Image, IntRangeSlider, Layout
)

from .utilities import (
    GenomeRange, GenomeLength, BUILT_IN_GENOMES, get_size, fig2bytes
)

from collections import OrderedDict

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class WidgetsBox(object):
    """
    Widgets panel base class.
    """
    pass


class SimpleWidgets(WidgetsBox):
    """
    Simplest widgets panel design.

    Layout:

    -----------------------------------------------------------------------------------------
    left_button | right_button | zoom_out_button | zoom_in_button | range_textbox | go_button
    -----------------------------------------------------------------------------------------
    range_min_label | range_slider | range_max_label
    --------------------------------------------------------------------
    auto_check_box | track_min_val_float_text | track_max_val_float_text
    --------------------------------------------------------------------
    frame
    -----

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
        widgets = OrderedDict([
            ("chromosomes_list",          Dropdown(options=self.browser.chrom_lengthes.keys())),
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

            ("frame",                     Image())
        ])
        return widgets

    def refresh_widgets(self, who=None):
        range_ = self.browser.current_range
        max_ = self.browser.chrom_lengthes[range_.chrom]
        # update range_textbox
        if who != "go_button":
            self.widgets_dict['range_textbox'].value = str(range_)
        # update range_slider
        if who != "range_slider":
            slider = self.widgets_dict['range_slider']
            slider.min = 1
            slider.max = max_
            slider.value = (range_.start, range_.end)
        # update chromosome list
        if who != "chromosomes_list":
            self.widgets_dict['chromosomes_list'].value = range_.chrom
        # update range max min label
        self.widgets_dict['range_min_label'].value = str(1)
        self.widgets_dict['range_max_label'].value = str(max_)

    def compose_panel(self, widgets_dict):
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

    def register_events_handler(self):

        # chromosomes_list value change
        def chrom_dropdown_val_change(change):
            new_chrom = change['new']
            current_range = self.browser.current_range
            # only change chromosome
            range_ = GenomeRange(new_chrom, current_range.start, current_range.end)
            range_ = self.browser.chrom_lengthes.bound_range(range_)
            self.browser.goto(range_)
            self.refresh_widgets(who="chromosomes_list")
            self.browser.refresh()
        self.widgets_dict['chromosomes_list'].observe(chrom_dropdown_val_change, names="value")

        # left_button click
        def left_button_click(b):
            self.browser.go_left()
            self.refresh_widgets()
            self.browser.refresh()
#            self.browser.preload_imgs('left')
        self.widgets_dict['left_button'].on_click(left_button_click)

        # right_button click
        def right_button_click(b):
            self.browser.go_right()
            self.refresh_widgets()
            self.browser.refresh()
#            self.browser.preload_imgs('right')
        self.widgets_dict['right_button'].on_click(right_button_click)

        # zoom_in_button click
        def zoom_in_button_click(b):
            self.browser.zoom_in()
            self.refresh_widgets()
            self.browser.refresh()
#            self.browser.preload_imgs('zoom-in')
        self.widgets_dict['zoom_in_button'].on_click(zoom_in_button_click)

        # zoom_out_button click
        def zoom_out_button_click(b):
            self.browser.zoom_out()
            self.refresh_widgets()
            self.browser.refresh()
#            self.browser.preload_imgs('zoom-out')
        self.widgets_dict['zoom_out_button'].on_click(zoom_out_button_click)

        # go_button click
        def go_button_click(b):
            range_str = self.widgets_dict['range_textbox'].value.strip("'")
            range_ = GenomeRange(range_str)
            self.browser.goto(range_)
            self.refresh_widgets(who="go_button")
            self.browser.refresh()
        self.widgets_dict['go_button'].on_click(go_button_click)

        # range_slider value change
        def range_slider_val_change(change):
            start, end = change['new']
            chrom = self.browser.current_range.chrom
            new_range = GenomeRange(chrom, start, end)
            self.browser.goto(new_range)
            self.refresh_widgets(who="range_slider")
            self.browser.refresh()
        self.widgets_dict['range_slider'].observe(range_slider_val_change, names="value")

        # auto_check_box value change
        def auto_check_box_val_change(change):
            if change['new'] == True:
                self.widgets_dict['track_min_val_float_text'].disabled = True
                self.widgets_dict['track_max_val_float_text'].disabled = True
                self.browser.frame.set_tracks_min_max('auto', 'auto')
            else:
                self.widgets_dict['track_min_val_float_text'].disabled = False
                self.widgets_dict['track_max_val_float_text'].disabled = False
                min_ = self.widgets_dict['track_min_val_float_text'].value
                max_ = self.widgets_dict['track_max_val_float_text'].value
                self.browser.frame.set_tracks_min_max(min_, max_)
            self.browser.clear_fig_cache()
            self.browser.refresh()
        self.widgets_dict['auto_check_box'].observe(auto_check_box_val_change, names="value")

        # track_max_value_float_text value change and
        #   track_min_value_float_text value change
        def track_float_text_val_change(change):
            min_ = self.widgets_dict['track_min_val_float_text'].value
            max_ = self.widgets_dict['track_max_val_float_text'].value
            self.browser.frame.set_tracks_min_max(min_, max_)
            self.browser.clear_fig_cache()
            self.browser.refresh()
        self.widgets_dict['track_min_val_float_text'].observe(track_float_text_val_change, names="value")
        self.widgets_dict['track_max_val_float_text'].observe(track_float_text_val_change, names="value")


class BrowserBase(object):
    """
    Browser base class.
    """

    MAX_CACHE_SIZE = 10_000_000 # max fig cache size.(bytes)

    def __init__(self, frame, reference_genome='hg19',
                 init_range=None, widgets_box=SimpleWidgets,
                 dpi=None, img_format='png'):
        """
        Parameters
        ----------
        frame : coolbox.core.Frame
            Browser's main frame.

        reference_genome : str, optional
            Reference genome,
            built-in references:('hg19', 'hg38', 'mm9', 'mm10')
            if you want use other genome, you can specify the "chromosome length file",
            that is a tab splited file, first column is the chromosomes,
            and second column is the length of correspond chromosome. ['hg19']

        init_range : str, optional
            Initial browser range.

        widgets_box : type, optional
            WidgetsBox sub class, default SimpleWidgets

        dpi : int, optional
            The dpi of frame's image.

        img_format : int, optional
            Frame image format, default png.
        """

        self.dpi = dpi
        self.img_format = img_format
        self.frame = frame

        if reference_genome in BUILT_IN_GENOMES:
            self.chrom_lengthes = BUILT_IN_GENOMES[reference_genome]
        else:
            self.chrom_lengthes = GenomeLength(reference_genome)
            if len(self.chrom_lengthes) == 0:
                raise IOError("chromosome lengthes file is not include any useful information."
                              "Please check file \"{}\".".format(reference_genome))

        if init_range is not None:
            self.current_range = GenomeRange(init_range)
        else:
            self.current_range = self.get_init_range()

        self.goto(self.current_range)
        self.widgets = widgets_box(self)
        self.fig = None

        # cache figs in dict, speed up the figure display process.
        #   key: genome range
        #   value: fig image bytes
        self.fig_cache = {}

    def get_init_range(self, chrom=None):
        """
        Generate an initial range within a chromosome.

        Args:
            chrom (str, optional): initial choromosome.

        Return:
            (:obj:`GenomeRange`)
        """

        if chrom is None:
            chrom = list(self.chrom_lengthes.keys())[0]

        default_length = 10_000_000

        if self.chrom_lengthes[chrom] > default_length:
            range_ = GenomeRange(chrom, 1, default_length)
        else:
            range_ = GenomeRange(chrom, 1, self.chrom_lengthes[chrom])

        return range_

    @property
    def window_size(self):
        size = self.current_range.end - self.current_range.start
        return size

    @property
    def center(self):
        center = (self.current_range.start + self.current_range.end) // 2
        return center

    def goto(self, genome_range):
        if isinstance(genome_range, str):
            genome_range = GenomeRange(genome_range)
        if not self.chrom_lengthes.check_range(genome_range):
            log.warning("The genome range {} is not valid.".format(genome_range))
            return
        self.current_range = genome_range
        frame_range = GenomeRange(genome_range.chrom,
                                  genome_range.start - 1, # NOTE: frame's start is zero based
                                  genome_range.end)
        self.frame.goto(frame_range)

    def go_left(self, step_ratio=0.5, dry_run=False):
        window_size = self.window_size
        step = int(window_size * step_ratio)
        start = self.current_range.start - step
        end = self.current_range.end - step
        genome_range = GenomeRange(self.current_range.chrom, start, end)
        genome_range = self.chrom_lengthes.bound_range(genome_range)
        if dry_run:
            return genome_range
        else:
            self.goto(genome_range)
    
    def go_right(self, step_ratio=0.5, dry_run=False):
        window_size = self.window_size
        step = int(window_size * step_ratio)
        start = self.current_range.start + step
        end = self.current_range.end + step
        genome_range = GenomeRange(self.current_range.chrom, start, end)
        genome_range = self.chrom_lengthes.bound_range(genome_range)
        self.goto(genome_range)
        if dry_run:
            return genome_range
        else:
            self.goto(genome_range)

    def zoom_in(self, zoom_ratio=2, dry_run=False):
        window_size = self.window_size
        window_size = window_size // zoom_ratio
        start = self.center - window_size // 2
        end = start + window_size
        genome_range = GenomeRange(self.current_range.chrom, start, end)
        genome_range = self.chrom_lengthes.bound_range(genome_range)
        self.goto(genome_range)
        if dry_run:
            return genome_range
        else:
            self.goto(genome_range)

    def zoom_out(self, zoom_ratio=2, dry_run=False):
        window_size = self.window_size
        window_size = window_size * zoom_ratio
        start = self.center - window_size // 2
        end = start + window_size
        genome_range = GenomeRange(self.current_range.chrom, start, end)
        genome_range = self.chrom_lengthes.bound_range(genome_range)
        self.goto(genome_range)
        if dry_run:
            return genome_range
        else:
            self.goto(genome_range)

    def show(self, *args, **kwargs):
        """
        Show widgets and frame.
        """
        display(self.widgets.panel)
        self.refresh()

    def refresh(self, hard=False):
        """
        Refresh the image display.
        """
        if (self.current_range in self.fig_cache) and (not hard):
            fig_bytes = self.fig_cache[self.current_range]
        else:
            fig_current = self.frame.show()
            fig_bytes = fig2bytes(fig_current, encode=self.img_format, dpi=self.dpi)
            self.fig_cache[self.current_range] = fig_bytes

        # auto clear fig cache for prevent memory leak
        if len(self.fig_cache) > 20 and \
           get_size(self.fig_cache) >= BrowserBase.MAX_CACHE_SIZE:
            self.clear_fig_cache()

        self.widgets.widgets_dict['frame'].value = fig_bytes

    def preload_imgs(self, directions):
        """
        Preloading images to self.fig_cache.

        Can load image in one of 4 directions:
            left, right, zoom-in, zoom-out
        or load all directions.
        """
        def preload(g_range):
            if g_range not in self.fig_cache:
                fig = self.frame.plot(g_range.chrom, g_range.start, g_range.end)
                self.fig_cache[g_range] = fig2bytes(fig, encode=self.img_format, dpi=self.dpi)

        possible_dires = ['left', 'right', 'zoom-in', 'zoom-out']

        if isinstance(directions, str):
            if directions == 'all':
                all_directions = possible_dires
            else:
                all_directions = [directions]
        else:
            assert isinstance(directions, list) or isinstance(directions, tuple), \
                "direction must str or list or tuple object."
            all_directions = directions

        for dire in all_directions:
            assert directions in possible_dires, \
                "direction must one of {}.".format(possible_dires)
            if dire == 'left':
                left_region = self.go_left(dry_run=True)
                preload(left_region)
            elif dire == 'right':
                right_region = self.go_right(dry_run=True)
                preload(right_region)
            elif dire == 'zoom-in':
                zoom_in_region = self.zoom_in(dry_run=True)
                preload(zoom_in_region)
            elif dire == 'zoom-out':
                zoom_out_region = self.zoom_out(dry_run=True)
                preload(zoom_out_region)

    def clear_fig_cache(self):
        """
        Clear the fig cache.
        """
        self.fig_cache = {}

    def save(self, path, dpi=None):
        """
        Save current frame's image to file.
        """
        c_fig = self.frame.show()
        dpi = dpi or self.dpi
        c_fig.savefig(path, dpi=dpi)
