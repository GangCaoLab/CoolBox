import logging

from IPython.display import display

from coolbox.utilities import (
    GenomeRange, GenomeLength, BUILT_IN_GENOMES, get_size, fig2bytes
)
from .widgets import SimpleWidgets, FullWidgets

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class Browser(object):
    """
    Browser base class.
    """

    MAX_CACHE_SIZE = 10 ** 7  # max fig cache size.(bytes)

    def __init__(self, frame, reference_genome='hg19',
                 init_range=None, widgets_box='simple',
                 dpi=None, img_format='svg'):
        """
        Parameters
        ----------
        frame : {coolbox.core.frame.Frame, coolbox.core.superframe.base.Superframe}
            Browser's main frame.

        reference_genome : str, optional
            Reference genome,
            built-in references:('hg19', 'hg38', 'mm9', 'mm10')
            if you want use other genome, you can specify the "chromosome length file",
            that is a tab splited file, first column is the chromosomes,
            and second column is the length of correspond chromosome. ['hg19']

        init_range : str, optional
            Initial browser range.

        widgets_box : {'simple', 'full'}, optional
            WidgetsBox sub class, default SimpleWidgets

        dpi : int, optional
            The dpi of frame's image.

        img_format : str, optional
            Frame image format, default svg.
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

        if widgets_box == 'simple':
            self.widgets = SimpleWidgets(self)
        elif widgets_box == 'full':
            self.widgets = FullWidgets(self)
        else:
            raise NotImplementedError("widgets type {} not support, please use 'simple' or 'full'".format(widgets_box))

        self.goto(self.current_range)
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

        default_length = 10 ** 7

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

    def goto(self, genome_range, who=None):
        if isinstance(genome_range, str):
            genome_range = GenomeRange(genome_range)
        if not self.chrom_lengthes.check_range(genome_range):
            log.warning("The genome range {} is not valid.".format(genome_range))
            return
        self.current_range = genome_range
        frame_range = GenomeRange(genome_range.chrom,
                                  genome_range.start - 1,  # NOTE: frame's start is zero based
                                  genome_range.end)
        self.frame.goto(frame_range)

        self.widgets.refresh_widgets(who=who)

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
            if self.img_format == 'svg':
                fig_bytes = fig_bytes.decode("utf-8")
            self.fig_cache[self.current_range] = fig_bytes

        # auto clear fig cache for prevent memory leak
        if len(self.fig_cache) > 20 and \
                get_size(self.fig_cache) >= Browser.MAX_CACHE_SIZE:
            self.clear_fig_cache()

        self.widgets.frame_widget.value = fig_bytes

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

    @property
    def tracks(self):
        return self.frame.tracks

    def fetch_data(self, genome_range=None):
        return self.frame.fetch_data(genome_range)
