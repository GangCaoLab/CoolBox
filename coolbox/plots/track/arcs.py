import numpy as np

from coolbox.utilities import (
    change_chrom_names,
    Interval, IntervalTree,
    get_logger
)

from coolbox.plots.track.base import TrackPlot


log = get_logger(__name__)


class PlotArcs(TrackPlot):

    def __init__(self, *args, **kwarg):
        TrackPlot.__init__(self, *args, **kwarg)
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is a score.

        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#') or self.__is_header(line):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2, score, *other = line.strip().split('\t')
                except Exception as detail:
                    msg = 'File not valid. The format is chrom1 start1, end1, ' \
                          'chrom2, start2, end2, score\nError: {}\n in line\n {}'.format(detail, line)
                    raise IOError(msg)

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    raise IOError(msg)

                assert start1 <= end1, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 larger than start2 in {}".format(line_number, line)
                try:
                    score = float(score)
                except ValueError as detail:
                    msg = "Error reading line: {}. The score is not valid {}. " \
                          "\nError message: {}".format(line_number, score, detail)
                    raise IOError(msg)

                if chrom1 != chrom2:
                    log.warning("Only links in same chromosome are used. Skipping line\n{}\n".format(line))
                    continue

                if chrom1 not in interval_tree:
                    interval_tree[chrom1] = IntervalTree()

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, score))
                valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

    def plot(self, ax, label_ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        """
        from matplotlib.patches import Arc
        height = 1
        max_diameter = 0
        count = 0
        if chrom_region not in list(self.interval_tree):
            chrom_region = change_chrom_names(chrom_region)
        arcs_in_region = sorted(self.interval_tree[chrom_region][region_start:region_end])

        for idx, interval in enumerate(arcs_in_region):
            # skip arcs whose start and end are outside the plotted region
            if interval.begin < region_start and interval.end > region_end:
                continue

            if 'line_width' in self.properties:
                line_width = float(self.properties['line_width'])
            else:
                line_width = 0.5 * np.sqrt(interval.data)

            diameter = (interval.end - interval.begin)
            center = (interval.begin + interval.end) / 2
            if diameter > max_diameter:
                max_diameter = diameter
            count += 1
            ax.plot([center], [diameter])
            ax.add_patch(Arc((center, 0), diameter,
                             height*2, 0, 0, 180, color=self.properties['color'], lw=line_width))

        # increase max_diameter slightly to avoid cropping of the arcs.
        #max_diameter += max_diameter * 0.05
        height += height * 0.05
        log.debug("{} were arcs plotted".format(count))
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(height, 0.001)
        else:
            ax.set_ylim(-0.001, height)

        ax.set_xlim(region_start, region_end)
        log.debug('title is {}'.format(self.properties['title']))
        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')

    def __is_header(self, line):
        fields = line.split()
        for idx in [1, 2, 4, 5]:
            if not fields[idx].isdigit():
                return True
        return False

