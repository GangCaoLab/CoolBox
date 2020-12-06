"""
format convert
"""

import collections
from .filetool import opener, to_string

_fields_rg = ("bin", "name", "chrom", "strand", "txStart", "txEnd",
              "cdsStart", "cdsEnd", "exonCount", "exonStart", "exonEnds",
              "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

_refGeneRec = collections.namedtuple("refGeneRec", _fields_rg)


class refGeneRec(_refGeneRec):

    def to_bed12_line(self):
        chrom = self.chrom
        start = self.txStart
        end = self.txEnd
        name = self.name2
        score = self.score
        strand = self.strand

        thick_start = start
        thick_end = end

        item_rgb = "0"

        block_count = self.exonCount
        block_starts = self.offset_zero(self.exonStart, start) + ","
        block_ends = self.exonEnds
        block_sizes = self.get_exons_size(self.exonStart, block_ends) + ","

        bed12_item = (chrom, start, end, name, score, strand, thick_start, thick_end, item_rgb,
                      block_count, block_sizes, block_starts)

        bed12_line = "\t".join(bed12_item)
        return bed12_line

    def to_line(self):
        return "\t".join(self)

    @staticmethod
    def offset_zero(block, offset):
        offset = int(offset)
        block = [int(i) for i in block.split(",") if i]
        block = [str(i - offset) for i in block]
        return ",".join(block)

    @staticmethod
    def get_exons_size(exons_start, exons_end):
        starts = exons_start.split(",")
        starts = [int(i) for i in starts if i]
        ends = exons_end.split(",")
        ends = [int(i) for i in ends if i]
        sizes = [ends[i] - starts[i] for i, _ in enumerate(starts)]
        exons_size = ",".join([str(i) for i in sizes])
        return exons_size


def refgene_txt_to_bed12(txt_file, bed_file):
    with opener(txt_file) as f_in, open(bed_file, 'w') as f_out:
        for line in f_in:
            line = to_string(line)
            items = line.strip().split("\t")
            refg_rec = refGeneRec._make(items)
            out_line = refg_rec.to_bed12_line() + "\n"
            f_out.write(out_line)
