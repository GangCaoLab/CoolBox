import os.path as osp

from coolbox.api import *

HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = "chr9:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')
empty_interval = "chr1:4000000-6000000"


def test_frame():
    frame = XAxis() + \
        Cool(f"{DATA_DIR}/cool_{test_itv}.mcool") + \
        Spacer(1) + \
        Arcs(f"{DATA_DIR}/arcs_{test_itv}.arcs") + Inverted() + \
        GTF(f"{DATA_DIR}/gtf_{test_itv}.gtf") + TrackHeight(7) + \
        Spacer(1) + \
        BAM(f"{DATA_DIR}/bam_{test_itv}.bam") +\
        BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw") + \
        Spacer(1) + \
        Bed(f"{DATA_DIR}/bed_{test_itv}.bed") + TrackHeight(10)
    fig = frame.plot(test_interval)
    tmp = "/tmp/test_coolbox_fig.pdf"
    fig.savefig(tmp)

    fig = frame.plot(empty_interval)
    tmp = "/tmp/test_coolbox_fig_empty.pdf"
    fig.savefig(tmp)

