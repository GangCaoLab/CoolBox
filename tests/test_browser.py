import os.path as osp

from coolbox.api import *

HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = "chr9:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')
empty_interval = "chr1:4000000-6000000"


def test_browser():
    frame = XAxis() + \
            Cool(f"{DATA_DIR}/cool_{test_itv}.mcool") + \
            Spacer(1) + \
            Arcs(f"{DATA_DIR}/bedpe_{test_itv}.bedpe") + Inverted() + \
            GTF(f"{DATA_DIR}/gtf_{test_itv}.gtf") + TrackHeight(7) + \
            Spacer(1) + \
            BAM(f"{DATA_DIR}/bam_{test_itv}.bam") + \
            BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw") + \
            Spacer(1) + \
            BED(f"{DATA_DIR}/bed_{test_itv}.bed") + TrackHeight(10)
    bsr = Browser(frame)
    bsr.goto(test_interval)
    tmp = "/tmp/test_coolbox_bsr.pdf"
    bsr.save(tmp)
    bsr.goto(empty_interval)
    tmp = "/tmp/test_coolbox_bsr_empty.pdf"
    bsr.save(tmp)


