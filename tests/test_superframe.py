import os.path as osp

from coolbox.api import *


HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = "chr9:4000000-6000000"
empty_interval = "chr10:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')


def test_joint_view():
    frame1 = XAxis() + GTF(f"{DATA_DIR}/gtf_{test_itv}.gtf") + Title("GTF")
    frame1 += BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw") + TrackHeight(2) + MinValue(0)
    frame2 = XAxis()
    frame2 += GTF(f"{DATA_DIR}/gtf_{test_itv}.gtf") + TrackHeight(5)

    cool1 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool", color_bar=False)

    jv = JointView(cool1, top=frame1, right=frame2, bottom=frame1, left=frame2, space=0, padding_left=0)
    fig = jv.plot("chr9:4750000-5100000", "chr9:5350000-5700000")
    fig.save("/tmp/test_coolbox_joint_view.svg")
