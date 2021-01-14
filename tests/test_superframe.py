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

    sub_frames = {
        "top": frame1,
        "right": frame2,
        "bottom": frame1,
        "left": frame2
    }

    jv = JointView(cool1, **sub_frames, space=0, padding_left=0)
    fig = jv.plot("chr9:4500000-5000000", "chr9:5200000-5850000")
    fig.save("/tmp/test_coolbox_joint_view.svg")

    tracks_data = jv.fetch_data()
    for pos, fr in sub_frames.items():
        assert len(tracks_data[pos]) == len(fr.fetch_data())

    import numpy as np
    gr, gr2 = jv.current_range
    assert np.all(cool1.fetch_data(gr, gr2=gr2) == tracks_data['center'])
