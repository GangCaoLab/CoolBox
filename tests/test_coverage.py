import os.path as osp
from coolbox.api import *

HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = "chr9:4000000-6000000"
empty_interval = "chr10:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')


def test_vlines():
    vlines = Vlines([("chr9", 4000000),
                     "chr9:5000000-6000000",
                     ("chr1", 4000000)])
    with vlines:
        frame = XAxis() + XAxis()
    fig = frame.plot(test_interval)
    fig.savefig("/tmp/test_coolbox_vline.pdf")


def test_highlights():
    highlights = HighLights([
        ("chr9", 4000000, 5000000),
        "chr9:5000000-6000000",
    ])
    with highlights:
        frame = XAxis() + XAxis()
    fig = frame.plot(test_interval)
    fig.savefig("/tmp/test_coolbox_highlights.pdf")


def test_bigwig_coverage():
    frame = XAxis() +\
        BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw", style="line:1") + \
        BigWigCoverage(f"{DATA_DIR}/bigwig_{test_itv}.bw", style="fill", color="red", alpha=0.3) +\
        BigWigCoverage(f"{DATA_DIR}/bigwig_{test_itv}.bw", style="fill", color="blue", alpha=0.3)
    fig = frame.plot(test_interval)
    fig.savefig("/tmp/test_coolbox_bwcov.pdf")


def test_arcs_coverage():
    frame = XAxis() + \
        BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw", style="fill", alpha=0.5) + \
        ArcsCoverage(f"{DATA_DIR}/bedpe_{test_itv}.bedpe")
    fig = frame.plot(test_interval)
    fig.savefig("/tmp/test_coolbox_arcscov.pdf")


def test_tad_coverage():
    frame = XAxis() + \
        Cool(f"{DATA_DIR}/cool_{test_itv}.mcool") + \
        TADCoverage(f"{DATA_DIR}/tad_{test_itv}.bed")
    fig = frame.plot(test_interval)
    fig.savefig("/tmp/test_coobox_tadcov.pdf")
