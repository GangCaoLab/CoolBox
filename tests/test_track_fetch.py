import os.path as osp
from coolbox.api import *


HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = "chr9:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')


def test_xaxis():
    x_axis = XAxis()
    assert x_axis


def test_bigwig():
    bw = BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw")
    assert bw.fetch_data(test_interval) is not None


def test_cool():
    cl = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    assert cl.fetch_data(test_interval) is not None
    assert cl.fetch_data(test_interval, test_interval) is not None


def test_bed():
    bed = Bed(f"{DATA_DIR}/bed_{test_itv}.bed")
    assert bed.fetch_data(test_interval) is not None


def test_arcs():
    arcs = Arcs(f"{DATA_DIR}/arcs_{test_itv}.arcs")
    assert arcs.fetch_data(test_interval) is not None
