import os
import os.path as osp
import pandas as pd

from coolbox.api import *
from coolbox.utilities import GenomeRange
import matplotlib.pyplot as plt

HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = GenomeRange("chr9:4000000-6000000")
empty_interval = GenomeRange("chr10:4000000-6000000")
sub_interval1 = GenomeRange("chr9:4500000-5000000")
sub_interval2 = GenomeRange("chr9:5200000-5850000")
test_itv = str(test_interval).replace(':', '_').replace('-', '_')


def test_xaxis():
    x_axis = XAxis()
    assert x_axis


def test_bigwig():
    bw = BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw")
    intervals = bw.fetch_data(test_interval)
    assert isinstance(intervals, pd.DataFrame)
    fig, ax = plt.subplots()
    bw.plot(ax, test_interval)
    bw.fetch_data(empty_interval)


def test_bed():
    path = f"{DATA_DIR}/bed_{test_itv}.bed"
    bed = BED(path)
    df = bed.fetch_data(test_interval)
    count = 0
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                count += 1
    assert df.shape[0] == count
    fig, ax = plt.subplots()
    bed.plot(ax, test_interval)
    bed.fetch_data(empty_interval)
    bed6 = BED(f"{DATA_DIR}/bed6_{test_itv}.bed")
    bed6.plot(ax, test_interval)
    assert bed6.properties['bed_type'] == 'bed6'
    bed9 = BED(f"{DATA_DIR}/bed9_{test_itv}.bed")
    bed9.plot(ax, test_interval)
    assert bed9.properties['bed_type'] == 'bed9'


def test_bedpe():
    arcs = BEDPE(f"{DATA_DIR}/bedpe_{test_itv}.bedpe")
    assert arcs.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    arcs.plot(ax, test_interval)
    arcs.fetch_data(empty_interval)


def test_pairs():
    arcs = Pairs(f"{DATA_DIR}/pairs_{test_itv}.pairs")
    assert arcs.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    arcs.plot(ax, test_interval)
    arcs.fetch_data(empty_interval)


def test_arcs():
    arcs = Arcs(f"{DATA_DIR}/bedpe_{test_itv}.bedpe")
    assert isinstance(arcs, BEDPE)
    arcs = Arcs(f"{DATA_DIR}/pairs_{test_itv}.pairs")
    assert isinstance(arcs, Pairs)


def test_gtf():
    gtf = GTF(f"{DATA_DIR}/gtf_{test_itv}.gtf")
    assert gtf.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    gtf.plot(ax, test_interval)
    gtf.fetch_data(empty_interval)
    gtf = GTF(f"{DATA_DIR}/gtf_{test_itv}_fake.gtf")
    assert gtf.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    gtf.plot(ax, test_interval)
    gtf.fetch_data(empty_interval)


def test_bam():
    bam_path = f"{DATA_DIR}/bam_{test_itv}.bam"
    bai_path = bam_path + ".bai"
    if osp.exists(bai_path):
        os.remove(bai_path)
    bam = BAM(bam_path, plot_type="alignment")
    assert bam.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    small_interval = GenomeRange("chr9:4998535-5006343")
    bam.plot(ax, small_interval)
    fig.savefig("/tmp/test_coolbox_bam.pdf")
    bam = BAM(bam_path, plot_type="coverage")
    fig, ax = plt.subplots()
    bam.plot(ax, test_interval)
    fig.savefig("/tmp/test_coolbox_bam_typecov.pdf")
    bam.fetch_data(empty_interval)


def test_bedgraph():
    bg_path = f"{DATA_DIR}/bedgraph_{test_itv}.bg"
    bg = BedGraph(bg_path)
    assert bg.fetch_data(test_interval) is not None
    bg.fetch_data(empty_interval)
    fig, ax = plt.subplots()
    bg.plot(ax, test_interval)
    fig.savefig("/tmp/test_coolbox_bg.pdf")


def test_tads():
    cool = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    tad = cool + TADCoverage(f"{DATA_DIR}/tad_{test_itv}.bed")
    assert tad.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    tad.plot(ax, test_interval)
    tad.fetch_data(empty_interval)
    tad = TAD(f"{DATA_DIR}/tad_{test_itv}.bed")
    tad.plot(ax, test_interval)


def test_hicdiff():
    cl1 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    cl2 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    diff = HiCDiff(cl1, cl2)
    assert diff.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    diff.plot(ax, test_interval)
    diff.fetch_data(empty_interval)


def test_selfish():
    cl1 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    cl2 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    sel = Selfish(cl1, cl2)
    assert sel.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    sel.plot(ax, test_interval)
    sel.fetch_data(empty_interval)


def test_cool():
    cl = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    exp_binsize = cl.infer_binsize(test_interval)
    assert cl.fetch_data(test_interval) is not None
    assert cl.fetch_data(test_interval, gr2=test_interval) is not None
    assert exp_binsize == cl.fetched_binsize
    fig, ax = plt.subplots()
    cl.plot(ax, test_interval)
    cl.fetch_data(empty_interval)
    exp_binsize = 10000
    cl2 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool", resolution=exp_binsize)
    cl2.fetch_data(test_interval)
    assert cl2.fetched_binsize == exp_binsize
    mat1 = cl.fetch_data(sub_interval1, gr2=sub_interval2)
    mat2 = cl.fetch_data(sub_interval2, gr2=sub_interval1)
    assert mat1.shape[0] != mat1.shape[1]
    assert mat2.shape[0] != mat2.shape[1]
    assert mat1.shape == (mat2.shape[1], mat2.shape[0])
    assert mat2[mat2>0].shape[0] > 0
    fig, ax = plt.subplots()
    cl.plot(ax, sub_interval1, gr2=sub_interval2)


def test_dothic():
    dothic_path = f"{DATA_DIR}/dothic_{test_itv}.hic"
    dot = DotHiC(dothic_path)
    exp_binsize = dot.infer_binsize(test_interval)
    assert dot.fetch_data(test_interval) is not None
    assert dot.fetch_data(test_interval, gr2=test_interval) is not None
    assert exp_binsize == dot.fetched_binsize
    fig, ax = plt.subplots()
    dot.plot(ax, test_interval)
    dot.fetch_data(empty_interval)
    exp_binsize = 10000
    dot = DotHiC(dothic_path, resolution=exp_binsize)
    dot.fetch_data(test_interval)
    assert dot.fetched_binsize == exp_binsize
    mat1 = dot.fetch_data(sub_interval1, gr2=sub_interval2)
    mat2 = dot.fetch_data(sub_interval2, gr2=sub_interval1)
    assert mat1.shape[0] != mat1.shape[1]
    assert mat2.shape[0] != mat2.shape[1]
    assert mat1.shape == (mat2.shape[1], mat2.shape[0])
    assert mat2[mat2>0].shape[0] > 0
    fig, ax = plt.subplots()
    dot.plot(ax, sub_interval1, gr2=sub_interval2)


def test_hicfeatures():
    cool_path = f"{DATA_DIR}/cool_{test_itv}.mcool"
    cl = Cool(cool_path)

    # test fot di_score and insu_score
    # compose from cool
    insu = InsuScore(cl)
    di = DiScore(cl)
    assert di.fetch_data(test_interval).shape == insu.fetch_data(test_interval).shape
    # compose from path
    insu = InsuScore(cool_path)
    di = DiScore(cool_path)
    assert di.fetch_data(test_interval).shape == insu.fetch_data(test_interval).shape
    # test for insu_score with mutlple window_size
    insu = InsuScore(cl, window_size="20-40")
    insus = insu.fetch_data(test_interval)
    assert len(insus.shape) == 2 and insus.shape[0] == 40 - 20

    # teset for virtual4c
    fig, ax = plt.subplots()
    cl.plot(ax, test_interval)
    chr_, _other = str(test_interval).split(":")
    s, e = _other.split("-")
    s, e = int(s), int(e)
    mid = (s + e) // 2
    mid_point = f"{chr_}:{mid}-{mid}"
    # compose from cool
    v4c = Virtual4C(cl, mid_point)
    assert v4c.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    v4c.plot(ax, test_interval)
    v4c.fetch_data(empty_interval)
    # compose from path
    v4c = Virtual4C(cool_path, mid_point)
    assert v4c.fetch_data(test_interval) is not None


if __name__ == "__main__":
    # test_cool()
    test_dothic()
    # test_xaxis()
    # test_gtf()
    # test_bam()
    # test_bedgraph()
    # test_arcs()
    # test_tads()
    # test_bed()
