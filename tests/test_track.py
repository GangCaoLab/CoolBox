import os
import os.path as osp
from coolbox.api import *
import matplotlib.pyplot as plt

HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = "chr9:4000000-6000000"
empty_interval = "chr10:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')


def test_xaxis():
    x_axis = XAxis()
    assert x_axis


def test_bigwig():
    bw = BigWig(f"{DATA_DIR}/bigwig_{test_itv}.bw")
    assert bw.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    bw.plot_genome_range(ax, test_interval)
    bw.fetch_data(empty_interval)


def test_bed():
    bed = BED(f"{DATA_DIR}/bed_{test_itv}.bed")
    assert bed.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    bed.plot_genome_range(ax, test_interval)
    bed.fetch_data(empty_interval)


def test_bedpe():
    arcs = BEDPE(f"{DATA_DIR}/bedpe_{test_itv}.bedpe")
    assert arcs.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    arcs.plot_genome_range(ax, test_interval)
    arcs.fetch_data(empty_interval)


def test_pairs():
    arcs = Pairs(f"{DATA_DIR}/pairs_{test_itv}.pairs")
    assert arcs.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    arcs.plot_genome_range(ax, test_interval)
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
    gtf.plot_genome_range(ax, test_interval)
    gtf.fetch_data(empty_interval)


def test_bam():
    bam_path = f"{DATA_DIR}/bam_{test_itv}.bam"
    bai_path = bam_path + ".bai"
    if osp.exists(bai_path):
        os.remove(bai_path)
    bam = BAM(bam_path, plot_type="alignment")
    assert bam.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    small_interval = "chr9:4998535-5006343"
    bam.plot_genome_range(ax, small_interval)
    fig.savefig("/tmp/test_coolbox_bam.pdf")
    bam = BAM(bam_path, plot_type="coverage")
    fig, ax = plt.subplots()
    bam.plot_genome_range(ax, test_interval)
    fig.savefig("/tmp/test_coolbox_bam_typecov.pdf")
    bam.fetch_data(empty_interval)


def test_bedgraph():
    bg_path = f"{DATA_DIR}/bedgraph_{test_itv}.bg"
    bg = BedGraph(bg_path)
    assert bg.fetch_data(test_interval) is not None
    bg.fetch_data(empty_interval)
    fig, ax = plt.subplots()
    bg.plot_genome_range(ax, test_interval)
    fig.savefig("/tmp/test_coolbox_bg.pdf")


def test_tads():
    tad = TADs(f"{DATA_DIR}/tad_{test_itv}.bed")
    assert tad.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    tad.plot_genome_range(ax, test_interval)
    tad.fetch_data(empty_interval)


def test_hicdiff():
    cl1 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    cl2 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    diff = HiCDiff(cl1, cl2)
    assert diff.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    diff.plot_genome_range(ax, test_interval)
    diff.fetch_data(empty_interval)


def test_selfish():
    cl1 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    cl2 = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    sel = Selfish(cl1, cl2)
    assert sel.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    sel.plot_genome_range(ax, test_interval)
    sel.fetch_data(empty_interval)


def test_cool():
    cl = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    exp_binsize = cl.infer_binsize(test_interval)
    assert cl.fetch_data(test_interval) is not None
    assert cl.fetch_data(test_interval, test_interval) is not None
    assert exp_binsize == cl.fetched_binsize
    fig, ax = plt.subplots()
    cl.plot_genome_range(ax, test_interval)
    cl.fetch_data(empty_interval)


def test_dothic():
    # .hic file is not easy to create test subset,
    # so check if symbol link is exists or not,
    # to decide whether test it.
    dothic_path = f"{DATA_DIR}/test.hic"
    if not osp.exists(dothic_path):
        return
    dot = DotHiC(dothic_path)
    exp_binsize = dot.infer_binsize(test_interval)
    assert dot.fetch_data(test_interval) is not None
    assert dot.fetch_data(test_interval, test_interval) is not None
    assert exp_binsize == dot.fetched_binsize
    fig, ax = plt.subplots()
    dot.plot_genome_range(ax, test_interval)
    dot.fetch_data(empty_interval)


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
    cl.plot_genome_range(ax, test_interval)
    chr_, _other = test_interval.split(":")
    s, e = _other.split("-")
    s, e = int(s), int(e)
    mid = (s + e) // 2
    mid_point = f"{chr_}:{mid}-{mid}"
    # compose from cool
    v4c = Virtual4C(cl, mid_point)
    assert v4c.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    v4c.plot_genome_range(ax, test_interval)
    v4c.fetch_data(empty_interval)
    # compose from path
    v4c = Virtual4C(cool_path, mid_point)
    assert v4c.fetch_data(test_interval) is not None


if __name__ == "__main__":
    # test_xaxis()
    # test_gtf()
    # test_bam()
    # test_bedgraph()
    # test_arcs()
    test_tads()
