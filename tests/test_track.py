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


def test_cool():
    cl = Cool(f"{DATA_DIR}/cool_{test_itv}.mcool")
    assert cl.fetch_data(test_interval) is not None
    assert cl.fetch_data(test_interval, test_interval) is not None
    fig, ax = plt.subplots()
    cl.plot_genome_range(ax, test_interval)
    cl.fetch_data(empty_interval)


def test_bed():
    bed = Bed(f"{DATA_DIR}/bed_{test_itv}.bed")
    assert bed.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    bed.plot_genome_range(ax, test_interval)
    bed.fetch_data(empty_interval)


def test_arcs():
    arcs = Arcs(f"{DATA_DIR}/arcs_{test_itv}.arcs")
    assert arcs.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    arcs.plot_genome_range(ax, test_interval)
    arcs.fetch_data(empty_interval)


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
    bam = BAM(bam_path, style="alignment")
    assert bam.fetch_data(test_interval) is not None
    fig, ax = plt.subplots()
    small_interval = "chr9:4998535-5006343"
    bam.plot_genome_range(ax, small_interval)
    fig.savefig("/tmp/test_coolbox_bam.pdf")
    bam = BAM(bam_path, style="coverage")
    fig, ax = plt.subplots()
    bam.plot_genome_range(ax, test_interval)
    fig.savefig("/tmp/test_coolbox_bam_stylecov.pdf")
    bam.fetch_data(empty_interval)


if __name__ == "__main__":
    #test_xaxis()
    #test_gtf()
    test_bam()

