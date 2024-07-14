from coolbox.api import *


def test_vlines(test_interval, tmp_dir):
    vlines = Vlines([("chr9", 4000000),
                     "chr9:5000000-6000000",
                     ("chr1", 4000000)])
    with vlines:
        frame = XAxis() + XAxis()
    fig = frame.plot(test_interval)
    fig.savefig(f"{tmp_dir}/test_coolbox_vline.pdf")


def test_highlights(test_interval, tmp_dir):
    highlights = HighLights([
        ("chr9", 4000000, 5000000),
        "chr9:5000000-6000000",
    ])
    with highlights:
        frame = XAxis() + XAxis()
    fig = frame.plot(test_interval, close_fig=False)
    fig.savefig(f"{tmp_dir}/test_coolbox_highlights.pdf")


def test_bigwig_coverage(data_dir, test_itv, test_interval, tmp_dir):
    frame = XAxis() +\
        BigWig(f"{data_dir}/bigwig_{test_itv}.bw", style="line:1") + \
        BigWigCoverage(f"{data_dir}/bigwig_{test_itv}.bw", style="fill", color="red", alpha=0.3) +\
        BigWigCoverage(f"{data_dir}/bigwig_{test_itv}.bw", style="fill", color="blue", alpha=0.3)
    fig = frame.plot(test_interval)
    fig.savefig(f"{tmp_dir}/test_coolbox_bwcov.pdf")


def test_arcs_coverage(data_dir, test_itv, test_interval, tmp_dir):
    frame = XAxis() + \
        BigWig(f"{data_dir}/bigwig_{test_itv}.bw", style="fill", alpha=0.5) + \
        ArcsCoverage(f"{data_dir}/bedpe_{test_itv}.bedpe")
    fig = frame.plot(test_interval)
    fig.savefig(f"{tmp_dir}/test_coolbox_arcscov.pdf")


def test_tad_coverage(data_dir, test_itv, test_interval, tmp_dir):
    frame = XAxis() + \
        Cool(f"{data_dir}/cool_{test_itv}.mcool") + \
        TADCoverage(f"{data_dir}/tad_{test_itv}.bed")
    fig = frame.plot(test_interval)
    fig.savefig(f"{tmp_dir}/test_coobox_tadcov.pdf")


if __name__ == '__main__':
    test_highlights(GenomeRange("chr9:3000000-7000000"), "/tmp")
    #test_vlines(GenomeRange("chr9:4000000-6000000"), "/tmp")
