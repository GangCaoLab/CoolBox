import os.path as osp
import subprocess as subp


HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"
test_interval = "chr9:4000000-6000000"
empty_interval = "chr10:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')


def test_cli_plot():
    cmd = [
        "python", "-m", "coolbox.cli",
        "add", "XAxis", "-",
        "add", "BigWig", f"{DATA_DIR}/bigwig_{test_itv}.bw", "-",
        "add", "BedGraph", f"{DATA_DIR}/bedgraph_{test_itv}.bg", "-",
        "add", "GTF", f"{DATA_DIR}/gtf_{test_itv}.gtf", "-",
        "goto", test_interval, "-",
        "plot", "/tmp/test_coolbox.pdf",
    ]
    subp.check_call(cmd)

    cmd = f"""
        python -m coolbox.cli
          joint_view top - 
            add XAxis - 
            add GTF {DATA_DIR}/gtf_{test_itv}.gtf - 
            add Title GTF - 
          joint_view right - 
            add XAxis - 
            add BigWig {DATA_DIR}/bigwig_{test_itv}.bw - 
            add TrackHeight 2 - 
            add MinValue 0 - 
          joint_view center - 
            add Cool {DATA_DIR}/cool_{test_itv}.mcool  - 
          goto 'chr9:4500000-5000000' 'chr9:5200000-5850000' - 
          plot /tmp/test_coolbox_joint_view.svg
        """

    subp.check_call(cmd.replace("\n", ""), shell=True)


def test_cli_gen_notebook():
    cmd = [
        "python", "-m", "coolbox.cli",
        "add", "XAxis", "-",
        "add", "BigWig", f"{DATA_DIR}/bigwig_{test_itv}.bw", "-",
        "add", "BedGraph", f"{DATA_DIR}/bedgraph_{test_itv}.bg", "-",
        "add", "GTF", f"{DATA_DIR}/gtf_{test_itv}.gtf", "-",
        "goto", test_interval, "-",
        "gen_notebook", "/tmp/test_coolbox.ipynb",
    ]
    subp.check_call(cmd)


def test_cli_import_custom():
    custom_path = osp.join(HERE, "custom_track.py")
    cmd = [
        "python", "-m", "coolbox.cli",
        "load_module", custom_path, "-",
        "add", "XAxis", "-",
        "add", "CustomTrack", "-",
        "goto", test_interval, "-",
        "plot", "/tmp/test_coolbox_custom.pdf"
    ]
    subp.check_call(cmd)
