import os.path as osp
import subprocess as subp

HERE = osp.dirname(osp.abspath(__file__))
test_interval = "chr9:4000000-6000000"
empty_interval = "chr10:4000000-6000000"
test_itv = test_interval.replace(':', '_').replace('-', '_')


def test_cli_plot(data_dir, tmp_dir):
    cmd = [
        "python", "-m", "coolbox.cli",
        "add", "XAxis", "-",
        "add", "BigWig", f"{data_dir}/bigwig_{test_itv}.bw", "-",
        "add", "BedGraph", f"{data_dir}/bedgraph_{test_itv}.bg", "-",
        "add", "GTF", f"{data_dir}/gtf_{test_itv}.gtf", "-",
        "goto", test_interval, "-",
        "plot", "/tmp/test_coolbox.pdf",
    ]
    subp.check_call(cmd)

    cmd = f"""
        python -m coolbox.cli
          joint_view top - 
            add XAxis - 
            add GTF {data_dir}/gtf_{test_itv}.gtf - 
            add Title GTF - 
          joint_view right - 
            add XAxis - 
            add BigWig {data_dir}/bigwig_{test_itv}.bw - 
            add TrackHeight 2 - 
            add MinValue 0 - 
          joint_view center - 
            add Cool {data_dir}/cool_{test_itv}.mcool  - 
          goto 'chr9:4500000-5000000' 'chr9:5200000-5850000' - 
          plot {tmp_dir}/test_coolbox_joint_view.svg
        """
    cmd = cmd.replace("\n", "")
    print(cmd)
    subp.check_call(cmd, shell=True)


def test_cli_gen_notebook(data_dir, tmp_dir):
    cmd = [
        "python", "-m", "coolbox.cli",
        "add", "XAxis", "-",
        "add", "BigWig", f"{data_dir}/bigwig_{test_itv}.bw", "-",
        "add", "BedGraph", f"{data_dir}/bedgraph_{test_itv}.bg", "-",
        "add", "GTF", f"{data_dir}/gtf_{test_itv}.gtf", "-",
        "goto", test_interval, "-",
        "gen_notebook", f"{tmp_dir}/test_coolbox.ipynb",
    ]
    cmd = " ".join(cmd)
    print(cmd)
    subp.check_call(cmd, shell=True)


def test_cli_import_custom(tmp_dir):
    custom_path = osp.join(HERE, "custom_track.py")
    cmd = [
        "python", "-m", "coolbox.cli",
        "load_module", custom_path, "-",
        "add", "XAxis", "-",
        "add", "CustomTrack", "-",
        "goto", test_interval, "-",
        "plot", f"{tmp_dir}/test_coolbox_custom.pdf"
    ]
    cmd = " ".join(cmd)
    print(cmd)
    subp.check_call(cmd, shell=True)
