import os
import os.path as osp
import subprocess as subp

import fire
import nbformat as nbf

from coolbox.api import *


class CLI(object):

    def __init__(self, genome="hg19", genome_range=None):
        self.source = "Frame()"
        self.genome = genome
        if genome_range and isinstance(genome_range, str):
            self.frame = Frame(genome_range=genome_range)
        else:
            self.frame = Frame()

    def set_genome(self, genome):
        self.genome = genome
        return self

    def goto(self, genome_range):
        self.frame.goto(genome_range)
        return self

    @property
    def current_range(self):
        return self.frame.current_range

    def add(self, elem_str, *args, **kwargs):
        elem_tp = eval(elem_str)
        elem = elem_tp(*args, **kwargs)
        self.frame = self.frame + elem
        source = " + " + elem_str + "("
        source += ", ".join(
            [repr(arg) for arg in args] +
            [f"{k} = {repr(v)}" for k, v in kwargs.items()]
        )
        source += ")"
        self.source += source
        return self

    def print_source(self):
        print(self.source)
        return self

    def gen_notebook(self, notebook_path, notes=True, figsave=True):
        code = nbf.v4.new_code_cell
        markdown = nbf.v4.new_markdown_cell
        nb = nbf.v4.new_notebook()
        cells = []
        if notes:
            cells.append(
                markdown("**Run the cell bellow to start the CoolBox browser.**")
            )
        cells.append(
            code(
                f"import os; os.chdir('{os.getcwd()}')\n"
                "from coolbox.api import *\n" +
                "frame = " + self.source + "\n" +
                f"bsr = Browser(frame, reference_genome='{self.genome}')\n" +
                (f"bsr.goto('{str(self.current_range)}')\n" if self.current_range else "") +
                "bsr.show()"
            ),
        )
        if figsave:
            if notes:
                cells.append(
                    markdown("**Run the cell bellow to save your figure.**"),
                )
            cells.append(
                code("bsr.save('test.pdf')"),
            )
        nb['cells'] = cells
        nbf.write(nb, notebook_path)
        return self

    def run_jupyter(self, jupyter_args="--ip=0.0.0.0"):
        i = 0
        tmp_notebook = lambda: f"/tmp/coolbox_tmp.{i}.ipynb"
        while osp.exists(tmp_notebook()):
            i += 1
        tmp = tmp_notebook()
        self.gen_notebook(tmp)
        subp.check_call(f"jupyter notebook {tmp} "+jupyter_args, shell=True)

    def plot(self, genome_range, fig_path):
        self.goto(genome_range)
        fig = self.frame.plot(str(self.current_range))
        fig.savefig(fig_path)

    def end(self):
        pass


fire.Fire(CLI)

