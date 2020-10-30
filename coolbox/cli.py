import os
import os.path as osp
import subprocess as subp

import fire
import nbformat as nbf

from coolbox.api import *
from coolbox.utilities import get_logger


log = get_logger("CoolBox CLI")


class CLI(object):
    """
    CoolBox Command Line Interface

    You can use this cli to create coolbox browser instance,
    visualize your data directly in shell.

    example:

    1. Draw tracks within a genome range, save figure to a pdf file:

        $ coolbox add XAxis - add BigWig test.bw - goto "chr1:5000000-6000000" - plot test.pdf

    2. Generate a notebook and run jupyter to open browser:

        $ coolbox add XAxis - add BigWig test.bw - goto "chr1:5000000-6000000" - run_jupyter

    3. Run a independent web application.

        $ coolbox add XAxis - add BigWig test.bw - goto "chr1:5000000-6000000" - run_webapp

    """

    def __init__(self, genome="hg19", genome_range=None):
        self.source = "Frame()"
        self.genome = genome
        if genome_range and isinstance(genome_range, str):
            self.frame = Frame(genome_range=genome_range)
        else:
            self.frame = Frame()

    def set_genome(self, genome):
        """Set reference genome for browser object.

        :param genome: Reference genome (hg19, hg38, mm9, mm10),
        or the path to chromosomes size file(tab splited chromname-size).
        """
        self.genome = genome
        return self

    def goto(self, genome_range):
        """Goto a genome range.

        :param genome_range: Genome range string, like "chr9:4000000-6000000".
        """
        log.info(f"Goto genome range '{genome_range}'")
        self.frame.goto(genome_range)
        return self

    @property
    def current_range(self):
        return self.frame.current_range

    def add(self, elem_str, *args, **kwargs):
        """Add a Element(Track, Coverage, Feature)

        :param elem_str: Element type string. Like BAM, BigWig, Cool ...
        :param args: Positional args for create elements.
        :param kwargs: Keyword args for create elements.
        """
        elem_tp = eval(elem_str)
        elem = elem_tp(*args, **kwargs)
        self.frame = self.frame + elem
        compose_code = elem_str + "("
        compose_code += ", ".join(
            [repr(arg) for arg in args] +
            [f"{k} = {repr(v)}" for k, v in kwargs.items()]
        )
        compose_code += ")"
        log.info(f"Create element, compose code: {compose_code}")
        self.source += " + " + compose_code
        return self

    def print_source(self):
        """Print the browser composing code."""
        print(self.source)
        return self

    def gen_notebook(self, notebook_path, notes=True, figsave=True):
        """Generate The notebook contain codes for run coolbox browser.

        :param notebook_path: The output notebook path.
        :param notes: Generate markdown notes in notebook or not.
        :param figsave: Generate codes for saving figure or not.
        """
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
        """Create a notebook according to command line, then start a jupyter process.

        :param jupyter_args: Arguments for run jupyter.
        """
        i = 0
        tmp_notebook = lambda: f"/tmp/coolbox_tmp.{i}.ipynb"
        while osp.exists(tmp_notebook()):
            i += 1
        tmp = tmp_notebook()
        self.gen_notebook(tmp)
        subp.check_call(f"jupyter notebook {tmp} "+jupyter_args, shell=True)

    def plot(self, fig_path, genome_range=None):
        """Draw a figure within a genome range and save to file

        :param fig_path: Figure save path
        :param genome_range: Genome range string, like "chr9:4000000-6000000".
        """
        if genome_range is None:
            if self.current_range is None:
                raise ValueError("Should specify the genome_range")
        else:
            self.goto(genome_range)
        fig = self.frame.plot(str(self.current_range))
        fig.savefig(fig_path)

    def run_webapp(self, voila_args=""):
        """Run a independent coolbox browser web app.
        (Create notebook and run voila)

        :param voila_args: Arguments for run jupyter.
        """
        i = 0
        tmp_notebook = lambda: f"/tmp/coolbox_tmp.{i}.ipynb"
        while osp.exists(tmp_notebook()):
            i += 1
        tmp = tmp_notebook()
        self.gen_notebook(tmp, notes=False, figsave=False)
        subp.check_call(f"voila {tmp} " + voila_args, shell=True)

    def end(self):
        """Terminate the CLI pipeline"""
        pass


if __name__ == "__main__":
    fire.Fire(CLI)

