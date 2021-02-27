import os
import os.path as osp
import subprocess as subp
from typing import Tuple
import runpy

import fire
import nbformat as nbf

import coolbox
from coolbox.utilities import get_logger
from coolbox.utilities.genome import GenomeRange
from coolbox.api import *

log = get_logger("CoolBox CLI")


def get_element_type_by_str(elem_str):
    import coolbox.api
    try:
        elem_tp = eval("coolbox.api." + elem_str)
    except NameError:
        log.error(
            f"No element type name as {elem_tp}, all elements type see: " +
            "https://gangcaolab.github.io/CoolBox/api.html"
        )
    return elem_tp


def get_compose_code(elem_str, args, kwargs):
    compose_code = elem_str + "("
    compose_code += ", ".join(
        [repr(arg) for arg in args] +
        [f"{k} = {repr(v)}" for k, v in kwargs.items()]
    )
    compose_code += ")"
    return compose_code


FRAME_POS = ("left", "right", "top", "bottom", "center")


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

    def __init__(self, genome="hg19", genome_range: str = None, genome_range2: str = None):
        self._indent = 0
        self.current_range = [None, None]
        self.goto(genome_range, genome_range2)
        self.genome = genome
        self.frame_pos = None
        self.frames = {}

    @staticmethod
    def version():
        """print coolbox version"""
        print(coolbox.__version__)

    def set_genome(self, genome):
        """Set reference genome for browser object.

        :param genome: Reference genome (hg19, hg38, mm9, mm10),
        or the path to chromosomes size file(tab splited chromname-size).
        """
        self.genome = genome
        return self

    def goto(self, genome_range: str = None, genome_range2: str = None):
        """Goto a genome range.

        :param genome_range: Genome range string, like "chr9:4000000-6000000".
        :param genome_range2: Genome range string, like "chr9:4000000-6000000. Only required in JointView mode".
        """
        if genome_range2 is None:
            genome_range2 = genome_range
        # validate genome range
        if genome_range:
            gr = GenomeRange(genome_range)
        if genome_range2:
            gr = GenomeRange(genome_range2)
        log.info(f"Goto genome range '{genome_range}' '{genome_range2}'")
        self.current_range = [genome_range, genome_range2]

        return self

    @staticmethod
    def show_doc(elem_str):
        """Print the document of specified Element type. For example: coolbox show_doc Cool"""
        elem_tp = get_element_type_by_str(elem_str)
        print(elem_tp.__doc__)

    def joint_view(self, frame_pos: str = "top"):
        """Start a new frame positioned at the specified frame_pos in the final joint view.
        The center frame should be a single Cool, HicMat, DotHic track.

        For example:
        coolbox - \
        joint_view top - \
            add BigWig BW_PATH - \
        joint_view right - \
            add GTF BTF_PATH - \
        joint_view center - \
            add Cool - \
        goto 'chr1:100000-200000' 'chr2:200000-30000' - \
        plot /tmp/test_joint_view.svg

        Parameters
        ----------
        frame_pos: str
        Add a frame in the given position in the joint view.
        Should be one of 'top', 'left', 'center', 'bottom', 'right'.

        Returns
        -------

        """
        if frame_pos not in FRAME_POS:
            raise ValueError(f"Frame position should be one of {FRAME_POS}")
        if frame_pos not in self.frames:
            self.frames[frame_pos] = "frame = Frame()\n"
        self.frame_pos = frame_pos

        return self

    def add(self, elem_str, *args, **kwargs):
        """Add a Element(Track, Coverage, Feature), for example: coolbox add XAxis

        :param elem_str: Element type string. Like BAM, BigWig, Cool ...
        Full list of Track types
        can be found here(https://gangcaolab.github.io/CoolBox/quick_start_API.html#Track-types).
        :param args: Positional args for create elements.
        :param kwargs: Keyword args for create elements.
        """
        if ("help" in args) or ("help" in kwargs):
            self.show_doc(elem_str)
            return
        if self.frame_pos is None:
            self.joint_view("top")

        compose_code = get_compose_code(elem_str, args, kwargs)
        log.info(f"Create element, compose code: {compose_code}")
        self.frames[self.frame_pos] += "    " * self._indent + "frame += " + compose_code + "\n"
        return self

    def start_with(self, elem_str, *args, **kwargs):
        """Start a 'with' block, apply the element to all elements within the block.

        :param elem_str: Element type string. Like VLines, Color, MinValue ...
        :param args: Positional args for create elements.
        :param kwargs: Keyword args for create elements.
        """
        if ("help" in args) or ("help" in kwargs):
            self.show_doc(elem_str)
            return
        if self.frame_pos is None:
            self.joint_view("top")

        compose_code = get_compose_code(elem_str, args, kwargs)
        log.info(f"Create a with block, compose code: {compose_code}")
        self.frames[self.frame_pos] += "    " * self._indent + f"with {compose_code}:\n"
        self._indent += 1
        return self

    def end_with(self):
        """End the with block"""
        if self._indent <= 0:
            raise ValueError("Expect a 'start_with' command before end_with clause.")
        self._indent -= 1
        return self

    @staticmethod
    def _fetch_frame_src(pos: str, source: str) -> Tuple[str, str]:
        """

        Parameters
        ----------
        pos
        source

        Returns
        -------
        Tuple[str, str].
        The first string is the frame's variable name, and the second string is the function's source code
        """
        # center frame return the first track. expected a Cool, DotHiC, HicMat
        if pos != 'center':
            source += "return frame\n"
        else:
            source += "return list(frame.tracks.values())[0]\n"
        source = "\n".join("    " + line for line in source.split("\n")) + "\n"
        # generate function and call
        frame_var = f"{pos}_frame"
        source = f"def fetch_{frame_var}():\n" + source
        source += f"{frame_var} = fetch_{frame_var}()\n"

        return frame_var, source

    def source(self) -> str:
        num_frames = len(self.frames)
        if num_frames == 0:
            raise RuntimeError("No frame has been added yet.")

        if num_frames > 1 and 'center' not in self.frames:
            raise RuntimeError("JointView mode needs a center frame. "
                               "Use `joint_view center ADD xxx` to add center track/frame.")
        gr1, gr2 = self.current_range
        if gr1 is None:
            raise RuntimeError("No genome range found."
                               "Use `goto chr1:5000000-6000000` to set the genome range.")

        frame_dict = {}
        source = ""
        for pos, src in self.frames.items():
            frame_var, frame_src = self._fetch_frame_src(pos, src)
            source += frame_src
            frame_dict[pos] = frame_var

        if 'center' in self.frames:
            source += f"frame = JointView({frame_dict['center']}, " \
                      f"left={frame_dict.get('left')}, " \
                      f"right={frame_dict.get('right')}, " \
                      f"bottom={frame_dict.get('bottom')}, " \
                      f"top={frame_dict.get('top')}" \
                      f")\n"
        else:
            source += f"frame = {list(frame_dict.values())[0]}\n"

        return source

    def print_source(self):
        """Print the browser composing code."""
        print(self.source())
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
                "import coolbox.api\n" +
                "from coolbox.api import *\n" +
                self.source() +
                f"bsr = Browser(frame, reference_genome='{self.genome}')\n" +
                (f"bsr.goto('{str(self.current_range[0])}')\n" if self.current_range[0] else "") +
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
        subp.check_call(f"jupyter notebook {tmp} " + jupyter_args, shell=True)

    def plot(self, fig_path, genome_range=None, genome_range2=None):
        """Draw a figure within a genome range and save to file

        :param fig_path: Figure save path
        :param genome_range: Genome range string, like "chr9:4000000-6000000".
        :param genome_range2: Genome range string, like "chr9:4000000-6000000".
        """
        if genome_range is None:
            if self.current_range[0] is None:
                raise ValueError("Should specify the gr")
        else:
            self.goto(genome_range, genome_range2)
        source = "from coolbox.api import *\n" + self.source() + "\n"
        gr1, gr2 = self.current_range
        if 'center' in self.frames:
            source += f"fig = frame.plot('{gr1}', '{gr2}')\n"
            # TODO: svgutils.compose.Figure can only save to svg, convert it?
            if not fig_path.endswith('.svg'):
                fig_path = fig_path[:-4] + '.svg'
                log.warning(f"The JointView only support save to svg now. Save fig to: {fig_path}.")
            source += f"fig.save('{fig_path}')\n"
        else:
            source += f"fig = frame.plot('{gr1}')\n"
            source += f"fig.savefig('{fig_path}')\n"
        try:
            code = compile(source, "coolbox_cli_source", "exec")
            eval(code)
        except Exception as e:
            log.error(
                "Error when execute the generated source code:\n\n" +
                "------------------------\n" +
                source + "\n" +
                "------------------------\n\n"
            )
            if type(e) == NameError:
                log.error(
                    f"All elements type see: " +
                    "https://gangcaolab.github.io/CoolBox/api.html"
                )
            raise e

    def run_webapp(self, voila_args="--Voila.ip=0.0.0.0"):
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

    def load_module(self, mod_str):
        """
        Import custom tracks from a module/package for example:

            $ coolbox - load_module ./my_custom.py - add XAxis - add CustomTrack - goto "chr1:5000000-6000000" - run_webapp

        :param mod_str: Path to the module.
        """
        globals().update(runpy.run_path(mod_str, init_globals=globals()))
        return self


if __name__ == "__main__":
    fire.Fire(CLI)

