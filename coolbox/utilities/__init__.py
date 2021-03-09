import typing as t
from os.path import abspath, dirname, join
from collections import deque

from .bed import *
from .figtools import *
from .filetool import *
from .fmtconvert import *
from .genome import *
from .logtools import *


def op_err_msg(a, b, op='+'):
    """
    Generate the error message of error operand type.
    """
    msg = "unsupported operand type(s) for {}: '{}' and '{}'".format(op, str(type(a)), str(type(b)))
    return msg


#
# read genome length file #

here = dirname(abspath(__file__))

HG19 = GenomeLength(join(here, "../genome/hg19.txt"), genome_name="hg19")
HG38 = GenomeLength(join(here, "../genome/hg38.txt"), genome_name="hg38")
MM9 = GenomeLength(join(here, "../genome/mm9.txt"), genome_name="mm9")
MM10 = GenomeLength(join(here, "../genome/mm10.txt"), genome_name="mm10")

BUILT_IN_GENOMES = {
    'hg19': HG19,
    'hg38': HG38,
    'mm9': MM9,
    'mm10': MM10,
}

FEATURES_STACK_NAME = "__COOLBOX_FEATURE_STACK__"
COVERAGE_STACK_NAME = "__COOLBOX_COVERAGE_STACK__"
IGNORE_TYPE_CHANGE = ['title','name']

def get_feature_stack():
    global_scope = globals()
    global_scope.setdefault(FEATURES_STACK_NAME, deque())
    return global_scope[FEATURES_STACK_NAME]


def get_coverage_stack():
    global_scope = globals()
    global_scope.setdefault(COVERAGE_STACK_NAME, deque())
    return global_scope[COVERAGE_STACK_NAME]


def format_properties(properties):
    """
    Convert bool value to 'yes' or 'no', for compatible with pyGenomeTracks
    Convert string to float if possible
    """
    for key, value in properties.items():
        if isinstance(value, bool):
            properties[key] = 'yes' if value else 'no'
        elif isinstance(value, str):
            if key in IGNORE_TYPE_CHANGE:
                continue
            try:
                float_val = float(value)
                properties[key] = float_val
            except Exception:
                pass
    return properties


def correspond_track(track_or_name: t.Union["Track", str]):
    """
    Search track by it's name.
    For mark corresponding track in a element definition, allow
    pass track name to find related object in CLI mode.
    """
    from coolbox.core.track.base import Track
    if isinstance(track_or_name, str):
        track = Track.tracks[track_or_name]
    else:
        track = track_or_name
    return track


if __name__ == "__main__":
    import doctest

    doctest.testmod()
