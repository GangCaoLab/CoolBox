import os
import os.path as osp
import pytest
from coolbox.utilities import GenomeRange


HERE = osp.dirname(osp.abspath(__file__))
DATA_DIR = f"{HERE}/test_data"


def clear_bgz_and_indexes():
    exts = ['bgz', 'tbi', 'px2', 'bai']
    files = os.listdir(DATA_DIR)
    for f in files:
        for ext in exts:
            if f.endswith('.'+ext):
                os.remove(f"{DATA_DIR}/{f}")


def pytest_sessionstart(session):
    clear_bgz_and_indexes()
    import matplotlib
    matplotlib.use("Agg")


def pytest_sessionfinish(session):
    clear_bgz_and_indexes()


@pytest.fixture
def data_dir():
    return DATA_DIR


@pytest.fixture
def test_interval():
    test_interval = GenomeRange("chr9:4000000-6000000")
    return test_interval


@pytest.fixture
def test_itv(test_interval):
    return str(test_interval).replace(':', '_').replace('-', '_')


@pytest.fixture
def empty_interval():
    empty_interval = GenomeRange("chr10:4000000-6000000")
    return empty_interval


@pytest.fixture
def sub_interval1():
    sub_interval1 = GenomeRange("chr9:4500000-5000000")
    return sub_interval1


@pytest.fixture
def sub_interval2():
    sub_interval2 = GenomeRange("chr9:5200000-5850000")
    return sub_interval2


@pytest.fixture
def tmp_dir():
    # yield a temporary directory
    # delete it after the test
    from tempfile import mkdtemp
    import shutil
    tmpdir = mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)
