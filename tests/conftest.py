import os
import os.path as osp


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

