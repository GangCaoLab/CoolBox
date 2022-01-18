import logging
import sys

LOG_LEVEL = logging.WARNING


def get_logger(name, file_=sys.stderr, level=LOG_LEVEL):
    FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
    formatter = logging.Formatter(fmt=FORMAT)
    if isinstance(file_, str):
        handler = logging.FileHandler(file_)
    else:
        handler = logging.StreamHandler(file_)
    handler.setFormatter(formatter)
    log = logging.getLogger(name)
    log.addHandler(handler)
    log.setLevel(level)
    return log
