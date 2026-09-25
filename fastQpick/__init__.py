import logging
import os
import sys

# fastQpick never calls BLAS, but the OpenBLAS bundled with numpy otherwise starts one
# busy-waiting worker thread per core in every process (including each pool worker), which
# inflates CPU usage without changing wall time and undercuts the `threads` budget. Pin the
# pool to one thread unless the caller has chosen a value. This only takes effect if numpy
# has not been imported yet, so importing numpy first leaves its configuration untouched.
if "numpy" not in sys.modules:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

# Default package logger. It is defined here, before importing the submodules, so that
# `main` and `utils` can do `from fastQpick import logger` without a circular import.
# A single StreamHandler is attached so that messages are emitted by default (mirroring
# the previous print-based output); guarding on existing handlers keeps repeated imports
# from stacking duplicate handlers.
logger = logging.getLogger("fastQpick")
if not logger.handlers:
    _handler = logging.StreamHandler(sys.stderr)
    _handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    logger.addHandler(_handler)
logger.setLevel(logging.INFO)

from .main import fastQpick

from ._version import __version__
__author__ = "Joseph Rich"
__email__ = "josephrich98@gmail.com"
