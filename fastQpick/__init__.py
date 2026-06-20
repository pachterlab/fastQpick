import logging
import sys

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
