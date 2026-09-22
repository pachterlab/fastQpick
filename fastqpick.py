"""Lowercase alias so that `import fastqpick` is equivalent to `import fastQpick`."""
import sys

import fastQpick

sys.modules[__name__] = fastQpick
