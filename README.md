# fastQpick

Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.

Install it with pip from pypi via:
pip install fastQpick

or with pip from the source code via:
pip install git+https://github.com/pachterlab/fastQpick.git

or by cloning the repository and running:
git clone https://github.com/pachterlab/fastQpick.git
cd fastQpick
python -m build
python -m pip install dist/fastQpick-x.x.x-py3-none-any.whl


bash
fastQpick --fraction FRACTION [OPTIONS] FASTQ_FILE1 FASTQ_FILE2 ...

python
from fastQpick import fastQpick
fastQpick(input_file_list=['FASTQ_FILE1', 'FASTQ_FILE2', ...], fraction=FRACTION, ...)

See fastQpick --help (shell) or help(fastQpick) (python) for more information.