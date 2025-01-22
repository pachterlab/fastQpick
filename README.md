# fastQpick

Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.

bash
fastQpick -i input.fastq -o output.fastq -n 1000

python
from fastQpick import fastQpick
fastQpick(input='input.fastq', output='output.fastq', n=1000)