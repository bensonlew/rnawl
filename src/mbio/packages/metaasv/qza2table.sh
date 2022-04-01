#! /bin/bash
source activate $1
qiime tools export --input-path $2 --output-path $3