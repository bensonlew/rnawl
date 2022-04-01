#! /bin/bash
source $8
conda activate $1
qiime feature-classifier classify-sklearn --i-classifier $2 --i-reads $3  --o-classification $4 --p-confidence $5 --p-n-jobs $6
qiime tools export --input-path $4 --output-path $7