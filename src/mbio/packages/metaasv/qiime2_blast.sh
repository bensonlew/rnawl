#! /bin/bash
source ${10}
conda activate $1
qiime feature-classifier classify-consensus-blast --i-query $2 --i-reference-reads $3 --i-reference-taxonomy $4 --p-maxaccepts $5 --p-perc-identity $6 --p-query-cov $7 --o-classification $8
qiime tools export --input-path $8 --output-path $9