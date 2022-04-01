#! /bin/bash
source ${11}
conda activate $1
qiime feature-classifier classify-consensus-vsearch --i-query $2  --i-reference-reads $3 --i-reference-taxonomy $4 --p-maxaccepts $5 --p-perc-identity $6 --p-query-cov $7 --p-threads $8 --o-classification $9 --p-top-hits-only
qiime tools export --input-path $9 --output-path ${10}