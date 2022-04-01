#! /bin/bash
source $5
conda activate $1
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads $2 --i-reference-taxonomy $3  --o-classifier $4