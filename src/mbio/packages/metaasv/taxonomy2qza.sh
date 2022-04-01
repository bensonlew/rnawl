#! /bin/bash
source activate $1
qiime tools import --type 'FeatureData[Taxonomy]' --input-path $2 --output-path $3