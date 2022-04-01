#! /bin/bash
source activate $1
qiime tools import --type 'FeatureData[Sequence]' --input-path $2 --output-path $3