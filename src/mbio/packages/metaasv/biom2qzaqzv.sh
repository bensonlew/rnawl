#! /bin/bash
source activate $1
qiime tools import --type 'FeatureTable[Frequency]' --input-format 'BIOMV210Format' --input-path $2 --output-path $3
qiime metadata tabulate --m-input-file $3 --o-visualization $4