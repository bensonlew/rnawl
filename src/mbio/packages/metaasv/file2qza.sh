#! /bin/bash
source activate $1
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-format SingleEndFastqManifestPhred33V2 --input-path $2 --output-path $3