#! /bin/bash
source activate $1
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-format SingleEndFastqManifestPhred33V2 --input-path $2 --output-path $3
qiime dada2 denoise-single --i-demultiplexed-seqs $3 --o-table $4 --p-trunc-len $5 --o-representative-sequences $6 --o-denoising-stats $7 --p-n-threads $8 --p-trunc-q ${12} --p-max-ee ${13}
qiime tools export --input-path $4 --output-path $9
qiime tools export --input-path $6 --output-path ${10}
qiime tools export --input-path $7 --output-path ${11}