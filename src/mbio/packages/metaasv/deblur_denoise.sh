#! /bin/bash
source activate $1
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-format SingleEndFastqManifestPhred33V2 --input-path ${13} --output-path $2
qiime deblur denoise-other --i-demultiplexed-seqs $2 --i-reference-seqs $3 --o-representative-sequences $4 --o-table $5  --o-stats $6 --p-jobs-to-start $7 --p-trim-length $8 --p-left-trim-len $9
qiime metadata tabulate --m-input-file $5 --o-visualization ${10}
qiime metadata tabulate --m-input-file $4 --o-visualization ${11}
qiime deblur visualize-stats --i-deblur-stats $6 --o-visualization ${12}