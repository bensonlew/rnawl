#! /bin/bash
source activate $1
qiime dada2 denoise-single --i-demultiplexed-seqs $2 --o-table $3 --p-trunc-len $4 --o-representative-sequences $5 --o-denoising-stats $6 --p-n-threads $7 --p-trunc-q 0
qiime metadata tabulate --m-input-file $3 --o-visualization $8
qiime metadata tabulate --m-input-file $5 --o-visualization $9
qiime metadata tabulate --m-input-file $6 --o-visualization ${10}