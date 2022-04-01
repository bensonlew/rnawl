#! /bin/bash
source activate $1

deblur workflow --seqs-fp $2 --output-dir $3 -O $4 --trim-length $5 --threads-per-sample $6 --pos-ref-fp $7 --pos-ref-db-fp $8 -w --min-size $9 --min-reads ${10}