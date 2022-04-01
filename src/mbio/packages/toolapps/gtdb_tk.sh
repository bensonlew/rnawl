#! /bin/bash
source $1 && conda activate gtdbtk && export GTDBTK_DATA_PATH=$2
gtdbtk identify --genome_dir $3 --out_dir $4 --cpus 2
gtdbtk align --identify_dir $4 --out_dir $5 --cpu 2
gtdbtk classify --genome_dir $3 --align_dir $5 --out_dir $6 --cpus 2