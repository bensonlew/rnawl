#! /bin/bash

source /mnt/ilustre/users/sanger-dev/app/bioinfo/Genomic/Sofware/smrtanalysis/smrtanalysis/current/etc/setup.sh &&$1 --distribute -D NPROC=4 -D CLUSTER=BASH -D MAX_THREADS=12 --output=$2 --params=settings.xml xml:input.xml 