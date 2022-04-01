#! /bin/bash
source $7
source $1 $2
roary -e --mafft -p 8 $3 -z -i $4 -iv $5
source $6 $2