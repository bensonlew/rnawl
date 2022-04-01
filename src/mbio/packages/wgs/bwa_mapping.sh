#! /bin/bash

echo $1 'mem -t $7 -M -R' $2 $3 $4 $5 '>' $6
$1 mem -t $7 -M -R $2 $3 $4 $5 > $6