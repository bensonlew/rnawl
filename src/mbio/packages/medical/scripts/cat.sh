#! /bin/bash

zcat $1 $2 >> $3
gzip $3
