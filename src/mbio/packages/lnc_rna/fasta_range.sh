#!/usr/bin/env bash
BIOAWK_PATH=$1
REF_RNA_PATH=$2
SAMPLE_NAME=$3
${BIOAWK_PATH}bioawk -c fastx '{print $name"\t1\t"length($seq)}' $REF_RNA_PATH > $SAMPLE_NAME