#!/usr/bin/env bash
BEDTOOLS_PATH=$1
RANGE_BED=$2
SAMPLE_GTF=$3
OUT_GTF=$4
${BEDTOOLS_PATH}bedtools intersect -a ${RANGE_BED} -b ${SAMPLE_GTF} -wb | awk -F "\t" '{print $4"\t"$5"\t"$6"\t"$2"\t"$3"\t"$9"\t"$10"\t"$11"\t"$12}'  > ${OUT_GTF}
