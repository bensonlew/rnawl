#!/bin/bash

GTF_FILE="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Mitochondrion/MT.GRCh38.96.gtf"
LIST_FILE="/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/wgbs/bsmap1/rsrs/output/dss.list.txt"
RESULT_FILE="/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/wgbs/bsmap1/rsrs/output/result.txt"
OUTPUT_DIR="/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/wgbs/bsmap1/rsrs/figure"

# set
set -e

mkdir -p ${OUTPUT_DIR}

# hist
echo -n "Begin hist "; date
/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python \
  /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/hist.py \
  ${RESULT_FILE} ${OUTPUT_DIR}
echo -n "Final hist "; date

# box
echo -n "Begin box "; date
/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python \
  /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/box.py \
  ${RESULT_FILE} ${OUTPUT_DIR}
echo -n "Final box "; date

# violin
echo -n "Begin violin "; date
/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python \
  /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/violin.py \
  ${RESULT_FILE} ${OUTPUT_DIR}
echo -n "Final violin "; date

# pie
echo -n "Begin pie "; date
/mnt/ilustre/users/sanger-dev/app/bioinfo/wgbs/miniconda3/bin/python3.7 \
  /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/pie.py \
  ${RESULT_FILE} ${OUTPUT_DIR}
echo -n "Final pie "; date

# bar
echo -n "Begin bar "; date
/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python \
  /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/bar.py \
  ${GTF_FILE} ${RESULT_FILE} ${OUTPUT_DIR}
echo -n "Final bar "; date

# heatmap
echo -n "Begin heatmap "; date
/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python \
  /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/difference/multi_contrast.py \
  ${LIST_FILE} ${RESULT_FILE} ${OUTPUT_DIR}

/mnt/ilustre/users/sanger-dev/app/bioinfo/wgbs/miniconda3/bin/python3.7 \
  /mnt/ilustre/users/sanger-dev/biocluster/src/mbio/packages/wgbs/visualization/heatmap.py \
  ${OUTPUT_DIR}/dml.result.txt ${RESULT_FILE} ${OUTPUT_DIR}
echo -n "Final heatmap "; date
