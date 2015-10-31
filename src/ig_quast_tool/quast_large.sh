#!/bin/bash

DIR="/johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/gzipped"
OUT_DIR="/home/ashlemov/barcodes_agressive/"

for i in `seq 1 9`; do
    ../barcode_metrics/barcode_quast.py --Bc ${OUT_DIR}/repertoire_large_${i}.fa \
      --Br ${OUT_DIR}/barcodes_large_${i}.rcm \
      -c ${OUT_DIR}/igrc_out_${i}/final_repertoire_large.fa \
      -r ${OUT_DIR}/igrc_out_${i}/final_repertoire.rcm \
      --out ${OUT_DIR}/quast_large_${i} \
      --rate 0.5 \
      --threads-num 2 &
done

wait
