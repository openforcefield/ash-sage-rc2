#!/usr/bin/env bash

QM_FFS=(
    -nf 'Sage 2.0.0' 'openff-2.0.0'
    -nf 'Sage 2.2.1' 'openff-2.2.1'
    -nf 'Sage 2.3.0rc1' 'openff-2.3.0rc1'
    -nf 'v1-k100'    'fb-fit-v1-single-mean-k100'
    -nf 'v3-k100'    'fb-fit-v3-single-mean-k100'
)
python combine-mm-qm-diff.py \
    "${QM_FFS[@]}" \
    -i  /Volumes/Nobbsy\ 1/ash-sage-rc2/04_benchmark/qm/topology-comparisons/pyarrow \
    -o  output/qm/mm-qm-diff/labelled-properties.csv > logs/combine-mm-qm-diff.log

