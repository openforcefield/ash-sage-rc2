#!/bin/bash

python compute-tfes.py          \
    -i output/freesolv.csv      \
    -i output/mnsol.csv         \
    -o output/tfes.csv > logs/compute-tfes.log

