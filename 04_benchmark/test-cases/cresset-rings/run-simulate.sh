#!/usr/bin/env bash

# python simulate.py \
#     -ff ../../forcefields/fb-fit-v3-single-mean-k100.offxml \
#     -o output/v3-k100

python simulate.py \
    -ff ../../forcefields/fb-fit-v1-single-mean-k100.offxml \
    -o output/v1-k100

python simulate.py \
    -ff openff-2.2.1.offxml \
    -o output/openff-2.2.1
