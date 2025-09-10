#!/usr/bin/env bash

python check-parameters.py \
    -ff ../../forcefields/fb-fit-v3-single-mean-k100.offxml  > check-parameters-v3.log

python check-parameters.py \
    -ff ../../forcefields/fb-fit-v1-single-mean-k100.offxml  > check-parameters-v1.log

