#!/usr/bin/env bash

mkdir -p logs/largest-differences

# python get-largest-differences.py \
#     -i output/phys-prop/density-all/stats/Sage-2.2.1.csv \
#     -i output/phys-prop/density-all/stats/v3-k100.csv \
#     -o output/phys-prop/density-all/largest-differences.csv \
#     -ref 'Sage 2.2.1' \
#     -tgt v3-k100 > logs/largest-differences/density-all.log

# python get-largest-differences.py \
#     -i output/phys-prop/density-validation/stats/Sage-2.2.1.csv \
#     -i output/phys-prop/density-validation/stats/v3-k100.csv \
#     -o output/phys-prop/density-validation/largest-differences.csv \
#     -ref 'Sage 2.2.1' \
#     -tgt v3-k100 > logs/largest-differences/density-validation.log

# python get-largest-differences.py \
#     -i output/phys-prop/dhmix-all/stats/Sage-2.2.1.csv \
#     -i output/phys-prop/dhmix-all/stats/v3-k100.csv \
#     -o output/phys-prop/dhmix-all/largest-differences.csv \
#     -ref 'Sage 2.2.1' \
#     -tgt v3-k100 > logs/largest-differences/dhmix-all.log


# python get-largest-differences.py \
#     -i output/phys-prop/dhmix-validation/stats/Sage-2.2.1.csv \
#     -i output/phys-prop/dhmix-validation/stats/v3-k100.csv \
#     -o output/phys-prop/dhmix-validation/largest-differences.csv \
#     -ref 'Sage 2.2.1' \
#     -tgt v3-k100 > logs/largest-differences/dhmix-validation.log

# python get-largest-differences.py \
#     -i output/sfes/freesolv/stats/Sage-2.2.1-+-ELF10.csv \
#     -i output/sfes/freesolv/stats/v3-k100.csv \
#     -o output/sfes/freesolv/largest-differences.csv \
#     -ref 'Sage 2.2.1 + ELF10' \
#     -tgt v3-k100 > logs/largest-differences/freesolv.log

# python get-largest-differences.py \
#     -i output/sfes/mnsol/stats/Sage-2.2.1-+-ELF10.csv \
#     -i output/sfes/mnsol/stats/v3-k100.csv \
#     -o output/sfes/mnsol/largest-differences.csv \
#     -ref 'Sage 2.2.1 + ELF10' \
#     -tgt v3-k100 > logs/largest-differences/mnsol.log

# python get-largest-differences.py \
#     -i output/sfes/tfes/stats/Sage-2.2.1-+-ELF10.csv \
#     -i output/sfes/tfes/stats/v3-k100.csv \
#     -o output/sfes/tfes/largest-differences.csv \
#     -ref 'Sage 2.2.1 + ELF10' \
#     -tgt v3-k100 > logs/largest-differences/tfes.log

python get-largest-differences.py \
    -i output/qm/mm-qm-diff/stats/Angles.csv \
    -o output/qm/mm-qm-diff/largest-differences-angles.csv \
    -ref 'Sage 2.2.1' \
    -tgt 'v3-k100' > logs/largest-differences/qm-angles.log

python get-largest-differences.py \
    -i output/qm/mm-qm-diff/stats/Bonds.csv \
    -o output/qm/mm-qm-diff/largest-differences-bonds.csv \
    -ref 'Sage 2.2.1' \
    -tgt 'v3-k100' > logs/largest-differences/qm-bonds.log
