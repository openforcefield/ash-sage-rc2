#!/usr/bin/env bash

mkdir -p logs/stats

# === sfes ===
SFE_KWARGS=(
    -cy 'Value (kcal / mol)'
    -cx 'Reference Value (kcal / mol)'
    -cye 'Uncertainty (kcal / mol)'
    -cxe 'Reference Uncertainty (kcal / mol)'
)

rm -f logs/stats/freesolv.log
rm -f logs/stats/mnsol.log
rm -f logs/stats/tfes.log

for FF in 'Sage 2.2.1 + AmberTools' \
          'Sage 2.2.1 + AshGC' \
          'Sage 2.2.1 + ELF10' \
          'Sage 2.3.0rc1' \
          'v1-k100' \
          'v3-k100'
do
    datetime=$(date +%Y-%m-%d_%H-%M-%S)
    echo "[$datetime] Starting statistics computation for SFEs with force field: $FF"
    # replace spaces with dashes for filenames
    FF_FN=$(echo $FF | tr ' ' '-')
    echo "Force field filename: $FF_FN"

    python compute-statistics.py                    \
        "${SFE_KWARGS[@]}"                          \
        -ff "$FF"                                   \
        -i output/sfes/freesolv/labelled-data.csv   \
        -o output/sfes/freesolv/stats/$FF_FN.csv    \
        >> logs/stats/freesolv.log

    python compute-statistics.py                    \
        "${SFE_KWARGS[@]}"                          \
        -ff "$FF"                                   \
        -i output/sfes/mnsol/labelled-data.csv      \
        -o output/sfes/mnsol/stats/$FF_FN.csv       \
        >> logs/stats/mnsol.log
    
    python compute-statistics.py                    \
        "${SFE_KWARGS[@]}"                          \
        -ff "$FF"                                   \
        -i output/sfes/tfes/labelled-data.csv       \
        -o output/sfes/tfes/stats/$FF_FN.csv        \
        >> logs/stats/tfes.log

    datetime=$(date +%Y-%m-%d_%H-%M-%S)
    echo "[$datetime] Completed statistics computation for SFEs with force field: $FF".
done



