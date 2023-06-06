#!/bin/bash
eval "$(conda shell.bash hook)"
not_aligned="$1"
reference="$2"
aligned="$3"

cat CNS_5/* > "$not_aligned"

conda activate nextstrain
augur align --sequences "$not_aligned" --reference-sequence "$reference" --output "$aligned"
conda deactivate