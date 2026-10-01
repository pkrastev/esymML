#!/bin/bash
#=====================================================================
# make_data.sh: regenerate the data the DNNs in models/ were trained on
#               EOS tables:    eos_paper_v1.h5 ... eos_paper_v4.h5
#               NS sequences:  ns_paper_v1.h5  ... ns_paper_v4.h5
#
#   Build eos_v3_cli.x and tov_ml_cli.x first (see README.md), then
#   run in this directory:   ./make_data.sh
#
#   Takes about 2-2.5 hours on one core (about 2 minutes per EOS file,
#   25-35 minutes per NS file). With the toolchain of README.md the
#   files are bit-for-bit identical to the ones used for training.
#=====================================================================
set -euo pipefail
cd "$(dirname "$0")"

for f in eos_v3_cli.x tov_ml_cli.x EOS_P_A18dvUIX.txt; do
    if [ ! -e "$f" ]; then
        echo "make_data.sh: $f not found (build the programs first, see README.md)"
        exit 1
    fi
done

# Number of accepted EOSs in the original files, for each seed
expected=(15070 14925 15181 14968)
status=0

# --- 1. EOS tables (seeds 101-104) ---
for v in 1 2 3 4; do
    seed=$((100 + v))
    echo "=== EOS set ${v} (seed ${seed})"
    out=$(./eos_v3_cli.x -n 100000 --paper -s ${seed} -o eos_paper_v${v}.h5 | tail -1)
    echo "${out}"
    n=$(echo "${out}" | awk '{print $2}')
    if [ "${n}" != "${expected[$((v-1))]}" ]; then
        echo "WARNING: expected ${expected[$((v-1))]} EOSs, got ${n} (different toolchain? see README.md)"
        status=1
    fi
done

# --- 2. NS sequences ---
for v in 1 2 3 4; do
    echo "=== NS sequences ${v}"
    ./tov_ml_cli.x -i eos_paper_v${v}.h5 -o ns_paper_v${v}.h5 | grep -v "^EOS "
done

if [ ${status} -eq 0 ]; then
    echo "make_data.sh: done, all EOS counts as expected"
else
    echo "make_data.sh: done, but some EOS counts differ from the original files"
fi
exit ${status}
