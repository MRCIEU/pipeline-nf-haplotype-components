#!/usr/bin/env bash
# Builds the HC-step test data in testdata/hc/: a combined chunk length matrix, its row sums and
# sample list, as written by the pipeline (--publish_matrix), for testing the HC step
# (calculate_haplotype_components) without painting.
#
# Unlike testdata/bcf/, this is painted at full HapMap 3 SNP density with the production
# settings (pbwt -paintSparse 100 2 500), so it has realistic sparsity and structure.
#
#   - 904 samples: 750 unrelated (150 per superpopulation), 50 parent-child trios (10 per
#     superpopulation, fewer where not available) and 2 sibling pairs (see select_samples.sh)
#   - the samples are painted against each other only (pipeline --keep_samples)
#
# Usage: testdata/make_hc_testdata.sh <dir containing chr{1..22}.hm3.snps.bcf>
# Needs: bcftools, GNU shuf (gshuf on macOS), curl, nextflow, docker (the pipeline's images)

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <dir containing chr{1..22}.hm3.snps.bcf>" >&2
    exit 1
fi
FULL_DIR=$(cd "$1" && pwd)
NPERPOP=150   # unrelated samples per superpopulation
NTRIOS=10     # trios per superpopulation
NSIBPAIRS=5   # sibling pairs in total (at most; few are available without their parents)

TESTDATA_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_DIR=$(dirname "${TESTDATA_DIR}")
OUT_DIR="${TESTDATA_DIR}/hc"
mkdir -p "${OUT_DIR}"

# --- samples (writes samples.txt, samples_superpop.txt, related_pairs.txt) ---
"${TESTDATA_DIR}/select_samples.sh" "${FULL_DIR}/chr1.hm3.snps.bcf" "${TESTDATA_DIR}/meta" \
    "${NPERPOP}" "${NTRIOS}" "${NSIBPAIRS}" "${OUT_DIR}/"

# --- paint the selected samples with the pipeline, and keep the matrix ---
TMP=$(mktemp -d)
trap 'rm -rf "${TMP}"' EXIT
for i in $(seq 1 22); do
    echo "${i} ${FULL_DIR}/chr${i}.hm3.snps.bcf"
done > "${TMP}/vcf_list.txt"

(
    cd "${REPO_DIR}"
    nextflow run main.nf -ansi-log false \
        --vcf_file_list "${TMP}/vcf_list.txt" \
        --map_file "${TESTDATA_DIR}/map_lengths.txt" \
        --keep_samples "${OUT_DIR}/samples.txt" \
        --nHCs 20 \
        --publish_matrix \
        --resultsdir "${TMP}/results"
)

# the pipeline's sample list must be the one we asked for, in the same order
cmp "${TMP}/results/matrix/samples.txt" "${OUT_DIR}/samples.txt"
cp "${TMP}/results/matrix/combined.chunklengths.txt.gz" "${TMP}/results/matrix/combined.chunklengths.sums" "${OUT_DIR}/"

echo "Wrote ${OUT_DIR}:" >&2
ls -l "${OUT_DIR}" >&2
