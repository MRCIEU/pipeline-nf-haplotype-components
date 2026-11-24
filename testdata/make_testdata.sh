#!/usr/bin/env bash
# Builds the small test dataset in testdata/ from 1000 Genomes high-coverage phased data
# (3202 samples), already subset to HapMap 3 SNPs (chr{1..22}.hm3.snps.bcf).
#
#   - 300 unrelated samples, 60 per superpopulation (AFR, AMR, EAS, EUR, SAS)
#   - plus 2 parent-child trios per superpopulation, and 2 sibling pairs (without their parents);
#     related_pairs.txt lists every parent-offspring and sibling pair (see select_samples.sh)
#   - sites polymorphic in those samples, INFO stripped, GT only
#   - every K-th site kept (whole chromosome, so whole-chromosome map lengths stay valid)
#
# Because of the thinning, the test profile also lowers the pbwt -paintSparse match-length
# cutoff (see nextflow.config), otherwise almost no matches survive.
#
# Usage: testdata/make_testdata.sh <dir containing chr{1..22}.hm3.snps.bcf>
# Needs: bcftools, GNU shuf (gshuf on macOS), GNU parallel, curl

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <dir containing chr{1..22}.hm3.snps.bcf>" >&2
    exit 1
fi
FULL_DIR=$1
K=10          # keep every K-th site
NPERPOP=60    # unrelated samples per superpopulation
NTRIOS=2      # trios per superpopulation
NSIBPAIRS=2   # sibling pairs in total

OUT_DIR=$(cd "$(dirname "$0")" && pwd)
mkdir -p "${OUT_DIR}/bcf"

# --- samples (writes samples.txt, samples_superpop.txt, related_pairs.txt) ---
"${OUT_DIR}/select_samples.sh" "${FULL_DIR}/chr1.hm3.snps.bcf" "${OUT_DIR}/meta" \
    "${NPERPOP}" "${NTRIOS}" "${NSIBPAIRS}" "${OUT_DIR}/"

# --- genotypes: subset samples, drop monomorphic sites, strip INFO, thin to every K-th site ---
export FULL_DIR OUT_DIR K
seq 1 22 | parallel -j8 '
    bcftools view -S "${OUT_DIR}/samples.txt" -Ou "${FULL_DIR}/chr{}.hm3.snps.bcf" |
    bcftools view --min-ac 1:minor -Ou |
    bcftools annotate -x INFO,^FORMAT/GT -Ov |
    awk -v k="${K}" "/^#/ || (++n % k == 1)" |
    bcftools view -Ob -o "${OUT_DIR}/bcf/chr$(printf %02d {}).bcf"
'

# --- pipeline inputs (paths relative to this file's directory) ---
for i in $(seq 1 22); do
    printf '%d bcf/chr%02d.bcf\n' "${i}" "${i}"
done > "${OUT_DIR}/vcf_list.txt"
