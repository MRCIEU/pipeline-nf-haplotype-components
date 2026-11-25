#!/usr/bin/env bash
# End-to-end tests of the pipeline on the test data (-profile test). Run from anywhere:
#
#   tests/test_pipeline.sh
#
# Needs nextflow and docker, with the pipeline's images built (./build.sh in each dependencies/<tool>/).
# Runs that share a painting reuse it (-resume). Exits non-zero if a test fails.

set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "${REPO}"

IMAGE_R=mrcieu/pipeline-nf-haplotype-components:calculate_HCs
IMAGE_BCFTOOLS=mrcieu/pipeline-nf-haplotype-components:bcftools-1.24

# temporary files go inside the repository, so docker can see them (e.g. Docker Desktop file sharing)
TMP="${REPO}/tests/.tmp"
rm -rf "${TMP}"
mkdir -p "${TMP}"

failures=0
check() {  # check <description> <command...>
    local what=$1; shift
    if "$@" > /dev/null 2>&1; then echo "  ok:   ${what}"; else echo "  FAIL: ${what}"; failures=$((failures + 1)); fi
}
nf() {  # nf <name> <args...>: run the pipeline on the test profile; results in ${TMP}/<name>, log in ${TMP}/<name>.log
    local name=$1; shift
    nextflow run main.nf -profile test -resume -ansi-log false --resultsdir "${TMP}/${name}" "$@" > "${TMP}/${name}.log" 2>&1
    echo $? > "${TMP}/${name}.status"
}
succeeded() { [ "$(cat "${TMP}/$1.status")" -eq 0 ]; }
failed_with() { ! succeeded "$1" && grep -q -- "$2" "${TMP}/$1.log"; }
check_hcs() {  # check_hcs <HCs.csv>: testdata/check_hcs.R, in the R image
    docker run --rm -v "${REPO}":"${REPO}" -w "${REPO}" "${IMAGE_R}" Rscript testdata/check_hcs.R "$1" 10
}
statuses() { tail -n +2 "$1" | awk -F'\t' '{print $NF}' | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}'; }

related_ids() { awk '{print $1; print $2}' testdata/related_pairs.txt | sort -u; }
reverse_lines() { awk '{a[NR]=$0} END {for (i=NR; i>=1; i--) print a[i]}' "$@"; }  # portable tac
grep -v -x -F -f <(related_ids) testdata/samples.txt > "${TMP}/unrelated.txt"
awk '$2=="AFR"{print $1}' testdata/samples_superpop.txt | reverse_lines > "${TMP}/afr_reversed.txt"

echo "Default test run"
nf default
check "succeeds" succeeded default
check "HCs pass testdata/check_hcs.R" check_hcs "${TMP}/default/HCs.csv"
check "diagnostics have one row per sample" test "$(($(wc -l < "${TMP}/default/HC_diagnostics.tsv") - 1))" -eq "$(wc -l < testdata/samples.txt)"
check "no matrix published without --publish_matrix" test ! -e "${TMP}/default/matrix"

echo "--publish_matrix"
nf publish --publish_matrix
check "succeeds" succeeded publish
check "matrix, sums and samples published" test -s "${TMP}/publish/matrix/combined.chunklengths.txt.gz" -a -s "${TMP}/publish/matrix/combined.chunklengths.sums" -a -s "${TMP}/publish/matrix/samples.txt"
check "published samples.txt is the BCF sample order" cmp "${TMP}/publish/matrix/samples.txt" testdata/samples.txt
check "HCs identical to the default run" cmp "${TMP}/publish/HCs.csv" "${TMP}/default/HCs.csv"

echo "--restrict + --remove_pairs + --reweight"
nf restrict --restrict 50000 --remove_pairs testdata/related_pairs.txt --reweight
check "succeeds" succeeded restrict
check "some rows lost entries to the cap (frac_kept < 1)" awk -F'\t' 'NR>1 && $5 < 0.99 {found=1} END {exit !found}' "${TMP}/restrict/HC_diagnostics.tsv"
check "relatives' rows lost entries (frac_removed > 0)" awk -F'\t' 'NR>1 && $7 > 0 {found=1} END {exit !found}' "${TMP}/restrict/HC_diagnostics.tsv"

echo "--restrict --restrict_mode row-fraction"
nf rowfrac --restrict 50000 --restrict_mode row-fraction
check "succeeds" succeeded rowfrac
# every row keeps the same share of its entries in the uncapped default run (global: 0.6% to 76% here)
kept_share_spread() {  # kept_share_spread <diagnostics>: min and max of n_entries / uncapped n_entries
    join -t $'\t' <(tail -n +2 "${TMP}/default/HC_diagnostics.tsv" | cut -f1,2 | sort) <(tail -n +2 "$1" | cut -f1,2 | sort) |
        awk -F'\t' '{r = $3 / $2; if (NR == 1 || r < mn) mn = r; if (r > mx) mx = r} END {print mn, mx}'
}
check "every row keeps the same share of its entries" awk -v s="$(kept_share_spread "${TMP}/rowfrac/HC_diagnostics.tsv")" \
    'BEGIN {split(s, r, " "); exit !(r[1] > 0 && r[1] >= 0.9 * r[2])}'
check "HCs pass testdata/check_hcs.R" check_hcs "${TMP}/rowfrac/HCs.csv"

echo "--fit_samples (unrelated) + --project + --reweight"
nf project --fit_samples "${TMP}/unrelated.txt" --project --reweight
check "succeeds" succeeded project
check "unrelated fitted, relatives projected" test "$(statuses "${TMP}/project/HC_diagnostics.tsv")" = "fitted=$(wc -l < "${TMP}/unrelated.txt" | tr -d ' ') projected=$(related_ids | wc -l | tr -d ' ') "
check "HCs pass testdata/check_hcs.R" check_hcs "${TMP}/project/HCs.csv"

echo "--keep_samples (AFR, reversed)"
nf keep --keep_samples "${TMP}/afr_reversed.txt"
check "succeeds" succeeded keep
check "rows are the kept samples, in the keep file's order" cmp <(tail -n +2 "${TMP}/keep/HCs.csv" | cut -d, -f1) "${TMP}/afr_reversed.txt"

echo "Errors"
nf e_project --project
check "--project without --fit_samples fails" failed_with e_project "needs --fit_samples"
nf e_mode --restrict_mode per-row
check "an unknown --restrict_mode fails" failed_with e_mode "must be global or row-fraction"
nf e_keep --keep_samples "${TMP}/does_not_exist.txt"
check "a missing --keep_samples file fails" failed_with e_keep "No such file"
printf 'NOT_A_SAMPLE\n' > "${TMP}/typo.txt"
nf e_typo --keep_samples "${TMP}/typo.txt"
check "an unknown ID in --keep_samples fails" failed_with e_typo "does not exist in header"
# chr22 with its samples in a different order
mkdir -p "${TMP}/reordered"
cp testdata/bcf/chr20.bcf testdata/bcf/chr21.bcf "${TMP}/reordered/"
reverse_lines testdata/samples.txt > "${TMP}/reordered/samples.txt"
docker run --rm --user "$(id -u):$(id -g)" -v "${REPO}":"${REPO}" -w "${REPO}" "${IMAGE_BCFTOOLS}" \
    bcftools view -S "${TMP}/reordered/samples.txt" -Ob -o "${TMP}/reordered/chr22.bcf" testdata/bcf/chr22.bcf
printf '20 chr20.bcf\n21 chr21.bcf\n22 chr22.bcf\n' > "${TMP}/reordered/vcf_list.txt"
nf e_order --vcf_file_list "${TMP}/reordered/vcf_list.txt"
check "input files with samples in a different order fail" failed_with e_order "does not have the same samples, in the same order"

if [ "${failures}" -eq 0 ]; then
    rm -rf "${TMP}"
    echo; echo "All pipeline tests passed"
else
    echo; echo "${failures} pipeline test(s) FAILED; logs in ${TMP}"
    exit 1
fi
