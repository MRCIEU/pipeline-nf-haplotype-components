#!/usr/bin/env bash
# Selects 1000 Genomes high-coverage samples for test data. Used by make_testdata.sh and
# make_hc_testdata.sh.
#
#   - NPERPOP unrelated samples per superpopulation (AFR, AMR, EAS, EUR, SAS), at random with a
#     fixed seed, from the 2504 samples not in the 698 added to complete trios
#   - NTRIOS parent-child trios per superpopulation (child + both parents; families with one child,
#     none of whom are among the unrelated samples), at random with a fixed seed
#   - up to NSIBPAIRS full-sibling pairs (the first two children of a family) whose parents and
#     siblings are not otherwise selected, first in sorted order of superpopulation and parents
#
# Writes, in the sample order of BCF:
#   <OUT_PREFIX>samples.txt           sample IDs
#   <OUT_PREFIX>samples_superpop.txt  sample ID, superpopulation
#   <OUT_PREFIX>related_pairs.txt     sample ID 1, sample ID 2, relationship (parent-offspring or
#                                     sibling): every such pair among the selected samples, from the
#                                     1000G pedigree
#
# Usage: select_samples.sh <BCF (for sample order)> <meta dir> <NPERPOP> <NTRIOS> <NSIBPAIRS> <OUT_PREFIX>
# Needs: bcftools, GNU shuf (gshuf on macOS), curl

set -euo pipefail

if [ $# -ne 6 ]; then
    echo "usage: $0 <BCF> <meta dir> <NPERPOP> <NTRIOS> <NSIBPAIRS> <OUT_PREFIX>" >&2
    exit 1
fi
BCF=$1; META_DIR=$2; NPERPOP=$3; NTRIOS=$4; NSIBPAIRS=$5; OUT_PREFIX=$6

SHUF=$(command -v gshuf || command -v shuf)
FTP=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
PED=20130606_g1k_3202_samples_ped_population.txt
RELATED=1000G_698_related_high_coverage.sequence.index

mkdir -p "${META_DIR}"
[ -s "${META_DIR}/${PED}" ]     || curl -sSf -o "${META_DIR}/${PED}" "${FTP}/${PED}"
[ -s "${META_DIR}/${RELATED}" ] || curl -sSf -o "${META_DIR}/${RELATED}" "${FTP}/${RELATED}"
PED="${META_DIR}/${PED}"

TMP=$(mktemp -d)
trap 'rm -rf "${TMP}"' EXIT

# --- unrelated: NPERPOP per superpopulation ---
grep -v '^#' "${META_DIR}/${RELATED}" | cut -f10 | sort -u > "${TMP}/related.txt"
for sp in AFR AMR EAS EUR SAS; do
    awk -v sp="${sp}" 'NR==FNR {rel[$1]=1; next}
                       FNR>1 && !($2 in rel) && $7==sp {print $2"\t"$7}' \
        "${TMP}/related.txt" "${PED}" |
        "${SHUF}" -n "${NPERPOP}" --random-source=<(yes)
done > "${TMP}/unrelated.tsv"

# --- trios: one-child families with both parents in the data, none already picked ---
# (one line per trio: superpop, child, father, mother)
awk 'NR==FNR {picked[$1]=1; next}
     FNR==1 {next}
     {sp[$2]=$7; father[$2]=$3; mother[$2]=$4; insample[$2]=1}
     $3!="0" && $4!="0" {nchild[$3"_"$4]++}
     END {
         for (c in father) {
             f = father[c]; m = mother[c]
             if (f=="0" || m=="0" || !(f in insample) || !(m in insample)) continue
             if (nchild[f"_"m] != 1) continue
             if ((c in picked) || (f in picked) || (m in picked)) continue
             print sp[c]"\t"c"\t"f"\t"m
         }
     }' "${TMP}/unrelated.tsv" "${PED}" | sort > "${TMP}/trio_candidates.tsv"
for sp in AFR AMR EAS EUR SAS; do
    awk -v sp="${sp}" '$1==sp' "${TMP}/trio_candidates.tsv" |
        "${SHUF}" -n "${NTRIOS}" --random-source=<(yes)
done > "${TMP}/trios.tsv"

# --- sibling pairs: the first two children of families with 2+ children; parents not included ---
awk 'NR==FNR {picked[$1]=1; next}
     FNR==1 {next}
     $3!="0" && $4!="0" {k=$3"_"$4; n[k]++; kids[k]=kids[k]" "$2; sp[k]=$7
                         if (($2 in picked) || ($3 in picked) || ($4 in picked)) bad[k]=1}
     END {for (k in n) if (n[k] >= 2 && !(k in bad)) {split(kids[k], a, " "); print sp[k]"\t"k"\t"a[1]"\t"a[2]}}' \
    <(cut -f2-4 "${TMP}/trios.tsv" | tr '\t' '\n'; cut -f1 "${TMP}/unrelated.tsv") "${PED}" |
    sort | head -n "${NSIBPAIRS}" \
    > "${TMP}/sibs.tsv"

# --- all selected samples, with superpopulation, in BCF order ---
{
    cat "${TMP}/unrelated.tsv"
    awk '{print $2"\t"$1; print $3"\t"$1; print $4"\t"$1}' "${TMP}/trios.tsv"
    awk '{print $3"\t"$1; print $4"\t"$1}' "${TMP}/sibs.tsv"
} > "${TMP}/selected.tsv"

if [ "$(cut -f1 "${TMP}/selected.tsv" | sort | uniq -d | wc -l)" -ne 0 ]; then
    echo "Error: a sample was selected twice" >&2
    exit 1
fi

bcftools query -l "${BCF}" |
    awk 'NR==FNR {p[$1]=$2; next} ($1 in p) {print $1"\t"p[$1]}' "${TMP}/selected.tsv" - \
    > "${OUT_PREFIX}samples_superpop.txt"
cut -f1 "${OUT_PREFIX}samples_superpop.txt" > "${OUT_PREFIX}samples.txt"

if [ "$(wc -l < "${OUT_PREFIX}samples.txt")" -ne "$(wc -l < "${TMP}/selected.tsv")" ]; then
    echo "Error: some selected samples are not in ${BCF}" >&2
    exit 1
fi

# --- related pairs among the selected samples, from the pedigree ---
awk 'NR==FNR {sel[$1]=1; next}
     FNR==1 {next}
     ($2 in sel) {
         if ($3 in sel) print $3"\t"$2"\tparent-offspring"
         if ($4 in sel) print $4"\t"$2"\tparent-offspring"
         if ($3!="0" && $4!="0") {k=$3"_"$4; kids[k]=kids[k]" "$2}
     }
     END {
         for (k in kids) {
             n = split(kids[k], a, " ")
             for (i = 1; i < n; i++) for (j = i + 1; j <= n; j++) print a[i]"\t"a[j]"\tsibling"
         }
     }' "${OUT_PREFIX}samples.txt" "${PED}" | sort > "${OUT_PREFIX}related_pairs.txt"

echo "Selected $(wc -l < "${OUT_PREFIX}samples.txt") samples:" \
     "$(wc -l < "${TMP}/unrelated.tsv" | tr -d ' ') unrelated," \
     "$(wc -l < "${TMP}/trios.tsv" | tr -d ' ') trios," \
     "$(wc -l < "${TMP}/sibs.tsv" | tr -d ' ') sibling pairs;" \
     "$(wc -l < "${OUT_PREFIX}related_pairs.txt" | tr -d ' ') related pairs" >&2
