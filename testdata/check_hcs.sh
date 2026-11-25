#!/usr/bin/env bash
# Sanity checks for the HCs produced by `nextflow run main.nf -profile test`.
#
# Usage: testdata/check_hcs.sh <HCs.csv> [nHCs (default 10)]
#
# Checks: a 'sample' column with the test samples in BCF order (as in samples_superpop.txt), nHCs HC
# columns, all values finite, and that the leading HCs separate the 1000G superpopulations (R² of an HC
# on superpopulation: between-superpopulation over total sum of squares). Exits non-zero if a check fails.

set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "usage: $0 <HCs.csv> [nHCs]" >&2
    exit 1
fi
HCS=$1
NHCS=${2:-10}
POPS="$(cd "$(dirname "$0")" && pwd)/samples_superpop.txt"

fail() { echo "FAIL: $*"; exit 1; }

[ "$(head -n 1 "${HCS}" | cut -d, -f1)" = "sample" ] || fail "the first column should be 'sample'"
cmp -s <(tail -n +2 "${HCS}" | cut -d, -f1) <(cut -f1 "${POPS}") || fail "sample IDs are not the test samples, in BCF order"
echo "ok: sample IDs match the test samples, in BCF order"

awk -F, -v n="${NHCS}" '
    NR == 1 { if (NF - 1 != n) { print "FAIL: expected " n " HC columns, got " NF - 1; exit 1 }; next }
    { for (i = 2; i <= NF; i++) if ($i !~ /^-?[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?$/) { print "FAIL: not a finite number: " $i; exit 1 } }
' "${HCS}"
echo "ok: $(($(wc -l < "${HCS}") - 1)) rows, ${NHCS} HC columns, all values finite"

# R² of HC1..HC5 on superpopulation, and mean HC1/HC2 by superpopulation
paste -d, <(cut -f2 "${POPS}") <(tail -n +2 "${HCS}" | cut -d, -f2-) | awk -F, -v k=5 '
    { g = $1; n[g]++; N++; for (i = 1; i <= k && i < NF; i++) { x = $(i + 1); s[g, i] += x; t[i] += x; tt[i] += x * x } }
    END {
        best = 0; line = "R^2 of HC1..HC" k " on superpopulation:"
        for (i = 1; i <= k; i++) {
            total = tt[i] - t[i] * t[i] / N; between = 0
            for (g in n) between += s[g, i] * s[g, i] / n[g]
            between -= t[i] * t[i] / N
            r2 = between / total; line = line sprintf(" %.2f", r2); if (r2 > best) best = r2
        }
        print line
        print "Mean HC1/HC2 by superpopulation:"
        for (g in n) printf "  %s %8.3f %8.3f\n", g, s[g, 1] / n[g], s[g, 2] / n[g]
        if (best < 0.5) { print "FAIL: no leading HC is strongly associated with superpopulation (max R^2 < 0.5)"; exit 1 }
        print "ok: leading HCs separate superpopulations"
    }'
