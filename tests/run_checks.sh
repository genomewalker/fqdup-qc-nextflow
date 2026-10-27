#!/usr/bin/env bash
# Regression checks. Stops at the first failed check with "FAIL[<case>]: ..." and exit 1.
# A broken harness (checker or workflow that cannot run, timeout) stops with "HARNESS: ..." and exit 2.
#   tests/run_checks.sh             trim-report properties and manifest routing (stub runs)
#   tests/run_checks.sh e2e FQDUP   also a real run with a single-lane single-end sample
# CHECKER and PIPELINE override the checker and workflow under test (tests/mutants.sh uses them).
# OUT keeps every case's files and cases.tsv (case, expected, observed) in that directory.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
checker=${CHECKER:-$repo/bin/check_trim_report.awk}
pipeline=${PIPELINE:-$repo}
seeds=${SEEDS:-500}
if [ -n "${OUT:-}" ]; then
    work=$OUT; mkdir -p "$work"
else
    work=$(mktemp -d "${TMPDIR:-/tmp}/fqdup-qc-checks.XXXXXX"); trap 'rm -rf "$work"' EXIT
fi
die() { echo "FAIL[$1]: $2" >&2; exit 1; }
harness() { echo "HARNESS: $*" >&2; exit 2; }
record() { printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$work/cases.tsv"; }

# Runs the checker on file $1; case $2 expects accept or reject ($3). The checker must exit 0 or 1
# and stay silent: gawk also exits 1 on a syntax error, which must not count as a rejection.
expect() {
    local rc=0 got
    "$checker" "$1" 2> "$work/checker.err" || rc=$?
    [ -s "$work/checker.err" ] && harness "checker stderr in $2: $(head -n 1 "$work/checker.err")"
    case $rc in 0) got=accept ;; 1) got=reject ;; *) harness "checker exit $rc in $2" ;; esac
    record "$2" "$3" "$got"
    [ "$got" = "$3" ] || die "$2" "expected $3, got $got"
}

# A raw_single_end_v1 report from seed $1, with corruption kind $2 (0 = none).
# Every corruption leaves the report invalid; kinds 8-11 and 13 keep the counts numerically
# balanced, so only the format checks can reject them.
report() {
    awk -v seed="$1" -v kind="$2" 'BEGIN {
        srand(seed)
        i = 1 + int(rand() * 2^40); if (seed % 7 == 0) i = 1 + int(rand() * 9)
        r = int(rand() * (i + 1)); if (seed % 5 == 0) r = i; if (seed % 11 == 0) r = 0
        d = i - r; e = int(rand() * (d + 1)); s = int(rand() * (d - e + 1)); l = d - e - s
        if (kind == 11) { r = i; d = 0; e = 0; s = 0; l = 0 }
        if (kind == 15) { i = 0; r = 0; d = 0; e = 0; s = 0; l = 0 }
        n = 0
        if (kind != 1) t[++n] = (kind == 2 ? "raw_single_end_v2" : "raw_single_end_v1")
        t[++n] = "quality_policy\tnone"; t[++n] = "min_phred\tunused"; t[++n] = "min_length\t16"
        t[++n] = "max_length\t0"; t[++n] = "clip5\t0"; t[++n] = "clip3\t0"; t[++n] = "adapter5\t"
        t[++n] = "adapter3\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"; t[++n] = "adapter_min_overlap\t5"
        t[++n] = "read_clip_categories\toverlap_allowed"
        if (kind != 3) t[++n] = "input\t" (kind == 8 ? "+" : "") sprintf("%d", i) (kind == 13 ? "\tx" : "")
        if (kind == 6) t[++n] = "input\t" sprintf("%d", i)
        if (kind != 4) t[++n] = "retained\t" sprintf("%d", r) (kind == 10 ? ".0" : "")
        t[++n] = "empty\t" e; t[++n] = "short_reads\t" s; t[++n] = "long_reads\t" l
        t[++n] = "input_bases\t" 150 * i; t[++n] = "retained_bases\t" 60 * r
        if (kind != 5) t[++n] = "dropped\t" (kind == 9 ? "0" : "") (kind == 11 ? "" : sprintf("%d", d + (kind == 12)))
        if (kind == 7) t[++n] = "dropped\t" sprintf("%d", d)
        for (k = 1; k <= n; k++) out = out t[k] "\n"
        if (kind == 14) { cut = 1 + int(rand() * (length(out) - 2)); out = substr(out, 1, cut) }
        printf "%s", out
    }'
}

# Counts beyond 2^53, where awk's doubles round: equal as doubles, unequal as integers.
counts() { printf 'raw_single_end_v1\ninput\t%s\nretained\t%s\ndropped\t%s\n' "$1" "$2" "$3" > "$work/r"; }

check_reports() {
    local seed kind
    counts 3 2 1; expect "$work/r" precision-small accept
    counts 9007199254740993 9007199254740992 1; expect "$work/r" precision-2p53-valid accept
    counts 9007199254740992 9007199254740993 0; expect "$work/r" precision-2p53 reject
    counts 100000000000000000000 99999999999999999999 1; expect "$work/r" precision-20digit-valid accept
    counts 100000000000000000001 99999999999999999999 1; expect "$work/r" precision-20digit reject
    report 1 0 > "$work/r"; expect "$work/r" report-valid accept
    for kind in $(seq 1 15); do
        report 3 "$kind" > "$work/r"; expect "$work/r" "report-kind-$kind" reject
    done
    for seed in $(seq 1 "$seeds"); do
        report "$seed" 0 > "$work/r"; expect "$work/r" report-valid accept
        kind=$((1 + seed % 15))
        report "$seed" "$kind" > "$work/r"; expect "$work/r" "report-kind-$kind" reject
    done
    echo "reports: 5 precision cases, $((seeds + 1)) valid accepted, $((seeds + 15)) corrupted rejected"
}

# Runs the workflow on manifest body $2 in a fresh directory $1, writing nf.log there.
# A stub run unless $3 names an fqdup binary. A run that cannot start or times out is a harness error.
nf() {
    mkdir "$work/$1"
    printf '%b' "$2" > "$work/$1/samples.tsv"
    local mode=(-stub) rc=0; [ -n "${3:-}" ] && mode=(--fqdup "$3")
    (cd "$work/$1" && timeout "${NF_TIMEOUT:-1800}" nextflow run "$pipeline" -profile test "${mode[@]}" \
        --manifest samples.tsv --outdir out > nf.log 2>&1) || rc=$?
    [ "$rc" = 124 ] && harness "$1: nextflow timed out"
    ! grep -qE '^Error .*:[0-9]+:[0-9]+:|compilation error|command not found' "$work/$1/nf.log" ||
        harness "$1: workflow did not start: $(grep -m1 -E '^Error|compilation error|command not found' "$work/$1/nf.log")"
    return 0
}
rejects() {
    nf "manifest-$1" "$2"
    local got=accept; grep -q SUCCESS "$work/manifest-$1/nf.log" || got=reject
    record "manifest-$1" "reject:$3" "$got:$(grep -m1 -o "$3" "$work/manifest-$1/nf.log" || true)"
    grep -q "$3" "$work/manifest-$1/nf.log" && [ "$got" = reject ] || die "manifest-$1" "expected rejection '$3'"
}
holds() { record "$1" true "$([ "$2" = 0 ] && echo true || echo false)"; [ "$2" = 0 ] || die "$1" "$3"; }

check_manifests() {
    local L=$repo/test h='sample\tr1_files\tr2_files\tlayout\tlibrary_type\n' o=$work/route-mixed/out
    rejects bad-layout "${h}a\t$L/lane1_R1.fq.gz\t\tse\tauto\n" 'Invalid layout for a: se'
    rejects pe-default-no-r2 "sample\tr1_files\tr2_files\na\t$L/lane1_R1.fq.gz\t\n" 'Unpaired raw lanes: a'
    rejects pe-no-r2 "${h}a\t$L/lane1_R1.fq.gz\t\tPE\tauto\n" 'Unpaired raw lanes: a'
    rejects pe-lane-mismatch "${h}a\t$L/lane1_R1.fq.gz,$L/lane2_R1.fq.gz\t$L/lane1_R2.fq.gz\tPE\tauto\n" 'Unpaired raw lanes: a'
    rejects se-with-r2 "${h}a\t$L/lane1_R1.fq.gz\t$L/lane1_R2.fq.gz\tSE\tauto\n" 'Single-end sample with r2_files: a'
    rejects bad-sample "${h}a b\t$L/lane1_R1.fq.gz\t\tSE\tauto\n" 'Invalid sample ID: a b'
    rejects bad-libtype "${h}a\t$L/lane1_R1.fq.gz\t\tSE\tdss\n" 'Invalid library_type for a: dss'
    rejects duplicate "${h}a\t$L/lane1_R1.fq.gz\t\tSE\tauto\na\t$L/lane2_R1.fq.gz\t\tSE\tauto\n" 'Duplicate sample IDs: a'
    nf route-mixed "${h}pe\t$L/lane1_R1.fq.gz,$L/lane2_R1.fq.gz\t$L/lane1_R2.fq.gz,$L/lane2_R2.fq.gz\tPE\tauto\nse\t$L/lane1_R1.fq.gz\t\tSE\tauto\n"
    local ok=0
    [ -s "$o/merge_qc/pe.merge-outcomes.tsv" ] && [ ! -e "$o/merge_qc/se.merge-outcomes.tsv" ] || ok=1
    holds route-merge "$ok" 'only the PE sample is merged'
    ok=0; [ -s "$o/trim_qc/se.trim-report.tsv" ] && [ ! -e "$o/trim_qc/pe.trim-report.tsv" ] || ok=1
    holds route-trim "$ok" 'only the SE sample is trimmed'
    ok=0; [ -s "$o/dedup/pe.dedup.fq.gz" ] && [ -s "$o/dedup/se.dedup.fq.gz" ] || ok=1
    holds route-dedup "$ok" 'both samples reach dedup'
    ok=0; grep -q 'SUCCESS.*completed=8' "$work/route-mixed/nf.log" || ok=1
    holds route-complete "$ok" 'mixed run completes 8 tasks'
    echo "manifests: 8 rejected with the expected error, mixed PE+SE routed"
}

check_e2e() {
    local L=$repo/test o=$work/e2e/out ok=0
    nf e2e "sample\tr1_files\tr2_files\tlayout\nse1\t$L/lane1_R1.fq.gz\t\tSE\n" "$1"
    grep -q 'SUCCESS.*completed=4' "$work/e2e/nf.log" || ok=1
    holds e2e-complete "$ok" 'single-lane SE run completes 4 tasks'
    expect "$o/trim_qc/se1.trim-report.tsv" e2e-report accept
    ok=0; grep -qx 'se_adapter_min_overlap=5' "$o/trim_qc/se1.trim.provenance.txt" || ok=1
    holds e2e-trim-provenance "$ok" 'trim provenance records the trim settings'
    ok=0; grep -qx 'layout=SE' "$o/dedup/se1.derep.provenance.txt" || ok=1
    holds e2e-derep-provenance "$ok" 'derep provenance records the layout'
    echo "e2e: single-lane SE completed; $(awk -F'\t' '$1 == "input" || $1 == "retained" { printf "%s=%s ", $1, $2 }' "$o/trim_qc/se1.trim-report.tsv")"
}

check_reports
check_manifests
if [ "${1:-}" = e2e ]; then check_e2e "${2:?fqdup binary}"; fi
echo PASS
