#!/usr/bin/env bash
# Plants one deliberate bug at a time in a copy of the workflow and runs tests/run_checks.sh on it,
# three times each. A mutant is killed only when the run fails the check named for it; a run that
# passes is a survivor, and any other outcome (another check, a harness error) fails the harness.
# An unmutated copy must pass first, in the same environment.
#   tests/mutants.sh             checker and routing mutants (stub runs)
#   tests/mutants.sh e2e FQDUP   the same mutants, with the control also running the real e2e check
# OUT (default: a new directory under TMPDIR, kept) holds every run's files and mutants.tsv.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
out=${OUT:-$(mktemp -d "${TMPDIR:-/tmp}/fqdup-qc-mutants.XXXXXX")}
mkdir -p "$out"
status=0

# id @ file @ sed expression @ check that must fail
mutants=(
    'header@bin/check_trim_report.awk@s/if (\$0 != "raw_single_end_v1") bad = 1; //@report-kind-1'
    'duplicate-key@bin/check_trim_report.awk@s/seen\[\$1\]++ || /(seen[$1] = 1) \&\& 0 || /@report-kind-6'
    'integer@bin/check_trim_report.awk@s/\^(0|\[1-9\]\[0-9\]\*)\$/^[-+0-9.]*$/@report-kind-9'
    'field-count@bin/check_trim_report.awk@s/^NF != 2/NF < 2/@report-kind-13'
    'conservation@bin/check_trim_report.awk@s/ \&\& v\["input"\] == add(v\["retained"\], v\["dropped"\])//@precision-2p53'
    'float-conservation@bin/check_trim_report.awk@s/v\["input"\] == add(v\["retained"\], v\["dropped"\])/v["input"] + 0 == v["retained"] + v["dropped"]/@precision-2p53'
    'input-positive@bin/check_trim_report.awk@s/v\["input"\] != "0" \&\& //@report-kind-15'
    'layout-default-se@main.nf@s/row.layout ?: '"'PE'"'/row.layout ?: '"'SE'"'/@manifest-pe-default-no-r2'
    'se-with-r2@main.nf@/Single-end sample with r2_files/d@manifest-se-with-r2'
    'route-all-paired@main.nf@s/paired: meta.layout == '"'PE'"'/paired: true/@route-merge'
    'duplicate-sample@main.nf@s/if (dup) error/if (false) error/@manifest-duplicate'
)

copy() {
    mkdir "$out/$1"
    cp -R "$repo/main.nf" "$repo/nextflow.config" "$repo/bin" "$repo/test" "$out/$1/"
}
run() {  # copy id, run number; prints the exit code
    local rc=0
    CHECKER=$out/$1/bin/check_trim_report.awk PIPELINE=$out/$1 OUT=$out/$1/run$2 SEEDS=${SEEDS:-100} \
        bash "$repo/tests/run_checks.sh" "${@:3}" > "$out/$1/run$2.log" 2>&1 || rc=$?
    echo "$rc"
}

printf 'mutant\trun\texit\tresult\tfirst_failure\n' > "$out/mutants.tsv"
copy control
rc=$(run control 1 "$@")
printf 'control\t1\t%s\t%s\t%s\n' "$rc" "$([ "$rc" = 0 ] && echo pass || echo HARNESS)" \
    "$(grep -m1 -E '^(FAIL|HARNESS)' "$out/control/run1.log" || true)" >> "$out/mutants.tsv"
[ "$rc" = 0 ] || { echo "HARNESS: unmutated control failed, see $out/control/run1.log" >&2; exit 2; }

for m in "${mutants[@]}"; do
    IFS='@' read -r id file expr want <<< "$m"
    copy "$id"
    sed -i.orig -e "$expr" "$out/$id/$file" || { echo "HARNESS: $id: sed failed" >&2; exit 2; }
    cmp -s "$out/$id/$file" "$out/$id/$file.orig" && { echo "HARNESS: $id: sed did not change $file" >&2; exit 2; }
    for n in 1 2 3; do
        rc=$(run "$id" "$n" "$@")
        first=$(grep -m1 -E '^(FAIL|HARNESS)' "$out/$id/run$n.log" || true)
        if [ "$rc" = 0 ]; then result=SURVIVED; status=1
        elif [ "$rc" = 1 ] && [[ $first == "FAIL[$want]:"* ]]; then result=killed
        else result=HARNESS; status=1
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$id" "$n" "$rc" "$result" "$first" >> "$out/mutants.tsv"
        echo "$id run $n: $result ($first)"
    done
done
[ "$status" = 0 ] && echo "PASS: every mutant killed 3/3 by its named check; logs in $out"
exit "$status"
