#!/usr/bin/awk -f
# Exits 0 only for a complete fqdup raw_single_end_v1 trim report whose read counts balance:
# the header line, then key<TAB>value lines with input, retained and dropped each exactly once,
# as unsigned integers without leading zeros, input > 0 and input = retained + dropped.
# Counts are added as decimal strings: awk numbers are doubles and round above 2^53.
function add(a, b,    i, j, carry, d, s) {
    s = ""; carry = 0; i = length(a); j = length(b)
    while (i > 0 || j > 0 || carry) {
        d = (i > 0 ? substr(a, i, 1) : 0) + (j > 0 ? substr(b, j, 1) : 0) + carry
        carry = d >= 10; s = (d % 10) s; i--; j--
    }
    return s
}
BEGIN { FS = "\t" }
NR == 1 { if ($0 != "raw_single_end_v1") bad = 1; next }
NF != 2 { bad = 1; next }
$1 == "input" || $1 == "retained" || $1 == "dropped" {
    if (seen[$1]++ || $2 !~ /^(0|[1-9][0-9]*)$/) bad = 1
    v[$1] = $2
}
END {
    ok = NR > 1 && !bad && seen["input"] == 1 && seen["retained"] == 1 && seen["dropped"] == 1
    exit !(ok && v["input"] != "0" && v["input"] == add(v["retained"], v["dropped"]))
}
