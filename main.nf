nextflow.enable.dsl = 2

params.manifest = null
params.outdir = 'results'
params.fqdup = 'fqdup'
params.fqdup_sha256 = null      // set to refuse any other binary
params.fqdup_revision = 'unrecorded'
params.merge_min_length = 16
params.library_type = 'auto'

// library_type is a per-sample property: a cohort may mix single- and double-stranded preps.
// A declared ds/ss is passed to every stage. For 'auto'/'unknown', profile and derep detect the
// type on their own input.
def libPlain(meta)  { meta.library_type in ['ds', 'ss'] ? "--library-type ${meta.library_type}" : '' }
// Every manifest column, verbatim, as manifest.<column>=<value> lines for the provenance files.
def manifestLines(meta) { meta.columns.collect { k, v -> "manifest.${k}=${v ?: ''}" }.join('\n') }
// Records the binary's SHA-256 in \$var and, when a pin is given, fails on any other binary.
def pinCheck(String var, String tool, String bin, pin) {
    "${var}=\$(sha256sum \"\$(command -v '${bin}')\" | cut -d' ' -f1)\n" +
        (pin ? "[ \"\$${var}\" = '${pin}' ] || { echo '${tool} SHA-256 mismatch' >&2; exit 1; }" : '')
}

process MERGE_RAW {
    tag "$meta.sample"
    cpus 8
    memory '16 GB'
    time '24h'
    publishDir "$params.outdir/merge_qc", mode: 'copy', overwrite: false,
        saveAs: { name -> name.endsWith('.merged.fq.gz') ? null : name }

    input:
    tuple val(meta), path(r1_lanes), path(r2_lanes)

    output:
    tuple val(meta), path('*.merged.fq.gz'), path('*.merge-outcomes.tsv'), path('*.merge.log'),
          path('*.merge.provenance.txt')

    script:
    def r1Array = r1_lanes.collect { "'" + it.toString() + "'" }.join(' ')
    def r2Array = r2_lanes.collect { "'" + it.toString() + "'" }.join(' ')
    """
    set -euo pipefail
    ${pinCheck('got', 'fqdup', params.fqdup, params.fqdup_sha256)}
    R1=( $r1Array )
    R2=( $r2Array )
    [ "\${#R1[@]}" -eq "\${#R2[@]}" ] && [ "\${#R1[@]}" -gt 0 ]
    r1_csv=\$(IFS=,; printf '%s' "\${R1[*]}")
    r2_csv=\$(IFS=,; printf '%s' "\${R2[*]}")
    "$params.fqdup" merge -1 "\$r1_csv" -2 "\$r2_csv" -o merge_out \
      ${libPlain(meta)} --min-length "$params.merge_min_length" \
      --outcomes summary -p ${task.cpus} 2> "${meta.sample}.merge.log" \
      || { cat "${meta.sample}.merge.log" >&2; exit 1; }
    awk -F'\t' 'NR > 1 && \$1 != "total" { s += \$2 } \$1 == "total" { t = \$2 } END { exit !(t > 0 && s == t) }' \
      merge_out/outcomes.summary.tsv || { echo 'Merge outcomes do not add up' >&2; exit 1; }
    cp merge_out/outcomes.summary.tsv "${meta.sample}.merge-outcomes.tsv"
    gzip -c merge_out/merged.fq > "${meta.sample}.merged.fq.gz"
    first=\$(gzip -cd "${meta.sample}.merged.fq.gz" | head -c 1 || true)
    [ "\$first" = '@' ] || { echo 'Merged FASTQ framing invalid' >&2; exit 1; }
    printf 'sample=%s\\nlibrary_type=%s\\nraw_lane_pairs=%s\\nfqdup_revision=%s\\nfqdup_sha256=%s\\nstamped=%s\\nhost=%s\\n' \
      "$meta.sample" "$meta.library_type" "$meta.lanes" \
      "$params.fqdup_revision" "\$got" "\$(date -Is)" "\$(hostname)" \
      > "${meta.sample}.merge.provenance.txt"
    cat >> "${meta.sample}.merge.provenance.txt" <<'MANIFEST_END'
${manifestLines(meta)}
MANIFEST_END
    rm -rf merge_out
    """

    stub:
    """
    printf '@stub\\nA\\n+\\n!\\n' | gzip -c > "${meta.sample}.merged.fq.gz"
    printf 'outcome\\tcount\\n' > "${meta.sample}.merge-outcomes.tsv"
    printf 'stub\\n' > "${meta.sample}.merge.log"
    printf 'stub\\n' > "${meta.sample}.merge.provenance.txt"
    """
}

process PROFILE_PRE {
    tag "$meta.sample"
    cpus 8
    memory '16 GB'
    time '12h'
    publishDir "$params.outdir/profiles/pre_derep", mode: 'copy', overwrite: false,
        saveAs: { name -> name ==~ /.*\.(profile-pre\.(json|html|tsv)|pre\.tdmg)/ ? name : null }

    input:
    tuple val(meta), path(merged), path(merge_outcomes), path(merge_log),
          path(merge_provenance)

    output:
    tuple val(meta), path(merged), path('*.profile-pre.json'), path('*.pre.tdmg'),
          path('*.profile-pre.html'), path('*.profile-pre.tsv')

    script:
    """
    set -euo pipefail
    ${pinCheck('got', 'fqdup', params.fqdup, params.fqdup_sha256)}
    "$params.fqdup" profile -i "$merged" \
      --json "${meta.sample}.profile-pre.json" \
      --damage-summary "${meta.sample}.pre.tdmg" \
      --html "${meta.sample}.profile-pre.html" --tsv "${meta.sample}.profile-pre.tsv" \
      ${libPlain(meta)} --length-bins auto -p ${task.cpus}
    python3 -c 'import json, sys; json.load(open(sys.argv[1], encoding="utf-8"))' "${meta.sample}.profile-pre.json"
    test -s "${meta.sample}.pre.tdmg" && test -s "${meta.sample}.profile-pre.html" && test -s "${meta.sample}.profile-pre.tsv"
    """

    stub:
    """
    printf '{"stub":true}\\n' > "${meta.sample}.profile-pre.json"
    printf 'stub\\n' > "${meta.sample}.pre.tdmg"
    printf 'stub\\n' > "${meta.sample}.profile-pre.html"
    printf 'stub\\n' > "${meta.sample}.profile-pre.tsv"
    """
}

process SORT_DEREP {
    tag "$meta.sample"
    cpus 32
    memory { Math.max(60, Math.min(400, 8 * (merged.size().intdiv(1000000000L) as int))).GB }
    time '72h'
    publishDir "$params.outdir/dedup", mode: 'copy', overwrite: false

    input:
    tuple val(meta), path(merged), path(pre_json), path(pre_tdmg), path(pre_html), path(pre_tsv)

    output:
    tuple val(meta), path('*.dedup.fq.gz'), path('*.clusters.tsv.gz'),
          path('*.derep.provenance.txt')

    script:
    def maxMemGb = Math.max(1, (task.memory.toGiga() * 0.9) as int)
    """
    set -euo pipefail
    ${pinCheck('got', 'fqdup', params.fqdup, params.fqdup_sha256)}
    mkdir -p tmp
    "$params.fqdup" sort -i "$merged" -o sorted.fq.gz \
      --max-memory ${maxMemGb}G -p ${task.cpus} -t tmp
    gzip -t sorted.fq.gz
    first=\$(gzip -cd sorted.fq.gz | head -c 1 || true)
    [ "\$first" = '@' ] || { echo 'Sorted FASTQ framing invalid' >&2; exit 1; }
    "$params.fqdup" derep -i sorted.fq.gz \
      -o "${meta.sample}.dedup.fq.gz" \
      -c "${meta.sample}.clusters.tsv.gz" \
      --damage report --library-type "$meta.library_type" \
      --max-mem "$maxMemGb" -p ${task.cpus}
    gzip -t "${meta.sample}.dedup.fq.gz" "${meta.sample}.clusters.tsv.gz"
    first=\$(gzip -cd "${meta.sample}.dedup.fq.gz" | head -c 1 || true)
    [ "\$first" = '@' ] || { echo 'Dedup FASTQ framing invalid' >&2; exit 1; }
    printf 'sample=%s\\ninput_merged=%s\\nlibrary_type=%s\\nfqdup_revision=%s\\nfqdup_sha256=%s\\nmax_mem_gb=%s\\nstamped=%s\\nhost=%s\\n' \
      "$meta.sample" "$merged" "$meta.library_type" "$params.fqdup_revision" "\$got" "$maxMemGb" "\$(date -Is)" "\$(hostname)" \
      > "${meta.sample}.derep.provenance.txt"
    rm -f sorted.fq.gz
    rm -rf tmp
    """

    stub:
    """
    printf '@stub\\nA\\n+\\n!\\n' | gzip -c > "${meta.sample}.dedup.fq.gz"
    printf 'stub\\n' | gzip -c > "${meta.sample}.clusters.tsv.gz"
    printf 'stub\\n' > "${meta.sample}.derep.provenance.txt"
    """
}

process PROFILE_POST {
    tag "$meta.sample"
    cpus 8
    memory '16 GB'
    time '12h'
    publishDir "$params.outdir/profiles/post_derep", mode: 'copy', overwrite: false,
        saveAs: { name -> name ==~ /.*\.(profile-post\.(json|html|tsv)|post\.tdmg)/ ? name : null }

    input:
    tuple val(meta), path(dedup), path(clusters), path(derep_provenance)

    output:
    tuple val(meta), path('*.profile-post.json'), path('*.post.tdmg'),
          path('*.profile-post.html'), path('*.profile-post.tsv')

    script:
    """
    set -euo pipefail
    ${pinCheck('got', 'fqdup', params.fqdup, params.fqdup_sha256)}
    "$params.fqdup" profile -i "$dedup" \
      --json "${meta.sample}.profile-post.json" \
      --damage-summary "${meta.sample}.post.tdmg" \
      --html "${meta.sample}.profile-post.html" --tsv "${meta.sample}.profile-post.tsv" \
      ${libPlain(meta)} --length-bins auto -p ${task.cpus}
    python3 -c 'import json, sys; json.load(open(sys.argv[1], encoding="utf-8"))' "${meta.sample}.profile-post.json"
    test -s "${meta.sample}.post.tdmg" && test -s "${meta.sample}.profile-post.html" && test -s "${meta.sample}.profile-post.tsv"
    """

    stub:
    """
    printf '{"stub":true}\\n' > "${meta.sample}.profile-post.json"
    printf 'stub\\n' > "${meta.sample}.post.tdmg"
    printf 'stub\\n' > "${meta.sample}.profile-post.html"
    printf 'stub\\n' > "${meta.sample}.profile-post.tsv"
    """
}

workflow {
    if (!params.manifest) error 'Provide --manifest: TSV with sample, r1_files, r2_files [, library_type, any other columns]'
    if (!(params.library_type in ['auto', 'ds', 'ss', 'unknown'])) error 'Invalid library_type'
    // Relative FASTQ paths resolve against the manifest's directory.
    def manifestDir = file(params.manifest).parent
    def lane = { String p -> file(p.startsWith('/') ? p : manifestDir.resolve(p).toString(), checkIfExists: true) }
    reads = Channel.fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def id = row.sample?.toString()
            if (!(id ==~ /[A-Za-z0-9_.-]+/)) error "Invalid sample ID: $id"
            def r1 = row.r1_files.split(',').collect(lane)
            def r2 = row.r2_files.split(',').collect(lane)
            if (!r1 || r1.size() != r2.size()) error "Unpaired raw lanes: $id"
            def lt = (row.library_type ?: params.library_type).toString().trim()
            if (!(lt in ['auto', 'ds', 'ss', 'unknown'])) error "Invalid library_type for $id: $lt"
            def meta = [sample: id, library_type: lt, lanes: r1.size(), columns: new LinkedHashMap(row)]
            [meta, r1, r2]
        }
    MERGE_RAW(reads)
    PROFILE_PRE(MERGE_RAW.out)
    SORT_DEREP(PROFILE_PRE.out)
    PROFILE_POST(SORT_DEREP.out)
}
