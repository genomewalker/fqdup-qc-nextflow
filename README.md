# fqdup-qc-nextflow

A Nextflow workflow that runs fqdup QC on paired-end ancient-DNA libraries. Each sample goes through four stages:

1. `MERGE_RAW`: optional fastp adapter and quality trimming per lane, then `fqdup merge` of all lanes into one merged FASTQ. Writes a merge report, a damage estimate and the substitution table.
2. `PROFILE_PRE`: `fqdup profile` damage profile of the merged reads.
3. `SORT_DEREP`: `fqdup sort` then `fqdup derep`. Writes the deduplicated FASTQ and a cluster table.
4. `PROFILE_POST`: `fqdup profile` damage profile of the deduplicated reads.

Every stage checks its outputs (gzip integrity, FASTQ framing, valid JSON) before it finishes. Merge and dedup write a provenance file that records the binary SHA-256, the host and a timestamp.

## Requirements

- Nextflow 25.10 or later.
- [fqdup](https://github.com/genomewalker/fqdup) on `PATH`, or pass `--fqdup /path/to/fqdup`.
- [fastp](https://github.com/OpenGene/fastp), only when `--trim_mode fastp` (the default).
- `python3`, `gzip` and `sha256sum`.

## Test

```bash
nextflow run genomewalker/fqdup-qc-nextflow -profile test --fqdup /path/to/fqdup
```

The test runs one synthetic sample: 5,000 read pairs in two lanes, from fqdup's simulated development set. It uses `--trim_mode fqdup`, so fastp is not needed.

## Usage

```bash
nextflow run genomewalker/fqdup-qc-nextflow --manifest samples.tsv --outdir results
```

The manifest is a tab-separated file with a header row:

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | sample ID: letters, digits, `_`, `.`, `-` |
| `r1_files` | yes | comma-separated R1 FASTQs, one per lane |
| `r2_files` | yes | comma-separated R2 FASTQs, in the same lane order |
| `library_type` | no | `ds`, `ss`, `auto` or `unknown`; defaults to `--library_type` |
| `batch`, `is_control` | no | copied into the merge provenance file |

Relative FASTQ paths resolve against the manifest's directory. Extra columns are ignored. See [test/samples.tsv](test/samples.tsv) for an example.

With `ds` or `ss`, merge runs with that library type in strict mode. With `auto` or `unknown`, fqdup detects the type during merge.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--manifest` | | sample manifest (required) |
| `--outdir` | `results` | output directory |
| `--fqdup` | `fqdup` | fqdup binary |
| `--fqdup_sha256` | | when set, every stage refuses any other fqdup binary |
| `--fqdup_revision` | `unrecorded` | label written to the provenance files |
| `--fastp` | `fastp` | fastp binary |
| `--fastp_sha256` | | when set, refuses any other fastp binary |
| `--trim_mode` | `fastp` | `fastp` trims before merge; `fqdup` uses only fqdup's merge filters |
| `--merge_min_length` | `16` | `fqdup merge --min-length` |
| `--merge_min_entropy` | `0.0` | `fqdup merge --min-entropy` |
| `--library_type` | `auto` | default when the manifest has no `library_type` column |

## Outputs

```
results/
  merge_qc/              <sample>.merge.json, .damage.json, .subst.tsv, .merge.provenance.txt, .fastp-reports.tar.gz
  profiles/pre_derep/    <sample>.profile-pre.json, .pre.tdmg
  dedup/                 <sample>.dedup.fq.gz, .clusters.tsv.gz, .derep.provenance.txt
  profiles/post_derep/   <sample>.profile-post.json, .post.tdmg
```

The merged FASTQ is an intermediate and is not published.

## Profiles

- `test`: the bundled synthetic sample on the local executor with small resources.
- `slurm`: submits each task to Slurm.

Put site settings in your own config file and load it with `-c site.config`: queue names, other executors such as HyperQueue, and binary paths and pins. Profiles with the same name merge with the ones here.

## License

MIT. See [LICENSE](LICENSE).
