# Tutorial

One run on the bundled test sample, then your own data. About five minutes on a laptop or a single node.

## 1. Install the tools

You need Nextflow 25.10 or later, fqdup 2.1.1 or later, `python3`, `gzip` and `sha256sum`.

```bash
nextflow -version
fqdup merge --help | head -3
```

If fqdup is not on your `PATH`, pass its location with `--fqdup /path/to/fqdup` in every command below.

## 2. Run the test

```bash
nextflow run genomewalker/fqdup-qc-nextflow -profile test --outdir test-results
```

The test has two samples built from the same synthetic double-stranded library with simulated damage:
- `synthetic_L1` is paired-end: 5,000 read pairs in two lanes, through `MERGE_RAW`;
- `synthetic_SE` is the R1 files alone, read as single-end, through `TRIM_RAW_SE`.

The run ends with:

```
[SUCCESS] completed=8 failed=0 cached=0
```

## 3. Look at the results

```bash
find test-results -type f | sort
```

Open them in this order:

1. **Merge yield.** Only merged reads go on to the later stages.

   ```bash
   cat test-results/merge_qc/synthetic_L1.merge-outcomes.tsv
   ```
   ```
   outcome	count
   merged	4409
   unmerged	591
   orphan_r1	0
   orphan_r2	0
   dropped	0
   total	5000
   ```

2. **Damage curve.** Open `test-results/profiles/post_derep/synthetic_L1.profile-post.html` in a browser. Both curves rise towards the read ends; that rise is the deamination signal:

   ![Damage curve of the test sample](img/damage-curve.svg)

3. **Damage verdict.**

   ```bash
   python3 -c '
   import json
   p = json.load(open("test-results/profiles/post_derep/synthetic_L1.profile-post.json"))
   v = p["characterization"]["verdict"]["deamination"]
   e = p["characterization"]["damage_context"]["evidence"]
   print(p["library_call"]["library_type"], v["call"], v["authenticated_by"], v["confidence"])
   print("d_max 5prime", e["d_max_5"], "3prime", e["d_max_3"], "lambda 5prime", e["lambda_5"])
   print("artefacts", p["qc"]["library_artifact_reasons"])'
   ```

   ```
   double-stranded present channel_b high
   d_max 5prime 0.200049 3prime 0.140327 lambda 5prime 0.223379
   artefacts ['adapter_remnant_3prime', 'position0_artifact_3prime']
   ```

   The damage call is confident, but the 3′ end also carries artefact flags. On real data, treat the 3′ values with care when these appear ([step 4](outputs.md#reading-the-profile-json)).

4. **Other damage types.** The profile also tests for oxidation, depurination and other end patterns, one detector per channel:

   ```bash
   python3 -c '
   import json
   p = json.load(open("test-results/profiles/post_derep/synthetic_L1.profile-post.json"))
   for r in p["characterization"]["damage_types"]:
       if not r["name"].endswith("_3prime"):
           print(r["channel"], r["name"], "detected" if r.get("detected") else "-")
   f = p["channels"]["fragmentation"]
   print("median length", f["median_length"], "under 50 bp", f["short_fraction_lt_50bp"])'
   ```

   ```
   A cytosine_deamination detected
   B stop_codon_conversion detected
   C 8_oxog_top_strand -
   D chargaff_gt_asymmetry -
   E purine_enrichment -
   F 8_oxog_complement -
   G cg_stop_enrichment -
   H at_stop_enrichment -
   median length 51.888283 under 50 bp 0.397343
   ```

   Here only deamination (A) and its stop-codon check (B) fire. [Damage types](outputs.md#damage-types-characterizationdamage_types) explains each channel.

5. **Duplication.** Clusters against merged reads:

   ```bash
   zcat test-results/dedup/synthetic_L1.clusters.tsv.gz |
     awk 'NR > 1 { n++; s += $3 } END { printf "%d reads, %d clusters, %.1f%% duplicates\n", s, n, 100 * (1 - n / s) }'
   ```

   ```
   4409 reads, 3312 clusters, 24.9% duplicates
   ```

6. **Single-end yield.** The single-end sample has no merge. Its trim report counts the reads and bases kept and clipped:

   ```bash
   grep -E '^(input|retained|dropped|adapter3_reads)\b' test-results/trim_qc/synthetic_SE.trim-report.tsv
   ```

   ```
   input	5000
   retained	5000
   adapter3_reads	4133
   dropped	0
   ```

   Its profiles and dedup files have the same layout as the paired-end sample's.

[Outputs](outputs.md) explains every file and field, and how to read the profile.

## 4. Run your own samples

Write a tab-separated manifest, one row per sample. List each lane's R1 and R2 files, comma-separated and in the same order. For a single-end sample, set `layout` to `SE`, list its lanes under `r1_files` and leave `r2_files` empty. A `PE` row (the default) with no R2 files stops the run.

```
sample	r1_files	r2_files	layout	library_type	source_id
lib01	/data/lib01_L001_R1.fq.gz,/data/lib01_L002_R1.fq.gz	/data/lib01_L001_R2.fq.gz,/data/lib01_L002_R2.fq.gz	PE	ds	S-0001
lib02	/data/lib02_L001_R1.fq.gz	/data/lib02_L001_R2.fq.gz	PE	auto	S-0002
lib03	/data/lib03_L001.fq.gz,/data/lib03_L002.fq.gz		SE	ss	S-0003
blank01	/data/blank01_L001_R1.fq.gz	/data/blank01_L001_R2.fq.gz	PE	ds	BLANK-01
```

- `sample` names the output files. Use letters, digits, `_`, `.` and `-`.
- Single-end reads are cut where the adapter starts. A read whose molecule was shorter than the cycle count then ends at the molecule's 3′ end; a read from a longer molecule does not. Damage at the 3′ end is measured only on the first kind. Set `--se_adapter3` if your libraries use an adapter other than TruSeq.
- `library_type` is `ds`, `ss`, `auto` or `unknown`. Declare it when you know it. A declared type is checked against fqdup's detection at merge and stops the run if the two clearly disagree. Single-end trimming does not use it; the profiles and derep do.
- Add any other columns you need, such as source IDs, library IDs or a control flag. The workflow ignores them for processing and copies them verbatim into `merge_qc/<sample>.merge.provenance.txt`, so every output can be traced back to your records. Single-end rows go to `trim_qc/<sample>.trim.provenance.txt`.

Then run:

```bash
nextflow run genomewalker/fqdup-qc-nextflow --manifest samples.tsv --outdir results -with-trace results/trace.tsv
```

Add `-resume` to continue an interrupted run without redoing finished tasks.

## 5. Resources and where to run

Each stage asks for the resources below. Each tool's thread count follows the CPUs given to its task.

| Stage | CPUs | Memory | Time |
|---|---|---|---|
| `MERGE_RAW` | 8 | 16 GB | 24 h |
| `TRIM_RAW_SE` | 8 | 4 GB | 24 h |
| `PROFILE_PRE`, `PROFILE_POST` | 8 | 16 GB | 12 h |
| `SORT_DEREP` | 32 | 8 GB per GB of merged FASTQ, between 60 and 400 GB | 72 h |

Sort keeps its buffers within 90% of the task's memory. Derep stops with an error, before its heaviest phase, if its estimated peak would exceed that limit. If that happens, give `SORT_DEREP` more memory.

To change these, or to run on a cluster, write a small config and pass it with `-c`:

```groovy
// site.config
params.fqdup = '/opt/fqdup/bin/fqdup'
params.fqdup_sha256 = '<sha256 of that binary>'   // refuse any other build
profiles {
    cluster {
        process.executor = 'slurm'
        process.queue = 'normal'
        process { withName: 'SORT_DEREP' { cpus = 16; memory = '120 GB' } }
    }
}
```

```bash
nextflow run genomewalker/fqdup-qc-nextflow -c site.config -profile cluster --manifest samples.tsv
```

Any Nextflow executor works the same way: Slurm, SGE, HyperQueue (`process.executor = 'hq'`), or cloud batch.

## 6. Pinning the fqdup build

For results you will publish, set `--fqdup_sha256`. Every stage then hashes the fqdup binary and stops if it is not the build you pinned. Each provenance file records the hash either way, and `--fqdup_revision` adds a readable label such as a git commit.
