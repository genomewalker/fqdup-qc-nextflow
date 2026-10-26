# Tutorial

One run on the bundled test sample, then your own data. About five minutes on a laptop or a single node.

## 1. Install the tools

You need Nextflow 25.10 or later, fqdup, `python3`, `gzip` and `sha256sum`. fastp is only needed if you trim with it (the default for real data).

```bash
nextflow -version
fqdup merge --help | head -3
```

If fqdup is not on your `PATH`, pass its location with `--fqdup /path/to/fqdup` in every command below.

## 2. Run the test

```bash
nextflow run genomewalker/fqdup-qc-nextflow -profile test --outdir test-results
```

The test sample is 5,000 synthetic read pairs in two lanes: a double-stranded library with simulated damage. The run ends with:

```
[SUCCESS] completed=4 failed=0 cached=0
```

## 3. Look at the results

```bash
find test-results -type f | sort
```

Open them in this order:

1. **Merge yield.** Only merged reads go on to the later stages.

   ```bash
   python3 -c 'import json; print(json.load(open("test-results/merge_qc/synthetic_L1.merge.json"))["counts"])'
   ```
   ```
   {'pairs_in': 5000, 'merged': 1567, 'unmerged': 516, 'orphan_r1': 18, 'orphan_r2': 17, 'dropped': 2882, 'balanced': True}
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
   print("d_max 5prime", e["d_max_5"], "3prime", e["d_max_3"], "lambda 5prime", e["lambda_5"])'
   ```

   ```
   double-stranded present channel_b high
   d_max 5prime 0.204246 3prime 0.163931 lambda 5prime 0.154537
   ```

4. **Duplication.** Clusters against merged reads:

   ```bash
   zcat test-results/dedup/synthetic_L1.clusters.tsv.gz |
     awk 'NR > 1 { n++; s += $3 } END { printf "%d reads, %d clusters, %.1f%% duplicates\n", s, n, 100 * (1 - n / s) }'
   ```

   ```
   1567 reads, 1179 clusters, 24.8% duplicates
   ```

[Outputs](outputs.md) explains every file and field, and how to read the profile.

## 4. Run your own samples

Write a tab-separated manifest, one row per sample. List each lane's R1 and R2 files, comma-separated and in the same order:

```
sample	r1_files	r2_files	library_type	source_id
lib01	/data/lib01_L001_R1.fq.gz,/data/lib01_L002_R1.fq.gz	/data/lib01_L001_R2.fq.gz,/data/lib01_L002_R2.fq.gz	ds	S-0001
lib02	/data/lib02_L001_R1.fq.gz	/data/lib02_L001_R2.fq.gz	auto	S-0002
blank01	/data/blank01_L001_R1.fq.gz	/data/blank01_L001_R2.fq.gz	ds	BLANK-01
```

- `sample` names the output files. Use letters, digits, `_`, `.` and `-`.
- `library_type` is `ds`, `ss`, `auto` or `unknown`. Declare it when you know it. A declared type is checked against fqdup's detection at merge and stops the run if the two clearly disagree.
- Add any other columns you need, such as source IDs, library IDs or a control flag. The workflow ignores them for processing and copies them verbatim into `merge_qc/<sample>.merge.provenance.txt`, so every output can be traced back to your records.

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
