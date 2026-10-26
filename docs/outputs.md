# Outputs and how to read them

Every file is named `<sample>.<kind>`. The example values come from the bundled test sample (`-profile test`), a synthetic double-stranded library with simulated ancient-DNA damage.

```
results/
  merge_qc/              merge.json  damage.json  subst.tsv  merge.provenance.txt  fastp-reports.tar.gz
  profiles/pre_derep/    profile-pre.json  profile-pre.html  profile-pre.tsv  pre.tdmg
  dedup/                 dedup.fq.gz  clusters.tsv.gz  derep.provenance.txt
  profiles/post_derep/   profile-post.json  profile-post.html  profile-post.tsv  post.tdmg
```

## 1. Merge (`merge_qc/`)

`fqdup merge` finds the overlap between R1 and the reverse complement of R2 and collapses each overlapping pair into one read. In ancient DNA most inserts are shorter than the read length, so most usable pairs overlap completely. All lanes of a sample go into one merged FASTQ. That FASTQ is passed to the next stages but is not published.

**Pairs that do not merge are not kept.** The workflow does not ask fqdup for unmerged or orphan reads, so everything downstream is merged reads only.

### `merge.json`: merge report

| Field | Meaning |
|---|---|
| `counts.pairs_in` | read pairs read from all lanes |
| `counts.merged` | pairs collapsed into one merged read; these go downstream |
| `counts.unmerged` | pairs that passed QC but had no acceptable overlap (discarded) |
| `counts.orphan_r1`, `orphan_r2` | pairs where only one mate survived QC (discarded) |
| `counts.dropped` | pairs where both mates failed QC, for example shorter than `--merge_min_length` after clipping, or adapter dimers |
| `counts.balanced` | `true` when the counts add up exactly; `false` means a bookkeeping error and the run should not be trusted |
| `library.type`, `type_source`, `type_confidence` | `ds` or `ss`; `declared` if you set it in the manifest, `detected` if fqdup inferred it, with the detector's confidence |
| `library.udg`, `skip_terminal` | a quick treatment call from overlap mismatches: `none` (both ends damaged), `half` (5′ end damaged, 3′ not), `full` (neither end). This is a coarse screen; use the profile for damage |
| `adapters` | adapter sequences fqdup used or detected |
| `complexity_gate` | data-derived entropy and dominant-base thresholds of the low-complexity filter |
| `params` | the merge settings used |
| `drops.by_reason` | reads dropped as technical artefacts, by reason |

Example: 5,000 pairs in, 1,567 merged, 516 unmerged, 35 orphans and 2,882 dropped. The test data contains many very short inserts, so a large `dropped` count is expected there. On real libraries, read `merged / pairs_in` as the merge yield. A sudden drop in yield across samples usually means adapter or length problems.

### `damage.json`: damage measured on the overlaps

This is a damage profile built from the overlapping region of each pair before consensus, in the same schema as the profile JSON below. It sees damage directly as disagreement between the two mates, independent of any reference. Read it the way you read the profile JSON in section 2. It is most useful as a second opinion, for example when the merged-read profile is ambiguous.

### `subst.tsv`: overlap substitution matrix

Comment lines start with `#`. They give the pair and base totals and explain the encoding. The table has the columns:

| Column | Meaning |
|---|---|
| `strand` | `fwd`: position counted from the 5′ end of the insert. `rev`: counted from the 3′ end |
| `pos` | distance from that end, 0 = terminal base |
| `r1` | base in R1 |
| `rc2` | base in the reverse complement of R2 at the same insert position |
| `count` | number of overlap positions with that base pair |

Where the two mates disagree, one of them carries damage or an error. `r1=T`, `rc2=C` at small `fwd` positions is 5′ deamination. `r1=G`, `rc2=A` at small `rev` positions is 3′ deamination in a double-stranded library. `r1=T`, `rc2=G` at any position is the 8-oxoguanine signature.

### `fastp-reports.tar.gz`

The per-lane fastp JSON and HTML reports, with `--trim_mode fastp`. With `--trim_mode fqdup` it holds only a README saying fastp was skipped.

### `merge.provenance.txt`

`key=value` lines: the sample, the library type as declared in the manifest, the lane count, the trim mode, the fqdup revision label, the SHA-256 of the fqdup and fastp binaries that ran, a timestamp and the host. It also contains every manifest column, verbatim, as `manifest.<column>=<value>`. Put your own identifiers in the manifest (source sample ID, library ID, control flag) and they are recorded here unchanged.

## 2. Damage profile before dedup (`profiles/pre_derep/`)

`fqdup profile` estimates ancient-DNA damage from the merged reads without a reference genome. It compares base composition at the read ends with the read interior. Deamination of cytosine in single-stranded overhangs makes C read as T near the 5′ end, and G read as A near the 3′ end of double-stranded libraries. The damage engine is [libtaph](https://github.com/genomewalker/libtaph). The profile runs before and after dedup; see [section 4](#4-damage-profile-after-dedup-profilespost_derep) for why.

### `profile-pre.html`: interactive report

Open this first. It is a self-contained page: the plotting library is embedded, so it opens offline. It plots the terminal substitution curves, the length-stratified fits and the QC panels, and it gives the headline numbers described below.

### `profile-pre.tsv`: per-position table

One row per position from each end (positions 0–14):

| Column | Meaning |
|---|---|
| `pos` | distance from the read end, 0 = terminal base |
| `freq5` | T / (T + C) at that position from the 5′ end |
| `freq3` | A / (A + G) at that position from the 3′ end (T / (T + C) for single-stranded libraries) |
| `freq5_gt` | T / (T + G) from the 5′ end, the 8-oxoguanine marker |
| `mask` | 1 if the position is damaged enough that downstream tools should mask it (excess deamination probability above 0.05) |
| `cov5`, `cov3`, `cov5_gt` | number of bases behind each frequency |

A value of `-1` means too little coverage to estimate.

These are composition fractions, not substitution rates, because there is no reference. Read them as a curve. **Damage shows as `freq5` and `freq3` rising towards position 0 and flattening to a plateau in the interior.** The size of the rise above the plateau matters; the absolute level does not, because it reflects the genome's base composition. In the test sample before dedup, `freq5` goes from about 0.50 at positions 13–14 to 0.58 at position 0, and positions 0–1 are masked.

### `profile-pre.json`: full profile

A large JSON. Read it from the top down and stop when you have your answer.

**1. The verdict: `characterization.verdict`.** This is the authoritative library-level call.

| Field | Meaning |
|---|---|
| `deamination.call` | `present` or `absent`. Is there authenticated age-related C→T damage? This is the field to use for "is this library ancient?" |
| `deamination.authenticated_by` | which independent evidence carried the call: `pi_estimate`, `channel_b` (stop-codon-context cross-check, independent of composition), `d5_decay_fit`, `pinned_shallow_ct`, or `none` |
| `deamination.confidence` | `high`, `shallow`, `low_shallow` or `none` |
| `deamination.d_max_5prime` | the 5′ terminal damage amplitude behind the call |
| `oxidation.call` | 8-oxoguanine (G→T) damage. Real but not ancient-specific: sample preparation also causes it |
| `preservation.tier` | an aggregate preservation index, separate from the deamination call |

**2. The magnitude: `characterization.damage_context.evidence`.**

- `d_max_5`, `d_max_3`: terminal damage above the interior baseline at each end. As a rough guide for untreated libraries, ≥ 0.20 is strong, 0.08–0.20 is typical of Holocene material, 0.02–0.08 is weak and needs the verdict's support, and < 0.02 is effectively absent. These are reference-free estimates and read somewhat below alignment-based tools such as metaDMG.
- `lambda_5`, `lambda_3`: how fast damage decays from the end. About 0.2–0.4 is typical. Smaller values mean damage reaches further into the read.

For a single number per library, use `characterization.damage_magnitude_recommended.absolute_estimate`. Check `absolute_estimate_valid`: `0.0` with `valid: true` means "no damage within noise", while `null` means "could not be determined".

**3. The library type: `library_call`.** `library_type` is `double-stranded`, `single-stranded` or unknown, decided by competing models whose BIC scores and posterior probabilities are in `bic`. `library_call_overridden: true` means you declared the type in the manifest.

**4. QC flags: `qc`.**
- `library_artifact_contaminated` with `library_artifact_reasons` flags adapter remnants or ligation artefacts at the read ends. When it is set, terminal damage numbers can be inflated or deflated; look at the HTML curves.
- `library_qc.position0_artifact_5prime` / `_3prime` mark a depleted terminal base with rising positions 1–4. That shape is a clipping or ligation artefact, not damage, and libtaph shifts its fit inward when it sees it.

**5. Detail: `channels` and `diagnostics`.** Per-mechanism detectors and the raw tables behind them. `channels.deamination` is a sensitive detector. It can read `AUTHENTICATED` while the verdict says `absent`, because the verdict is deliberately stricter. Use the verdict for decisions and the channels when you want to understand why.

Test sample, before dedup: 1,567 reads; double-stranded; deamination `present`, authenticated by `channel_b`, confidence `high`; `d_max_5` 0.18, `d_max_3` 0.14, `lambda_5` 0.16; oxidation `none`.

Do not judge a library on one number. A high `d_max` with a confident verdict, consistent curves in the HTML and no artefact flags is a solid call. A high `d_max` with `library_artifact_contaminated` set, or with the verdict `absent`, needs a closer look.

### `pre.tdmg`: binary damage summary

A small binary sidecar (`TDM1` header) holding the fitted damage model. [DART](https://github.com/genomewalker/dart) reads it with `--library-profile`, so it reuses this fit instead of re-estimating damage. It is not meant to be read by hand.

## 3. Dedup (`dedup/`)

`fqdup sort` orders the merged reads, then `fqdup derep` collapses PCR duplicates: reads with the same sequence, including reverse-complement copies, become one cluster. One representative per cluster is written. Dedup runs per sample: lanes are pooled (they were merged together), but samples are never pooled. Damage is fitted and reported (`--damage report`) but not used for clustering.

### `dedup.fq.gz`

The deduplicated merged reads, one per cluster. Read names are the original R1 names.

### `clusters.tsv.gz`

One row per cluster:

| Column | Meaning |
|---|---|
| `hash` | 128-bit fingerprint of the cluster's sequence, in hexadecimal |
| `seq_len` | length of the representative read |
| `count` | number of input reads in the cluster |

The `count` column sums to the number of merged reads. The duplication rate is `1 − clusters / merged reads`. The test sample has 1,567 reads in 1,179 clusters, so 24.8% were duplicates, and 297 clusters had more than one read. The distribution of `count` shows library complexity: many large clusters mean the library was sequenced close to exhaustion.

### `derep.provenance.txt`

The sample, the merged input, the library type passed to derep, the fqdup SHA-256 and revision label, the memory limit, a timestamp and the host.

## 4. Damage profile after dedup (`profiles/post_derep/`)

The same four files as section 2, computed on the deduplicated reads. **Use the post-dedup profile when reporting damage for a library.** PCR duplicates give extra weight to whichever molecules were over-amplified, and dedup removes that bias. Compare it with the pre-dedup profile: the two should agree closely. A large shift means duplicates were concentrated in a damage-biased part of the library.

Test sample, after dedup: 1,179 reads; deamination `present` (`channel_b`, `high`); `d_max_5` 0.20, `d_max_3` 0.16, close to the pre-dedup values.

## Library type through the stages

- **Declared** (`ds` or `ss` in the manifest): merge runs in strict mode and stops if the declared type contradicts a confident detection. The type is then passed to both profiles and to derep.
- **`auto` or `unknown`**: each stage detects the type on its own input. Merge's call is recorded in `merge.json` but is not passed on. Check that `merge.json` `library.type` and `library_call.library_type` in both profiles agree. A disagreement is worth investigating, and declaring the type in the manifest settles it.

## What the built-in checks cover

Each stage fails the task, rather than publishing, if:
- a gzip file fails `gzip -t`;
- a FASTQ does not start with `@`;
- a JSON does not parse;
- an expected output is empty;
- a binary's SHA-256 does not match a pin you set.

These checks catch truncated or empty files. They are not a full validation of every FASTQ record.
