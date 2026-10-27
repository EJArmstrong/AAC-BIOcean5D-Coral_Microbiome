# BIOcean5D coral microbiome — methods and results walkthrough

*Working document, begun 26 October 2026. Each section summarises one step of the analysis pipeline as
it currently stands, with methods, recorded results, and an explicit note of what the unbalanced design
does and does not support. Numbers are quoted from recorded output (knitted `.md` files, saved tables,
`Master_SampleKey.xlsx`); anything not yet recorded is marked as such rather than estimated.*

---

## 0. Experimental design as sampled

### 0.1 Intended design

Ten source colonies per host species — five at Opunohu Bay (turbid, impacted) and five at Tiahura
(pristine) — were fragmented and reciprocally transplanted between bays, with sampling at 2 weeks
(short-term shifts), 4 months (longer-term shifts) and immediately after an acute heat challenge
(stress-driven shifts).

### 0.2 Sample-name convention

Sample IDs read `<origin bay><transplant type><colony><fragment><species>`: `TBC01F001` is fragment 001
of colony C01, originating at **Tiahura** and back-transplanted to Tiahura. The second character is `B`
for back transplant (returned to its own bay), `T` for out transplant (moved to the other bay) and `G`
for a day-0 genotype sample. **The `Site` field therefore records the colony's bay of origin, not where
the fragment was when it was sampled.** An out-transplant's actual location is the opposite bay, and no
destination field exists in the metadata — it has to be derived. This distinction governs the
interpretation of every analysis downstream and is handled incorrectly in step 05 (section 5).

### 0.3 Timepoints as recorded

| Node | Day | Date sampled | n (coral) |
|---|---|---|---|
| Genotype (baseline) | 0 | 30–31 Oct 2023 | 31 |
| 2-Weeks | 14 | 17 Nov 2023 | 39 |
| 4-Months | 120 | 13 Mar 2024 | 30 |
| Pre-HS | 121 | 15 Mar 2024 | 29 |
| Post-HS | 122 | 15 Mar 2024 | 29 |

Two points of correction against the design as remembered. The second field timepoint is **4 months
(day 120)**, not one month. And `Pre-HS` is a *separate sampling two days after* the 4-month field
collection, taken after a short period in common-garden aquaria — so pooling the two, as step 05 does,
mixes field-collected with aquarium-held fragments and places two fragments from each colony in a single
timeline node.

### 0.4 What was actually obtained

Two separate kinds of loss have to be kept apart: fragments that **died** and were never collected, and
fragments that were collected but yielded **no usable library** and were dropped at the empty-well /
low-read filter in step 01. They fall on different species.

![Design grid: fragments present in the analysed phyloseq object, per origin-to-destination path and timepoint. Hatched cells were lost to mortality; stippled cells were collected but produced no usable library. Pocillopora has no 2-week data at all; Porites loses both paths ending at Tiahura after 2 weeks.](AAC-Walkthrough_fig0_DesignGrid.png)

**All 18 *Pocillopora* 2-week libraries failed.** The fragments were collected and appear in
`Master_SampleKey.xlsx`, but none survived the step-01 read-count filter, so *Pocillopora* has no
short-term timepoint in any analysed object. From 4 months onward it is otherwise complete: all four
origin → destination paths at 4 months, Pre-HS and Post-HS, five fragments per cell except
Opunohu → Opunohu at 4 months (three). Day-0 baselines exist for five colonies per bay, each sampled
twice (`Genotype` and `Genotype2`, both retained in the analysed object).

***Porites* lost every fragment that was sitting at Tiahura after the 2-week sampling** — both the
Tiahura → Tiahura back-transplants and the Opunohu → Tiahura out-transplants. Its 2-week timepoint is,
by contrast, the most complete cell block in the whole experiment: all four paths, 21 fragments. From
4 months onward only fragments located at **Opunohu** remain: four to five Opunohu → Opunohu
back-transplants and four to five Tiahura → Opunohu out-transplants.

The two species therefore have **complementary gaps**: *Pocillopora* is missing the early timepoint,
*Porites* is missing the later Tiahura environment. Neither gap can be filled by the other species.

### 0.5 Which contrasts the design supports

| Contrast | *Pocillopora* | *Porites* |
|---|---|---|
| Origin bay (colony provenance) | Estimable from 4 months onward | Estimable, but from 4 months onward perfectly confounded with transplant status |
| Destination bay (environment experienced) | Estimable at 4 months, Pre-HS and Post-HS; fully crossed with origin | **Only at 2 weeks.** Constant (Opunohu) from 4 months onward |
| Origin × destination (the transplant effect proper) | Estimable from 4 months onward | Only at 2 weeks |
| Short-term response (2 weeks) | **Not estimable — no libraries** | Estimable, all four paths |
| Longer-term response (4 months) | Estimable | Estimable only for fragments at Opunohu |
| Acute heat shock (Pre-HS vs Post-HS) | Estimable within colony, all four paths | Estimable within the two surviving paths |
| Cross-species comparison | Day 0, and Opunohu-destination paths from 4 months onward | — |

Four consequences follow, and they bear on every result reported so far.

1. **A `Site` term is an origin term.** Because origin and destination are fully crossed in
   *Pocillopora* from 4 months onward, a model containing `Site` alone tests colony provenance while
   pooling fragments from both bays — it does not test environmental exposure. The exposure contrast
   requires a destination term, which no model in step 05 contains.
2. **For *Porites* after 2 weeks, origin, transplant status and destination collapse into one
   another.** The two surviving paths are "stayed home at Opunohu" and "moved to Opunohu", so a
   significant `Site` effect could equally be a provenance difference or the act of being transplanted,
   and a null result cannot be read as "no site effect" in the environmental sense.
3. **Dropping the 2-week samples costs *Porites* its only destination contrast and costs
   *Pocillopora* nothing.** Step 05 discards them to impose a symmetric three-node timeline across
   species. Because *Pocillopora* has no 2-week libraries anyway, that symmetry is achieved trivially,
   and the entire cost of the decision falls on *Porites* — removing the one timepoint at which it has
   all four paths.
4. **The original objectives map onto different species.** The short-term question can only be
   answered in *Porites*; the longer-term and heat-shock questions are answered fully only in
   *Pocillopora*, and in *Porites* only for the Opunohu environment.


---

## 1. Step 01 — sequence variants and QC

*Source: `Analysis/01-AAC-AssignSeqVariants/01-AAC-AssignSeqVariants.md` (knitted 13 May 2026) and
`Data/01-AAC-AssignSeqVariants/01-AAC-AssignSeqVariants_track_read_counts.RData`.*

### 1.1 Methods

Primers were removed with `cutadapt` in an orientation-aware two-pass scheme, targeting the 515F
(`GTGYCAGCMGCCGCGGTAA`) / 926R (`CCGYCAATTYMTTTRAGTTT`) pair, i.e. the V4–V5 region. Host reads were
then depleted with `bowtie2` against host reference assemblies (for *P. rus*,
`GCA_964035705.1_jaPorRusx1.1`), retaining only pairs in which both mates failed to map
(`samtools view -f 12 -F 256`).

ASV inference used DADA2 on the **forward reads only** — the script records that reverse-read quality
was too poor to merge without heavy loss. Reads were filtered with `filterAndTrim(truncQ = 3,
maxEE = 2)` and no `truncLen`; error rates were learned per run, sequences dereplicated, and ASVs
inferred, followed by a length filter (`200:250` bp) and `removeBimeraDenovo`. Taxonomy was assigned
against SILVA nr99 v138.2 with the RDP naive Bayesian classifier and refined with `addSpecies`. Outputs
are the ASV count tables (numbered and sequence-named), the taxonomy table, an ASV fasta and the raw
phyloseq object.

### 1.2 Results

Of 187 collected libraries, **157 entered DADA2**; 30 were dropped beforehand at the empty-well /
low-read-count filter. Read retention through the pipeline was consistent:

| Stage | Reads | Per-sample median |
|---|---|---|
| Input | 2,362,715 | 13,273 |
| Filtered | 2,289,279 | 12,846 |
| Denoised (forward) | 1,858,525 | 10,214 |
| Non-chimeric | 1,795,938 | 10,214 |

Overall retention was **76.0%**, with a per-sample median of 77.3% and a narrow range (60.3–86.4%) —
no sample behaved anomalously. Chimera removal cost 3.4% of reads and took the ASV count from 15,650 to
**15,067**. Depth is modest: median 10,214 reads per sample, minimum 1,133 (`58-OBC03F40POC`), with five
samples below 2,000 reads and 24 below 5,000. Inferred ASVs are 206–216 bp.

### 1.3 What this step does to the design

1. **The 30 dropped libraries are not a random sample of the design.** They comprise all 18
   *Pocillopora* 2-week libraries, three *Pocillopora* 4-month fragments, one *Porites* Pre-HS fragment,
   six sediment samples, one seawater sample and one extraction blank. The loss of an entire
   species × timepoint cell is a structural change to the experiment, not attrition — and it is the
   reason *Pocillopora* has no short-term data (section 0.4).
2. **Only one extraction blank survived into DADA2**, which weakens the contamination assessment
   available to step 02.
3. **The amplicon is forward-read only.** ASVs are ~215 bp of V4 rather than a merged V4–V5 product,
   so genus- and species-level assignment is correspondingly less certain, and the comment in the
   script about trimming ASVs "to 370 bp" is stale: the applied `200:250` filter retained every ASV
   (15,650 before and after, read ratio exactly 1).
4. **Six of twelve sediment samples were lost**, which matters because sediment is one of the three
   environmental reference media used in the step-05 tracking analysis.


---

## 2. Step 02 — decontamination and pruning

*Source: `Analysis/02-AAC-PreProcessing/02-AAC-PreProcessing.md` (knitted 15 May 2026), verified against
the ASV counts table and the step-03 fasta.*

### 2.1 Methods

The raw phyloseq object (15,067 ASVs × 157 samples) was pruned in three moves: removal of ASVs assigned
to chloroplasts (`Order != "Chloroplast"`), removal of mitochondrial ASVs (`Family != "Mitochondria"`),
and an attempt to remove every ASV detected in a negative control, followed by removal of the control
samples themselves. Both taxonomic filters retain `NA` assignments. Contamination was assessed by
comparing raw abundances of control-derived ASVs between blanks and true samples, with a significance
test plotted via `ggpubr`; no prevalence-based model (`decontam` or equivalent) was fitted.

### 2.2 Results

| Filter | ASVs removed | ASVs remaining |
|---|---|---|
| Chloroplast | 510 | 14,557 |
| Mitochondria | 98 | 14,459 |
| Control-derived | 0 (see below) | 14,459 |

Organelle removal cost 608 ASVs, retaining 95.96%. Four blanks survived step 01 — one field blank
(`60-NC`), one PCR blank and two extraction/soil-kit blanks — and 168 ASVs were detected in at least one
of them, dominated by recognisable reagent and skin contaminants (*Cutibacterium*, *Staphylococcus
aureus*, *Sphingobium*, *Sediminibacterium*). The four control samples were removed as intended,
taking the object from 157 to **153 samples**, which is exactly the sample count carried into step 04.

### 2.3 A latent bug, without consequence

**The control-ASV removal never executed.** The call is
`subset_taxa(., !(ASV %in% ASVs_control))`, but `ASVs_control` is a data frame, not the character vector
the earlier code used (`ASVs_control$ASV`); `%in%` against a data frame matches nothing, so the negation
is `TRUE` for every taxon. This is confirmed directly: the fasta handed to step 03
(`03-AAC-PhylogeneticTree_PhyseqASVs.fasta`) contains **14,459 sequences**, exactly
15,067 − 510 − 98. Not one control ASV was dropped.

The effect on results is nevertheless negligible. Reconstructing the control-ASV set from the counts
table gives the same 168 ASVs, and in the 153 true samples those ASVs carry **0.0% of non-organelle
reads** — median 0.0% per sample, maximum 0.1%, with the single largest contributing 52 reads across
the entire dataset, and *Cutibacterium* and *S. aureus* contributing none. The blanks and the samples
simply do not share abundant taxa, so the intended filter would have removed almost nothing.

Two things follow. The bug should be fixed for correctness and reported honestly, but **no downstream
result needs revisiting on its account** — which is worth stating explicitly, because the alternative
assumption (that reagent contaminants pervade the ordinations) is wrong. Separately, the blanket
"remove any ASV seen in a blank" strategy is more aggressive than a prevalence-based approach would be;
had the blanks shared genuine marine taxa with the samples, executing it would have deleted real
signal. If it is fixed, it is worth switching to `decontam` at the same time.


---

## 3. Step 03 — phylogenetic tree

*Source: `Analysis/03-AAC-PhylogeneticTree/03-AAC-PhylogeneticTree.md` (knitted 15 May 2026), script
re-run 17 Aug 2026.*

### 3.1 Methods

The 14,459 ASV sequences carried out of step 02 were written to fasta, aligned with `MAFFT --auto`, and
a tree inferred with `FastTree2` (`-nt -gtr -fastest`). The resulting unrooted tree was merged into the
phyloseq object after an `ntaxa` identity check, eukaryotic ASVs were dropped with a `pop_taxa` helper,
and the tree was midpoint-rooted (`phytools::midpoint.root`) before being merged back to produce
`physeq_midpointTree` — the object every later analysis loads.

### 3.2 Results

The unrooted tree has **14,459 tips and 14,457 internal nodes**. Removing ASVs assigned to Eukaryota
drops **384 ASVs**, leaving **14,075** in the final tree-bearing object. Those eukaryotic ASVs were not
trivial in abundance: they carried **4.85% of non-organelle reads** in the true samples, so this is a
substantive filter rather than a formality — appropriate for a 16S study, but worth reporting as a
quantity rather than a footnote.

**206 archaeal ASVs are retained.** Every downstream analysis is described as addressing "bacterial
communities", but the object is bacterial *and archaeal*. Either the wording should be relaxed to
"prokaryotic" or the archaea should be dropped explicitly; at present the text and the data disagree.

### 3.3 Provenance note

The alignment, tree and `FastTree.log` are dated 15 May 2026, while the input fasta and
`physeq_midpointTree.RData` were regenerated on 17 August 2026 — the MAFFT and FastTree steps are bash
chunks that were not re-executed on the later run. The `stopifnot(ntaxa(...) == ntaxa(unrooted_tree))`
guard in the script is what makes this safe: it would have failed had the regenerated ASV set differed
from the one that was aligned. The check passed at 14,459, so tree and table correspond.

### 3.4 Consequence for the distance metrics

All four metrics used in step 05 rest on this object, and two of them (unweighted and weighted UniFrac)
rest on this tree. The distinction matters for interpreting step 05: unweighted UniFrac is the one
metric in which the site term is non-significant globally (p = 0.105) while the three
abundance-sensitive metrics detect it. A `--fastest` FastTree topology from a ~215 bp forward-read
alignment is a coarse phylogeny, so presence/absence phylogenetic distances are the most error-prone of
the four, and that null should be read with the tree's resolution in mind rather than as a clean
biological statement about lineage pools.


---

## 4. Step 04 — diversity analyses

*Source: `Analysis/04-AAC-DiversityAnalyses/04-AAC-DiversityAnalyses.Rmd` (**never knitted**) and
`Data/04-AAC-DiversityAnalyses/04-AAC-DiversityAnalyses_iNEXT_data.RData`. Results below were computed
here directly from the saved iNEXT object, because no recorded output exists.*

### 4.1 Methods

Hill numbers of order q = 0 (richness), 1 (Shannon) and 2 (Simpson) were computed with `iNEXT` on
abundance data for all 153 samples, with rarefaction and extrapolation curves plotted via `ggiNEXT` and
re-plotted manually. The `iNEXT()` call itself is commented out in the script and the object is loaded
from the saved `.RData`. Diversity was then plotted as boxplots by `FullTreatment` for each host
(`04-AAC-DiversityAnalyses_plot_Hill_Both.png/.pdf`), and richness was tested with a cumulative link
mixed model, `clmm(factor(qD) ~ FullTreatment + (1 | Individual))`.

### 4.2 Results

Across the 153 samples, library sizes are a median of 8,639 reads (minimum 1,068) and observed richness
a median of 135 ASVs (range 21–663). Median observed Shannon diversity, in effective numbers of ASVs,
laid out by where each fragment came from and where it was:

| Species | Origin → destination | Day 0 | 2 wk | 4 mo | Pre-HS | Post-HS |
|---|---|---|---|---|---|---|
| *Pocillopora* | Opunohu → Opunohu | 31.3 | — | 31.6 | 60.7 | 35.6 |
| | Opunohu → Tiahura | — | — | 13.5 | 48.0 | 71.4 |
| | Tiahura → Opunohu | — | — | 59.3 | 103.4 | 80.3 |
| | Tiahura → Tiahura | 33.7 | — | 38.1 | 42.1 | 48.3 |
| *Porites* | Opunohu → Opunohu | 49.0 | 82.5 | 117.7 | 158.3 | 87.0 |
| | Opunohu → Tiahura | — | 52.8 | died | died | died |
| | Tiahura → Opunohu | — | 63.0 | 156.3 | 114.3 | 102.3 |
| | Tiahura → Tiahura | 42.9 | 44.3 | died | died | died |

![Shannon diversity of individual fragments, plotted against the bay each fragment was actually in when sampled. Points are fragments, bars are medians, y-axis logarithmic. Pocillopora has no 2-week libraries; Porites has no fragments at Tiahura after 2 weeks.](AAC-Walkthrough_fig4_AlphaByDestination.png)

Diversity is higher in *Porites* than in *Pocillopora* throughout, and in both hosts the later timepoints
sit above the day-0 baseline. **That apparent increase is not tested here** and sits awkwardly beside
step 05's finding of no timeline effect on Shannon diversity (section 5); the two differ in that step 05
worked on scaled libraries, merged 4-Months with Pre-HS, and used three timeline nodes rather than five.
Reconciling them is on the re-analysis list (section 6).

### 4.3 Why this step has never been knitted

The final third of the script is **verbatim residue from an unrelated tutorial**. It selects
`water_tempC:DO_mgL`, plots diversity against `salinity_psu`, colours by `station` using
`station_colors`, and interprets a freshwater-to-marine salinity gradient ("Freshwater = 0–0.5,
Brackish = 0.5–30, Marine = >30"). None of those variables exist in this project's metadata — there is
no salinity column in `Master_SampleKey.xlsx` — and there is no salinity gradient in a two-bay coral
transplant. These chunks will error, which is the most likely reason the script has no `.md` or `.html`
output. Deleting or disabling them should be enough to make it knit.

### 4.4 Two further problems with the step as written

1. **The iNEXT estimation framework is degenerate on DADA2 output.** DADA2 does not emit singleton
   ASVs, and indeed **f1 = 0 in all 153 samples**. Good's coverage is therefore exactly 1.0 for every
   sample — the reported "complete coverage" is an artefact of the denoiser, not evidence of saturated
   sampling — and the Chao-type asymptotic richness estimator reduces to the observed count: `Observed`
   equals `Estimator` in all 153 samples. Extrapolated richness and the asymptotic estimates carry no
   information beyond observed richness, so only the observed Hill values should be reported, and the
   extrapolation half of the rarefaction figure should be dropped or labelled as non-identifiable.
2. **The statistical test covers one cell of eight.** The CLMM is fitted to `df_richness_POC_Opu` —
   *Pocillopora* of Opunohu origin only — so it speaks to neither host generally nor to any
   destination contrast. It also treats richness as an ordered factor with as many levels as distinct
   values, which is not a natural model for a count; a GLMM with a Poisson or negative-binomial family,
   or a linear model on Hill numbers, would be the conventional choice.


---

## 5. Step 05 — community analyses, re-read against the real design

*Source: `Analysis/05-AAC-CommunityAnalyses/05-AAC-CommunityAnalyses.md` (knitted 26 Oct 2026). Every
statistic below is quoted from that knit; the verdicts apply the design semantics of section 0.*

### 5.1 Methods

The curated object drops the 2-week samples, keeps only `Genotype2` as the *Pocillopora* baseline, and
merges `4-Months` with `Pre-HS` into a single `4-Month-Pre-HS` node. Samples below a read-depth
threshold and those whose rarefaction curves had not flattened (terminal slope > 0.05 new ASVs per read)
were removed, retaining **111 of 122 samples (91.0%)**. Libraries were scaled to a within-group minimum
by proportional scaling — 3,072 reads for *Pocillopora*, 2,066 for *Porites* — and community structure
tested by PERMANOVA on Sørensen, Bray–Curtis, unweighted and weighted UniFrac, with permutations
blocked by mother colony, plus `betadisper`, `pairwise.adonis2`, environmental tracking, alpha-diversity
mixed models, heat-shock models, a three-method differential-abundance suite, Baselga turnover
partitioning, Hill profiles and baseline-drift plots.

### 5.2 Verdicts

| Analysis | Verdict |
|---|---|
| Global PERMANOVA, `Subtype` main effect | Survives, with a dispersion caveat |
| Global PERMANOVA, `Site` main effect and `Subtype × Site` | **Withdraw** — `Site` is not one factor |
| `betadisper` | Survives |
| Pairwise adonis across compartments | Survives, with a dispersion caveat |
| Host-specific PERMANOVA, `Treatment` | Survives |
| Host-specific PERMANOVA, `Site` | Survives **relabelled as origin**; exposure untested |
| Host-specific post-hocs | Survives relabelled |
| Environmental tracking | **Withdraw and recompute** — sign inverted |
| Alpha diversity LMEs and Tukey contrasts | Survives as a null for origin; destination untested |
| Heat-shock ordinations and PERMANOVAs | **Withdraw** — no Pre-HS samples in the subset |
| Differential abundance, transplant tables | Survives relabelled as an origin contrast; misses the real effect |
| Differential abundance, heat-shock tables | **Withdraw** — Post-HS only, origin contrast |
| Baselga turnover partitioning | Survives (descriptive) |
| Hill-number LMEs | Survives as a null for origin |
| Distance-to-baseline drift | Survives relabelled as origin |

### 5.3 What survives

**Compartment dominates community structure.** `Subtype` explains 10.8% (Sørensen) to 23.6% (weighted
UniFrac) of variation, p = 0.001 in all four metrics, and every pairwise compartment contrast is
significant — with *Porites* more distinct from the environmental media than *Pocillopora* is
(weighted UniFrac R² = 0.287 for *Porites* vs 0.2 µm seawater, against 0.129 for *Pocillopora*). The
caveat is that `betadisper` on `Subtype` is significant (F = 10.74, p = 0.001), so some of this is
unequal within-group spread rather than centroid separation. Dispersion across `Site` is homogeneous
(F = 0.197, p = 0.653).

**Time in the field matters more than provenance.** In the host-specific models the `Treatment` term is
significant in all four metrics for both hosts — *Pocillopora* 6.6–8.1% (p ≤ 0.003), *Porites*
10.2–19.4% (p = 0.001) — while the origin term explains 1.9–3.8% in *Pocillopora* (significant in three
of four metrics) and nothing detectable in *Porites* (p = 0.066–0.746). Post-hocs place the shift in the
interval between the day-0 baseline and later sampling rather than across the heat challenge.

**Turnover, not nestedness.** Baselga partitioning gives total Sørensen dissimilarity of 0.943
(*Pocillopora*) and 0.952 (*Porites*), of which 98.0% and 98.5% is pure species turnover. Communities
are replaced rather than depleted — consistent with the flat alpha-diversity results.

### 5.4 What must be withdrawn

**The `Site` term in the global model is not a single factor.** For coral samples `Site` is the colony's
bay of origin; for seawater and sediment it is where the sample was physically collected. The main
effect therefore averages a provenance contrast over corals with a genuine spatial contrast over
environmental media, and the `Subtype × Site` interaction — presented as the mathematical justification
for splitting the dataset by host — is an interaction between compartment and a factor whose meaning
changes across compartments. The conclusion that hosts should be analysed separately is reasonable on
other grounds; this particular test does not support it.

**The environmental-tracking result has its sign inverted.** Line 1092 reads
`t_site <- coral_samples$Site[i] # Current Destination Site` and derives the origin as the opposite bay,
but `Site` is the origin. The fragment's "local" references are therefore drawn from the bay it came
from and its "source" references from the bay it was in, so `Individual_Shift = To_Source − To_Local` is
the negative of the intended quantity and the plot's "Destination Site" axis is the origin site.
Corrected, and relabelled by where each fragment actually was:

| Fragments | Free-living (0.2 µm) | Particle-bound (3 µm) | Sediment |
|---|---|---|---|
| *Pocillopora*, Opunohu → **Tiahura** (n = 13) | +0.014 (FDR 0.004) | +0.008 (n.s.) | +0.030 (FDR 0.0007) |
| *Pocillopora*, Tiahura → **Opunohu** (n = 15) | −0.018 (FDR 0.0003) | −0.027 (FDR 0.008) | −0.028 (FDR 0.0003) |
| *Porites*, Tiahura → **Opunohu** (n = 12) | −0.018 (FDR 0.0009) | −0.041 (FDR 0.008) | −0.023 (FDR 0.0009) |

Positive means the fragment resembles the environment it now sits in. Read this way, **every group is
closer to the Tiahura environmental community regardless of where it was** — fragments at Tiahura look
locally acquired, fragments at Opunohu retain a Tiahura signature. That is equally consistent with coral
microbiomes simply sitting nearer Tiahura's environmental pool than Opunohu's, and the analysis cannot
separate the two because it uses out-transplants only. The discriminating group — Opunohu-origin
fragments held at Opunohu — is in the data but excluded. A further limitation: six of the twelve
sediment samples were lost at step 01 (section 1.3), so the sediment reference is half-sampled.

**The heat-shock analyses contain no heat-shock contrast.** The curation recodes `Pre-HS` into
`4-Month-Pre-HS` at line 147, and the heat-shock subsets then ask for
`Treatment %in% c("Pre-HS", "Post-HS")` — which, after the recode, matches **only Post-HS**. The sample
counts confirm it: `ps_poc_hs` has 20 samples, exactly the *Pocillopora* Post-HS total, and `ps_por_hs`
has 8. This is why `Treatment` and the interaction drop out of both models, leaving `Site` alone
(*Pocillopora* R² = 9.4%, p = 0.077; *Porites* p = 0.5 on seven possible permutations). The two
`HeatShock_DiffAbunTaxa` tables run on the same subsets — 20 and 8 samples — so the 12 and 173 ASVs they
report are origin contrasts among Post-HS fragments, not responses to heat. **No statement about the
acute heat shock is supported anywhere in step 05.** The contrast itself is estimable in both hosts
(section 0.5); it simply was not run.

### 5.5 What was missed

**Differential abundance rests on one method and the wrong factor.** LEfSe and MaAsLin2 returned zero
significant features in all four contrasts, so the "three-method consensus" reduces to DESeq2 alone;
fold changes are unshrunken, which is why magnitudes reach log₂FC ≈ 23 and why they should be read as
direction only. More importantly the contrast is coded on `Site` — origin — and each origin group
contains fragments from both bays, which splits any destination effect down the middle.

![Endozoicomonas relative abundance per colony, by the bay each fragment was in when sampled. Points are colonies (fragments averaged within colony); bars are medians. At day 0 the genus separates by origin; by 4 months it is present in every fragment at Tiahura and absent from every fragment at Opunohu, irrespective of origin.](AAC-Walkthrough_fig5_Endozoicomonas.png)

Taking all 20 *Endozoicomonas* ASVs in the taxonomy rather than the three DESeq2 flagged, and collapsing
to one value per colony, the genus is governed by **destination**. At day 0 it separates by origin —
Tiahura colonies 10/10 positive at a mean of 18.2%, Opunohu colonies 3/10 at a mean of 0.35%. By
4 months every Opunohu-origin colony moved to Tiahura has **gained** it (6.0–23.7%, from near zero), and
every Tiahura-origin colony moved to Opunohu has **lost it entirely** (zero reads). Within each origin
separately: Mann–Whitney p = 0.0075 and Fisher on presence p = 0.0079 (5 vs 5 colonies — 0.0075 is the
smallest value attainable at that sample size). Across all timepoints, 35 of 35 fragments sampled at
Tiahura carry the genus against 2 of 33 at Opunohu, both at trace level. This is a colonisation result,
not merely a loss of symbiont, and it is the strongest effect in the dataset — unestimated by any model
in the step because no model contains a destination term.


---

## 6. Re-analysis specification for the unbalanced design

What follows is the minimum set of changes that makes the analysis answer the questions the experiment
was designed to ask. Code is written against the existing objects and follows the project's explicit
`dplyr::` convention.

### 6.1 Derive the destination, and stop using `Site` as if it were one

Every model in step 05 uses `Site`, which is origin. The exposure factor has to be constructed:

``` r
# Derive where each fragment actually was, and keep origin as a separate factor
meta <- data.frame(phyloseq::sample_data(physeq_midpointTree)) %>%
  dplyr::mutate(
    Origin      = as.character(Site),
    Destination = dplyr::case_when(
      Type != "Coral"                  ~ as.character(Site),   # water/sediment: collected there
      Transplant == "Out Transplant"   ~ dplyr::if_else(Site == "Tiahura", "Opunohu", "Tiahura"),
      TRUE                             ~ as.character(Site)    # back transplants and genotypes stay
    ),
    Moved       = Transplant == "Out Transplant",
    Colony_ID   = paste0(Host, "_", Origin, "_", Colony)       # unique across bays
  )
phyloseq::sample_data(physeq_midpointTree) <- phyloseq::sample_data(meta)
```

`Colony_ID` matters: colony labels `C01`–`C05` repeat across bays, so blocking on `Colony` alone pools
unrelated colonies. The step-05 script builds a `Unique_Colony` for this purpose — keep that, and use it
consistently as the permutation block and the random effect.

### 6.2 Keep the 2-week samples

Dropping them is what costs *Porites* its only destination contrast (section 0.5). Because
*Pocillopora* has no 2-week libraries at all, nothing is gained in symmetry by discarding them. Analyse
the two species separately with their own timelines rather than forcing a shared one.

### 6.3 Un-merge the timeline, or declare the pseudo-replication

`4-Months` (day 120, field) and `Pre-HS` (day 121, after a short common-garden period) are currently
pooled into one node, which puts two fragments from each colony in the same level and mixes field with
aquarium-held material. Either keep them as separate levels, or keep the merge and state that the node
contains two fragments per colony — the colony random effect absorbs this, but it should be explicit.

### 6.4 Models the design supports

``` r
# --- Pocillopora: origin and destination are fully crossed from 4 months onward ---
perm_poc <- permute::how(nperm = 999)
permute::setBlocks(perm_poc) <- poc_meta$Unique_Colony

vegan::adonis2(
  dist_poc ~ Origin * Destination + Treatment,
  data = poc_meta, permutations = perm_poc, by = "terms"
)

# --- Porites, stratum 1: 2 weeks only, where all four paths exist ---
vegan::adonis2(
  dist_por_2wk ~ Origin * Destination,
  data = por_2wk_meta, permutations = perm_por_2wk, by = "terms"
)

# --- Porites, stratum 2: 4 months onward, all survivors are at Opunohu ---
# Destination is constant here, so the only available term confounds origin with having been moved.
vegan::adonis2(
  dist_por_late ~ Moved + Treatment,
  data = por_late_meta, permutations = perm_por_late, by = "terms"
)
```

The *Porites* late-stratum term must be reported as "origin / transplant status among fragments held at
Opunohu", not as a site effect. A null there is not evidence that environment has no effect.

### 6.5 Run the heat-shock test that was never run

``` r
# Subset BEFORE the Pre-HS recode, or the subset matches Post-HS only
ps_hs <- physeq_midpointTree %>%
  phyloseq::subset_samples(Type == "Coral" & Treatment %in% c("Pre-HS", "Post-HS")) %>%
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0, .)

phyloseq::sample_data(ps_hs)$Treatment <-
  factor(phyloseq::sample_data(ps_hs)$Treatment, levels = c("Pre-HS", "Post-HS"))

# Pre-HS and Post-HS are different fragments from the same colony: block on colony
perm_hs <- permute::how(nperm = 999)
permute::setBlocks(perm_hs) <- hs_meta$Unique_Colony
vegan::adonis2(dist_hs ~ Destination * Treatment, data = hs_meta,
               permutations = perm_hs, by = "terms")
```

Check `table(hs_meta$Treatment)` before fitting — it should show both levels, roughly 20 Pre-HS and
20 Post-HS for *Pocillopora*. If `Treatment` still drops out of the model, the subset is wrong, not the
design.

### 6.6 Differential abundance on destination, with shrinkage

``` r
ds <- phyloseq::phyloseq_to_deseq2(ps_filtered, ~ Origin + Destination)
ds <- DESeq2::DESeq(ds, test = "Wald", fitType = "parametric", sfType = "poscounts")

# Shrink before quoting any effect size
res <- DESeq2::lfcShrink(ds, coef = "Destination_Tiahura_vs_Opunohu", type = "apeglm")
```

Including `Origin` alongside `Destination` separates a provenance difference present at day 0 from a
change acquired after transplantation — the distinction that matters for *Endozoicomonas* (section 5.5).
Report direction and prevalence rather than unshrunken fold change. Given LEfSe and MaAsLin2 returned
nothing in any contrast, either drop the consensus framing or report the zero counts explicitly.

### 6.7 Recompute environmental tracking, with the control it needs

Swap `local` and `source` (section 5.4), and add the back-transplants. Out-transplants alone cannot
distinguish "tracks the local environment" from "is generally closer to the Tiahura pool"; fragments
that stayed home supply the baseline that separates them. Note that only six of twelve sediment samples
survived step 01, so the sediment arm is half-sampled and should carry that caveat.

### 6.8 Smaller fixes, in priority order

1. **Step 05** — fix the heat-shock subset (section 6.5). Highest impact: it restores an entire
   experimental objective.
2. **Step 05** — recompute tracking with corrected polarity (section 6.7).
3. **Step 04** — delete the salinity / `station` / `water_tempC:DO_mgL` chunks so the script knits, and
   report observed Hill numbers only, not the extrapolated or asymptotic values (section 4.4).
4. **Step 02** — repair the control-ASV filter (`ASVs_control$ASV`, or switch to `decontam`). Correctness
   only: the retained ASVs carry 0.0% of reads, so no result changes.
5. **Wording** — the object contains 206 archaeal ASVs, so "bacterial communities" should read
   "prokaryotic", or the archaea should be dropped explicitly.
6. **Step 01** — note the forward-read-only amplicon (~215 bp of V4) in the methods, and remove the
   stale "trim to 370 bp" comment, which describes a filter that retained every ASV.


---

## 7. The heat-shock test, actually run

*Computed here on 26 October 2026, independently of the step-05 pipeline. The scripts have been patched
(section 7.3) but not re-knitted; the UniFrac versions of these models will come from that knit.*

### 7.1 What was run

The step-01 ASV table was filtered to drop chloroplast, mitochondrial and eukaryotic ASVs by taxonomy,
then each host's heat-shock window (`Pre-HS` and `Post-HS` fragments, selected on the original timepoint
rather than the recoded `Treatment`) was scaled proportionally to its own minimum depth — 1,334 reads for
*Pocillopora*, 1,819 for *Porites*. PERMANOVA was run with `vegan::adonis2(..., by = "terms")` on
Bray–Curtis and Sørensen distances, with permutations blocked on `Colony_ID` (host × origin bay ×
colony). Note this differs from step 05 in two respects: no low-depth or rarefaction-slope sample
removal was applied, and only the two non-phylogenetic metrics could be computed without the phyloseq
tree pipeline.

The resulting design is balanced for *Pocillopora*: **40 fragments, exactly 10 in each of the four
destination × timepoint cells, across 10 colony blocks**. For *Porites*, all 17 surviving fragments were
at Opunohu (8 Pre-HS, 9 Post-HS, 9 blocks), so destination is constant and only the heat term is
estimable; the block structure admits just 383 distinct permutations, which were enumerated completely.

### 7.2 Results

| Host | Metric | Term | R² | F | p |
|---|---|---|---|---|---|
| *Pocillopora* | Bray–Curtis | **Destination** | **0.073** | 2.99 | **0.001** |
| | | Timepoint (heat shock) | 0.026 | 1.05 | 0.290 |
| | | Destination × Timepoint | 0.021 | 0.86 | 0.712 |
| *Pocillopora* | Sørensen | **Destination** | **0.046** | 1.85 | **0.001** |
| | | Timepoint (heat shock) | 0.027 | 1.07 | 0.221 |
| | | Destination × Timepoint | 0.023 | 0.92 | 0.728 |
| *Porites* | Bray–Curtis | Timepoint (heat shock) | 0.069 | 1.10 | 0.125 |
| *Porites* | Sørensen | Timepoint (heat shock) | 0.066 | 1.06 | 0.250 |

**The acute heat shock did not detectably shift community composition in either host.** In
*Pocillopora* the heat term explains 2.6–2.7% of variation (p = 0.22–0.29) and does not interact with
destination (p ≥ 0.71); in *Porites* it explains 6.6–6.9% but is not significant (p = 0.125–0.25) on a
design that permits only 383 permutations. This is now a tested null rather than an absence of analysis —
and with ten fragments per cell in a balanced *Pocillopora* design, it is a reasonably informative one.

**Destination is significant, and roughly twice the size of the origin effect it was conflated with.**
The bay a *Pocillopora* fragment was actually in explains 4.6% (Sørensen) to 7.3% (Bray–Curtis) of
variation at p = 0.001, against the 1.9–3.8% that step 05 attributed to `Site` — which is origin. The
abundance-weighted metric gives the larger effect, consistent with the compositional pattern in
section 5.5. Dispersion is homogeneous for both terms (`betadisper` on destination F = 0.26, p = 0.609;
on timepoint F = 1.33, p = 0.267), so this is a centroid shift rather than a spread artefact.

Read together with section 5.5, the picture is coherent: within the heat-shock window, where every
fragment has been in its destination bay for four months, community composition is organised by
**where the fragment is**, and a 1–2 day acute heat challenge does not measurably perturb it.

### 7.3 Patches applied to the repository

`Analysis/05-AAC-CommunityAnalyses/05-AAC-CommunityAnalyses.Rmd`:

* A `Timepoint` column now preserves the original sampling timepoint before the `case_when` recode folds
  `Pre-HS` into `4-Month-Pre-HS`; both heat-shock ordination subsets and both heat-shock
  differential-abundance subsets select on it.
* `Colony_ID` (host × origin × colony) added alongside the existing `Unique_Colony`, which omits the bay
  and therefore merges colony `C01` of Tiahura with colony `C01` of Opunohu into one permutation block.
  The heat-shock models now block on `Colony_ID`; **the earlier analyses still use `Unique_Colony` and
  should be switched over**, which will change previously reported numbers.
* `Origin` and `Destination` columns added; the heat-shock models now test
  `Destination * Timepoint`, and section 7.2's breakdown tests `Origin * Destination * Timepoint`.
* `stopifnot` guards and a cross-tabulation print were added after the subsets so an empty factor level
  fails loudly instead of silently dropping a term.

`Analysis/04-AAC-DiversityAnalyses/04-AAC-DiversityAnalyses.Rmd`: the four tutorial-residue chunks
(`environmental-pairs`, `div-vs-salinity`, `richness-vs-salinity`, `plot-richness-sal`) are set to
`eval=FALSE, include=FALSE` with a note explaining why. Nothing was deleted, pending review.


### 7.4 Figure outputs added to the scripts

Neither script wrote its figures anywhere except into the knit's `_files/figure-html/` directory under
chunk-generated names. Both now save to the convention folders with the `pref_out` prefix:

| Script | File written to `Figures/NN-AAC-Name/` |
|---|---|
| 05 | `05-AAC-CommunityAnalyses_plot_EnvTracking.png/.pdf` |
| 05 | `05-AAC-CommunityAnalyses_plot_AlphaDiversity_Shannon.png/.pdf` |
| 05 | `05-AAC-CommunityAnalyses_plot_HeatShockOrdination.png/.pdf` |
| 05 | `05-AAC-CommunityAnalyses_plot_DiffAbund_Phylum.png/.pdf` |
| 05 | `05-AAC-CommunityAnalyses_plot_DiffAbund_Genus.png/.pdf` |
| 05 | `05-AAC-CommunityAnalyses_plot_BaselineDrift.png/.pdf` |
| 04 | `04-AAC-DiversityAnalyses_plot_RarefactionCurves.png/.pdf` (the Hill-profile figure was already saved) |

Each plot that was previously built as a bare expression is now assigned to a named object
(`plot_env_tracking`, `plot_alpha_shannon`, `plot_heatshock_ordination`, `plot_diffabund_phylum`,
`plot_diffabund_genus`, `plot_baseline_drift`, `plot_rarefaction_iNEXT`), printed, then saved in both
formats. `dir_out_fig` and `pref_out` were already defined in each script's setup chunk and are reused
unchanged.

### 7.5 The tracking analysis was also repaired

Because a full re-run would otherwise regenerate the inverted figure, the polarity fix from section 5.4
is now applied in the script: `t_origin` takes the `Site` value, `t_dest` is its opposite,
`local_env_ids` are drawn from `t_dest` and `source_env_ids` from `t_origin`, the output frame carries
both `Origin` and `Destination`, and the summaries, the join and the plot's x-axis all key on
`Destination`. The axis label "Destination Site" is now accurate rather than misleading.

### 7.6 Why these models were not re-run here

The sandbox cannot load phyloseq. Compiled R packages installed into the workspace library are refused
at `dyn.load` because the path is writable, and the conda environment's trusted library is read-only, so
there is nowhere to put Biostrings and the rest of the Bioconductor chain; bioconda has no Windows
builds either. The heat-shock PERMANOVAs in section 7.2 were therefore computed with `vegan` directly
from the ASV table, which is why only Bray–Curtis and Sørensen appear. Knitting the patched script in
RStudio will reproduce them on all four metrics and generate every figure listed above.


### 7.7 Permutation blocking corrected throughout

`Unique_Colony = paste0(Host, "_", Colony)` includes the host but not the bay, and colony labels
`C01`–`C05` restart in each bay. Among the coral samples in the analysed object it therefore yields
**10 blocks where there are 20 real colonies, and every single one of the 10 merged two unrelated
colonies** — one from each bay. All active uses now key on `Colony_ID = Host_Site_Colony`:

* the host-specific PERMANOVA permutation blocks and plot strata;
* the random effect in both Shannon mixed models and both Hill-number mixed models;
* the `geom_path` grouping in the ordination and heat-shock trajectory plots, which was drawing a path
  between two different colonies as though it were one fragment lineage;
* the baseline-drift colony lookup.

The old definition is retained but marked deprecated so it is not reused.

**The drift lookup was the most damaging of these.** It locates a fragment's day-0 reference with
`filter(Colony_ID == current_colony & Treatment == "Genotype") %>% pull(SampleID)`. Under the merged
identifier this returned *two* baseline samples — one per bay — so
`wunifrac_mat[current_sample, baseline_id]` evaluated to a length-2 vector and the surrounding
`data.frame()` recycled the scalar columns into **two rows per fragment, one of them measuring drift
against the wrong bay's colony**. A guard now averages across baselines, which also handles the one
Opunohu *Porites* colony that has two day-0 samples.

**How much does the blocking change move the numbers?** Only p-values can move: `R²` and pseudo-`F` are
computed from the data alone, and the permutation structure affects nothing but the reference
distribution. Re-running the corrected *Pocillopora* heat-shock model both ways confirms this — with the
merged blocking the heat term gives p = 0.315, with the corrected blocking p = 0.290, and `R²` (0.0731
destination, 0.0257 timepoint) and `F` are identical to four decimal places. The implication for
previously reported results is reassuring but should not be over-generalised: **the effect sizes in
section 5.3 stand unchanged**, and only their p-values are open to revision, by an amount this one test
suggests is small. The models that could not be re-run here may behave differently, and the knit will
show.


### 7.8 The knit, and a correction the blocking fix made necessary

Script 05 was committed, knitted and recommitted on 26 October. Three things came through as intended.

* **The heat-shock contrast is live.** The `stopifnot` guards passed and the cross-tabulation confirms
  both levels: *Pocillopora* 10 Pre-HS / 10 Post-HS at Opunohu and 9 / 10 at Tiahura (39 samples — the
  step-05 depth filter removed one that my independent run in section 7.2 retained), *Porites* 8 / 8.
* **Weighted UniFrac reproduces the section 7.2 result.** Destination R² = 8.05%, F = 3.19, p = 0.011;
  heat shock R² = 2.47%, F = 0.98, p = 0.420; interaction p = 0.936. For *Porites*, heat shock
  R² = 5.24%, F = 0.77, p = 0.734. The acute heat shock is a null on the phylogenetic metric as well as
  on Bray–Curtis and Sørensen.
* **Tracking is corrected.** The output table now carries a `Destination` column with Opunohu negative
  (−0.018, −0.027, −0.028 for *Pocillopora*; −0.018, −0.041, −0.023 for *Porites*) and Tiahura positive
  (+0.014, +0.008, +0.030), matching the corrected values in section 5.4 exactly.

**But the blocking fix broke the origin term, and that is my error to own.** Because `Colony_ID`
contains the origin bay, origin is constant within every block — 0 of 10 blocks contain more than one
origin, while destination varies within all 10. Permuting samples within colonies cannot generate a null
distribution for a factor that does not vary within a colony, so the host-specific `Site` p-values in
that knit (all four metrics reading p = 0.001, with `R²` unchanged, including unweighted UniFrac
apparently flipping from p = 0.205) are artefacts of the scheme and must not be used.

The design is **split-plot**: origin is a between-colony (whole-plot) factor; destination and timepoint
are within-colony (split-plot) factors. Each requires its own permutation scheme, demonstrated on the
heat-shock window:

| Term | Permutation scheme | R² | F | p |
|---|---|---|---|---|
| Origin | whole colonies (`Plots(type="free")`, `Within(type="none")`) | 0.046 | 1.83 | **0.009** |
| Origin | within colonies — **invalid for this term** | 0.046 | 1.83 | — |
| Destination | within colonies | 0.073 | 2.99 | **0.001** |
| Timepoint | within colonies | 0.026 | 1.05 | 0.290 |

Origin does have a real effect; it simply has to be tested by permuting whole colonies. Note too that
`R²` and `F` are identical across schemes, which is the general point from section 7.7: permutation
structure moves p-values only.

`run_host_permanova_table()` has therefore been rewritten to build both schemes, test
`Destination * Treatment` within colonies and `Origin` across whole colonies, and report a `Scheme`
column alongside each term. Section 7.2's breakdown is split the same way. Two `stopifnot` guards assert
that origin does *not* vary within colony and that the timeline does, so if the design assumption is
ever violated the script fails loudly rather than returning a spurious p-value; a cross-tabulation of
origin against destination is printed for each host so the tested cells are visible in the knit.

Script 04 was not re-knitted, so `04-AAC-DiversityAnalyses_plot_RarefactionCurves.*` is still pending;
the six step-05 figures were written as expected.
