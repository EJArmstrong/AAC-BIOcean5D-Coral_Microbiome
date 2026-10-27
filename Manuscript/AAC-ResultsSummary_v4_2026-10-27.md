# BIOcean5D — bacterial and Symbiodiniaceae responses to reciprocal transplantation between an impacted and a pristine bay

Results summary v4, 27 October 2026 · *Pocillopora meandrina* (POC) and *Porites rus* (POR) · Opunohu Bay (turbid, impacted) and Tiahura (pristine), Mo'orea.

Supersedes v3 of 26 October, which compressed the 16S results to make room for ITS2 and in doing so omitted the compartment effect, the timeline effect and the environmental-tracking analysis. All three are restored here, so this version is the complete record of everything analysed to date. **The *Pocillopora* heat-shock RNAseq is deliberately excluded** — it has not been quantified. Statistics from `Analysis/05-AAC-CommunityAnalyses.Rmd` (16S) and `Analysis/06-AAC-ITS2Symbionts.Rmd` (ITS2); full audit, including the corrections applied to the original pipeline, in `Manuscript/AAC-MethodsResults_Walkthrough.md`.

## Methods

Colonies of both hosts (5 per bay per species) were reciprocally transplanted between Opunohu and Tiahura and sampled as mother-colony genotypes (day 0), at 2 weeks, at 4 months, immediately before an acute heat challenge and immediately after it, alongside free-living (0.2 µm) and particle-bound (3 µm) seawater and sediment references.

Bacterial V4–V5 was amplified with 515F/926R, primers were removed with cutadapt, host reads depleted against the *Porites rus* assembly GCA_964035705.1 with bowtie2, and ASVs inferred with DADA2 from forward reads only (reverse quality was insufficient; `truncQ = 3`, `maxEE = 2`), with taxonomy from SILVA nr99 v138.2. Chimera removal left **1,795,938 reads (76.0% of input) in 15,067 ASVs**; after removal of chloroplast, mitochondrial and eukaryotic sequences, **14,075 prokaryotic ASVs** remained at a median depth of 10,214 reads. Sequences were aligned with MAFFT and a midpoint-rooted tree built with FastTree2. Community analyses used 111 of 122 coral samples (91.0%), scaled to a within-group minimum by proportional scaling — 3,072 reads for *Pocillopora*, 2,066 for *Porites* — and were tested on Sørensen, Bray–Curtis, and unweighted and weighted UniFrac distances.

Symbiodiniaceae ITS2 amplicons from the same fragments (135 samples) were processed externally with SymPortal. Because each sample carries essentially one ITS2 type profile (mean 1.1–1.4 detected), the profile table is close to categorical, so all symbiont statistics use the **post-MED sequence** layer — SymPortal's clade-C UniFrac distances and Bray–Curtis on post-MED sequence abundances — with the type profiles used for description only.

Community structure was tested by PERMANOVA (`vegan::adonis2`, 999 permutations) under a split-plot permutation design: whole colonies permuted freely for the between-colony term (bay of origin), samples permuted within colony for the within-colony terms (destination, timepoint). R² and pseudo-F are computed from the data alone and so are unaffected by the permutation scheme; only p-values depend on it. Dispersion homogeneity was checked with `betadisper` wherever a PERMANOVA term was significant.

## Design, and which contrasts it supports

Sample identifiers encode the **bay of origin**, so a fragment's actual location is the opposite bay whenever it was out-transplanted. This distinction is the one on which the whole interpretation turns, and the original analysis had it inverted.

![**The design as sampled.** Values are fragments in the analysed phyloseq object. The two hosts failed in complementary ways: *Pocillopora* lost every 2-week library at sequencing, and every *Porites* fragment held at Tiahura died after the 2-week timepoint.]({{artifact:446109d1-c804-48aa-870f-4fbdec0ed9d8}})

The consequence is that **origin and destination are fully crossed only in *Pocillopora*, from 4 months onward** — which is where the environmental and heat-shock contrasts are estimable. For *Porites*, destination is testable only at 2 weeks; from 4 months on every surviving fragment sat in Opunohu, so destination is constant and origin is perfectly confounded with having been moved. Neither gap can be filled by the other species.

## Bacterial community

**1. Compartment dominates overall community structure.** Host or medium explains 10.8% (Sørensen) to 23.6% (weighted UniFrac) of variation across the full dataset, p = 0.001 in all four metrics, and every pairwise compartment contrast is significant. *Porites* is more distinct from the environmental media than *Pocillopora* is (weighted UniFrac R² = 0.287 against 0.2 µm seawater, versus 0.129 for *Pocillopora*). One caveat attaches: `betadisper` on compartment is significant (F = 10.74, p = 0.001), so part of this reflects unequal within-group spread rather than centroid separation alone. Dispersion across bay is homogeneous (F = 0.197, p = 0.653).

**2. Destination, not origin, governs coral composition — and *Endozoicomonas* is colonised from the surrounding water.** In the balanced *Pocillopora* heat-shock window (n = 40, 10 colonies), the bay a fragment sat in explains **4.6% (Sørensen) to 7.3% (Bray–Curtis)** of variation at p = 0.001, and **8.05% of weighted UniFrac** (F = 3.19, p = 0.011). Bay of origin, tested correctly by permuting whole colonies, explains 4.6% (F = 1.83, p = 0.009). The abundance-weighted metric gives the larger effect, consistent with the compositional pattern below.

The taxon-level pattern is far stronger than any fitted model. At day 0 *Endozoicomonas* separates by provenance — Tiahura colonies 10/10 positive at a mean of 18.2% of reads, Opunohu colonies 3/10 at 0.35%. By 4 months every Opunohu-origin colony moved to Tiahura had **gained** the genus (6.0–23.7%, from near zero) and every Tiahura-origin colony moved to Opunohu had **lost it entirely** (zero reads); within each origin separately, Mann–Whitney p = 0.0075 and Fisher on presence p = 0.0079 (5 vs 5 colonies, where 0.0075 is the smallest value attainable). Across all timepoints, 35 of 35 fragments at Tiahura carry it against 2 of 33 at Opunohu. This is colonisation, not merely loss: corals from the impacted bay **acquire** the genus when placed in the pristine one.

![***Endozoicomonas* is governed by destination, not origin.** Each point is a colony (fragments averaged within colony); bars are medians. Colour marks the bay of origin, position the bay the fragments were in when sampled. All 20 *Endozoicomonas* ASVs are included.]({{artifact:0fb82a6b-1509-468e-bfd9-a5399acc1e58}})

**3. Time in the field outweighs colony provenance.** In host-specific models blocked by mother colony, the timeline term is significant in all four metrics for both hosts — *Pocillopora* 6.6–8.1% (p ≤ 0.003), *Porites* 10.2–19.4% (p = 0.001) — while origin explains 1.9–3.8% in *Pocillopora* (significant in three of four metrics) and nothing detectable in *Porites* (p = 0.066–0.746). Post-hoc tests place the shift in the interval between the day-0 baseline and later sampling, not across the heat challenge.

**4. Transplanted fragments resemble the Tiahura environmental pool wherever they sit.** Comparing each out-transplanted fragment's phylogenetic distance to local versus origin-bay reference communities, with positive meaning the fragment resembles the environment it now occupies:

| Fragments | Free-living (0.2 µm) | Particle-bound (3 µm) | Sediment |
|---|---|---|---|
| *Pocillopora*, Opunohu → **Tiahura** (n = 13) | +0.014 (FDR 0.004) | +0.008 (n.s.) | +0.030 (FDR 0.0007) |
| *Pocillopora*, Tiahura → **Opunohu** (n = 15) | −0.018 (FDR 0.0003) | −0.027 (FDR 0.008) | −0.028 (FDR 0.0003) |
| *Porites*, Tiahura → **Opunohu** (n = 12) | −0.018 (FDR 0.0009) | −0.041 (FDR 0.008) | −0.023 (FDR 0.0009) |

Read this way, **every group is closer to the Tiahura environmental community regardless of where it was**: fragments at Tiahura look locally acquired, fragments at Opunohu retain a Tiahura signature. That is equally consistent with coral microbiomes simply sitting nearer Tiahura's environmental pool than Opunohu's, and this analysis cannot separate the two, because it uses out-transplants only — the discriminating group, Opunohu-origin fragments held at Opunohu, is in the dataset but excluded. Six of the twelve sediment samples were also lost at sequencing, so the sediment reference is half-sampled. **Treat this as suggestive, not as demonstrated environmental tracking.**

**5. Communities are replaced, not depleted.** Baselga partitioning gives total Sørensen dissimilarity of 0.943 (*Pocillopora*) and 0.952 (*Porites*), of which **98.0% and 98.5% is pure species turnover** rather than nestedness. Consistently, no alpha-diversity measure responds to provenance or timeline in either host (mixed models, all p ≥ 0.27).

**6. The acute heat shock is a tested null for bacteria.** With the heat-shock subsets taken before the timepoint recode that had previously removed the contrast entirely, Pre-HS versus Post-HS is non-significant in *Pocillopora* on every metric (Bray–Curtis R² = 0.026, p = 0.290; Sørensen 0.027, p = 0.221; weighted UniFrac 2.47%, p = 0.420), and the destination × heat-shock interaction is flat (p ≥ 0.712; weighted UniFrac p = 0.936). *Porites* shows nothing either (Bray–Curtis p = 0.125; weighted UniFrac R² = 5.24%, p = 0.734). Dispersion is homogeneous for both terms (destination F = 0.26, p = 0.609; timepoint F = 1.33, p = 0.267), so this is a genuine absence of response rather than a power artefact of unequal spread.

## Symbiodiniaceae community

**7. Symbionts are host-controlled and site-insensitive — the mirror image of the bacteria.** Host identity explains **72.5% of clade-C UniFrac variation** (p = 0.001; 65.9% on Bray–Curtis): *Pocillopora* carries C1/C42-type profiles, *Porites* the C15 lineage, the canonical host-specific split within *Cladocopium*. Both hosts are essentially pure *Cladocopium* (97.3% and 98.9% of post-MED reads) and **Durusdinium is absent** (0.000% and 0.004%), so these corals hold no thermotolerant symbiont to shuffle towards. Unlike the bacteria, neither origin (p = 0.587) nor destination (p = 0.208) shifted the symbiont community of *Pocillopora* over four months in the opposite bay.

**8. Heat shock destabilised the symbiont community instead of shifting it — an Anna Karenina response.** The mean symbiont community did not move across the heat challenge (timepoint p = 0.568 UniFrac, 0.065 Bray–Curtis), but **among-fragment variability rose about 1.75-fold** in *Pocillopora* (UniFrac F = 8.19, p = 0.005; Bray–Curtis F = 8.61, p = 0.003). Variance rising while the mean holds still is the diagnostic signature of the Anna Karenina principle (Zaneveld, McMinds & Vega Thurber 2017), and it survives the obvious alternatives: the median dispersion ratio (2.30) exceeds the mean, so it is not outlier-driven; permutations blocked by colony give p = 0.003, and a paired test across the 20 colony × destination pairs carrying both timepoints gives 16 of 20 increasing (p = 0.0042); and although Post-HS libraries are shallower (median 39,213 versus 48,336 reads, p = 0.017), depth does not predict distance-to-centroid (ρ = −0.106, p = 0.515) and sequence richness is unchanged.

The destination × heat-shock interaction is significant on both layers (R² = 7.5%, p = 0.021; 7.2%, p = 0.010) while neither main effect is, so the responses are idiosyncratic with a site-dependent direction rather than purely stochastic (Opunohu R² = 7.8%, p = 0.006; Tiahura R² = 11.2%, p = 0.043). *Porites* showed no such increase (ratio 1.24, p = 0.482).

![**Symbiodiniaceae response to acute heat shock.** (A) PCoA of clade-C UniFrac distances for the *Pocillopora* heat-shock window with 95% ellipses by timepoint; the Post-HS ellipse is visibly wider. (B) Distance to group centroid rises after heat shock in fragments held in either bay.]({{artifact:37304051-a86b-44b6-ae71-1c5a49f91eda}})

## Synthesis

The two symbiotic partners responded on orthogonal axes. **Bacteria tracked the environment and ignored the heat shock**: destination governs composition, *Endozoicomonas* presence is set by the bay a fragment occupies, and Pre-HS versus Post-HS is a clean null. **Symbionts ignored the environment and came apart under the heat shock**: four months in the opposite bay moved nothing, while the acute challenge raised among-fragment dispersion 1.75-fold without shifting the mean.

The direction of that contrast is worth stating plainly, because it inverts the usual framing. The Anna Karenina principle was formulated for the bacterial microbiome, yet here the bacterial community of the same fragments did **not** destabilise (betadisper on the heat-shock term F = 1.33, p = 0.267) while the obligate, host-selected symbiont did. Whatever the heat challenge did to these corals, it was not registered by the environmentally-acquired community that is usually expected to report it first.

## Limits

1. **The Anna Karenina reading rests on treatment, not measured stress.** No bleaching score, photochemical yield or symbiont density has been brought into the analysis, and the assemblage being destabilised is one to two closely related *Cladocopium* type profiles, not a diverse community.
2. **Environmental tracking is not demonstrated** — see result 4. The back-transplant controls that would discriminate are in the dataset but were not included, and the sediment reference is half-sampled.
3. **Differential abundance rests on DESeq2 alone.** LEfSe and MaAsLin2 returned zero features in all four contrasts, and fold changes are unshrunken, so direction and prevalence are reportable but magnitudes are not. Those tables still contrast origin and need refitting on destination with origin as a covariate and `lfcShrink` applied.
4. **The *Porites* heat-shock window (n = 18) lies entirely at Opunohu**, so no destination contrast exists for that host and its null results are low-powered.
5. **The compartment effect carries a dispersion caveat** (result 1), and SymPortal's clade-C phylogeny is taken as published rather than re-derived.

## Next steps

Refit the differential-abundance models on destination with origin as a covariate and shrunken fold changes; recompute environmental tracking with the back-transplant controls included, which converts result 4 from suggestive to testable; and quantify the *Pocillopora* heat-shock RNAseq (34 single-end libraries, each sequenced twice), which is the one dataset that can establish whether these fragments were physiologically stressed at all, and whether the host transcriptional response aligns with the symbiont destabilisation reported here.
