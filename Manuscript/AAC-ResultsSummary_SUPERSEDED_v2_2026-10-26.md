# BIOcean5D — coral microbiome responses to reciprocal transplantation between an impacted and a pristine bay

Results summary v2, 26 October 2026 · *Pocillopora meandrina* (POC) and *Porites rus* (POR) · Opunohu Bay (turbid, impacted) and Tiahura (pristine), Mo'orea · Supersedes the v1 summary of the same date, which misread the site factor. Statistics from the knitted Analysis/05-AAC-CommunityAnalyses.Rmd, with the design audit and taxon-level analysis in Manuscript/AAC-MethodsResults_Walkthrough.md.


## Design, and what it can support

Colonies of both hosts were reciprocally transplanted between the two bays and sampled at day 0 (mother-colony baseline), 2 weeks, 4 months, immediately before an acute heat challenge, and immediately after it. Sample identifiers encode the **bay of origin**, so a fragment's actual location is the opposite bay whenever it was out-transplanted — the distinction on which the interpretation below turns.

Losses were substantial and fell differently on the two species. **All 18 *Pocillopora* 2-week libraries failed** at sequencing, so that host has no short-term timepoint. **Every *Porites* fragment located at Tiahura died after 2 weeks**, so from 4 months onward all surviving *Porites* sat in Opunohu. The gaps are complementary: the short-term question is answerable only in *Porites*, the longer-term and heat-shock questions fully only in *Pocillopora*. For *Pocillopora* from 4 months onward, origin and destination are fully crossed, so environmental exposure is estimable; for *Porites* after 2 weeks it is not, because destination is constant and origin is confounded with having been moved.


## Main results

**1. *Endozoicomonas* presence is determined by the bay a fragment is in, not where it came from.** At day 0 the genus separates by provenance — Tiahura colonies 10/10 positive at a mean of 18.2% of reads, Opunohu colonies 3/10 at 0.35%. By 4 months, every Opunohu-origin colony moved to Tiahura had **gained** it (6.0–23.7%, from near zero) and every Tiahura-origin colony moved to Opunohu had **lost it entirely** (zero reads). Within each origin separately, Mann–Whitney p = 0.0075 and Fisher on presence p = 0.0079 (5 vs 5 colonies; 0.0075 is the smallest value attainable at that size). Across all timepoints, 35 of 35 fragments at Tiahura carry the genus against 2 of 33 at Opunohu. This is *colonisation*, not only loss of a symbiont: corals from the impacted bay acquire the genus when placed in the pristine one.


![***Endozoicomonas* is governed by destination, not origin.** Each point is a colony (fragments averaged within colony); bars are medians. Colour marks the bay the colony came from; position marks the bay its fragments were in when sampled. At day 0 the genus tracks provenance; by 4 months it is present in every fragment at Tiahura and absent from every fragment at Opunohu. All 20 *Endozoicomonas* ASVs in the taxonomy are included.](AAC-Walkthrough_fig5_Endozoicomonas.png)

**2. Compartment dominates overall community structure.** Host or medium explains 10.8% (Sørensen) to 23.6% (weighted UniFrac) of variation across the dataset, p = 0.001 in all four distance metrics, and every pairwise compartment contrast is significant. *Porites* is more distinct from the environmental media than *Pocillopora* (weighted UniFrac R² = 0.287 vs 0.2 µm seawater, against 0.129 for *Pocillopora*). Dispersion is heterogeneous across compartments (betadisper F = 10.74, p = 0.001), so part of this reflects unequal spread rather than centroid separation.

**3. Time in the field outweighs colony provenance.** In host-specific models blocked by mother colony, the timeline term is significant in all four metrics for both hosts — *Pocillopora* 6.6–8.1% (p ≤ 0.003), *Porites* 10.2–19.4% (p = 0.001) — while provenance explains 1.9–3.8% in *Pocillopora* and nothing detectable in *Porites*. Post-hoc tests place the shift between the day-0 baseline and later sampling.

**4. Communities are replaced, not depleted.** Total Sørensen dissimilarity is 0.943 (*Pocillopora*) and 0.952 (*Porites*), of which **98.0% and 98.5% is pure species turnover** rather than nestedness. Consistently, no alpha-diversity measure responds to provenance or timeline in either host (mixed models, all p ≥ 0.27).

**5. Transplanted fragments resemble the Tiahura environmental pool wherever they sit.** Comparing each out-transplanted fragment's phylogenetic distance to local versus origin-bay reference communities, fragments moved to Tiahura resemble their new surroundings (sediment +0.030, FDR = 0.0007; free-living water +0.014, FDR = 0.004) while fragments moved to Opunohu retain a Tiahura signature (−0.018 to −0.041, FDR ≤ 0.008). Because the test uses out-transplants only, it cannot yet distinguish genuine environmental tracking from coral microbiomes simply sitting nearer the Tiahura pool; the back-transplant controls needed to separate these are in the dataset but were not included.


## What the current analysis cannot support

1. **Nothing about the acute heat shock.** The heat-shock subsets request `Treatment %in% c("Pre-HS", "Post-HS")` after a recode has already folded `Pre-HS` into the 4-month node, so they match Post-HS only — 20 *Pocillopora* and 8 *Porites* samples. The term therefore drops out of both models and both "heat shock" differential-abundance tables are provenance contrasts among post-heat fragments. The contrast is estimable in both hosts; it has simply not been run.

2. **No site effect in the environmental sense.** Every `Site` term in the fitted models is bay of origin. In the whole-dataset model the factor is not even internally consistent — origin for corals, collection location for water and sediment — so that main effect and the compartment × site interaction should be set aside.

3. **Effect sizes from differential abundance.** LEfSe and MaAsLin2 returned zero features in all four contrasts, leaving DESeq2 alone, and fold changes are unshrunken (log₂FC reaching 23 on sparse data). Direction and prevalence are reportable; magnitudes are not.


## Immediate next steps

Re-run the heat-shock models with the subset taken before the recode; recompute environmental tracking with corrected polarity and the back-transplant controls included; refit differential abundance on destination with origin as a covariate and `lfcShrink` applied; and retain the 2-week samples, which carry the only *Porites* environmental contrast in the experiment.
