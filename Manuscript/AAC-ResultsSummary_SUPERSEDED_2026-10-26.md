# BIOcean5D — coral microbiome responses to reciprocal transplantation and acute heat shock

Results summary, 26 October 2026 · *Pocillopora meandrina* (POC) and *Porites rus* (POR) · Opunohu Bay (turbid, impacted) versus Tiahura (pristine), Mo'orea · All values from the knitted run of Analysis/05-AAC-CommunityAnalyses.Rmd (26 Oct 2026).


## Design and dataset

Fragments of both hosts were reciprocally transplanted between the two bays and sampled on a three-point timeline — mother-colony baseline (*Genotype*), four months in situ (*4-Month-Pre-HS*), and immediately after an acute heat challenge (*Post-HS*) — alongside environmental references from free-living seawater (0.2 µm), particle-bound seawater (3 µm) and benthic sediment. After removal of low-depth samples and non-asymptotic rarefaction curves, **111 of 122 samples (91.0%) were retained**; libraries were scaled to the within-group minimum, 3,072 reads for *Pocillopora* and 2,066 for *Porites*. Community structure was tested by PERMANOVA on four distance metrics (Sørensen, Bray–Curtis, unweighted and weighted UniFrac) with permutations blocked by mother colony.


## Main results

**1. Compartment, not site, dominates.** Host/medium (*Subtype*) explains 10.8% (Sørensen) to 23.6% (weighted UniFrac) of variation (all p = 0.001); site explains 0.95–1.55%. Site is detectable in three of four metrics (p = 0.002–0.013) but **not** in unweighted UniFrac (p = 0.105): the pool of lineages present is shared between bays, and site differences are carried by relative abundances. A significant *Subtype × Site* interaction in every metric (4.3–4.8%, p ≤ 0.015) justifies analysing hosts separately.

**2. Site sensitivity is host-specific.** With permutations blocked by colony, *Pocillopora* differs by site in three metrics — Bray–Curtis R² = 3.2% (p = 0.001), weighted UniFrac R² = 3.8% (p = 0.008), Sørensen R² = 1.9% (p = 0.014) — but not unweighted UniFrac (p = 0.205). *Porites* shows **no site effect in any metric** (p = 0.066–0.746).


| Term | Pocillopora R² (p) | Porites R² (p) |
|---|---|---|
| Site | 1.9–3.8% (0.001–0.205) | 2.1–4.0% (0.066–0.746) |
| Treatment (timeline) | 6.6–8.1% (0.001–0.003) | 10.2–19.4% (all 0.001) |
| Site × Treatment | 2.7–4.1% (0.018–0.788) | 6.2–8.0% (0.153–0.508) |

*Ranges span the four distance metrics; weighted UniFrac gives the largest timeline effect in both hosts.*

**3. Time in the field outweighs site.** The timeline term is significant in all four metrics in both hosts and is roughly twice as large in *Porites*. Post-hoc tests place the shift in the four-month interval, not the heat challenge: in *Pocillopora*, *Genotype* vs *Post-HS* is significant at both sites (p = 0.001–0.021) while *4-Month-Pre-HS* vs *Post-HS* is significant in none (p = 0.056–0.245); in Opunohu *Porites* only *Genotype* vs *4-Month-Pre-HS* survives correction (p = 0.041–0.044).

**4. Fragments track local sources in Tiahura but resist in Opunohu.** Comparing each out-transplanted fragment's unweighted UniFrac distance to local versus source-site references (paired Wilcoxon, FDR-corrected), every fragment moved *into* Tiahura converged on local sources (+0.018 to +0.041, FDR ≤ 0.008 in all six host × medium combinations). *Pocillopora* moved into Opunohu did the opposite, staying closer to its Tiahura source signature for free-living water (−0.014, FDR = 0.004) and sediment (−0.030, FDR = 0.0007); the particle-bound fraction was not significant (−0.008, FDR = 0.34). This asymmetry is the most distinctive result in the dataset.


![**Directional environmental tracking of transplanted fragments.** Mean ± SE of the per-fragment difference in unweighted UniFrac distance (distance to source-site reference minus distance to local reference); values above zero indicate convergence on the destination bay's community. Asterisks mark FDR-corrected paired Wilcoxon p < 0.05. n = 13 and 15 fragments for *Pocillopora* in Opunohu and Tiahura, 12 for *Porites* in Tiahura; no estimable values exist for *Porites* out-transplants into Opunohu (n.d.).](AAC-ResultsSummary_EnvTracking.png)

**5. Diversity is flat while membership turns over almost completely.** Neither Shannon diversity nor Hill q1 responds to site, timeline or their interaction in either host (mixed models with colony as random effect, Type III tests: all p ≥ 0.27; Tukey contrasts p ≥ 0.43). Yet Baselga partitioning gives total Sørensen dissimilarity of 0.943 (*Pocillopora*) and 0.952 (*Porites*), of which **98.0% and 98.5% is pure species turnover** rather than nestedness — communities are replaced, not depleted.

**6. *Endozoicomonas* loss marks the impacted bay in *Pocillopora*.** The *Pocillopora* transplant contrast yields seven differentially abundant ASVs, of which **three are *Endozoicomonas*, all depleted in Opunohu** (log₂FC −23.2 to −23.9, padj ≤ 2.6 × 10⁻¹³), against enrichment of *Pseudoalteromonas* and a Xenococcaceae cyanobacterium. Given this genus's status as a core *Pocillopora* symbiont, it is the clearest candidate dysbiosis signal here. Other contrasts: 19 ASVs (*Porites* transplant, 15 enriched in Opunohu), 12 (*Pocillopora* heat-shock window), 173 (*Porites* heat-shock window, 138 depleted).


## Limitations that bear on interpretation

1. **Dispersion confounds result 1.** Betadisper is significant for *Subtype* (F = 10.74, p = 0.001), so part of that result may be unequal within-group spread rather than centroid separation. Site dispersion is homogeneous (F = 0.197, p = 0.653), so results 2 and 4 are unaffected.

2. **The acute heat shock was not testable as specified.** In both heat-shock PERMANOVAs the *Treatment* term and interaction drop out, leaving site alone (*Pocillopora* R² = 9.4%, p = 0.077, n = 20; *Porites* p = 0.5 on only **7 possible permutations**, n = 8). No heat-shock effect should be claimed from these models, and the "heat shock" differential-abundance tables are site contrasts *within* that window, not pre- versus post-heat comparisons.

3. **Differential abundance rests on one method.** LEfSe and MaAsLin2 returned zero significant features in all four contrasts, so the three-method consensus reduces to DESeq2 alone. Fold changes are unshrunken (no `lfcShrink`), which is why magnitudes reach log₂FC ≈ 23; report direction, not magnitude.

4. **Sequencing depth is shallow.** Analysis depths of 2,066–3,072 reads per sample under-sample rare taxa.

5. **The reciprocal design is incomplete for *Porites*.** No tracking values are estimable for *Porites* fragments transplanted into Opunohu, so result 4 is a one-host statement about the impacted bay.


## Next steps

Recode the heat term so it is estimable within colony blocks; apply `lfcShrink` before quoting any effect size; and either recover the Opunohu *Porites* out-transplants or declare the asymmetry explicitly in the manuscript.
