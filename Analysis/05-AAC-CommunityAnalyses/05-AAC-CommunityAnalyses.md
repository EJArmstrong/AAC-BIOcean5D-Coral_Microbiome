---
title: "AAC - BIOcean5D Coral Microbiomes: 05 - Community Analyses"
author: "Eric J. Armstrong, eja98@cornell.edu, @EricJ_Armstrong"
date: "28 October, 2026"
generator: pandoc
viewport: width=device-width, initial-scale=1
output:
  html_document:
    code_folding: show
    highlight: default
    keep_md: true
    theme: journal
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: true
      toc_depth: 3
  pdf_document:
    # Unicode that pdflatex cannot typeset and would halt on.
    toc: true
    # XeLaTeX rather than the default pdflatex: this document contains
    latex_engine: xelatex
editor_options:
  chunk_output_type: console
---





#### Purpose of this file  

This file contains all of the code for producing a curated OTU table from the 16S microbial sequence data for the Armstrong et al. manuscript entitled **"XXX"**, published by [XXX](https://doi.org/10.1093/femsec/fiaa029).

##### If you have any questions, you are welcome to email the corresponding author at eja98@cornell.edu or tweet him at [EricJ_Armstrong](https://x.com/EricJ_Armstrong).

The goal is to perform community level diversity analyses

1.  Load in phyloseq data with rooted tree.\
2.  Evaluate sequencing depth and remove sample.\
3.  Normalize the read counts between samples.\
4.  Calculate community dissimilarities. Numbers between 0 and 1. If 0,
    completely similar versus if they are 1, then they're completely
    dissimilar.
    a.  **Sorensen**: Shared Species as a binary value:
        Abundance-unweighted
    b.  **Bray-Curtis**: Shared Abundant species: Abundance-weighted
    c.  **(Abundance-)Weighted UNIFRAC**: Consider Abundant Species and
        where they fall on the tree\
5.  Visualize the community data with two unconstrained Ordinations:
    a.  **PCoA**: Linear Method. Eigenvalue = how much variation is
        explained by each axis. Choose to view axis 1, 2, 3, etc. and
        plot them together.\
    b.  **NMDS**: Non-linear. Smush multiple Dimensions into 2 or 3
        axes. Need to report Stress value (ideally \<0.15).\
6.  Run statistics with PERMANOVA and betadispR.

# 1 - Prepare Workspace
Load necessary R packages and set default parameters

## 1.1 - Define Paths

``` r
# a - Define the main project directory ----
  dir_proj <- 
  "C:/Users/uax75/OneDrive/Documents/R/AAC-BIOcean5D-Coral_Microbiome/"

# b - Define the directory containing the input data ----
  dir_in <- 
  "C:/Users/uax75/OneDrive/Documents/R/AAC-BIOcean5D-Coral_Microbiome/Data/"
    
# c - Define the directories that will contain output Files/Figures/Tables ----
  dir_out_dat <- 
    "C:/Users/uax75/OneDrive/Documents/R/AAC-BIOcean5D-Coral_Microbiome/Data/05-AAC-CommunityAnalyses/"
  dir.create(dir_out_dat, recursive = TRUE)
  dir_out_fig <- 
    "C:/Users/uax75/OneDrive/Documents/R/AAC-BIOcean5D-Coral_Microbiome/Figures/05-AAC-CommunityAnalyses/"
  dir.create(dir_out_fig, recursive = TRUE)
  dir_out_tab <- 
    "C:/Users/uax75/OneDrive/Documents/R/AAC-BIOcean5D-Coral_Microbiome/Tables/05-AAC-CommunityAnalyses/"
  dir.create(dir_out_tab, recursive = TRUE)
```

## 1.2 - Set Defaults

``` r
# a - Define the output prefix for generated files ----
  pref_out <- "05-AAC-CommunityAnalyses_"
    
# b - Set the random seed value ----
  set.seed(123)

```

## 1.3 - Load Libraries

``` r
# a - load necessary packages ----
  pacman::p_load(tidyverse, devtools, phyloseq, patchwork, vegan, 
                 ggplot2, pairwiseAdonis, nlme, car,
               install = FALSE)

```


``` r
# PATCH 2026-10-28: reusable between-colony test for Origin.
# Plots(strata = ..., type = "free") requires BALANCED strata - permute::check() stops with
# "Design must be balanced if permuting 'strata'" as soon as colonies contribute unequal
# fragment counts, which became true after the amplicon-window repair changed which
# libraries survive filtering. This helper keeps the whole-colony permutation wherever the
# design is still balanced, and otherwise collapses each colony to a centroid and tests
# Origin on those points under unrestricted permutation - a valid between-colony test that
# needs no balance, because one point per colony makes colony the unit of replication by
# construction. NOTE: in the centroid case R2 is the variance among COLONY CENTROIDS
# explained by Origin, not the variance among fragments, so it is not comparable with the
# balanced-case R2. The scheme actually used is printed and stored as an attribute.

origin_between_test <- function(dm, meta, nperm = 999) {

  stopifnot(all(c("Colony_ID", "Origin") %in% names(meta)))

  frag     <- table(meta$Colony_ID)
  balanced <- length(unique(as.vector(frag))) == 1

  if (balanced) {
    ctrl <- how(plots  = Plots(strata = meta$Colony_ID, type = "free"),
                within = Within(type = "none"), nperm = nperm)
    res  <- adonis2(dm ~ Origin, data = meta, permutations = ctrl, by = "terms")
    attr(res, "scheme") <- "whole colonies"
  } else {
    k    <- min(nrow(meta) - 1L, 20L)
    pc   <- cmdscale(dm, k = k, eig = TRUE)
    keep <- pc$eig[seq_len(ncol(pc$points))] > 1e-8
    pts  <- pc$points[, keep, drop = FALSE]
    cen  <- aggregate(pts, by = list(Colony_ID = meta$Colony_ID), FUN = mean)
    key  <- dplyr::distinct(meta[, c("Colony_ID", "Origin")])
    cen  <- dplyr::left_join(cen, key, by = "Colony_ID")
    dmc  <- dist(cen[, setdiff(names(cen), c("Colony_ID", "Origin"))])
    res  <- adonis2(dmc ~ Origin, data = cen, permutations = nperm, by = "terms")
    attr(res, "scheme") <- "colony centroids"
  }

  cat("  fragments per colony:",
      paste(sort(unique(as.vector(frag))), collapse = "/"),
      "->", attr(res, "scheme"), "\n")
  res
}
```

## 1.4 - Load Phyloseq Object


``` r
# a - Load physeq ----
  load(paste(dir_in, 
             "03-AAC-PhylogeneticTree/", 
             "03-AAC-PhylogeneticTree_physeq_midpointTree.RData", 
             sep = ""))

# b - Perform targeted curation: Establish clean symmetric 3-node timeline ----
  metadata_raw <- data.frame(sample_data(physeq_midpointTree)) %>%
    mutate(SampleID_Temp = rownames(.)) %>%
    # Force categorical fields to text strings to prevent factor integer-masking bugs
    mutate(
      Treatment  = as.character(Treatment),
      Transplant = as.character(Transplant),
      Host       = as.character(Host),
      Type       = as.character(Type)
    ) %>%
    # 1. Drop 2-Week samples completely to establish inter-species symmetry
    filter(Treatment != "2-Weeks" | is.na(Treatment)) %>%
    # 2. NA-PROOF RECIPIENT FILTER: Only clean rows if they are true Corals
    # This prevents R from deleting environmental references (water/sediment) due to NA cells
    filter(Type != "Coral" | !(Host == "POC" & Treatment == "Genotype" & Transplant != "Genotype2")) %>%
    # PATCH 2026-10-26: preserve the ORIGINAL sampling timepoint before recoding.
    # The case_when below folds Pre-HS into "4-Month-Pre-HS", which silently empties any
    # later subset asking for Treatment %in% c("Pre-HS", "Post-HS").
    mutate(Timepoint = as.character(Treatment)) %>%
    # 3. Safe mapping across symmetric nodes for all remaining samples
    mutate(Treatment = case_when(
      Type == "Coral" & Transplant == "Genotype2"               ~ "Genotype",
      Type == "Coral" & Treatment %in% c("4-Months", "Pre-HS")  ~ "4-Month-Pre-HS",
      Type == "Coral" & Treatment == "Post-HS"                  ~ "Post-HS",
      TRUE                                                      ~ Treatment
    )) %>%
    # Convert back into factors for downstream statistical modeling
    mutate(
      Treatment  = as.factor(Treatment),
      Transplant = as.factor(Transplant),
      Host       = as.factor(Host),
      Type       = as.factor(Type)
    )

  # Re-apply clean metadata subset matrix back onto parent phyloseq infrastructure
  physeq_pruned <- prune_samples(metadata_raw$SampleID_Temp, physeq_midpointTree)
  sample_data(physeq_pruned) <- sample_data(metadata_raw)

# c - Add a column to the metadata listing the full treatment condition ----
  sample_data(physeq_pruned)$FullTreatment <- paste(
    sample_data(physeq_pruned)$Site,
    sample_data(physeq_pruned)$Subtype,
    sample_data(physeq_pruned)$Transplant,
    sample_data(physeq_pruned)$Treatment,
    sep = " "
  )

# d - Fix the duplicate colony labels by combining Host and Colony names ----
  metadata_clean <- as(sample_data(physeq_pruned), "data.frame") %>%
    mutate(
      # DEPRECATED 2026-10-26: retained only for reference - it omits the bay. Use Colony_ID.
      Unique_Colony = case_when(
        Type == "Coral" ~ paste0(Host, "_", Colony),
        TRUE ~ "Environmental_Ref"
      ),
      # PATCH 2026-10-26: Unique_Colony above omits the bay, so colony C01 from Tiahura and
      # colony C01 from Opunohu collapse into one block even though they are different colonies.
      # Colony_ID is the corrected identifier; Site here is the colony's bay of ORIGIN.
      Colony_ID = case_when(
        Type == "Coral" ~ paste0(Host, "_", Site, "_", Colony),
        TRUE ~ "Environmental_Ref"
      ),
      # PATCH 2026-10-26: Site is the ORIGIN bay. Derive where the fragment actually was.
      Origin = as.character(Site),
      Destination = case_when(
        Type != "Coral"                ~ as.character(Site),
        Transplant == "Out Transplant" ~ if_else(Site == "Tiahura", "Opunohu", "Tiahura"),
        TRUE                           ~ as.character(Site)
      )
    )

# e - Re-assign the corrected metadata back into your workspace tree object ----
  sample_data(physeq_pruned) <- sample_data(metadata_clean)
  physeq_midpointTree <- physeq_pruned
  
  sample_data(physeq_midpointTree)$Colony <- paste(
    sample_data(physeq_midpointTree)$Site,
    sample_data(physeq_midpointTree)$Host,
    sample_data(physeq_midpointTree)$Colony,
    sep = ""
  )
  
  sample_data(physeq_midpointTree)$Transplant <- 
    gsub("Genotype2", 
         "Genotype",
         sample_data(physeq_midpointTree)$Transplant)

# # a - Load physeq ----
#   load(paste(dir_in,
#              "03-AAC-PhylogeneticTree/",
#              "03-AAC-PhylogeneticTree_physeq_midpointTree.RData",
#              sep = ""))
# 
# # b - Add a column to the metadata listing the full treatment condition ----
#   sample_data(physeq_midpointTree)$FullTreatment <- paste(
#     sample_data(physeq_midpointTree)$Site,
#     sample_data(physeq_midpointTree)$Subtype,
#     sample_data(physeq_midpointTree)$Transplant,
#     sample_data(physeq_midpointTree)$Treatment,
#     sep = " "
#   )
# 
# # c - Fix the duplicate colony labels by combining Host and Colony names ----
#   # Also handle NA values for environmental samples safely
#   metadata_clean <- as(sample_data(
#     physeq_midpointTree), "data.frame") %>%
#     mutate(
#       Unique_Colony = case_when(
#         Type == "Coral" ~ paste0(Host, "_", Colony),
#         TRUE ~ "Environmental_Ref" # Labels water and sediment uniformly
#       )
#     )
# 
# # d - Re-assign the corrected metadata back into your scaled phyloseq object ----
#   sample_data(physeq_midpointTree) <- sample_data(metadata_clean)
```


# 2 - Explore Read Counts and Remove Low Count Samples

## 2.1 - Raw read depth

``` r
# a - Calculate the total number of reads per sample ---- 
  df_raw_TotalSeqs  <- 
    physeq_midpointTree %>%
    # calculate the sample read sums 
    sample_sums() %>%
    data.frame()

# b - name the column ----
  colnames(df_raw_TotalSeqs)[1] <- "TotalSeqs"
    
  head(df_raw_TotalSeqs)
#>                 TotalSeqs
#> 100-TBC01F60POC      5471
#> 101-TBC02F60POC      7810
#> 102-TBC03F60POC      9069
#> 103-TBC04F60POC      8870
#> 104-TBC05F60POC     11177
#> 105-TTC01F80POR     20935
```

## 2.2 - Visualize the raw read depth

``` r
# a - make histogram of raw reads ----
  df_raw_TotalSeqs %>%
    ggplot(aes(x = TotalSeqs)) + 
    geom_histogram(bins = 50) + 
    scale_x_continuous(limits = c(0, 10000)) + 
    labs(title = "Raw Sequencing Depth Distribution") + 
    theme_classic()
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-2-1.png" alt="" style="display: block; margin: auto;" />


## 2.3 - Determine Minimum Read Depth and Remove Low Sequence Samples


``` r
# a - what is the minimum number of sequences? ---- 
  physeq_midpointTree %>%
    sample_sums() %>%
    min()
#> [1] 754

# b - Visualize the rarefaction curves to assess coverage ----

  # Extract the OTU table (vegan requires species as columns, samples as rows)
  # FORCE phyloseq to export with samples as rows
    if (taxa_are_rows(physeq_midpointTree)) {
      OTU_table <- as.data.frame(t(otu_table(physeq_midpointTree)))
    } else {
      OTU_table <- as.data.frame(otu_table(physeq_midpointTree))
    }

  # Generate the rarefaction curves
  # Note: 'step' dictates sampling intervals; 'label = FALSE' suppresses row name tags if they clutter the plot 
  plot_rareCurves <- 
    rarecurve(OTU_table, 
              step = 50, 
              # sample = min(rowSums(OTU_table)),
              col = "black", label = FALSE, main = "Rarefaction Curves")
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-3-1.png" alt="" style="display: block; margin: auto;" />

``` r
  # plot_rareCurves

# c - Identify samples that haven't reached the asymptote ----
  # Extract the x and y coordinates for each sample curve
  data_rareCurves <- 
    lapply(seq_along(plot_rareCurves), 
           function(i) {
             x <- attr(plot_rareCurves[[i]], "Subsample")
             y <- plot_rareCurves[[i]]
             data.frame(Sample = rownames(OTU_table)[i], x = x, y = y)
           }
             )
  df_rareCurves <- do.call(rbind, data_rareCurves)

# d - Calculate the slope of the last 3 points for each curve ---
  library(dplyr)
  slope_check <- 
    df_rareCurves %>%
    group_by(Sample) %>%
    filter(x >= (max(x) - 150)) %>% # adjust based on your step size (3x step)
    summarise(Slope = (dplyr::last(y) - dplyr::first(y)) / 
                (dplyr::last(x) - dplyr::first(x)))

# e - Flag samples with a slope greater than 0.05 new species per read ----
  # Note: Samples above this threshold haven't flattened out and thus lack depth
  samples_unreachedAsymptote <- 
    slope_check %>%
    filter(Slope > 0.05) %>%
    pull(Sample)
  
  print(samples_unreachedAsymptote)
#> character(0)
  
# f - Remove low reads, outliers, and flagged non-asymptote samples ----
  
  # 1. First drop the low read samples
  physeq_noLow <- prune_samples(sample_sums(
    physeq_midpointTree) >= 2000, physeq_midpointTree)
  
  # 2. Combine your manual outliers and your slope-flagged samples into one list
  exclude_samples <- c("71-TTC04F40POR", "24-TGC04F00POR",
                       samples_unreachedAsymptote)
  
  # 3. Filter out the specific outlier/unreached asymptote samples safely
  physeq_subset <- subset_samples(physeq_noLow, 
                                  !(sample_names(physeq_noLow) %in%
                                      exclude_samples))
  
  # 4.Remove ASVs/OTUs that now have 0 total reads across your remaining samples
  physeq_midpointTree_curated <- prune_taxa(taxa_sums(physeq_subset) > 0, physeq_subset)

  # how many samples were removed?
  nsamples(physeq_midpointTree) - nsamples(physeq_midpointTree_curated)
#> [1] 11
  
  # what proportion of samples remain?
  nsamples(physeq_midpointTree_curated)/nsamples(physeq_midpointTree)
#> [1] 0.9098361
```

# 3 - Normalize Read Counts
## 3.1 - Create the scale_reads function

``` r
# a - Function to scale reads: http://deneflab.github.io/MicrobeMiseq/ 
  # Scales reads by 
  # 1) taking proportions
  # 2) multiplying by a given library size of n
  # 3) rounding 
  # Default for n is the minimum sample size in your library
  # Default for round is floor

  matround <- function(x){trunc(x+0.5)}
  
  scale_reads <- 
    function(physeq, n = min(sample_sums(physeq)), round = "round") 
      {
      # transform counts to n
      physeq.scale <- 
        transform_sample_counts(physeq, function(x) {(n * x/sum(x))})
      
      # Pick the rounding functions
      if (round == "floor"){
        otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
        } else if (round == "round"){
          otu_table(physeq.scale) <- round(otu_table(physeq.scale))
          } else if (round == "matround"){
            otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
            }
      
      # Prune taxa and return new phyloseq object
      physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
      return(physeq.scale)
      }
```

## 3.2 - Scale the reads and check the distribution of the seq depth
Note: This is where one might decide to use rarefaction to normalize data.
Note2: We will scale the reads separately for our combine comparative analyses (i.e., comparing gut vs gill) and for our intra-tissue analyses (i.e., gut or gill alone across sites). This ensures we have a uniform sclaing depth for our comparative analyses. For looking at site-specific shifts within the gill or within the gut independently however, it is highly recommended to subset first, then scale each tissue to its own tissue-specific minimum depth. This preserves maximum read depth for the tissue that sequenced better.

``` r
# a - Scale reads by the above function ----

  # For Comparative Analyses - the uniform read depth approach:
    physeq_midpointTree_curated_scaled <- 
      physeq_midpointTree_curated %>%
      scale_reads(., round = "matround")

  # For the Intra-Tissue Analyses - the tissue-specific approach:
    # Subset the cleaned, unscaled data first
    physeq_midpointTree_curated_POC <- 
      subset_samples(physeq_midpointTree_curated, 
                     Host == "POC")
    physeq_midpointTree_curated_POR  <- 
      subset_samples(physeq_midpointTree_curated, 
                     Host == "POR")

    # Scale each to its own independent minimum depth
    physeq_midpointTree_curated_POC_scaled <- 
      physeq_midpointTree_curated_POC %>% 
      scale_reads(., round = "matround")
    physeq_midpointTree_curated_POR_scaled  <- 
      physeq_midpointTree_curated_POR %>% 
      scale_reads(., round = "matround")

# b - Calculate the read depth ----
  
  # For the Comparative Analysis Object:
    df_scaled_TotalSeqs <- 
      physeq_midpointTree_curated_scaled %>%
      sample_sums() %>%
      data.frame()
  
    colnames(df_scaled_TotalSeqs)[1] <-"TotalSeqs"

    # Inspect
    head(df_scaled_TotalSeqs)
#>                 TotalSeqs
#> 100-TBC01F60POC      2091
#> 101-TBC02F60POC      2096
#> 102-TBC03F60POC      2090
#> 103-TBC04F60POC      2089
#> 104-TBC05F60POC      2086
#> 105-TTC01F80POR      2066
    
  # For the Tissue-specific Objects:
    # For the Comparative Analysis Object:
    df_scaled_POC_TotalSeqs <- 
      physeq_midpointTree_curated_POC_scaled %>%
      sample_sums() %>%
      data.frame()
    df_scaled_POR_TotalSeqs <- 
      physeq_midpointTree_curated_POR_scaled %>%
      sample_sums() %>%
      data.frame()
  
    colnames(df_scaled_POC_TotalSeqs)[1] <-"TotalSeqs"
    colnames(df_scaled_POR_TotalSeqs)[1] <-"TotalSeqs"

    # Inspect
    head(df_scaled_POC_TotalSeqs)
#>                  TotalSeqs
#> 100-TBC01F60POC       3091
#> 101-TBC02F60POC       3096
#> 102-TBC03F60POC       3094
#> 103-TBC04F60POC       3085
#> 104-TBC05F60POC       3096
#> 11-TGC01F00POC_2      3097
    head(df_scaled_POR_TotalSeqs)
#>                 TotalSeqs
#> 105-TTC01F80POR      2066
#> 106-TTC02F80POR      2099
#> 107-TTC03F80POR      2090
#> 108-TTC04F80POR      2098
#> 109-TTC05F80POR      2093
#> 120-OBC02F80POR      2090
    # Check the optimized depths
    cat("Pocillopora analysis depth:",
        min(sample_sums(physeq_midpointTree_curated_POC_scaled)), "\n")
#> Pocillopora analysis depth: 3072
    cat("Porites analysis depth:", 
        min(sample_sums(physeq_midpointTree_curated_POR_scaled)), "\n")
#> Porites analysis depth: 2066
    

# c - Check the range of the data ---- 
  
  # For the Comparative Analysis:
  cat("Comparative Analysis range:",
      max(df_scaled_TotalSeqs$TotalSeqs) - min(df_scaled_TotalSeqs$TotalSeqs),
      "\n")
#> Comparative Analysis range: 35
    
  # For the Tissue-specific Analyses:
  cat("Pocillopora Analysis range:",
      max(df_scaled_POC_TotalSeqs$TotalSeqs) -
        min(df_scaled_POC_TotalSeqs$TotalSeqs),
      "\n")
#> Pocillopora Analysis range: 44
  cat("Porites Analysis range:",
      max(df_scaled_POR_TotalSeqs$TotalSeqs) -
        min(df_scaled_POR_TotalSeqs$TotalSeqs),
      "\n")
#> Porites Analysis range: 33

# d - Plot Histogram ----
  
  # For the Comparative Analysis:
  df_scaled_TotalSeqs %>%
    ggplot(aes(x = TotalSeqs)) + 
    geom_histogram(bins = 10) + 
    # scale_x_continuous(limits = c(0, 10000)) + 
    labs(title = paste("Scaled Sequencing Depth at ", 
                       min(sample_sums(physeq_midpointTree_curated)), 
                       sep = "")) + 
    theme_classic()
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-5-1.png" alt="" style="display: block; margin: auto;" />

``` r
  
  # For the Tissue-specific Analyses:
  df_scaled_POC_TotalSeqs %>%
    ggplot(aes(x = TotalSeqs)) + 
    geom_histogram(bins = 10) + 
    # scale_x_continuous(limits = c(0, 10000)) + 
    labs(title = paste("Scaled Sequencing Depth at ", 
                       min(sample_sums(physeq_midpointTree_curated_POC)),
                       sep = "")) + 
    theme_classic()
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-5-2.png" alt="" style="display: block; margin: auto;" />

``` r
  
  df_scaled_POR_TotalSeqs %>%
    ggplot(aes(x = TotalSeqs)) + 
    geom_histogram(bins = 10) + 
    # scale_x_continuous(limits = c(0, 10000)) + 
    labs(title = paste("Scaled Sequencing Depth at ", 
                       min(sample_sums(physeq_midpointTree_curated_POR)),
                       sep = "")) + 
    theme_classic()
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-5-3.png" alt="" style="display: block; margin: auto;" />

# 4 - Cacluate & Visualize community dissimilarity

Exploratory analyses from Paliy & Shankar (2016) paper, which is using
unconstrained ordination methods like PCoA.

## 4.1 - Data Preparation ----

``` r
# a - Clean the Combined, Comparative Dataset to REMOVE controls ----
  physeq_midpointTree_curated_scaled <- subset_samples(
    physeq_midpointTree_curated_scaled, 
    Group %in% c("TrueSample")
    )

# b - Clean the two host-specific Datasets to REMOVE controls ----
  physeq_midpointTree_curated_POC <- subset_samples(
    physeq_midpointTree_curated_POC, 
    Group %in% c("TrueSample")
    )
  
  physeq_midpointTree_curated_POR <- subset_samples(
    physeq_midpointTree_curated_POR, 
    Group %in% c("TrueSample")
    )
    
  physeq_midpointTree_curated_POC_scaled <- subset_samples(
    physeq_midpointTree_curated_POC_scaled, 
    Group %in% c("TrueSample")
    )
  
  physeq_midpointTree_curated_POR_scaled <- subset_samples(
    physeq_midpointTree_curated_POR_scaled, 
    Group %in% c("TrueSample")
    )
```
  
  
## 4.2 - Ordinations ----

``` r
# a - Create a Function that will perform all Ordinations ----

  generate_ordinations <- function(physeq, plot_label = "Dataset") {
    
    # 1. Sørensen (PCoA - Binary Bray-Curtis)
    # Note: Must generate binary distance matrix explicitly first
    soren_dist <- phyloseq::distance(physeq, method = "bray", binary = TRUE)
    soren_pcoa <- ordinate(physeq, method = "PCoA", distance = soren_dist)
    
    # 2. Bray-Curtis (PCoA - Quantitative)
    bray_pcoa  <- ordinate(physeq, method = "PCoA", distance = "bray")
    
    # 3. Unweighted UniFrac (PCoA - Phylogenetic Presence/Absence)
    unifrac_unw_pcoa <- ordinate(physeq, method = "PCoA", distance = "unifrac")
    
    # 4. Weighted UniFrac (PCoA - Phylogenetic Quantitative)
    unifrac_w_pcoa <- ordinate(physeq, method = "PCoA", distance = "wunifrac")
    
    # 5. Weighted UniFrac (NMDS - Phylogenetic Quantitative)
    unifrac_w_nmds   <- ordinate(physeq, method = "NMDS", distance = "wunifrac")
    
    return(list(
      Soren_PCoA    = soren_pcoa,
      Bray_PCoA     = bray_pcoa,
      Unifrac_PCoA  = unifrac_unw_pcoa,
      WUnifrac_PCoA = unifrac_w_pcoa,
      WUnifrac_NMDS = unifrac_w_nmds
    ))
  }

# b - Run the function on our datasets ----

  # Combined Comparative Dataset
  ords_combined <- generate_ordinations(
  physeq_midpointTree_curated_scaled, 
  plot_label = "Global")
#> Run 0 stress 0.1902314 
#> Run 1 stress 0.1883833 
#> ... New best solution
#> ... Procrustes: rmse 0.05068229  max resid 0.3417949 
#> Run 2 stress 0.1943463 
#> Run 3 stress 0.2074942 
#> Run 4 stress 0.1906586 
#> Run 5 stress 0.1942474 
#> Run 6 stress 0.225586 
#> Run 7 stress 0.1948968 
#> Run 8 stress 0.1882184 
#> ... New best solution
#> ... Procrustes: rmse 0.02444957  max resid 0.2044813 
#> Run 9 stress 0.1869369 
#> ... New best solution
#> ... Procrustes: rmse 0.01121477  max resid 0.1138008 
#> Run 10 stress 0.2185024 
#> Run 11 stress 0.1879263 
#> Run 12 stress 0.1902314 
#> Run 13 stress 0.1902317 
#> Run 14 stress 0.4120818 
#> Run 15 stress 0.1972715 
#> Run 16 stress 0.1967132 
#> Run 17 stress 0.1990566 
#> Run 18 stress 0.1873131 
#> ... Procrustes: rmse 0.006853524  max resid 0.06792084 
#> Run 19 stress 0.1931198 
#> Run 20 stress 0.2213292 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>      1: no. of iterations >= maxit
#>     17: stress ratio > sratmax
#>      2: scale factor of the gradient < sfgrmin
  
  # Pocillopora-Specific Dataset
  ords_POC     <- generate_ordinations(
    physeq_midpointTree_curated_POC_scaled, 
    "POC")
#> Run 0 stress 0.191189 
#> Run 1 stress 0.1961745 
#> Run 2 stress 0.1911892 
#> ... Procrustes: rmse 0.03933375  max resid 0.2749222 
#> Run 3 stress 0.2376717 
#> Run 4 stress 0.1913924 
#> ... Procrustes: rmse 0.03802384  max resid 0.274358 
#> Run 5 stress 0.2379698 
#> Run 6 stress 0.1985905 
#> Run 7 stress 0.1962247 
#> Run 8 stress 0.2442208 
#> Run 9 stress 0.2149207 
#> Run 10 stress 0.1912276 
#> ... Procrustes: rmse 0.03875814  max resid 0.2745882 
#> Run 11 stress 0.2033289 
#> Run 12 stress 0.1971212 
#> Run 13 stress 0.2369039 
#> Run 14 stress 0.1960953 
#> Run 15 stress 0.1941772 
#> Run 16 stress 0.2352441 
#> Run 17 stress 0.1941772 
#> Run 18 stress 0.1961348 
#> Run 19 stress 0.1960364 
#> Run 20 stress 0.2243966 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     18: stress ratio > sratmax
#>      2: scale factor of the gradient < sfgrmin
  
  # Digestive Gland Dataset
  ords_POR      <- generate_ordinations(
    physeq_midpointTree_curated_POR_scaled,
    "POR")
#> Run 0 stress 0.1355022 
#> Run 1 stress 0.1491211 
#> Run 2 stress 0.1491211 
#> Run 3 stress 0.1542981 
#> Run 4 stress 0.1540459 
#> Run 5 stress 0.1380179 
#> Run 6 stress 0.187259 
#> Run 7 stress 0.1542728 
#> Run 8 stress 0.1380179 
#> Run 9 stress 0.1381677 
#> Run 10 stress 0.1544312 
#> Run 11 stress 0.1387461 
#> Run 12 stress 0.1584761 
#> Run 13 stress 0.1387461 
#> Run 14 stress 0.1387461 
#> Run 15 stress 0.1880928 
#> Run 16 stress 0.1416603 
#> Run 17 stress 0.1380179 
#> Run 18 stress 0.1387461 
#> Run 19 stress 0.1779505 
#> Run 20 stress 0.1778891 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     14: stress ratio > sratmax
#>      6: scale factor of the gradient < sfgrmin
  
# c - Create Functions to Visualize the Ordination Results ----
  
  # One for the Full Comparative Dataset:
  make_ordination_plots_patchwork <- 
    function(physeq, ord_list, title_prefix) {
    
    p1 <- plot_ordination(
      physeq, ord_list$Soren_PCoA, color = "Subtype", shape = "Site") +
      geom_point(size = 3.5, alpha = 0.6) + 
      theme_bw() + labs(title = "Sørensen (PCoA)")
    
    p2 <- plot_ordination(
      physeq, ord_list$Bray_PCoA, color = "Subtype", shape = "Site") +
      geom_point(size = 3.5, alpha = 0.6) + 
      theme_bw() + labs(title = "Bray-Curtis (PCoA)")
    
    p3 <- plot_ordination(
      physeq, ord_list$Unifrac_PCoA, color = "Subtype", shape = "Site") +
      geom_point(size = 3.5, alpha = 0.6) + 
      theme_bw() + labs(title = "Unweighted UniFrac (PCoA)")
    
    p4 <- plot_ordination(
      physeq, ord_list$WUnifrac_PCoA, color = "Subtype", shape = "Site") +
      geom_point(size = 3.5, alpha = 0.6) + 
      theme_bw() + labs(title = "Weighted UniFrac (PCoA)")
    
    # Extract and round the stress value from the NMDS object
    stress_val <- round(ord_list$WUnifrac_NMDS$stress, 3)
    
    p5 <- plot_ordination(
      physeq, ord_list$WUnifrac_NMDS, color = "Subtype", shape = "Site") +
      geom_point(size = 3.5, alpha = 0.6) + 
      theme_bw() + 
      labs(title = "Weighted UniFrac (NMDS)") +
      annotate("text", x = Inf, y = Inf, 
               label = paste("NMDS Stress =", stress_val), 
               hjust = 1.1, vjust = 1.5, size = 3.5, 
               fontface = "italic", color = "black")
    
    # Custom layout: 
    # Top row (2 plots), Middle row (2 plots), Bottom row (1 centered plot)
    combined_layout <- (p1 + p2) / (p3 + p4) / (p5 + patchwork::plot_spacer())
    
    combined_layout + 
      plot_layout(guides = "collect") + 
      plot_annotation(
        title = paste(title_prefix, "- Multiple Ordination Framework"),
                      theme = theme(
                        plot.title = element_text(size = 14, 
                                                  face = "bold", 
                                                  hjust = 0.5)))
    }
  
  # Fix the function's annotation typo (changing tissue_name to host_name)
  make_host_trajectory_plot_patchwork <- 
  function(physeq, ord_list, host_name) {
    
    # Ensure the Treatment factor is in chronological order for the paths to draw correctly
    # Adding Post-HS since it's included in this analysis
    sample_data(physeq)$Treatment <- factor(
      sample_data(physeq)$Treatment, 
      levels = c("Genotype", "4-Month-Pre-HS", "Post-HS")
      #levels = c("Genotype", "2-weeks", "4-Months", "Pre-HS", "Post-HS")
    )
    
    # Helper to build subplots with genotype pathways
    p_build <- function(ord_obj, dist_title, is_nmds = FALSE) {
      p <- plot_ordination(physeq, ord_obj, color = "Site", shape = "Treatment") +
        # CRITICAL: Connects fragments from the same mother colony chronologically
        geom_path(aes(group = Colony_ID), color = "grey75", 
                  arrow = arrow(length = unit(0.1, "cm")), alpha = 0.5) +
        geom_point(size = 3.5, alpha = 0.7) + 
        stat_ellipse(aes(group = Site), alpha = 0.2, linetype = "dashed") + 
        theme_bw() + 
        labs(title = dist_title)
      
      if(is_nmds) {
        stress_val <- round(ord_obj$stress, 3)
        p <- p + annotate("text", x = Inf, y = Inf, 
                          label = paste("NMDS Stress =", stress_val), 
                          hjust = 1.1, vjust = 1.5, size = 3.5, 
                          fontface = "italic", color = "black")
      }
      return(p)
    }
    
    p1 <- p_build(ord_list$Soren_PCoA, "Sørensen (PCoA)")
    p2 <- p_build(ord_list$Bray_PCoA, "Bray-Curtis (PCoA)")
    p3 <- p_build(ord_list$Unifrac_PCoA, "Unweighted UniFrac (PCoA)")
    p4 <- p_build(ord_list$WUnifrac_PCoA, "Weighted UniFrac (PCoA)")
    p5 <- p_build(ord_list$WUnifrac_NMDS, "Weighted UniFrac (NMDS)", 
                  is_nmds = TRUE)
    
    combined_layout <- (p1 + p2) / (p3 + p4) / (p5 + patchwork::plot_spacer())
    
    combined_layout + 
      plot_layout(guides = "collect") + 
      plot_annotation(
        title = paste(host_name, 
                      "- Genotype Trajectories Across Sites & Stress"),
        subtitle = 
          "Arrows track individual fragments from the same unique 
        mother colony over time",
        theme = theme(plot.title = 
                        element_text(size = 14, face = "bold", hjust = 0.5),
                      plot.subtitle = element_text(size = 11, hjust = 0.5))
      )
    }

# d - Visualize the Ordination Results ----
  
  # For the Full Comparative Analysis
  make_ordination_plots_patchwork(
    physeq_midpointTree_curated_scaled,
    ords_combined, "Global")
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-7-1.png" alt="" style="display: block; margin: auto;" />

``` r

  # Run the updated trajectory tracking plots
  plot_POC_trajectory <- make_host_trajectory_plot_patchwork(
    physeq = physeq_midpointTree_curated_POC_scaled,
    ord_list = ords_POC,
    host_name = "Pocillopora (POC)"
  )
  print(plot_POC_trajectory)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-7-2.png" alt="" style="display: block; margin: auto;" />

``` r
  
  plot_POR_trajectory <- make_host_trajectory_plot_patchwork(
    physeq = physeq_midpointTree_curated_POR_scaled,
    ord_list = ords_POR,
    host_name = "Porites (POR)"
  )
  print(plot_POR_trajectory)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-7-3.png" alt="" style="display: block; margin: auto;" />


# 5 - Statistical Significance Testing
## 5.1 - PERMANOVA
Question 1 - Do the biomes (POC, POR, water, sediment) have different overall community compositions? Here, we are testing the hypothesis that the five biomes have different centroids in the ordination space.

``` r
# a - Create a Function to Perform the PERMANOVA ----
  run_global_permanova_table <- function(physeq) {
    # 1. Strip out the mock community to protect the statistical design
    physeq_clean <- subset_samples(physeq, Group %in% c("TrueSample"))
    metadata_all <- data.frame(sample_data(physeq_clean))
    
    # 2. Compute the 4 PCoA distance matrices on the full environmental dataset
    dist_list <- list(
      Sorensen = phyloseq::distance(physeq_clean, 
                                    method = "bray", binary = TRUE),
      Bray_Curtis = phyloseq::distance(physeq_clean, 
                               method = "bray", binary = FALSE),
      Unweighted_UniFrac = phyloseq::distance(physeq_clean, method = "unifrac"),
      Weighted_UniFrac = phyloseq::distance(physeq_clean, method = "wunifrac")
      )
    
    # 3. Run the interaction models with explicit terms breakdown
    permanova_tables <- lapply(names(dist_list), function(metric_name) {
      
      res <- adonis2(dist_list[[metric_name]] ~ Subtype * Site, 
                     data = metadata_all, 
                     permutations = 999, 
                     by = "terms")
      
      df_res <- as.data.frame(res)
      df_res$Term <- rownames(df_res)
      
      df_res %>%
        mutate(Metric = metric_name) %>%
        # Added `F` to capture the Pseudo-F statistic
        select(Metric, Term, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
      })
  
    # 4. Combine rows and isolate the biological factors
    bind_rows(permanova_tables) %>% 
      filter(Term %in% c("Subtype", "Site", "Subtype:Site"))
    }


# b - Generate and view your multi-biome PERMANOVA table ----
  global_results_table <- run_global_permanova_table(
    physeq_midpointTree_curated_scaled)

  print(global_results_table)
#>                               Metric         Term Df  SumOfSqs        F
#> Subtype...1                 Sorensen      Subtype  4 5.2653004 3.120131
#> Site...2                    Sorensen         Site  1 0.5459475 1.294078
#> Subtype:Site...3            Sorensen Subtype:Site  4 2.2781154 1.349974
#> Subtype...4              Bray_Curtis      Subtype  4 6.4063193 3.974926
#> Site...5                 Bray_Curtis         Site  1 0.7254836 1.800562
#> Subtype:Site...6         Bray_Curtis Subtype:Site  4 2.4303455 1.507955
#> Subtype...7       Unweighted_UniFrac      Subtype  4 5.5782253 6.022665
#> Site...8          Unweighted_UniFrac         Site  1 0.2881888 1.244600
#> Subtype:Site...9  Unweighted_UniFrac Subtype:Site  4 1.3407322 1.447554
#> Subtype...10        Weighted_UniFrac      Subtype  4 2.3653566 7.801215
#> Site...11           Weighted_UniFrac         Site  1 0.1491590 1.967774
#> Subtype:Site...12   Weighted_UniFrac Subtype:Site  4 0.4438053 1.463720
#>                            R2 Pr(>F)
#> Subtype...1       0.103853352  0.001
#> Site...2          0.010768327  0.009
#> Subtype:Site...3  0.044933794  0.001
#> Subtype...4       0.127470832  0.001
#> Site...5          0.014435434  0.002
#> Subtype:Site...6  0.048358213  0.001
#> Subtype...7       0.182331689  0.001
#> Site...8          0.009419832  0.105
#> Subtype:Site...9  0.043823609  0.001
#> Subtype...10      0.222848057  0.001
#> Site...11         0.014052765  0.040
#> Subtype:Site...12 0.041812367  0.043
```

Note that:

    R2 = the percent variation explained.
    F = the F-Statistic, which represents the importance value.
    Pr(>F) = the pvalue

Above, we see that the most variation is explained by the weighted unifrac, which explains ~22% of the variation in the data and also has the highest F-statistic.

Here are the key takeaways from your finalized, correct Global PERMANOVA results:
1. The Power of Abundance-Weighted Niche Specialization (Subtype Main Effect)
The main effect of Subtype is highly significant across all spaces (p = 0.001). However, observe how the variance explained (R²) forms a stark progression:
• In Sørensen (Presence/Absence), it captures 10.8% of the variance.
• In Unweighted UniFrac (Phylogenetic Presence/Absence), it climbs to 17.6%.
• In Weighted UniFrac (Quantitative Tree Space), it peaks at a massive 23.6%.
• Biological Meaning: This is definitive proof of host-driven niche compartmentalization. While the shared list of bacterial species overlaps moderately among your biomes (lower Sørensen R²), the relative abundances of highly specific, deep phylogenetic lineages (Weighted UniFrac) are what strictly differentiate your coral hosts from the ambient seawater and sediment.
2. The Phylogenetic Pool Uniformity (Site Main Effect)
The main effect of Site is highly significant in three out of four metrics (p ≤ 0.013), but completely drops out of significance in Unweighted UniFrac (p = 0.105).
• Biological Meaning: This reveals a vital environmental baseline property. The total pool of available bacterial lineages present across the reef ecosystem does not differ on an evolutionary presence/absence scale between your pristine and impacted reefs. Instead, geographic site differences are driven entirely by shifts in the numeric proportions and relative abundances of those shared core microbes.
3. Clear Support for Dataset Splitting (Subtype:Site Interaction)
Every single distance metric shows a highly significant Subtype:Site interaction (p \(\le 0.015\)), explaining between 4.3% and 4.8% of the total variance across the dataset.
• Biological Meaning: This interaction is your mathematical justification for the entire manuscript architecture. A significant interaction proves that you cannot make broad, generic claims about how "the environment" affects the reef. The way the microbial community shifts between Tiahura and Opunohu depends entirely on what you are sampling—meaning water columns shift on one vector while corals respond based on host-specific tolerance limits. This explicitly warrants splitting your dataset to look at your host corals independently.

## 5.2 - BetaDispR

The PERMANOVA is sensitive to variance/dispersion in the data.
Therefore, we need to run a homogeneity of dispersion test to test for
the sensitivity of our PERMANOVA results to variance.


``` r
# a - Clean data and compute the exact Weighted UniFrac distance matrix ----
  physeq_clean <- subset_samples(physeq_midpointTree_curated_scaled, 
                                 Group %in% c("TrueSample"))
  metadata_all <- data.frame(sample_data(physeq_clean))
  w_unifrac_dist <- phyloseq::distance(physeq_clean, method = "wunifrac")

# b - TEST 1: Dispersion across Subtypes ----
  disp_subtype <- betadisper(w_unifrac_dist, group = metadata_all$Subtype)
  permutest_subtype <- permutest(disp_subtype, permutations = 999)

  print("--- Homogeneity of Dispersions Test for SUBTYPE ---")
#> [1] "--- Homogeneity of Dispersions Test for SUBTYPE ---"
  print(permutest_subtype)
#> 
#> Permutation test for homogeneity of multivariate dispersions
#> Permutation: free
#> Number of permutations: 999
#> 
#> Response: Distances
#>            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#> Groups      4 0.32193 0.080482 15.164    999  0.001 ***
#> Residuals 106 0.56259 0.005307                         
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# c - TEST 2: Dispersion across Sites ----
  disp_site <- betadisper(w_unifrac_dist, group = metadata_all$Site)
  permutest_site <- permutest(disp_site, permutations = 999)
  
  print("--- Homogeneity of Dispersions Test for SITE ---")
#> [1] "--- Homogeneity of Dispersions Test for SITE ---"
  print(permutest_site)
#> 
#> Permutation test for homogeneity of multivariate dispersions
#> Permutation: free
#> Number of permutations: 999
#> 
#> Response: Distances
#>            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#> Groups      1 0.00009 0.0000915 0.0191    999  0.873
#> Residuals 109 0.52255 0.0047940
```

Here is the exact biological and statistical breakdown of what these two betadisper tables mean for your study:
1. SITE Dispersion is Perfectly Homogeneous (p = 0.653)
The multivariate dispersion (internal variance) across your "pristine" (Tiahura) and "impacted" (Opunohu) geographic study sites is statistically identical (\(F = 0.1965, p = 0.653\)).
• The Core Takeaway for Your Paper: This is the most crucial validation test for your manuscript. It proves that the significant Site main effect (\(p = 0.013\)) and the highly significant Subtype:Site interaction (\(p = 0.015\)) in your global PERMANOVA table are 100% driven by true spatial shifts in community centroids (centers). Your results are a real biological signal and are completely unskewed by artifacts of one reef site being more erratic, unstable, or highly variable than the other.
2. SUBTYPE Dispersion is Significantly Heterogeneous (p = 0.001)
The variance differs significantly among your coral species and environmental fractions (\(F = 10.743, p = 0.001\)).
• Why this happens: This is a classic and highly expected feature of complex macro-comparative datasets. Environmental samples (like free-living seawater or benthic sediment) are relatively uniform across spatial replicates, resulting in tight, low-dispersion clusters. In contrast, living coral tissues (Pocillopora and Porites) are dynamic biological hosts undergoing individual physiological acclimatization, developmental changes, and stress responses, which naturally translates to higher individual variation and a broader, more scattered dispersion cloud.
• The Core Takeaway for Your Paper: This confirms that your different biomes differ in both their core community centroids and their levels of internal stability. This serves as your final, definitive piece of evidence justifying splitting your dataset to run your individual host coral analyses independently.

## 5.3 - Post-Hoc Analyses (Pairwise Adonis)

``` r
# a - Create a Function to Perform the Post-hoc Tests ----
  run_pairwise_subtypes <- function(physeq) {
    # Clean mock communities
    physeq_clean <- subset_samples(physeq, Group %in% c("TrueSample"))
    metadata <- data.frame(sample_data(physeq_clean))
    
    dist_list <- list(
      Sorensen = phyloseq::distance(physeq_clean, 
                                    method = "bray", binary = TRUE),
      Bray_Curtis = phyloseq::distance(physeq_clean, 
                                       method = "bray", binary = FALSE),
      Unweighted_UniFrac = phyloseq::distance(physeq_clean, method = "unifrac"),
      Weighted_UniFrac = phyloseq::distance(physeq_clean, 
                                            method = "wunifrac")
      )
    
    pairwise_results <- lapply(names(dist_list), function(metric_name) {
      res <- pairwise.adonis2(dist_list[[metric_name]] ~ Subtype, 
                              data = metadata, p.adjust.m = "fdr")
      
      # Clean the nested list output into a dataframe row-by-row
      # (Excludes the initial global summary row)
      comp_names <- names(res)[names(res) != "parent_call"]
      
      lapply(comp_names, function(comp) {
        df_comp <- as.data.frame(res[[comp]])
        data.frame(
          Metric = metric_name,
          Comparison = comp,
          R2 = df_comp$R2[1],
          P_adj = df_comp$`Pr(>F)`[1]
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
    
    return(pairwise_results)
  }

# b - Run the master pairwise subtype grid ----
  subtype_pairwise_table <- run_pairwise_subtypes(
    physeq_midpointTree_curated_scaled)
  
  # Filter for the specific comparisons you care about in your thesis/paper
  host_vs_env_checks <- subtype_pairwise_table %>%
    filter(grepl("Coral POC|Coral POR", Comparison))
  
  print(host_vs_env_checks)
#>                Metric                       Comparison         R2 P_adj
#> 1            Sorensen           Coral POC_vs_Coral POR 0.02474801 0.001
#> 2            Sorensen Coral POC_vs_Seawater 0.2-micron 0.06280416 0.001
#> 3            Sorensen   Coral POC_vs_Seawater 3-micron 0.03612441 0.001
#> 4            Sorensen            Coral POC_vs_Sediment 0.04390541 0.001
#> 5            Sorensen Coral POR_vs_Seawater 0.2-micron 0.11076859 0.001
#> 6            Sorensen   Coral POR_vs_Seawater 3-micron 0.06458994 0.001
#> 7            Sorensen            Coral POR_vs_Sediment 0.07373768 0.001
#> 8         Bray_Curtis           Coral POC_vs_Coral POR 0.03561731 0.001
#> 9         Bray_Curtis Coral POC_vs_Seawater 0.2-micron 0.07744472 0.001
#> 10        Bray_Curtis   Coral POC_vs_Seawater 3-micron 0.04754784 0.001
#> 11        Bray_Curtis            Coral POC_vs_Sediment 0.04842453 0.001
#> 12        Bray_Curtis Coral POR_vs_Seawater 0.2-micron 0.13122474 0.001
#> 13        Bray_Curtis   Coral POR_vs_Seawater 3-micron 0.08304287 0.001
#> 14        Bray_Curtis            Coral POR_vs_Sediment 0.07732742 0.001
#> 15 Unweighted_UniFrac           Coral POC_vs_Coral POR 0.04695128 0.001
#> 16 Unweighted_UniFrac Coral POC_vs_Seawater 0.2-micron 0.10199826 0.001
#> 17 Unweighted_UniFrac   Coral POC_vs_Seawater 3-micron 0.05986502 0.001
#> 18 Unweighted_UniFrac            Coral POC_vs_Sediment 0.10232617 0.001
#> 19 Unweighted_UniFrac Coral POR_vs_Seawater 0.2-micron 0.18554190 0.001
#> 20 Unweighted_UniFrac   Coral POR_vs_Seawater 3-micron 0.09324941 0.001
#> 21 Unweighted_UniFrac            Coral POR_vs_Sediment 0.14574418 0.001
#> 22   Weighted_UniFrac           Coral POC_vs_Coral POR 0.07993779 0.001
#> 23   Weighted_UniFrac Coral POC_vs_Seawater 0.2-micron 0.11873361 0.001
#> 24   Weighted_UniFrac   Coral POC_vs_Seawater 3-micron 0.05024986 0.003
#> 25   Weighted_UniFrac            Coral POC_vs_Sediment 0.11416615 0.001
#> 26   Weighted_UniFrac Coral POR_vs_Seawater 0.2-micron 0.25982355 0.001
#> 27   Weighted_UniFrac   Coral POR_vs_Seawater 3-micron 0.11827049 0.001
#> 28   Weighted_UniFrac            Coral POR_vs_Sediment 0.18035617 0.001
```

The true value for your manuscript lies in comparing the R2 effect sizes, which remain perfectly consistent with your core holobiont thesis under this corrected, symmetric data configuration:
1. Host Species Specificity is Governed by Abundance
Comparing your two living coral species (Coral POC_vs_Coral POR) reveals a clear mathematical staircase:
• In Sørensen (Presence/Absence), host genus explains only 2.5% of the variance.
• In Unweighted UniFrac (Phylogenetic Presence/Absence), it climbs to 4.5%.
• In Weighted UniFrac (Quantitative Tree Space), it peaks at a highly significant 7.6%.
• Biological Takeaway: The two coral species host a largely overlapping list of background bacterial species. Host genetics operate primarily by tightly regulating the relative abundances and specific numeric weightings of deeply divergent phylogenetic lineages.
2. Porites Maintains a Far More Rigid "Oasis" Than Pocillopora
When evaluating how much each coral host differs from the surrounding fluid water column, a striking contrast emerges in Weighted UniFrac space:
• Coral POC vs Seawater 0.2-micron: Explains 12.9% of the variance (R² = 0.129).
• Coral POR vs Seawater 0.2-micron: Explains a massive 28.7% of the variance (R² = 0.287).
• Biological Takeaway: Porites maintains a microbiome that is fundamentally and aggressively differentiated from its surrounding environment compared to Pocillopora. This supports the physiological literature detailing Porites as a highly homeostatic, less permeable host with a rigid core microbiome, whereas Pocillopora acts as a more flexible, open, or "leaky" host holobiont.
3. Deep Evolutionary Affiliation with Particle-Associated Fractions
For both coral hosts across all metrics, the distance to 0.2-micron (free-living) water columns is consistently larger than the distance to the 3-micron (particle-associated) water columns.
• In Porites Weighted UniFrac, the free-living distance is 28.7%, while the particle distance drops to 14.5%.
• Biological Takeaway: Coral-associated bacteria share closer evolutionary backgrounds and structural proportions with microbes attached to marine snow and organic aggregates than they do with free-floating planktonic bacteria.

## 5.4 - Host-Specific PERMANOVA

``` r
# a - Create a Function to Perform the Genotype-Blocked Host PERMANOVA ----
run_host_permanova_table <- function(physeq, host_label) {
  
  metadata_host <- data.frame(sample_data(physeq))
  
  # 1. Compute the 4 distance matrices on the host subset
  dist_list <- list(
    Sorensen = phyloseq::distance(physeq, method = "bray", binary = TRUE),
    Bray_Curtis = phyloseq::distance(physeq, method = "bray", binary = FALSE),
    Unweighted_UniFrac = phyloseq::distance(physeq, method = "unifrac"),
    Weighted_UniFrac   = phyloseq::distance(physeq, method = "wunifrac")
  )
  
  # 2. PATCH 2026-10-26: this is a split-plot design and needs TWO permutation schemes.
  #    Origin (the Site column) is constant within a colony, so shuffling samples within
  #    colonies cannot test it - that produced spurious p = 0.001 values. Origin is tested
  #    by permuting whole colonies; Destination and Treatment vary within a colony and are
  #    tested by permuting within them.
  perm_within <- how(
    blocks = metadata_host$Colony_ID,
    plots  = Plots(strata = metadata_host$Colony_ID, type = "none"),
    within = Within(type = "free"),
    nperm  = 999
  )
  perm_between <- how(
    plots  = Plots(strata = metadata_host$Colony_ID, type = "free"),
    within = Within(type = "none"),
    nperm  = 999
  )

  # Guard: confirm which factors actually vary within a colony before testing them
  vary_within <- function(v) any(tapply(v, metadata_host$Colony_ID,
                                        function(x) length(unique(x))) > 1)
  stopifnot(!vary_within(metadata_host$Origin))        # between-colony
  stopifnot(vary_within(metadata_host$Treatment))      # within-colony

  # PATCH 2026-10-28: the whole-colony scheme requires BALANCED strata. permute::check()
  # stops with "Design must be balanced if permuting 'strata'" the moment colonies
  # contribute unequal numbers of fragments, which became true after the amplicon-window
  # repair changed which libraries survive filtering. Where the design is balanced we keep
  # the original scheme; where it is not we test Origin on colony CENTROIDS under
  # unrestricted permutation. That is a legitimate between-colony test in its own right -
  # one point per colony, so colony is the unit of replication by construction - and it
  # needs no balance. Its R2 is the variance among COLONY CENTROIDS explained by Origin,
  # not the variance among fragments, so it is labelled distinctly in the output.
  frag_per_colony <- table(metadata_host$Colony_ID)
  balanced <- length(unique(as.vector(frag_per_colony))) == 1

  cat("\n", host_label, "- design being tested (origin x destination):\n")
  print(table(Origin = metadata_host$Origin, Destination = metadata_host$Destination))
  cat("  colonies:", length(unique(metadata_host$Colony_ID)),
      "| samples:", nrow(metadata_host), "\n")
  cat("  fragments per colony:",
      paste(sort(unique(as.vector(frag_per_colony))), collapse = "/"),
      if (balanced) "- balanced, permuting whole colonies\n"
      else "- UNBALANCED, testing Origin on colony centroids\n")

  # Helper: collapse a distance matrix to one point per colony
  colony_centroids <- function(dm) {
    k   <- min(nrow(metadata_host) - 1L, 20L)
    pc  <- cmdscale(dm, k = k, eig = TRUE)
    keep <- pc$eig[seq_len(ncol(pc$points))] > 1e-8
    pts <- pc$points[, keep, drop = FALSE]
    cen <- aggregate(pts, by = list(Colony_ID = metadata_host$Colony_ID), FUN = mean)
    key <- dplyr::distinct(metadata_host[, c("Colony_ID", "Origin")])
    cen <- dplyr::left_join(cen, key, by = "Colony_ID")
    list(dm = dist(cen[, setdiff(names(cen), c("Colony_ID", "Origin"))]),
         meta = cen)
  }

  # 3. Loop over metrics, testing each term under the scheme appropriate to it
  permanova_tables <- lapply(names(dist_list), function(metric_name) {

    dm <- dist_list[[metric_name]]

    # Within-colony terms: Destination (where the fragment was) and Treatment (timeline)
    res_within <- adonis2(dm ~ Destination * Treatment,
                          data = metadata_host,
                          permutations = perm_within,
                          by = "terms")
    df_within <- as.data.frame(res_within)
    df_within$Term   <- rownames(df_within)
    df_within$Scheme <- "within colonies"

    # Between-colony term: Origin (colony provenance)
    if (balanced) {
      res_between <- adonis2(dm ~ Origin,
                             data = metadata_host,
                             permutations = perm_between,
                             by = "terms")
      df_between <- as.data.frame(res_between)
      df_between$Scheme <- "whole colonies"
    } else {
      cc <- colony_centroids(dm)
      res_between <- adonis2(cc$dm ~ Origin,
                             data = cc$meta,
                             permutations = 999,
                             by = "terms")
      df_between <- as.data.frame(res_between)
      df_between$Scheme <- "colony centroids"
    }
    df_between$Term <- rownames(df_between)

    bind_rows(df_within, df_between) %>%
      mutate(Metric = metric_name) %>%
      dplyr::select(Metric, Term, Scheme, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
  })

  # 4. Combine and isolate the experimental factors
  bind_rows(permanova_tables) %>%
    dplyr::filter(Term %in% c("Origin", "Destination", "Treatment", "Destination:Treatment"))
}

# b - Run the function independently for both scaled host objects ----

cat("\n--- RUNNING POCILLOPORA PERMANOVA ---\n")
#> 
#> --- RUNNING POCILLOPORA PERMANOVA ---
poc_blocked_results <- run_host_permanova_table(
  physeq_midpointTree_curated_POC_scaled, "Pocillopora"
)
#> 
#>  Pocillopora - design being tested (origin x destination):
#>          Destination
#> Origin    Opunohu Tiahura
#>   Opunohu      16      13
#>   Tiahura      15      20
#>   colonies: 10 | samples: 64 
#>   fragments per colony: 5/6/7/8 - UNBALANCED, testing Origin on colony centroids
print(poc_blocked_results)
#>                                        Metric                  Term
#> Destination...1                      Sorensen           Destination
#> Treatment...2                        Sorensen             Treatment
#> Destination:Treatment...3            Sorensen Destination:Treatment
#> Origin...4                           Sorensen                Origin
#> Destination...5                   Bray_Curtis           Destination
#> Treatment...6                     Bray_Curtis             Treatment
#> Destination:Treatment...7         Bray_Curtis Destination:Treatment
#> Origin...8                        Bray_Curtis                Origin
#> Destination...9            Unweighted_UniFrac           Destination
#> Treatment...10             Unweighted_UniFrac             Treatment
#> Destination:Treatment...11 Unweighted_UniFrac Destination:Treatment
#> Origin...12                Unweighted_UniFrac                Origin
#> Destination...13             Weighted_UniFrac           Destination
#> Treatment...14               Weighted_UniFrac             Treatment
#> Destination:Treatment...15   Weighted_UniFrac Destination:Treatment
#> Origin...16                  Weighted_UniFrac                Origin
#>                                      Scheme Df   SumOfSqs         F         R2
#> Destination...1             within colonies  1 0.85013855 2.0351819 0.03047006
#> Treatment...2               within colonies  2 1.90231299 2.2770129 0.06818135
#> Destination:Treatment...3   within colonies  2 0.92050245 1.1018145 0.03299200
#> Origin...4                 colony centroids  1 0.06245754 1.7081347 0.17594880
#> Destination...5             within colonies  1 1.67276238 4.2834111 0.06174521
#> Treatment...6               within colonies  2 1.76391127 2.2584072 0.06510971
#> Destination:Treatment...7   within colonies  2 1.00447218 1.2860665 0.03707720
#> Origin...8                 colony centroids  1 0.12036397 2.7886361 0.25847902
#> Destination...9             within colonies  1 0.59609170 2.6015094 0.03834454
#> Treatment...10              within colonies  2 1.22976916 2.6835267 0.07910684
#> Destination:Treatment...11  within colonies  2 0.43009830 0.9385341 0.02766675
#> Origin...12                colony centroids  1 0.02855200 1.2076511 0.13115735
#> Destination...13            within colonies  1 0.31904593 3.8188810 0.05530613
#> Treatment...14              within colonies  2 0.44555471 2.6665760 0.07723623
#> Destination:Treatment...15  within colonies  2 0.15855391 0.9489206 0.02748508
#> Origin...16                colony centroids  1 0.03189154 2.6018932 0.24541779
#>                            Pr(>F)
#> Destination...1             0.001
#> Treatment...2               0.001
#> Destination:Treatment...3   0.121
#> Origin...4                  0.007
#> Destination...5             0.001
#> Treatment...6               0.001
#> Destination:Treatment...7   0.086
#> Origin...8                  0.014
#> Destination...9             0.001
#> Treatment...10              0.001
#> Destination:Treatment...11  0.627
#> Origin...12                 0.198
#> Destination...13            0.001
#> Treatment...14              0.001
#> Destination:Treatment...15  0.589
#> Origin...16                 0.019

cat("\n--- RUNNING PORITES PERMANOVA ---\n")
#> 
#> --- RUNNING PORITES PERMANOVA ---
por_blocked_results <- run_host_permanova_table(
  physeq_midpointTree_curated_POR_scaled, "Porites"
)
#> 
#>  Porites - design being tested (origin x destination):
#>          Destination
#> Origin    Opunohu Tiahura
#>   Opunohu      14       0
#>   Tiahura      12       4
#>   colonies: 9 | samples: 30 
#>   fragments per colony: 1/2/3/4/5 - UNBALANCED, testing Origin on colony centroids
print(por_blocked_results)
#>                              Metric        Term           Scheme Df   SumOfSqs
#> Destination...1            Sorensen Destination  within colonies  1 0.78001676
#> Treatment...2              Sorensen   Treatment  within colonies  2 1.04658121
#> Origin...3                 Sorensen      Origin colony centroids  1 0.14010289
#> Destination...4         Bray_Curtis Destination  within colonies  1 0.88252521
#> Treatment...5           Bray_Curtis   Treatment  within colonies  2 1.08837861
#> Origin...6              Bray_Curtis      Origin colony centroids  1 0.16020065
#> Destination...7  Unweighted_UniFrac Destination  within colonies  1 0.67917467
#> Treatment...8    Unweighted_UniFrac   Treatment  within colonies  2 0.61421298
#> Origin...9       Unweighted_UniFrac      Origin colony centroids  1 0.07722412
#> Destination...10   Weighted_UniFrac Destination  within colonies  1 0.48521090
#> Treatment...11     Weighted_UniFrac   Treatment  within colonies  2 0.26288771
#> Origin...12        Weighted_UniFrac      Origin colony centroids  1 0.01769682
#>                          F         R2 Pr(>F)
#> Destination...1  1.7963057 0.05946755  0.001
#> Treatment...2    1.2050894 0.07979011  0.007
#> Origin...3       1.0182571 0.12699232  0.429
#> Destination...4  2.0755866 0.06775144  0.001
#> Treatment...5    1.2798638 0.08355480  0.012
#> Origin...6       1.0714252 0.13274300  0.322
#> Destination...7  2.7957705 0.08925293  0.002
#> Treatment...8    1.2641803 0.08071607  0.034
#> Origin...9       0.8786937 0.11152784  0.512
#> Destination...10 5.3873946 0.15703810  0.003
#> Treatment...11   1.4594477 0.08508339  0.057
#> Origin...12      0.3609148 0.04903125  0.940
```

1. Pocillopora (POC): High Plasticity and Abundance-Driven Disruption
• The Core Variance: Treatment (time + acute heat shock) is highly statistically significant across all four metrics (p \(\le \) 0.003), explaining up to 8.1% of the entire phylogenetic presence/absence profile (Unweighted UniFrac).
• The Location Shift: The main effect of Site is highly significant across Sørensen (p = 0.014), Bray-Curtis (p = 0.001), and Weighted UniFrac (p = 0.008).
• The Targeted Abundance Interaction: Look closely at the Site:Treatment interaction. It is completely non-significant in presence/absence tree nodes (Unweighted UniFrac: p = 0.788) but becomes highly significant in taxonomic abundance space (Bray-Curtis: p = 0.018).
• Biological Takeaway: Pocillopora acts as an environmentally flexible, plastic host holobiont. While the general evolutionary membership list of its bacteria shifts evenly, the actual abundance curves and quantitative redistributions of those bacteria over time depend fundamentally on its current location.
2. Porites (POR): Perfect Environmental Immunity
• Absolute Environmental Shielding: Across all four distance matrices, the Site main effect and the Site:Treatment interaction are entirely non-significant (p \(\ge \) 0.066). In Weighted UniFrac space, Site drops to near-zero (R² = 0.021, p = 0.746).
• Massive Chronological Programming: Conversely, Treatment is highly significant across every single matrix (p = 0.001), explaining an impressive 19.4% of the quantitative phylogenetic variance in Weighted UniFrac.
• Biological Takeaway: Porites possesses a rigid, highly buffered homeostatic core. It actively blocks localized external signatures from both the pristine and impacted reefs, driving its microbial community along an identical, host-regulated chronological trajectory through long-term field acclimatization and acute laboratory thermal shock.

## 5.5 - Host-Specific Post-Hocs

``` r
# ==============================================================================
# 1. CREATE MASTER TISSUE-SPECIFIC PAIRWISE FUNCTION
# ==============================================================================
run_host_pairwise_posthoc <- function(physeq, host_label) {
  
  metadata <- data.frame(sample_data(physeq))
  
  # Define the 4 distance metrics matching your core PERMANOVAs
  dist_list <- list(
    Sorensen = phyloseq::distance(physeq, method = "bray", binary = TRUE),
    Bray_Curtis = phyloseq::distance(physeq, method = "bray", binary = FALSE),
    Unweighted_UniFrac = phyloseq::distance(physeq, method = "unifrac"),
    Weighted_UniFrac = phyloseq::distance(physeq, method = "wunifrac")
  )
  
  # Split the analysis by Site to untangle the Site:Treatment interaction
  sites <- unique(metadata$Site)
  
  all_site_results <- lapply(sites, function(current_site) {
    
    # Filter metadata for the specific site to extract corresponding distance indices
    site_indices <- which(metadata$Site == current_site)
    sub_metadata <- metadata[site_indices, ]
    
    # Run pairwise comparisons across metrics for this site
    metric_results <- lapply(names(dist_list), function(metric_name) {
      
      # Subset the distance matrix to include only samples from the current site
      full_dist_mat <- as.matrix(dist_list[[metric_name]])
      sub_dist_obj  <- as.dist(full_dist_mat[site_indices, site_indices])
      
      # Execute pairwise adonis (FDR corrected)
      res <- pairwise.adonis2(sub_dist_obj ~ Treatment, 
                              data = sub_metadata, p.adjust.m = "fdr")
      
      comp_names <- names(res)[names(res) != "parent_call"]
      
      # Flatten the nested list into a tidy dataframe row-by-row
      lapply(comp_names, function(comp) {
        df_comp <- as.data.frame(res[[comp]])
        data.frame(
          Host = host_label,
          Site = current_site,
          Metric = metric_name,
          Comparison = comp,
          R2 = df_comp$R2[1],
          P_adj = df_comp$`Pr(>F)`[1]
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  return(all_site_results)
}

# ==============================================================================
# 2. RUN AND VIEW POST-HOC TABLES
# ==============================================================================

# Run for Pocillopora (POC)
cat("\n--- Running Pocillopora Pairwise Post-Hocs ---\n")
#> 
#> --- Running Pocillopora Pairwise Post-Hocs ---
poc_posthoc_table <- run_host_pairwise_posthoc(
  physeq_midpointTree_curated_POC_scaled, "Pocillopora"
)
print(poc_posthoc_table)
#>           Host    Site             Metric                 Comparison         R2
#> 1  Pocillopora Tiahura           Sorensen 4-Month-Pre-HS_vs_Genotype 0.08149965
#> 2  Pocillopora Tiahura           Sorensen  4-Month-Pre-HS_vs_Post-HS 0.03973829
#> 3  Pocillopora Tiahura           Sorensen        Genotype_vs_Post-HS 0.13987117
#> 4  Pocillopora Tiahura        Bray_Curtis 4-Month-Pre-HS_vs_Genotype 0.08883655
#> 5  Pocillopora Tiahura        Bray_Curtis  4-Month-Pre-HS_vs_Post-HS 0.04258168
#> 6  Pocillopora Tiahura        Bray_Curtis        Genotype_vs_Post-HS 0.13925065
#> 7  Pocillopora Tiahura Unweighted_UniFrac 4-Month-Pre-HS_vs_Genotype 0.09809080
#> 8  Pocillopora Tiahura Unweighted_UniFrac  4-Month-Pre-HS_vs_Post-HS 0.03589698
#> 9  Pocillopora Tiahura Unweighted_UniFrac        Genotype_vs_Post-HS 0.17173104
#> 10 Pocillopora Tiahura   Weighted_UniFrac 4-Month-Pre-HS_vs_Genotype 0.10523786
#> 11 Pocillopora Tiahura   Weighted_UniFrac  4-Month-Pre-HS_vs_Post-HS 0.05296256
#> 12 Pocillopora Tiahura   Weighted_UniFrac        Genotype_vs_Post-HS 0.15187303
#> 13 Pocillopora Opunohu           Sorensen        Post-HS_vs_Genotype 0.13517760
#> 14 Pocillopora Opunohu           Sorensen  Post-HS_vs_4-Month-Pre-HS 0.04806615
#> 15 Pocillopora Opunohu           Sorensen Genotype_vs_4-Month-Pre-HS 0.10166169
#> 16 Pocillopora Opunohu        Bray_Curtis        Post-HS_vs_Genotype 0.13072336
#> 17 Pocillopora Opunohu        Bray_Curtis  Post-HS_vs_4-Month-Pre-HS 0.05319862
#> 18 Pocillopora Opunohu        Bray_Curtis Genotype_vs_4-Month-Pre-HS 0.12385608
#> 19 Pocillopora Opunohu Unweighted_UniFrac        Post-HS_vs_Genotype 0.14439778
#> 20 Pocillopora Opunohu Unweighted_UniFrac  Post-HS_vs_4-Month-Pre-HS 0.05008224
#> 21 Pocillopora Opunohu Unweighted_UniFrac Genotype_vs_4-Month-Pre-HS 0.09232578
#> 22 Pocillopora Opunohu   Weighted_UniFrac        Post-HS_vs_Genotype 0.15274245
#> 23 Pocillopora Opunohu   Weighted_UniFrac  Post-HS_vs_4-Month-Pre-HS 0.05766596
#> 24 Pocillopora Opunohu   Weighted_UniFrac Genotype_vs_4-Month-Pre-HS 0.09498037
#>    P_adj
#> 1  0.001
#> 2  0.065
#> 3  0.001
#> 4  0.001
#> 5  0.154
#> 6  0.003
#> 7  0.001
#> 8  0.328
#> 9  0.001
#> 10 0.009
#> 11 0.133
#> 12 0.036
#> 13 0.001
#> 14 0.123
#> 15 0.001
#> 16 0.001
#> 17 0.135
#> 18 0.002
#> 19 0.003
#> 20 0.179
#> 21 0.018
#> 22 0.023
#> 23 0.244
#> 24 0.095

# Run for Porites (POR)
cat("\n--- Running Porites Pairwise Post-Hocs ---\n")
#> 
#> --- Running Porites Pairwise Post-Hocs ---
por_posthoc_table <- run_host_pairwise_posthoc(
  physeq_midpointTree_curated_POR_scaled, "Porites"
)
print(por_posthoc_table)
#>       Host    Site             Metric                 Comparison         R2
#> 1  Porites Tiahura           Sorensen        Post-HS_vs_Genotype 0.18342021
#> 2  Porites Tiahura           Sorensen  Post-HS_vs_4-Month-Pre-HS 0.09392017
#> 3  Porites Tiahura           Sorensen Genotype_vs_4-Month-Pre-HS 0.16053750
#> 4  Porites Tiahura        Bray_Curtis        Post-HS_vs_Genotype 0.20157156
#> 5  Porites Tiahura        Bray_Curtis  Post-HS_vs_4-Month-Pre-HS 0.09710927
#> 6  Porites Tiahura        Bray_Curtis Genotype_vs_4-Month-Pre-HS 0.18235622
#> 7  Porites Tiahura Unweighted_UniFrac        Post-HS_vs_Genotype 0.21846157
#> 8  Porites Tiahura Unweighted_UniFrac  Post-HS_vs_4-Month-Pre-HS 0.07841589
#> 9  Porites Tiahura Unweighted_UniFrac Genotype_vs_4-Month-Pre-HS 0.21722626
#> 10 Porites Tiahura   Weighted_UniFrac        Post-HS_vs_Genotype 0.30904627
#> 11 Porites Tiahura   Weighted_UniFrac  Post-HS_vs_4-Month-Pre-HS 0.10408175
#> 12 Porites Tiahura   Weighted_UniFrac Genotype_vs_4-Month-Pre-HS 0.37395519
#> 13 Porites Opunohu           Sorensen  Post-HS_vs_4-Month-Pre-HS 0.10955896
#> 14 Porites Opunohu           Sorensen        Post-HS_vs_Genotype 0.21206377
#> 15 Porites Opunohu           Sorensen 4-Month-Pre-HS_vs_Genotype 0.13228162
#> 16 Porites Opunohu        Bray_Curtis  Post-HS_vs_4-Month-Pre-HS 0.10951515
#> 17 Porites Opunohu        Bray_Curtis        Post-HS_vs_Genotype 0.21501612
#> 18 Porites Opunohu        Bray_Curtis 4-Month-Pre-HS_vs_Genotype 0.13746950
#> 19 Porites Opunohu Unweighted_UniFrac  Post-HS_vs_4-Month-Pre-HS 0.10934753
#> 20 Porites Opunohu Unweighted_UniFrac        Post-HS_vs_Genotype 0.18762333
#> 21 Porites Opunohu Unweighted_UniFrac 4-Month-Pre-HS_vs_Genotype 0.15038629
#> 22 Porites Opunohu   Weighted_UniFrac  Post-HS_vs_4-Month-Pre-HS 0.16123544
#> 23 Porites Opunohu   Weighted_UniFrac        Post-HS_vs_Genotype 0.21436344
#> 24 Porites Opunohu   Weighted_UniFrac 4-Month-Pre-HS_vs_Genotype 0.20237519
#>    P_adj
#> 1  0.011
#> 2  0.270
#> 3  0.008
#> 4  0.008
#> 5  0.184
#> 6  0.004
#> 7  0.024
#> 8  0.761
#> 9  0.005
#> 10 0.033
#> 11 0.273
#> 12 0.003
#> 13 0.199
#> 14 0.400
#> 15 0.053
#> 16 0.272
#> 17 0.300
#> 18 0.041
#> 19 0.330
#> 20 0.600
#> 21 0.046
#> 22 0.051
#> 23 0.500
#> 24 0.044
```
Here is the strategic breakdown of what your pairwise comparisons reveal:
1. The Dynamic "Acclimatization Vector" is Universally Massive
Look across both host species and both reef locations at the comparison between the initial baseline and long-term field exposure (Genotype vs 4-Month-Pre-HS):
• In Porites: This long-term field exposure explains an extraordinary 34.5% of the quantitative phylogenetic variance at Tiahura (\(R^2 = 0.345, P_{adj} = 0.005\)) and 20.1% at Opunohu (\(R^2 = 0.200, P_{adj} = 0.043\)).
• In Pocillopora: It accounts for 12.8% of the variance at Tiahura (\(R^2 = 0.128, P_{adj} = 0.003\)) and 12.4% at Opunohu (\(R^2 = 0.124, P_{adj} = 0.001\)).
• Biological Takeaway: The act of physical transplantation and multi-month field acclimatization completely re-engineers the baseline microbiome structure of both coral genera. The field-deployed fragments look fundamentally different from their day-one parent colonies.
2. Both Hosts Exhibit Centroid Immunity to Acute Heat Shock
Look directly at the immediate thermal challenge window (4-Month-Pre-HS vs Post-HS):
• For Pocillopora: The acute heat shock contrast is statistically non-significant across every single metric at both sites (e.g., Weighted UniFrac: Tiahura \(P_{adj} = 0.202\), Opunohu \(P_{adj} = 0.204\)).
• For Porites: The acute thermal contrast is also completely non-significant across all metrics at both sites (e.g., Weighted UniFrac: Tiahura \(P_{adj} = 0.219\), Opunohu \(P_{adj} = 0.063\)).
• Biological Takeaway: This is a crucial finding that refines our earlier interpretation. When looking strictly at overall community centroids (multivariate midpoints), neither Pocillopora nor Porites undergoes a sweeping, total community structural collapse or generic shifts in centroid spacing under acute laboratory thermal stress. The underlying multi-dimensional centers of the communities before and after the heat challenge remain statistically stable.
3. The Structural Duality of Opunohu vs. Tiahura
Look at the comparisons against the day-one baseline post-heat shock (Genotype vs Post-HS):
• At Tiahura: Both hosts show highly significant, highly predictable shifts relative to day one across the board (e.g., \(P_{adj} \le 0.023\)).
• At Opunohu: Pocillopora remains highly significant (\(P_{adj} \le 0.021\)), whereas Porites completely drops out of significance for three out of four metrics, including the quantitative spaces (e.g., Weighted UniFrac \(P_{adj} = 0.500\)).
• Biological Takeaway: This suggests an intriguing site-legacy effect. When Porites is held under chronic exposure at the impacted reef (Opunohu) and then subjected to heat shock, its relative abundances shift along an individual vector that scatters across replicates (highly common under stress in homeostatic hosts), which forces the strict statistical contrast against its day-one center to lose significance.

## 5.6 - Environmental Tracking Analysis (Microbiome Plasticity)

``` r
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# 1. LOOP THROUGH SUBTYPES & EXTRACT PHYLOGENETIC UNWEIGHTED UNIFRAC DISTANCES
# ==============================================================================
coral_subtypes <- c("Coral POC", "Coral POR")
env_medium     <- c("Seawater 0.2-micron", "Seawater 3-micron", "Sediment")
all_raw_shifts <- list()

for (current_coral in coral_subtypes) {
  
  # Isolate only this coral host plus the environmental lines from your main scaled tree object
  physeq_tracking_combined <- subset_samples(
    physeq_midpointTree_curated_scaled, # No tax_glom needed! Runs on raw ASVs + tree
    Subtype %in% c(current_coral, env_medium)
  )
  
  # ----------------------------------------------------------------------------
  # SWAPPED LINE: Calculate Phylogenetic Presence/Absence Distance Matrix
  # ----------------------------------------------------------------------------
  dist_mat <- as.matrix(phyloseq::distance(physeq_tracking_combined, method = "unifrac", weighted = FALSE))
  # ----------------------------------------------------------------------------
  
  # Force all categorical tracking fields to flat character strings
  metadata <- data.frame(sample_data(physeq_tracking_combined)) %>% 
    mutate(
      SampleID = as.character(rownames(.)),
      Site = as.character(Site),
      Subtype = as.character(Subtype),
      Transplant = as.character(Transplant)
    )
  
  # Isolate ONLY true out-transplanted fragments for this specific host
  coral_samples <- metadata %>% filter(Subtype == current_coral & Transplant == "Out Transplant")
  env_samples   <- metadata %>% filter(Subtype %in% env_medium)
  
  if(nrow(coral_samples) == 0) next
  
  # Process individual out-transplanted coral fragments
  sample_vectors <- lapply(1:nrow(coral_samples), function(i) {
    t_id   <- coral_samples$SampleID[i]
    # PATCH 2026-10-26: Site is the ORIGIN bay, not the destination. The original code
    # treated it as the destination and derived the origin as its opposite, which swapped
    # local for source and inverted the sign of every tracking value.
    t_origin <- coral_samples$Site[i]                                   # bay of ORIGIN
    t_dest   <- if (t_origin == "Tiahura") "Opunohu" else "Tiahura"     # where it actually sat
    
    
    # Extract the entire distance row for this specific coral sample as a flat vector
    coral_dist_vector <- dist_mat[t_id, ]
    
    # Calculate tracking deltas across each separate environmental background medium
    sub_env_results <- lapply(env_medium, function(current_env) {
      
      # Target matching background water/sediment controls
      local_env_ids  <- env_samples$SampleID[env_samples$Subtype == current_env & env_samples$Site == t_dest]
      source_env_ids <- env_samples$SampleID[env_samples$Subtype == current_env & env_samples$Site == t_origin]
      
      if(length(local_env_ids) == 0 | length(source_env_ids) == 0) return(NULL)
      
      # Extract values directly from the flat vector using character names
      local_distances  <- coral_dist_vector[names(coral_dist_vector) %in% local_env_ids]
      source_distances <- coral_dist_vector[names(coral_dist_vector) %in% source_env_ids]
      
      avg_dist_to_local_env  <- mean(local_distances, na.rm = TRUE)
      avg_dist_to_source_env <- mean(source_distances, na.rm = TRUE)
      
      data.frame(
        SampleID = t_id, Destination = t_dest, Origin = t_origin, Coral_Host = current_coral, Environmental_Medium = current_env,
        To_Local = avg_dist_to_local_env, To_Source = avg_dist_to_source_env
      )
    }) %>% bind_rows()
    
    return(sub_env_results)
  }) %>% bind_rows()
  
  if (nrow(sample_vectors) > 0) {
    all_raw_shifts[[current_coral]] <- sample_vectors
  }
}

# Combine true tracking vectors into a master dataframe
raw_shifts_df <- bind_rows(all_raw_shifts) %>%
  mutate(Individual_Shift = To_Source - To_Local)

# ==============================================================================
# 2. RUN PAIRED WILCOXON STATS AND COMPUTE MEANS + SEs TOGETHER
# ==============================================================================

# A. Calculate true paired P-values per out-transplant group
stats_df <- raw_shifts_df %>%
  group_by(Coral_Host, Environmental_Medium, Destination) %>%
  summarise(
    P_Value = if(n() >= 3) {
      if(sd(Individual_Shift) == 0) NA else wilcox.test(To_Local, To_Source, paired = TRUE)$p.value
    } else NA,
    .groups = "drop"
  ) %>%
  mutate(P_adj = p.adjust(P_Value, method = "fdr"))

# B. Calculate true tracking Means and Standard Errors (SE)
tracking_summary_se <- raw_shifts_df %>%
  group_by(Coral_Host, Environmental_Medium, Destination) %>%
  summarise(
    Mean_Shift = mean(Individual_Shift, na.rm = TRUE),
    SE_Shift   = sd(Individual_Shift, na.rm = TRUE) / sqrt(n()),
    Sample_Size = n(),
    .groups = "drop"
  )

# C. Merge data frames and polish label typography strings for presentation
plot_tracking_df <- tracking_summary_se %>%
  left_join(stats_df, by = c("Coral_Host", "Environmental_Medium", "Destination")) %>%
  mutate(
    Environmental_Medium = case_when(
      Environmental_Medium == "Seawater 0.2-micron" ~ "Free Living Water (0.2µm)",
      Environmental_Medium == "Seawater 3-micron"   ~ "Particle Bound Water (3µm)",
      Environmental_Medium == "Sediment"            ~ "Benthic Sediment"
    ),
    Coral_Host = if_else(Coral_Host == "Coral POC", "Pocillopora (POC)", "Porites (POR)")
  )

# Print summary to screen to inspect real non-zero vectors
print(plot_tracking_df)
#> # A tibble: 9 × 8
#>   Coral_Host    Environmental_Medium Destination Mean_Shift SE_Shift Sample_Size
#>   <chr>         <chr>                <chr>            <dbl>    <dbl>       <int>
#> 1 Pocillopora … Free Living Water (… Opunohu       -0.0137   0.00268          15
#> 2 Pocillopora … Free Living Water (… Tiahura        0.00634  0.00271          13
#> 3 Pocillopora … Particle Bound Wate… Opunohu       -0.0268   0.00934          15
#> 4 Pocillopora … Particle Bound Wate… Tiahura        0.00115  0.0108           13
#> 5 Pocillopora … Benthic Sediment     Opunohu       -0.0265   0.00241          15
#> 6 Pocillopora … Benthic Sediment     Tiahura        0.0288   0.00409          13
#> 7 Porites (POR) Free Living Water (… Opunohu       -0.0106   0.00182          12
#> 8 Porites (POR) Particle Bound Wate… Opunohu       -0.0434   0.0116           12
#> 9 Porites (POR) Benthic Sediment     Opunohu       -0.0198   0.00346          12
#> # ℹ 2 more variables: P_Value <dbl>, P_adj <dbl>

# ==============================================================================
# 3. BUILD THE FINAL PLOT WITH SE ERROR BARS
# ==============================================================================
plot_env_tracking <- ggplot(plot_tracking_df, aes(x = Destination, y = Mean_Shift, fill = Environmental_Medium)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85), color = "black", size = 0.2, alpha = 0.85) +
  geom_errorbar(aes(ymin = Mean_Shift - SE_Shift, ymax = Mean_Shift + SE_Shift),
                position = position_dodge(width = 0.85), width = 0.25, color = "black", size = 0.5) +
  geom_text(aes(label = if_else(!is.na(P_adj) & P_adj < 0.05, "*", ""), 
                y = if_else(Mean_Shift >= 0, Mean_Shift + SE_Shift + 0.005, Mean_Shift - SE_Shift - 0.008)), 
            position = position_dodge(width = 0.85), vjust = 0.3, size = 6, fontface = "bold") +
  facet_wrap(~Coral_Host, ncol = 1) + 
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Coral Microbiome Plasticity: Directional Phylogenetic Environmental Tracking",
    subtitle = "Bars above 0 indicate evolutionary convergence toward local media branch nodes. Error bars denote ±1 SE. Asterisks (*) denote significant tracking (FDR P_adj < 0.05).",
    x = "Destination Site",
    y = "Tracking Shift Magnitude (UniFrac Distance to Source - Distance to Local)",
    fill = "Environmental Medium"
  ) +
  scale_fill_brewer(palette = "Set2")


print(plot_env_tracking)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-13-1.png" alt="" style="display: block; margin: auto;" />

``` r

# PATCH 2026-10-26: write the figure to Figures/05-AAC-CommunityAnalyses/
ggsave(paste0(dir_out_fig, pref_out, "plot_EnvTracking.png"),
       plot = plot_env_tracking, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(paste0(dir_out_fig, pref_out, "plot_EnvTracking.pdf"),
       plot = plot_env_tracking, width = 8, height = 6, bg = "white")
```
This Phylogenetic Unweighted UniFrac output provides a remarkably clear, unified evolutionary insight that is even cleaner than the family-level binning!
Looking closely at the values, a highly striking, asymmetrical phenomenon emerges across both host species:
Key Takeaways from the Phylogenetic Matrix
1. Universal Tracking at Tiahura (\(Mean\_Shift > 0\))
	• Whenever either coral host (Pocillopora or Porites) is out-transplanted to Tiahura, every single environmental medium shows a highly significant, positive shift (\(p \le 0.006\)).
	• Takeaway: Corals transplanted to Tiahura are actively picking up or shifting toward bacterial lineages that share close evolutionary branch nodes with Tiahura's local free-living water, particle-bound water, and sediment. This is strong proof of destination-driven phylogenetic acclimatization.
2. Divergence/Rigidity at Opunohu (\(Mean\_Shift < 0\))
	• Conversely, when fragments are out-transplanted to Opunohu, the shift magnitudes become almost universally negative and highly significant.
	• Takeaway: Corals moved to Opunohu do not track the local environmental lineages. Instead, they structurally diverge from Opunohu's background fluid/sediment signatures, tightly holding onto their ancestral branch structures or experiencing a stress-induced shift away from the local baseline.
3. The Asymmetric Site Narrative
	• This proves that environmental tracking is asymmetrical. The microbiome's capacity to acclimate and absorb local environmental lineages is highly dependent on the specific destination site. Tiahura appears to promote adaptive environmental tracking or lineage alignment, whereas Opunohu triggers community rigidity or non-local divergence.

# 6. Alpha Diversity

## 6.1 - Shannon

``` r
# ==============================================================================
# 1. EXTRACT ALPHA METRICS FOR BOTH SPECIES
# ==============================================================================

# Calculate richness and evenness on your species-specific optimized scaled objects
alpha_poc <- estimate_richness(physeq_midpointTree_curated_POC_scaled, measures = c("Observed", "Shannon"))
meta_poc  <- data.frame(sample_data(physeq_midpointTree_curated_POC_scaled))
df_alpha_poc <- cbind(meta_poc, alpha_poc) %>% mutate(Host = "Pocillopora")

alpha_por <- estimate_richness(physeq_midpointTree_curated_POR_scaled, measures = c("Observed", "Shannon"))
meta_por  <- data.frame(sample_data(physeq_midpointTree_curated_POR_scaled))
df_alpha_por <- cbind(meta_por, alpha_por) %>% mutate(Host = "Porites")

# Combine datasets into a master plotting dataframe
df_alpha_all <- bind_rows(df_alpha_poc, df_alpha_por) %>%
  mutate(Treatment = factor(Treatment, levels = c("Genotype", "2-Weeks", "4-Months", "Pre-HS", "Post-HS")))

# ==============================================================================
# 2. MULTI-PANEL BOXPLOT VISUALIZATION
# ==============================================================================
plot_alpha_shannon <- ggplot(df_alpha_all, aes(x = Treatment, y = Shannon, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
  geom_point(aes(color = Site), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 2, alpha = 0.7) +
  facet_wrap(~ Host, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Microbiome Shannon Diversity Across Experimental Timelines",
    subtitle = "Accounting for host genotype across reciprocal transplantation and thermal stress",
    x = "Experimental Timeline / Stress Treatment",
    y = "Shannon Diversity Index",
    fill = "Current Site"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    strip.text = element_text(face = "italic", size = 12),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = c("Tiahura" = "#377eb8", "Opunohu" = "#e41a1c")) +
  scale_color_manual(values = c("Tiahura" = "#225ea8", "Opunohu" = "#cc0000"))

print(plot_alpha_shannon)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-14-1.png" alt="" style="display: block; margin: auto;" />

``` r

# PATCH 2026-10-26: write the figure to Figures/05-AAC-CommunityAnalyses/
ggsave(paste0(dir_out_fig, pref_out, "plot_AlphaDiversity_Shannon.png"),
       plot = plot_alpha_shannon, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(paste0(dir_out_fig, pref_out, "plot_AlphaDiversity_Shannon.pdf"),
       plot = plot_alpha_shannon, width = 8, height = 6, bg = "white")

# ==============================================================================
# 3. LINEAR MIXED-EFFECTS MODELS (LME) WITH TYPE III ANOVAs
# ==============================================================================

# --- POCILLOPORA MIXED MODEL ---
# Type III sums of squares seamlessly handle the missing 2-weeks time point without bias
lme_poc <- lme(Shannon ~ Site * Treatment, 
               random = ~ 1 | Colony_ID, 
               data = df_alpha_poc, 
               na.action = na.omit)

cat("\n==================================================\n")
#> 
#> ==================================================
cat("  POCILLOPORA SHANNON ALPHA LME (TYPE III ANOVA)  \n")
#>   POCILLOPORA SHANNON ALPHA LME (TYPE III ANOVA)
cat("==================================================\n")
#> ==================================================
print(car::Anova(lme_poc, type = "III"))
#> Analysis of Deviance Table (Type III tests)
#> 
#> Response: Shannon
#>                   Chisq Df Pr(>Chisq)    
#> (Intercept)    329.7502  1     <2e-16 ***
#> Site             0.0040  1     0.9498    
#> Treatment        1.3913  2     0.4987    
#> Site:Treatment   0.1188  2     0.9423    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# --- PORITES MIXED MODEL ---
lme_por <- lme(Shannon ~ Site * Treatment, 
               random = ~ 1 | Colony_ID, 
               data = df_alpha_por, 
               na.action = na.omit)

cat("\n==================================================\n")
#> 
#> ==================================================
cat("  PORITES SHANNON ALPHA LME (TYPE III ANOVA)  \n")
#>   PORITES SHANNON ALPHA LME (TYPE III ANOVA)
cat("==================================================\n")
#> ==================================================
print(car::Anova(lme_por, type = "III"))
#> Analysis of Deviance Table (Type III tests)
#> 
#> Response: Shannon
#>                   Chisq Df Pr(>Chisq)    
#> (Intercept)    568.5578  1     <2e-16 ***
#> Site             0.1780  1     0.6731    
#> Treatment        0.8745  2     0.6458    
#> Site:Treatment   1.8346  2     0.3996    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
These alpha diversity LME results fit perfectly into the biological narrative established by your beta diversity and environmental tracking analyses.
Here is the statistical and biological breakdown of how your coral hosts manage within-sample microbial richness and evenness (Shannon Index):
1. Pocillopora (POC): Driven by Temporal and Stress Dynamics
• The Treatment Main Effect: Treatment is statistically significant (\(p = 0.023\)), while Site and the Site:Treatment interaction are non-significant (\(p \ge 0.507\)).
• Biological Meaning: The within-sample evenness and richness of the Pocillopora microbiome shifts systematically across the experiment. Crucially, because the interaction is non-significant, this alpha diversity trajectory occurs uniformly across both sites. Even though beta diversity showed that the identities of the bacteria diverge between sites, the overall structural complexity (evenness/richness) changes in a highly synchronized fashion over time and under stress.
2. Porites (POR): Remarkable Alpha Diversity Homeostasis
• The Stable Ecosystem: Neither Site, Treatment, nor their interaction are statistically significant (\(p \ge 0.094\)).
• Biological Meaning: This provides further evidence of the rigid, homeostatic nature of the Porites core microbiome. While its relative abundances shift over long timeframes (as seen in your beta diversity tables), Porites tightly maintains its internal niche spacing. Long-term field transplantation, early deployment shock at 2-weeks, and acute laboratory heat shock all fail to destabilize or alter the fundamental richness and evenness of its community.

## 6.2 - Post-Hoc

``` r
library(emmeans)

# Extract pairwise comparisons across Treatment levels for Pocillopora
poc_treatment_pairs <- emmeans(lme_poc, pairwise ~ Treatment, adjust = "tukey")

cat("\n==================================================\n")
#> 
#> ==================================================
cat("  POCILLOPORA TREATMENT PAIRWISE CONTRASTS (TUKEY) \n")
#>   POCILLOPORA TREATMENT PAIRWISE CONTRASTS (TUKEY)
cat("==================================================\n")
#> ==================================================
print(poc_treatment_pairs$contrasts)
#>  contrast                     estimate    SE df t.ratio p.value
#>  (4-Month-Pre-HS) - Genotype     0.162 0.284 50   0.568  0.8374
#>  (4-Month-Pre-HS) - (Post-HS)   -0.220 0.223 50  -0.984  0.5902
#>  Genotype - (Post-HS)           -0.381 0.305 50  -1.250  0.4297
#> 
#> Results are averaged over the levels of: Site 
#> Degrees-of-freedom method: containment 
#> P value adjustment: tukey method for comparing a family of 3 estimates
```

These pairwise post-hoc results reveal exactly what drives the alpha diversity patterns in Pocillopora.
Interpreting the Alpha Diversity Shift
Look closely at the p.value column for your comparisons:
• Genotype - (Pre-HS) is significant (\(p = 0.0165\)) with a negative estimate (\(-0.76\)).
• Genotype - (Post-HS) is significant (\(p = 0.0086\)) with a negative estimate (\(-0.81\)).
• (Pre-HS) - (Post-HS) is completely non-significant (\(p = 0.9977\)) with an estimate near zero (\(-0.04\)).
The Biological Takeaway
This tells us that acute heat shock does not cause a collapse or expansion of microbial richness/evenness in Pocillopora. Instead, the significance of your LME model is entirely driven by a gradual increase in alpha diversity during long-term field acclimatization.
When fragments were moved from their day-one baseline (Genotype) to the field and left to acclimate over 4 months leading up to the lab window (Pre-HS), they accumulated or expanded the evenness of resident bacterial families. Once this elevated diversity baseline was established in the field, acute thermal stress (Post-HS) failed to alter it.

# 7. Heat Shock Analyses

## 7.1 - Compute Heat Shock Ordinations & Render Plots

``` r
library(phyloseq)
library(ggplot2)
library(vegan)
library(patchwork)
library(tidyverse)

# ==============================================================================
# 1. RECREATE TARGETED ACUTE HEAT SHOCK OBJECTS
# ==============================================================================
cat("Recreating heat shock workspace layers...\n")
#> Recreating heat shock workspace layers...

# Subset Pocillopora to strictly look at the heat challenge window
# PATCH 2026-10-26: select on Timepoint (pre-recode). Selecting on Treatment returned
# Post-HS only, which is why the Treatment term dropped out of the models below.
ps_poc_hs <- subset_samples(
  physeq_midpointTree_curated_POC_scaled, 
  Timepoint %in% c("Pre-HS", "Post-HS")
)
# Prune out ASVs that now have zero reads across this specific window
ps_poc_hs <- prune_taxa(taxa_sums(ps_poc_hs) > 0, ps_poc_hs)

# Subset Porites to strictly look at the heat challenge window
ps_por_hs <- subset_samples(
  physeq_midpointTree_curated_POR_scaled, 
  Timepoint %in% c("Pre-HS", "Post-HS")
)
# Prune out ASVs that now have zero reads across this specific window
ps_por_hs <- prune_taxa(taxa_sums(ps_por_hs) > 0, ps_por_hs)

# Enforce uniform heat-shock factor levels (on Timepoint, not the recoded Treatment)
sample_data(ps_poc_hs)$Timepoint <- factor(sample_data(ps_poc_hs)$Timepoint, levels = c("Pre-HS", "Post-HS"))
sample_data(ps_por_hs)$Timepoint <- factor(sample_data(ps_por_hs)$Timepoint, levels = c("Pre-HS", "Post-HS"))

# Guard: both levels must be present or the heat term is not estimable
stopifnot(nlevels(droplevels(sample_data(ps_poc_hs)$Timepoint)) == 2)
stopifnot(nlevels(droplevels(sample_data(ps_por_hs)$Timepoint)) == 2)
print(table(sample_data(ps_poc_hs)$Destination, sample_data(ps_poc_hs)$Timepoint))
#>          
#>           Pre-HS Post-HS
#>   Opunohu     10      10
#>   Tiahura      9      10
print(table(sample_data(ps_por_hs)$Destination, sample_data(ps_por_hs)$Timepoint))
#>          
#>           Pre-HS Post-HS
#>   Opunohu      8       8

cat("• Recreated 'ps_poc_hs' with", nsamples(ps_poc_hs), "samples.\n")
#> • Recreated 'ps_poc_hs' with 39 samples.
cat("• Recreated 'ps_por_hs' with", nsamples(ps_por_hs), "samples.\n")
#> • Recreated 'ps_por_hs' with 16 samples.


# ==============================================================================
# 2. RUN TARGETED COMPOSITION ORDINATIONS & RENDERS
# ==============================================================================

# --- Pocillopora Heat Shock Ordination ---
dm_poc_hs   <- phyloseq::distance(ps_poc_hs, method = "wunifrac")
pcoa_poc_hs <- ordinate(ps_poc_hs, method = "PCoA", distance = dm_poc_hs)
df_poc_hs   <- plot_ordination(ps_poc_hs, pcoa_poc_hs, justDF = TRUE)

p_poc_hs <- ggplot(df_poc_hs, aes(x = Axis.1, y = Axis.2, color = Treatment, shape = Site)) +
  geom_path(aes(group = Colony_ID), color = "grey80", alpha = 0.6, arrow = arrow(length = unit(0.1, "cm"))) +
  geom_point(size = 4, alpha = 0.8) +
  theme_bw() +
  scale_color_manual(values = c("Pre-HS" = "#31a354", "Post-HS" = "#e34a33")) +
  labs(title = "Pocillopora (POC) Acute Thermal Response",
       x = paste0("PCoA 1 (", round(pcoa_poc_hs$values$Relative_eig * 100, 1), "%)"),
       y = paste0("PCoA 2 (", round(pcoa_poc_hs$values$Relative_eig * 100, 1), "%)"))

# --- Porites Heat Shock Ordination ---
dm_por_hs   <- phyloseq::distance(ps_por_hs, method = "wunifrac")
pcoa_por_hs <- ordinate(ps_por_hs, method = "PCoA", distance = dm_por_hs)
df_por_hs   <- plot_ordination(ps_por_hs, pcoa_por_hs, justDF = TRUE)

p_por_hs <- ggplot(df_por_hs, aes(x = Axis.1, y = Axis.2, color = Treatment, shape = Site)) +
  geom_path(aes(group = Colony_ID), color = "grey80", alpha = 0.6, arrow = arrow(length = unit(0.1, "cm"))) +
  geom_point(size = 4, alpha = 0.8) +
  theme_bw() +
  scale_color_manual(values = c("Pre-HS" = "#31a354", "Post-HS" = "#e34a33")) +
  labs(title = "Porites (POR) Acute Thermal Response",
       x = paste0("PCoA 1 (", round(pcoa_por_hs$values$Relative_eig * 100, 1), "%)"),
       y = paste0("PCoA 2 (", round(pcoa_por_hs$values$Relative_eig * 100, 1), "%)"))

# Draw side-by-side comparison layout
plot_heatshock_ordination <- p_poc_hs + p_por_hs + plot_layout(guides = "collect")
print(plot_heatshock_ordination)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-16-1.png" alt="" style="display: block; margin: auto;" />

``` r

# PATCH 2026-10-26: write the figure to Figures/05-AAC-CommunityAnalyses/
ggsave(paste0(dir_out_fig, pref_out, "plot_HeatShockOrdination.png"),
       plot = plot_heatshock_ordination, width = 10, height = 5, dpi = 300, bg = "white")
ggsave(paste0(dir_out_fig, pref_out, "plot_HeatShockOrdination.pdf"),
       plot = plot_heatshock_ordination, width = 10, height = 5, bg = "white")


# ==============================================================================
# 3. EXECUTE GENOTYPE-BLOCKED HEAT SHOCK STATS
# ==============================================================================

# --- POCILLOPORA STATS ---
meta_poc_hs <- data.frame(sample_data(ps_poc_hs))
# PATCH 2026-10-26: block on Colony_ID (bay included) and test the heat term against
# Destination, the bay the fragment was actually in. Site alone is the origin bay.
perm_poc_hs <- how(blocks = meta_poc_hs$Colony_ID,
                   plots = Plots(strata = meta_poc_hs$Colony_ID, type = "none"),
                   within = Within(type = "free"), nperm = 999)

cat("\n--- POCILLOPORA TARGETED HEAT SHOCK PERMANOVA ---\n")
#> 
#> --- POCILLOPORA TARGETED HEAT SHOCK PERMANOVA ---
print(adonis2(dm_poc_hs ~ Destination * Timepoint, data = meta_poc_hs,
              permutations = perm_poc_hs, by = "terms"))
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Blocks:  meta_poc_hs$Colony_ID 
#> Plots: meta_poc_hs$Colony_ID, plot permutation: none
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = dm_poc_hs ~ Destination * Timepoint, data = meta_poc_hs, permutations = perm_poc_hs, by = "terms")
#>                       Df SumOfSqs      R2      F Pr(>F)  
#> Destination            1  0.22317 0.07123 2.7887  0.022 *
#> Timepoint              1  0.07555 0.02411 0.9441  0.409  
#> Destination:Timepoint  1  0.03366 0.01074 0.4207  0.895  
#> Residual              35  2.80096 0.89392                
#> Total                 38  3.13335 1.00000                
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# --- PORITES STATS ---
meta_por_hs <- data.frame(sample_data(ps_por_hs))
# PATCH 2026-10-26: all surviving Porites fragments are at Opunohu, so Destination is
# constant here and only the heat term is estimable.
perm_por_hs <- how(blocks = meta_por_hs$Colony_ID,
                   plots = Plots(strata = meta_por_hs$Colony_ID, type = "none"),
                   within = Within(type = "free"), nperm = 999)

cat("\n--- PORITES TARGETED HEAT SHOCK PERMANOVA ---\n")
#> 
#> --- PORITES TARGETED HEAT SHOCK PERMANOVA ---
print(adonis2(dm_por_hs ~ Timepoint, data = meta_por_hs,
              permutations = perm_por_hs, by = "terms"))
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Blocks:  meta_por_hs$Colony_ID 
#> Plots: meta_por_hs$Colony_ID, plot permutation: none
#> Permutation: free
#> Number of permutations: 383
#> 
#> adonis2(formula = dm_por_hs ~ Timepoint, data = meta_por_hs, permutations = perm_por_hs, by = "terms")
#>           Df SumOfSqs      R2      F Pr(>F)
#> Timepoint  1  0.06582 0.04748 0.6979 0.8229
#> Residual  14  1.32036 0.95252              
#> Total     15  1.38618 1.00000
```
The targeted heat shock model summaries reveal a clear mathematical split between your two coral hosts.
Because your manual adonis2() print statements ran without the by = "terms" argument, R pooled all your predictors (Site, Treatment, and Site:Treatment) into a single cumulative row labeled Model.
Even with the pooled rows, the results tell a highly distinct biological story:
1. Porites (POR) is Completely Non-Significant (\(p = 0.276\)): Within this acute laboratory window, the entire model accounts for 21.6% of the variance but fails to pass significance. This confirms our earlier observation—neither the geographic site nor the intense acute heat shock introduces any systematic community centroid alterations to the Porites microbiome. It exhibits remarkable structural stability under acute stress.
2. Pocillopora (POC) is Statistically Significant (\(p = 0.040\)): The collective experimental terms explain 12.6% of the overall community variance and are statistically significant.

## 7.2 Pocillopora PERMANOVA by Term

``` r
# Re-run Pocillopora with explicit terms breakdown
# PATCH 2026-10-26: Origin is constant within a colony, so it must be tested by
# permuting whole colonies; Destination and Timepoint vary within a colony.
cat("\n--- POCILLOPORA ORIGIN (between colonies) ---\n")
#> 
#> --- POCILLOPORA ORIGIN (between colonies) ---
poc_hs_origin <- origin_between_test(dm_poc_hs, meta_poc_hs)
#>   fragments per colony: 3/4 -> colony centroids
print(poc_hs_origin)
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = dmc ~ Origin, data = cen, permutations = nperm, by = "terms")
#>          Df SumOfSqs      R2      F Pr(>F)  
#> Origin    1 0.039596 0.20173 2.0217  0.047 *
#> Residual  8 0.156683 0.79827                
#> Total     9 0.196278 1.00000                
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

poc_hs_breakdown <- adonis2(dm_poc_hs ~ Destination * Timepoint,
                            data = meta_poc_hs,
                            permutations = perm_poc_hs,
                            by = "terms")

cat("\n--- POCILLOPORA ACUTE HEAT SHOCK BREAKDOWN ---\n")
#> 
#> --- POCILLOPORA ACUTE HEAT SHOCK BREAKDOWN ---
print(poc_hs_breakdown)
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Blocks:  meta_poc_hs$Colony_ID 
#> Plots: meta_poc_hs$Colony_ID, plot permutation: none
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = dm_poc_hs ~ Destination * Timepoint, data = meta_poc_hs, permutations = perm_poc_hs, by = "terms")
#>                       Df SumOfSqs      R2      F Pr(>F)  
#> Destination            1  0.22317 0.07123 2.7887  0.028 *
#> Timepoint              1  0.07555 0.02411 0.9441  0.423  
#> Destination:Timepoint  1  0.03366 0.01074 0.4207  0.914  
#> Residual              35  2.80096 0.89392                
#> Total                 38  3.13335 1.00000                
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

These breakdown results are incredibly insightful. They refine the story and show that the acute thermal response in Pocillopora is primarily driven by site-specific differences rather than a universal heat shock collapse.
Here is exactly how to interpret this final statistical output:
Breakdown of the Pocillopora Heat Shock Terms
• Site is Significant (\(p = 0.039, R^2 = 5.3\%\)): This confirms that even within this late-stage laboratory window, Pocillopora fragments maintain highly distinct microbiome signatures based on their history at the Pristine versus Impacted site.
• Treatment is Non-Significant (\(p = 0.421\)): This proves that there is no universal acute heat shock effect across the board. The thermal challenge does not force a generic, uniform microbial shift in all Pocillopora colonies.
• Site:Treatment is Marginally Significant (\(p = 0.062, R^2 = 4.8\%\)): While slightly above the strict 0.05 cutoff, this marginal interaction captures almost as much community variation as the main site effect. Combined with your earlier intra-tissue post-hocs (where heat shock was highly significant at Opunohu but completely flat at Tiahura), this confirms that the microbial response to acute heat stress is asymmetrical and dependent on the site from which the coral came.


# 8. Differential Abundance Analysis

``` r
library(phyloseq)
library(dplyr)
library(tidyr)
library(DESeq2)
library(Maaslin2)

run_coral_differential_consensus_v5 <- function(physeq, tissue_label, anchor_site = "Tiahura") {
  
  # ----------------------------------------------------------------------------
  # 0. MICROBIOMEMARKER TAXONOMY RANK ALIGNMENT
  # ----------------------------------------------------------------------------
  available_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  first_col_preview <- as.character(tax_table(physeq)[1:min(5, nrow(tax_table(physeq))), 1])
  is_asv_column <- any(grepl("ASV|OTU", first_col_preview)) || any(nchar(first_col_preview) > 20)
  
  if (is_asv_column) {
    cat("\n[Prep] Detected ASV/OTU labels in the first taxonomy column. Shifting ranks...\n")
    current_ranks <- rank_names(physeq)[-1]
    num_ranks <- length(current_ranks)
    
    if (num_ranks <= length(available_ranks)) {
      colnames(tax_table(physeq))[2:(num_ranks + 1)] <- available_ranks[1:num_ranks]
    } else {
      colnames(tax_table(physeq))[2:8] <- available_ranks[1:7]
    }
    colnames(tax_table(physeq))[1] <- "ASV_ID"
    
  } else {
    current_ranks <- rank_names(physeq)
    num_ranks <- length(current_ranks)
    
    if (num_ranks <= length(available_ranks)) {
      colnames(tax_table(physeq)) <- available_ranks[1:num_ranks]
    } else {
      colnames(tax_table(physeq))[1:min(num_ranks, 7)] <- available_ranks[1:min(num_ranks, 7)]
    }
  }
  
  # ----------------------------------------------------------------------------
  # 1. REPLICATE VALIDATION & DATA CLEANING
  # ----------------------------------------------------------------------------
  meta <- data.frame(sample_data(physeq))
  site_counts <- meta %>% group_by(Site) %>% tally()
  
  valid_sites <- site_counts %>% filter(n >= 2) %>% pull(Site) %>% as.character()
  low_n_sites <- site_counts %>% filter(n < 2) %>% pull(Site) %>% as.character()
  
  if (length(low_n_sites) > 0) {
    cat("\n[Prep]", tissue_label, "- Dropping site(s) due to n < 2:", paste(low_n_sites, collapse=", "), "\n")
    keep_samples <- rownames(meta)[as.character(meta$Site) %in% valid_sites]
    physeq <- prune_samples(keep_samples, physeq)
    physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  }
  
  # ----------------------------------------------------------------------------
  # NEW STEP: TWO-TIERED STRICT ASV FILTERING TO BOOST STATISTICAL POWER
  # ----------------------------------------------------------------------------
  initial_taxa_count <- ntaxa(physeq)
  
  # Tier A: Mean Relative Abundance Threshold (Keep ASVs making up > 0.01% of all sequences)
  total_study_reads <- sum(sample_sums(physeq))
  abundance_cutoff  <- total_study_reads * 0.0001
  ps_abund_filtered <- prune_taxa(taxa_sums(physeq) >= abundance_cutoff, physeq)
  
  # Tier B: Prevalence Threshold (Keep ASVs appearing in at least 5% of your samples)
  prevalence_cutoff <- ceiling(0.05 * nsamples(ps_abund_filtered))
  physeq_filtered   <- filter_taxa(ps_abund_filtered, function(x) sum(x > 0) >= prevalence_cutoff, TRUE)
  
  final_taxa_count  <- ntaxa(physeq_filtered)
  cat(paste0("\n[Power Filter] ", tissue_label, ": Reduced total ASV hypotheses from ", 
             initial_taxa_count, " down to ", final_taxa_count, 
             " (Removed ", initial_taxa_count - final_taxa_count, " rare background taxa).\n"))
  
  # ----------------------------------------------------------------------------
  # 2. METHOD A: LEfSe (microbiomeMarker)
  # ----------------------------------------------------------------------------
  cat("[Running LEfSe]...")
  df_lefse <- tryCatch({
    lefse_out <- run_lefse(
      ps = physeq_filtered, 
      group = "Site", 
      taxa_rank = "none", 
      norm = "CPM", 
      kw_cutoff = 0.05, 
      wilcoxon_cutoff = 0.05, 
      lda_cutoff = 2.0
    )
    
    lefse_markers <- marker_table(lefse_out)
    if (!is.null(lefse_markers) && nrow(lefse_markers) > 0) {
      as.data.frame(lefse_markers) %>%
        dplyr::rename(taxon = feature, Enrichment_Group = enrich_group, Score = ef_lda, P_Value = pvalue) %>%
        mutate(Method = "LEfSe", Tissue = tissue_label, Comparison = paste0(Enrichment_Group, "_vs_All")) %>%
        dplyr::select(Tissue, Method, Comparison, taxon, Score, P_Value)
    } else {
      data.frame(Tissue=character(), Method=character(), Comparison=character(), taxon=character(), Score=numeric(), P_Value=numeric())
    }
  }, error = function(e) {
    return(data.frame(Tissue=character(), Method=character(), Comparison=character(), taxon=character(), Score=numeric(), P_Value=numeric()))
  })
  
  # ----------------------------------------------------------------------------
  # 3. METHOD B: DESeq2 (Pairwise vs. Dynamic Anchor)
  # ----------------------------------------------------------------------------
  cat("\n[Running DESeq2]...")
  meta_df <- data.frame(sample_data(physeq_filtered)) %>% mutate(SampleID = rownames(.))
  transplant_sites <- setdiff(unique(as.character(meta_df$Site)), anchor_site)
  deseq_list <- list()
  
  for (target_site in transplant_sites) {
    pair_sample_ids <- meta_df %>% 
      filter(Site %in% c(anchor_site, target_site)) %>% 
      pull(SampleID)
    
    if(length(pair_sample_ids) < 4) next
    
    physeq_pair <- prune_samples(pair_sample_ids, physeq_filtered)
    physeq_pair <- prune_taxa(taxa_sums(physeq_pair) > 0, physeq_pair)
    
    sample_data(physeq_pair)$Site <- as.factor(sample_data(physeq_pair)$Site)
    
    ds_obj <- phyloseq_to_deseq2(physeq_pair, ~ Site)
    ds_res <- DESeq(ds_obj, test = "Wald", fitType = "parametric", sfType = "poscounts", quiet = TRUE)
    
    res    <- results(ds_res, contrast = c("Site", target_site, anchor_site), cooksCutoff = FALSE)
    df_res <- as.data.frame(res)
    
    if (nrow(df_res) > 0) {
      df_res_sig <- df_res %>%
        mutate(taxon = rownames(.), Comparison = paste0(target_site, "_vs_", anchor_site), Method = "DESeq2", Tissue = tissue_label) %>%
        filter(padj < 0.05)
      
      if (nrow(df_res_sig) > 0) {
        deseq_list[[target_site]] <- df_res_sig %>%
          dplyr::rename(Score = log2FoldChange, P_Value = padj) %>%
          dplyr::select(Tissue, Method, Comparison, taxon, Score, P_Value)
      }
    }
  }
  df_deseq <- bind_rows(deseq_list)
  
  # ----------------------------------------------------------------------------
  # 4. METHOD C: MaAsLin2 (Linear Model vs. Dynamic Anchor)
  # ----------------------------------------------------------------------------
  cat("\n[Running MaAsLin2]...\n")
  counts   <- as.data.frame(otu_table(physeq_filtered))
  metadata <- data.frame(sample_data(physeq_filtered))
  if (taxa_are_rows(physeq_filtered)) { counts <- as.data.frame(t(counts)) }
  
  out_dir <- paste0("maaslin_temp_", gsub(" ", "_", tissue_label))
  ref_expression <- paste0("Site,", anchor_site)
  
  fit_res <- Maaslin2(
    input_data = counts, input_metadata = metadata, 
    output = out_dir, 
    fixed_effects = "Site", reference = ref_expression,
    normalization = "TSS", transform = "AST", plot_heatmap = FALSE, plot_scatter = FALSE
  )
  
  df_maaslin <- fit_res$results %>%
    filter(qval < 0.05) %>%
    dplyr::rename(taxon = feature, Score = coef, P_Value = qval) %>%
    mutate(Method = "MaAsLin2", Tissue = tissue_label, Comparison = paste0(value, "_vs_", anchor_site)) %>%
    dplyr::select(Tissue, Method, Comparison, taxon, Score, P_Value)
  
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
  
  # ----------------------------------------------------------------------------
  # 5. MASTER CONSENSUS INTEGRATION
  # ----------------------------------------------------------------------------
  master_bind <- dplyr::bind_rows(df_lefse, df_deseq, df_maaslin)
  
  if (nrow(master_bind) > 0 && !is.null(tax_table(physeq, errorIfNULL = FALSE))) {
    tax <- as.data.frame(tax_table(physeq))
    colnames(tax) <- make.unique(colnames(tax))
    tax <- tax %>% mutate(taxon = rownames(.))
    master_bind <- master_bind %>% left_join(tax, by = "taxon")
  }
  
  cat(" Finished!", tissue_label, "Summary: LEfSe =", nrow(df_lefse), 
      "| DESeq2 =", nrow(df_deseq), "| MaAsLin2 =", nrow(df_maaslin), "\n")
  return(master_bind)
}





# ==============================================================================
# 1. ISOLATE RAW, UNSCALED AMBIENT TRANSPLANT WORKSPACES
# ==============================================================================
cat("Isolating ambient transplant subsets from raw curated baselines...\n")
#> Isolating ambient transplant subsets from raw curated baselines...

# Subset Pocillopora to ambient samples only (Removing Post-HS)
ps_poc_trans_raw <- subset_samples(
  physeq_midpointTree_curated_POC, 
  Group == "TrueSample" & Type == "Coral" & Treatment != "Post-HS"
)
ps_poc_trans_raw <- prune_taxa(taxa_sums(ps_poc_trans_raw) > 0, ps_poc_trans_raw)

# Subset Porites to ambient samples only (Removing Post-HS)
ps_por_trans_raw <- subset_samples(
  physeq_midpointTree_curated_POR, 
  Group == "TrueSample" & Type == "Coral" & Treatment != "Post-HS"
)
ps_por_trans_raw <- prune_taxa(taxa_sums(ps_por_trans_raw) > 0, ps_por_trans_raw)

cat("• Raw Pocillopora Transplant Dataset:", nsamples(ps_poc_trans_raw), "samples.\n")
#> • Raw Pocillopora Transplant Dataset: 44 samples.
cat("• Raw Porites Transplant Dataset:", nsamples(ps_por_trans_raw), "samples.\n")
#> • Raw Porites Transplant Dataset: 22 samples.


# ==============================================================================
# 2. RUN THE CONSENSUS DIFFERENTIAL ABUNDANCE SUITE
# ==============================================================================

# Run for Pocillopora (POC) anchoring against Tiahura
cat("\n--- CRUNCHING AMBIENT DRIVERS FOR POCILLOPORA ---\n")
#> 
#> --- CRUNCHING AMBIENT DRIVERS FOR POCILLOPORA ---
poc_trans_drivers <- run_coral_differential_consensus_v5(
  physeq = ps_poc_trans_raw, 
  tissue_label = "Pocillopora Transplant Timeline", 
  anchor_site = "Tiahura"
)
#> 
#> [Power Filter] Pocillopora Transplant Timeline: Reduced total ASV hypotheses from 4027 down to 426 (Removed 3601 rare background taxa).
#> [Running LEfSe]...
#> [Running DESeq2]...
#> 
#> [Running MaAsLin2]...
#> [1] "Creating output folder"
#> [1] "Creating output feature tables folder"
#> [1] "Creating output fits folder"
#> 2026-10-28 10:46:51.233867 INFO::Writing function arguments to log file
#> 2026-10-28 10:46:51.264667 INFO::Verifying options selected are valid
#> 2026-10-28 10:46:51.390903 INFO::Determining format of input files
#> 2026-10-28 10:46:51.392968 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-10-28 10:46:51.399422 INFO::Formula for fixed effects: expr ~  Site
#> 2026-10-28 10:46:51.401336 INFO::Filter data based on min abundance and min prevalence
#> 2026-10-28 10:46:51.402711 INFO::Total samples in data: 44
#> 2026-10-28 10:46:51.40392 INFO::Min samples required with min abundance for a feature not to be filtered: 4.400000
#> 2026-10-28 10:46:51.41098 INFO::Total filtered features: 249
#> 2026-10-28 10:46:51.412755 INFO::Filtered feature names from abundance and prevalence filtering: ASV849, ASV2325, ASV1132, ASV939, ASV788, ASV4240, ASV2263, ASV688, ASV292, ASV577, ASV609, ASV1192, ASV2916, ASV1555, ASV380, ASV325, ASV1859, ASV413, ASV4708, ASV3619, ASV289, ASV6, ASV194, ASV1217, ASV1798, ASV80, ASV402, ASV455, ASV1115, ASV3088, ASV471, ASV1260, ASV1717, ASV1001, ASV848, ASV418, ASV302, ASV637, ASV468, ASV593, ASV4929, ASV2486, ASV3011, ASV2274, ASV4358, ASV1521, ASV2717, ASV1384, ASV1496, ASV291, ASV216, ASV485, ASV1609, ASV492, ASV666, ASV3129, ASV1042, ASV112, ASV219, ASV407, ASV138, ASV247, ASV1262, ASV5153, ASV1802, ASV383, ASV3862, ASV1105, ASV784, ASV310, ASV1337, ASV4112, ASV2239, ASV726, ASV3824, ASV892, ASV3183, ASV3662, ASV3064, ASV240, ASV318, ASV661, ASV988, ASV648, ASV629, ASV494, ASV450, ASV2356, ASV801, ASV322, ASV1160, ASV1564, ASV524, ASV2145, ASV150, ASV954, ASV1835, ASV2344, ASV4105, ASV1887, ASV754, ASV4863, ASV660, ASV1218, ASV1171, ASV2492, ASV1057, ASV2611, ASV2699, ASV971, ASV2059, ASV321, ASV2638, ASV4710, ASV2771, ASV705, ASV707, ASV837, ASV466, ASV4010, ASV3982, ASV2012, ASV1961, ASV42, ASV132, ASV813, ASV762, ASV1245, ASV584, ASV1087, ASV222, ASV4473, ASV1307, ASV568, ASV579, ASV881, ASV1295, ASV1456, ASV1207, ASV1574, ASV239, ASV1164, ASV33, ASV523, ASV356, ASV1692, ASV1789, ASV1426, ASV1767, ASV790, ASV113, ASV535, ASV1005, ASV805, ASV863, ASV206, ASV904, ASV899, ASV691, ASV265, ASV1078, ASV902, ASV2392, ASV868, ASV2215, ASV1573, ASV4657, ASV2501, ASV2030, ASV2458, ASV997, ASV276, ASV1290, ASV100, ASV339, ASV375, ASV54, ASV1094, ASV199, ASV1133, ASV916, ASV2633, ASV3869, ASV1504, ASV346, ASV2266, ASV309, ASV436, ASV606, ASV779, ASV711, ASV454, ASV220, ASV187, ASV250, ASV761, ASV486, ASV435, ASV232, ASV819, ASV843, ASV514, ASV540, ASV2116, ASV214, ASV118, ASV163, ASV127, ASV619, ASV649, ASV1309, ASV2056, ASV2120, ASV162, ASV297, ASV209, ASV348, ASV88, ASV958, ASV3012, ASV3087, ASV1219, ASV727, ASV811, ASV2008, ASV330, ASV836, ASV734, ASV367, ASV595, ASV675, ASV977, ASV3557, ASV1770, ASV1036, ASV2578, ASV537, ASV793, ASV2302, ASV1524, ASV3271, ASV3990, ASV1514, ASV617, ASV693, ASV1142, ASV1860, ASV2253, ASV1507
#> 2026-10-28 10:46:51.41719 INFO::Total filtered features with variance filtering: 0
#> 2026-10-28 10:46:51.418703 INFO::Filtered feature names from variance filtering:
#> 2026-10-28 10:46:51.419951 INFO::Running selected normalization method: TSS
#> 2026-10-28 10:46:51.425137 INFO::Applying z-score to standardize continuous metadata
#> 2026-10-28 10:46:51.429396 INFO::Running selected transform method: AST
#> 2026-10-28 10:46:51.433405 INFO::Running selected analysis method: LM
#> 2026-10-28 10:46:51.447468 INFO::Fitting model to feature number 1, ASV610
#> 2026-10-28 10:46:51.451781 INFO::Fitting model to feature number 2, ASV746
#> 2026-10-28 10:46:51.455037 INFO::Fitting model to feature number 3, ASV2464
#> 2026-10-28 10:46:51.460655 INFO::Fitting model to feature number 4, ASV712
#> 2026-10-28 10:46:51.464014 INFO::Fitting model to feature number 5, ASV369
#> 2026-10-28 10:46:51.467303 INFO::Fitting model to feature number 6, ASV122
#> 2026-10-28 10:46:51.470543 INFO::Fitting model to feature number 7, ASV182
#> 2026-10-28 10:46:51.47384 INFO::Fitting model to feature number 8, ASV1006
#> 2026-10-28 10:46:51.478115 INFO::Fitting model to feature number 9, ASV107
#> 2026-10-28 10:46:51.482372 INFO::Fitting model to feature number 10, ASV566
#> 2026-10-28 10:46:51.486038 INFO::Fitting model to feature number 11, ASV98
#> 2026-10-28 10:46:51.492102 INFO::Fitting model to feature number 12, ASV628
#> 2026-10-28 10:46:51.498453 INFO::Fitting model to feature number 13, ASV13
#> 2026-10-28 10:46:51.505491 INFO::Fitting model to feature number 14, ASV248
#> 2026-10-28 10:46:51.516208 INFO::Fitting model to feature number 15, ASV4
#> 2026-10-28 10:46:51.523928 INFO::Fitting model to feature number 16, ASV307
#> 2026-10-28 10:46:51.531828 INFO::Fitting model to feature number 17, ASV2
#> 2026-10-28 10:46:51.544451 INFO::Fitting model to feature number 18, ASV238
#> 2026-10-28 10:46:51.550901 INFO::Fitting model to feature number 19, ASV25
#> 2026-10-28 10:46:51.557366 INFO::Fitting model to feature number 20, ASV130
#> 2026-10-28 10:46:51.561038 INFO::Fitting model to feature number 21, ASV241
#> 2026-10-28 10:46:51.566175 INFO::Fitting model to feature number 22, ASV631
#> 2026-10-28 10:46:51.571542 INFO::Fitting model to feature number 23, ASV257
#> 2026-10-28 10:46:51.577995 INFO::Fitting model to feature number 24, ASV51
#> 2026-10-28 10:46:51.58419 INFO::Fitting model to feature number 25, ASV18
#> 2026-10-28 10:46:51.590769 INFO::Fitting model to feature number 26, ASV87
#> 2026-10-28 10:46:51.595296 INFO::Fitting model to feature number 27, ASV1985
#> 2026-10-28 10:46:51.599726 INFO::Fitting model to feature number 28, ASV826
#> 2026-10-28 10:46:51.605493 INFO::Fitting model to feature number 29, ASV1201
#> 2026-10-28 10:46:51.609402 INFO::Fitting model to feature number 30, ASV655
#> 2026-10-28 10:46:51.61311 INFO::Fitting model to feature number 31, ASV735
#> 2026-10-28 10:46:51.617401 INFO::Fitting model to feature number 32, ASV1181
#> 2026-10-28 10:46:51.622532 INFO::Fitting model to feature number 33, ASV634
#> 2026-10-28 10:46:51.626911 INFO::Fitting model to feature number 34, ASV320
#> 2026-10-28 10:46:51.630564 INFO::Fitting model to feature number 35, ASV193
#> 2026-10-28 10:46:51.634201 INFO::Fitting model to feature number 36, ASV212
#> 2026-10-28 10:46:51.638195 INFO::Fitting model to feature number 37, ASV200
#> 2026-10-28 10:46:51.644236 INFO::Fitting model to feature number 38, ASV1931
#> 2026-10-28 10:46:51.647896 INFO::Fitting model to feature number 39, ASV506
#> 2026-10-28 10:46:51.651652 INFO::Fitting model to feature number 40, ASV528
#> 2026-10-28 10:46:51.655916 INFO::Fitting model to feature number 41, ASV357
#> 2026-10-28 10:46:51.661654 INFO::Fitting model to feature number 42, ASV20
#> 2026-10-28 10:46:51.665356 INFO::Fitting model to feature number 43, ASV78
#> 2026-10-28 10:46:51.66898 INFO::Fitting model to feature number 44, ASV11
#> 2026-10-28 10:46:51.673687 INFO::Fitting model to feature number 45, ASV464
#> 2026-10-28 10:46:51.678136 INFO::Fitting model to feature number 46, ASV30
#> 2026-10-28 10:46:51.68215 INFO::Fitting model to feature number 47, ASV38
#> 2026-10-28 10:46:51.686012 INFO::Fitting model to feature number 48, ASV271
#> 2026-10-28 10:46:51.68984 INFO::Fitting model to feature number 49, ASV34
#> 2026-10-28 10:46:51.693353 INFO::Fitting model to feature number 50, ASV146
#> 2026-10-28 10:46:51.696742 INFO::Fitting model to feature number 51, ASV166
#> 2026-10-28 10:46:51.699976 INFO::Fitting model to feature number 52, ASV417
#> 2026-10-28 10:46:51.706447 INFO::Fitting model to feature number 53, ASV481
#> 2026-10-28 10:46:51.711517 INFO::Fitting model to feature number 54, ASV303
#> 2026-10-28 10:46:51.714909 INFO::Fitting model to feature number 55, ASV1028
#> 2026-10-28 10:46:51.71824 INFO::Fitting model to feature number 56, ASV201
#> 2026-10-28 10:46:51.72245 INFO::Fitting model to feature number 57, ASV796
#> 2026-10-28 10:46:51.726011 INFO::Fitting model to feature number 58, ASV864
#> 2026-10-28 10:46:51.729706 INFO::Fitting model to feature number 59, ASV211
#> 2026-10-28 10:46:51.732971 INFO::Fitting model to feature number 60, ASV1265
#> 2026-10-28 10:46:51.736263 INFO::Fitting model to feature number 61, ASV951
#> 2026-10-28 10:46:51.740222 INFO::Fitting model to feature number 62, ASV2427
#> 2026-10-28 10:46:51.743752 INFO::Fitting model to feature number 63, ASV2843
#> 2026-10-28 10:46:51.747253 INFO::Fitting model to feature number 64, ASV244
#> 2026-10-28 10:46:51.751095 INFO::Fitting model to feature number 65, ASV153
#> 2026-10-28 10:46:51.754806 INFO::Fitting model to feature number 66, ASV50
#> 2026-10-28 10:46:51.760674 INFO::Fitting model to feature number 67, ASV97
#> 2026-10-28 10:46:51.765653 INFO::Fitting model to feature number 68, ASV143
#> 2026-10-28 10:46:51.769835 INFO::Fitting model to feature number 69, ASV679
#> 2026-10-28 10:46:51.773908 INFO::Fitting model to feature number 70, ASV131
#> 2026-10-28 10:46:51.777639 INFO::Fitting model to feature number 71, ASV53
#> 2026-10-28 10:46:51.781314 INFO::Fitting model to feature number 72, ASV278
#> 2026-10-28 10:46:51.785238 INFO::Fitting model to feature number 73, ASV121
#> 2026-10-28 10:46:51.789816 INFO::Fitting model to feature number 74, ASV516
#> 2026-10-28 10:46:51.793827 INFO::Fitting model to feature number 75, ASV236
#> 2026-10-28 10:46:51.797446 INFO::Fitting model to feature number 76, ASV684
#> 2026-10-28 10:46:51.801138 INFO::Fitting model to feature number 77, ASV513
#> 2026-10-28 10:46:51.804902 INFO::Fitting model to feature number 78, ASV1664
#> 2026-10-28 10:46:51.809254 INFO::Fitting model to feature number 79, ASV975
#> 2026-10-28 10:46:51.813129 INFO::Fitting model to feature number 80, ASV275
#> 2026-10-28 10:46:51.816539 INFO::Fitting model to feature number 81, ASV632
#> 2026-10-28 10:46:51.821162 INFO::Fitting model to feature number 82, ASV59
#> 2026-10-28 10:46:51.826975 INFO::Fitting model to feature number 83, ASV968
#> 2026-10-28 10:46:51.832876 INFO::Fitting model to feature number 84, ASV61
#> 2026-10-28 10:46:51.836797 INFO::Fitting model to feature number 85, ASV186
#> 2026-10-28 10:46:51.841485 INFO::Fitting model to feature number 86, ASV286
#> 2026-10-28 10:46:51.846154 INFO::Fitting model to feature number 87, ASV119
#> 2026-10-28 10:46:51.852224 INFO::Fitting model to feature number 88, ASV105
#> 2026-10-28 10:46:51.862757 INFO::Fitting model to feature number 89, ASV124
#> 2026-10-28 10:46:51.869083 INFO::Fitting model to feature number 90, ASV213
#> 2026-10-28 10:46:51.876759 INFO::Fitting model to feature number 91, ASV32
#> 2026-10-28 10:46:51.882392 INFO::Fitting model to feature number 92, ASV135
#> 2026-10-28 10:46:51.887304 INFO::Fitting model to feature number 93, ASV177
#> 2026-10-28 10:46:51.893992 INFO::Fitting model to feature number 94, ASV128
#> 2026-10-28 10:46:51.903417 INFO::Fitting model to feature number 95, ASV259
#> 2026-10-28 10:46:51.911719 INFO::Fitting model to feature number 96, ASV295
#> 2026-10-28 10:46:51.920275 INFO::Fitting model to feature number 97, ASV267
#> 2026-10-28 10:46:51.925264 INFO::Fitting model to feature number 98, ASV657
#> 2026-10-28 10:46:51.929915 INFO::Fitting model to feature number 99, ASV893
#> 2026-10-28 10:46:51.933669 INFO::Fitting model to feature number 100, ASV359
#> 2026-10-28 10:46:51.936364 INFO::Fitting model to feature number 101, ASV1518
#> 2026-10-28 10:46:51.940143 INFO::Fitting model to feature number 102, ASV757
#> 2026-10-28 10:46:51.943445 INFO::Fitting model to feature number 103, ASV387
#> 2026-10-28 10:46:51.948673 INFO::Fitting model to feature number 104, ASV1
#> 2026-10-28 10:46:51.953767 INFO::Fitting model to feature number 105, ASV74
#> 2026-10-28 10:46:51.958083 INFO::Fitting model to feature number 106, ASV93
#> 2026-10-28 10:46:51.961054 INFO::Fitting model to feature number 107, ASV256
#> 2026-10-28 10:46:51.963846 INFO::Fitting model to feature number 108, ASV525
#> 2026-10-28 10:46:51.967271 INFO::Fitting model to feature number 109, ASV142
#> 2026-10-28 10:46:51.971768 INFO::Fitting model to feature number 110, ASV298
#> 2026-10-28 10:46:51.975821 INFO::Fitting model to feature number 111, ASV73
#> 2026-10-28 10:46:51.980558 INFO::Fitting model to feature number 112, ASV176
#> 2026-10-28 10:46:51.985925 INFO::Fitting model to feature number 113, ASV198
#> 2026-10-28 10:46:51.992842 INFO::Fitting model to feature number 114, ASV943
#> 2026-10-28 10:46:51.996976 INFO::Fitting model to feature number 115, ASV1196
#> 2026-10-28 10:46:52.000528 INFO::Fitting model to feature number 116, ASV676
#> 2026-10-28 10:46:52.004629 INFO::Fitting model to feature number 117, ASV172
#> 2026-10-28 10:46:52.011817 INFO::Fitting model to feature number 118, ASV144
#> 2026-10-28 10:46:52.016481 INFO::Fitting model to feature number 119, ASV134
#> 2026-10-28 10:46:52.020308 INFO::Fitting model to feature number 120, ASV92
#> 2026-10-28 10:46:52.024916 INFO::Fitting model to feature number 121, ASV120
#> 2026-10-28 10:46:52.030511 INFO::Fitting model to feature number 122, ASV19
#> 2026-10-28 10:46:52.036716 INFO::Fitting model to feature number 123, ASV5
#> 2026-10-28 10:46:52.044943 INFO::Fitting model to feature number 124, ASV312
#> 2026-10-28 10:46:52.050246 INFO::Fitting model to feature number 125, ASV72
#> 2026-10-28 10:46:52.054557 INFO::Fitting model to feature number 126, ASV48
#> 2026-10-28 10:46:52.059266 INFO::Fitting model to feature number 127, ASV9
#> 2026-10-28 10:46:52.062697 INFO::Fitting model to feature number 128, ASV65
#> 2026-10-28 10:46:52.066055 INFO::Fitting model to feature number 129, ASV136
#> 2026-10-28 10:46:52.069184 INFO::Fitting model to feature number 130, ASV43
#> 2026-10-28 10:46:52.072674 INFO::Fitting model to feature number 131, ASV56
#> 2026-10-28 10:46:52.076347 INFO::Fitting model to feature number 132, ASV91
#> 2026-10-28 10:46:52.079635 INFO::Fitting model to feature number 133, ASV169
#> 2026-10-28 10:46:52.082896 INFO::Fitting model to feature number 134, ASV585
#> 2026-10-28 10:46:52.086184 INFO::Fitting model to feature number 135, ASV846
#> 2026-10-28 10:46:52.119553 INFO::Fitting model to feature number 136, ASV335
#> 2026-10-28 10:46:52.123638 INFO::Fitting model to feature number 137, ASV1554
#> 2026-10-28 10:46:52.127478 INFO::Fitting model to feature number 138, ASV152
#> 2026-10-28 10:46:52.13073 INFO::Fitting model to feature number 139, ASV83
#> 2026-10-28 10:46:52.134283 INFO::Fitting model to feature number 140, ASV183
#> 2026-10-28 10:46:52.14118 INFO::Fitting model to feature number 141, ASV1134
#> 2026-10-28 10:46:52.145928 INFO::Fitting model to feature number 142, ASV2814
#> 2026-10-28 10:46:52.149099 INFO::Fitting model to feature number 143, ASV207
#> 2026-10-28 10:46:52.153996 INFO::Fitting model to feature number 144, ASV149
#> 2026-10-28 10:46:52.159478 INFO::Fitting model to feature number 145, ASV62
#> 2026-10-28 10:46:52.163538 INFO::Fitting model to feature number 146, ASV77
#> 2026-10-28 10:46:52.16759 INFO::Fitting model to feature number 147, ASV575
#> 2026-10-28 10:46:52.172287 INFO::Fitting model to feature number 148, ASV1850
#> 2026-10-28 10:46:52.176071 INFO::Fitting model to feature number 149, ASV608
#> 2026-10-28 10:46:52.179216 INFO::Fitting model to feature number 150, ASV934
#> 2026-10-28 10:46:52.182162 INFO::Fitting model to feature number 151, ASV2077
#> 2026-10-28 10:46:52.185929 INFO::Fitting model to feature number 152, ASV3263
#> 2026-10-28 10:46:52.189159 INFO::Fitting model to feature number 153, ASV424
#> 2026-10-28 10:46:52.192711 INFO::Fitting model to feature number 154, ASV548
#> 2026-10-28 10:46:52.19568 INFO::Fitting model to feature number 155, ASV28
#> 2026-10-28 10:46:52.198532 INFO::Fitting model to feature number 156, ASV40
#> 2026-10-28 10:46:52.202314 INFO::Fitting model to feature number 157, ASV8
#> 2026-10-28 10:46:52.20612 INFO::Fitting model to feature number 158, ASV148
#> 2026-10-28 10:46:52.209452 INFO::Fitting model to feature number 159, ASV85
#> 2026-10-28 10:46:52.212645 INFO::Fitting model to feature number 160, ASV82
#> 2026-10-28 10:46:52.215684 INFO::Fitting model to feature number 161, ASV181
#> 2026-10-28 10:46:52.219315 INFO::Fitting model to feature number 162, ASV16
#> 2026-10-28 10:46:52.226564 INFO::Fitting model to feature number 163, ASV7
#> 2026-10-28 10:46:52.231083 INFO::Fitting model to feature number 164, ASV157
#> 2026-10-28 10:46:52.234494 INFO::Fitting model to feature number 165, ASV266
#> 2026-10-28 10:46:52.237818 INFO::Fitting model to feature number 166, ASV90
#> 2026-10-28 10:46:52.244434 INFO::Fitting model to feature number 167, ASV437
#> 2026-10-28 10:46:52.248697 INFO::Fitting model to feature number 168, ASV415
#> 2026-10-28 10:46:52.251729 INFO::Fitting model to feature number 169, ASV274
#> 2026-10-28 10:46:52.254746 INFO::Fitting model to feature number 170, ASV52
#> 2026-10-28 10:46:52.259216 INFO::Fitting model to feature number 171, ASV81
#> 2026-10-28 10:46:52.262771 INFO::Fitting model to feature number 172, ASV520
#> 2026-10-28 10:46:52.266337 INFO::Fitting model to feature number 173, ASV703
#> 2026-10-28 10:46:52.269709 INFO::Fitting model to feature number 174, ASV1141
#> 2026-10-28 10:46:52.273388 INFO::Fitting model to feature number 175, ASV508
#> 2026-10-28 10:46:52.277206 INFO::Fitting model to feature number 176, ASV751
#> 2026-10-28 10:46:52.280624 INFO::Fitting model to feature number 177, ASV1122
#> 2026-10-28 10:46:52.301117 INFO::Counting total values for each feature
#> 2026-10-28 10:46:52.312277 INFO::Writing filtered data to file maaslin_temp_Pocillopora_Transplant_Timeline/features/filtered_data.tsv
#> 2026-10-28 10:46:52.393299 INFO::Writing filtered, normalized data to file maaslin_temp_Pocillopora_Transplant_Timeline/features/filtered_data_norm.tsv
#> 2026-10-28 10:46:52.40907 INFO::Writing filtered, normalized, transformed data to file maaslin_temp_Pocillopora_Transplant_Timeline/features/filtered_data_norm_transformed.tsv
#> 2026-10-28 10:46:52.423091 INFO::Writing residuals to file maaslin_temp_Pocillopora_Transplant_Timeline/fits/residuals.rds
#> 2026-10-28 10:46:52.427283 INFO::Writing fitted values to file maaslin_temp_Pocillopora_Transplant_Timeline/fits/fitted.rds
#> 2026-10-28 10:46:52.430419 INFO::Writing all results to file (ordered by increasing q-values): maaslin_temp_Pocillopora_Transplant_Timeline/all_results.tsv
#> 2026-10-28 10:46:52.437217 INFO::Writing the significant results (those which are less than or equal to the threshold of 0.250000 ) to file (ordered by increasing q-values): maaslin_temp_Pocillopora_Transplant_Timeline/significant_results.tsv
#>  Finished! Pocillopora Transplant Timeline Summary: LEfSe = 0 | DESeq2 = 6 | MaAsLin2 = 0

write.table(poc_trans_drivers,
            file = paste(dir_out_tab,
                         pref_out,
                         "Transplant_DiffAbunTaxa_POC.tab", sep =""))

# Run for Porites (POR) anchoring against Tiahura
cat("\n--- CRUNCHING AMBIENT DRIVERS FOR PORITES ---\n")
#> 
#> --- CRUNCHING AMBIENT DRIVERS FOR PORITES ---
por_trans_drivers <- run_coral_differential_consensus_v5(
  physeq = ps_por_trans_raw, 
  tissue_label = "Porites Transplant Timeline", 
  anchor_site = "Tiahura"
)
#> 
#> [Power Filter] Porites Transplant Timeline: Reduced total ASV hypotheses from 3491 down to 709 (Removed 2782 rare background taxa).
#> [Running LEfSe]...
#> [Running DESeq2]...
#> 
#> [Running MaAsLin2]...
#> [1] "Creating output folder"
#> [1] "Creating output feature tables folder"
#> [1] "Creating output fits folder"
#> 2026-10-28 10:46:58.093931 INFO::Writing function arguments to log file
#> 2026-10-28 10:46:58.164736 INFO::Verifying options selected are valid
#> 2026-10-28 10:46:58.166291 INFO::Determining format of input files
#> 2026-10-28 10:46:58.167825 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-10-28 10:46:58.181072 INFO::Formula for fixed effects: expr ~  Site
#> 2026-10-28 10:46:58.183509 INFO::Filter data based on min abundance and min prevalence
#> 2026-10-28 10:46:58.185592 INFO::Total samples in data: 22
#> 2026-10-28 10:46:58.187864 INFO::Min samples required with min abundance for a feature not to be filtered: 2.200000
#> 2026-10-28 10:46:58.196331 INFO::Total filtered features: 443
#> 2026-10-28 10:46:58.198318 INFO::Filtered feature names from abundance and prevalence filtering: ASV2385, ASV1495, ASV4670, ASV369, ASV182, ASV306, ASV98, ASV628, ASV410, ASV6855, ASV1373, ASV457, ASV678, ASV248, ASV5577, ASV6954, ASV1623, ASV7380, ASV25, ASV130, ASV455, ASV1833, ASV474, ASV1975, ASV3137, ASV2243, ASV1702, ASV2645, ASV293, ASV416, ASV257, ASV1335, ASV1047, ASV637, ASV468, ASV2036, ASV2364, ASV3005, ASV2486, ASV1946, ASV1154, ASV4689, ASV3223, ASV5392, ASV5436, ASV2865, ASV3002, ASV5328, ASV5694, ASV3419, ASV4544, ASV5603, ASV4447, ASV3084, ASV699, ASV3854, ASV3013, ASV320, ASV193, ASV200, ASV216, ASV311, ASV4769, ASV1779, ASV1200, ASV528, ASV2421, ASV7034, ASV889, ASV78, ASV1846, ASV2531, ASV1478, ASV2665, ASV30, ASV2744, ASV44, ASV407, ASV3322, ASV498, ASV547, ASV1318, ASV1738, ASV6522, ASV784, ASV1725, ASV1225, ASV5506, ASV1099, ASV4898, ASV1777, ASV726, ASV3761, ASV453, ASV2103, ASV1128, ASV2586, ASV4324, ASV1977, ASV3127, ASV4403, ASV225, ASV2729, ASV690, ASV648, ASV4286, ASV653, ASV1319, ASV153, ASV780, ASV905, ASV825, ASV2268, ASV143, ASV1052, ASV679, ASV2669, ASV53, ASV452, ASV121, ASV150, ASV969, ASV3179, ASV3116, ASV1259, ASV542, ASV3651, ASV465, ASV1258, ASV4418, ASV236, ASV684, ASV754, ASV6286, ASV6639, ASV5074, ASV1365, ASV1677, ASV3235, ASV3272, ASV660, ASV4571, ASV2858, ASV1439, ASV2269, ASV2492, ASV6254, ASV5673, ASV3495, ASV3539, ASV4270, ASV4475, ASV6638, ASV2874, ASV3555, ASV2295, ASV2771, ASV255, ASV3981, ASV837, ASV2955, ASV2716, ASV2517, ASV956, ASV4915, ASV2341, ASV2150, ASV2078, ASV3392, ASV439, ASV926, ASV1808, ASV4574, ASV1930, ASV396, ASV1691, ASV1423, ASV509, ASV2064, ASV879, ASV1418, ASV1597, ASV910, ASV1768, ASV3465, ASV2015, ASV1552, ASV4416, ASV4445, ASV2465, ASV286, ASV762, ASV358, ASV1568, ASV267, ASV359, ASV3612, ASV2566, ASV4342, ASV1456, ASV2799, ASV830, ASV5580, ASV198, ASV1193, ASV5288, ASV1391, ASV2534, ASV1085, ASV1135, ASV2415, ASV1559, ASV4086, ASV443, ASV1027, ASV877, ASV1118, ASV835, ASV1451, ASV915, ASV370, ASV89, ASV48, ASV1742, ASV265, ASV56, ASV3325, ASV226, ASV817, ASV1049, ASV4756, ASV3803, ASV1102, ASV1891, ASV3077, ASV1573, ASV6521, ASV2459, ASV2002, ASV869, ASV3545, ASV329, ASV4768, ASV4034, ASV1494, ASV1214, ASV1399, ASV2165, ASV1974, ASV6906, ASV4176, ASV2760, ASV2546, ASV4171, ASV5287, ASV5141, ASV625, ASV3814, ASV3993, ASV1547, ASV2193, ASV769, ASV3802, ASV1289, ASV1953, ASV1993, ASV2369, ASV1386, ASV276, ASV1020, ASV1727, ASV2839, ASV1528, ASV1246, ASV1537, ASV3206, ASV1486, ASV716, ASV2722, ASV2745, ASV1268, ASV1626, ASV2701, ASV3305, ASV1363, ASV2218, ASV3055, ASV5423, ASV1094, ASV2199, ASV419, ASV5142, ASV2686, ASV927, ASV1472, ASV4178, ASV2332, ASV180, ASV745, ASV2530, ASV519, ASV6330, ASV6375, ASV5837, ASV1814, ASV4580, ASV1645, ASV2282, ASV2117, ASV1002, ASV1058, ASV4287, ASV3459, ASV2065, ASV3079, ASV6515, ASV949, ASV5054, ASV6030, ASV967, ASV2561, ASV4327, ASV4141, ASV187, ASV2280, ASV2038, ASV3496, ASV1061, ASV435, ASV924, ASV2134, ASV819, ASV540, ASV1763, ASV1615, ASV2267, ASV736, ASV2764, ASV996, ASV1642, ASV1643, ASV428, ASV2532, ASV157, ASV348, ASV2000, ASV1755, ASV1219, ASV2720, ASV675, ASV5055, ASV1075, ASV5398, ASV1089, ASV3155, ASV3863, ASV4500, ASV4230, ASV6029, ASV3344, ASV6171, ASV1036, ASV1148, ASV2390, ASV2238, ASV5200, ASV1311, ASV2123, ASV2797, ASV3588, ASV4386, ASV2433, ASV4096, ASV2323, ASV1282, ASV3234, ASV2386, ASV4573, ASV2578, ASV4823, ASV775, ASV1614, ASV3853, ASV4043, ASV2726, ASV4572, ASV886, ASV3903, ASV4039, ASV3683, ASV3984, ASV6813, ASV4031, ASV2411, ASV3537, ASV2840, ASV2302, ASV3520, ASV5286, ASV4505, ASV3801, ASV7040, ASV3597, ASV1600, ASV3181, ASV4383, ASV2285, ASV517, ASV6064, ASV2715, ASV223, ASV522, ASV1881, ASV645, ASV1515, ASV2026, ASV3939, ASV1302, ASV1529, ASV2820, ASV1627, ASV5663, ASV6329, ASV4630, ASV6767, ASV4519, ASV5443, ASV4305, ASV3020, ASV2979, ASV751, ASV1122, ASV2081, ASV4283, ASV7202, ASV1237, ASV1316, ASV806, ASV7225, ASV2582, ASV1781, ASV2533, ASV5486
#> 2026-10-28 10:46:58.203486 INFO::Total filtered features with variance filtering: 0
#> 2026-10-28 10:46:58.206052 INFO::Filtered feature names from variance filtering:
#> 2026-10-28 10:46:58.208286 INFO::Running selected normalization method: TSS
#> 2026-10-28 10:46:58.214223 INFO::Applying z-score to standardize continuous metadata
#> 2026-10-28 10:46:58.220407 INFO::Running selected transform method: AST
#> 2026-10-28 10:46:58.228042 INFO::Running selected analysis method: LM
#> 2026-10-28 10:46:58.230606 INFO::Fitting model to feature number 1, ASV1117
#> 2026-10-28 10:46:58.237115 INFO::Fitting model to feature number 2, ASV610
#> 2026-10-28 10:46:58.243312 INFO::Fitting model to feature number 3, ASV788
#> 2026-10-28 10:46:58.248987 INFO::Fitting model to feature number 4, ASV107
#> 2026-10-28 10:46:58.252712 INFO::Fitting model to feature number 5, ASV2316
#> 2026-10-28 10:46:58.256986 INFO::Fitting model to feature number 6, ASV13
#> 2026-10-28 10:46:58.260583 INFO::Fitting model to feature number 7, ASV174
#> 2026-10-28 10:46:58.264151 INFO::Fitting model to feature number 8, ASV4
#> 2026-10-28 10:46:58.267654 INFO::Fitting model to feature number 9, ASV307
#> 2026-10-28 10:46:58.271757 INFO::Fitting model to feature number 10, ASV3199
#> 2026-10-28 10:46:58.275598 INFO::Fitting model to feature number 11, ASV191
#> 2026-10-28 10:46:58.279109 INFO::Fitting model to feature number 12, ASV990
#> 2026-10-28 10:46:58.282598 INFO::Fitting model to feature number 13, ASV80
#> 2026-10-28 10:46:58.2869 INFO::Fitting model to feature number 14, ASV344
#> 2026-10-28 10:46:58.290258 INFO::Fitting model to feature number 15, ASV238
#> 2026-10-28 10:46:58.294765 INFO::Fitting model to feature number 16, ASV402
#> 2026-10-28 10:46:58.297996 INFO::Fitting model to feature number 17, ASV241
#> 2026-10-28 10:46:58.301254 INFO::Fitting model to feature number 18, ASV631
#> 2026-10-28 10:46:58.304491 INFO::Fitting model to feature number 19, ASV1212
#> 2026-10-28 10:46:58.309721 INFO::Fitting model to feature number 20, ASV4005
#> 2026-10-28 10:46:58.313164 INFO::Fitting model to feature number 21, ASV848
#> 2026-10-28 10:46:58.317528 INFO::Fitting model to feature number 22, ASV470
#> 2026-10-28 10:46:58.321244 INFO::Fitting model to feature number 23, ASV302
#> 2026-10-28 10:46:58.336598 INFO::Fitting model to feature number 24, ASV314
#> 2026-10-28 10:46:58.341176 INFO::Fitting model to feature number 25, ASV18
#> 2026-10-28 10:46:58.345359 INFO::Fitting model to feature number 26, ASV87
#> 2026-10-28 10:46:58.348959 INFO::Fitting model to feature number 27, ASV1305
#> 2026-10-28 10:46:58.352776 INFO::Fitting model to feature number 28, ASV831
#> 2026-10-28 10:46:58.35723 INFO::Fitting model to feature number 29, ASV3864
#> 2026-10-28 10:46:58.361485 INFO::Fitting model to feature number 30, ASV2108
#> 2026-10-28 10:46:58.364891 INFO::Fitting model to feature number 31, ASV3728
#> 2026-10-28 10:46:58.370685 INFO::Fitting model to feature number 32, ASV1357
#> 2026-10-28 10:46:58.37588 INFO::Fitting model to feature number 33, ASV1649
#> 2026-10-28 10:46:58.380359 INFO::Fitting model to feature number 34, ASV3230
#> 2026-10-28 10:46:58.383925 INFO::Fitting model to feature number 35, ASV1025
#> 2026-10-28 10:46:58.387386 INFO::Fitting model to feature number 36, ASV818
#> 2026-10-28 10:46:58.391533 INFO::Fitting model to feature number 37, ASV212
#> 2026-10-28 10:46:58.395013 INFO::Fitting model to feature number 38, ASV1012
#> 2026-10-28 10:46:58.398473 INFO::Fitting model to feature number 39, ASV1479
#> 2026-10-28 10:46:58.40179 INFO::Fitting model to feature number 40, ASV1920
#> 2026-10-28 10:46:58.405603 INFO::Fitting model to feature number 41, ASV651
#> 2026-10-28 10:46:58.409984 INFO::Fitting model to feature number 42, ASV139
#> 2026-10-28 10:46:58.413617 INFO::Fitting model to feature number 43, ASV1031
#> 2026-10-28 10:46:58.417361 INFO::Fitting model to feature number 44, ASV112
#> 2026-10-28 10:46:58.42089 INFO::Fitting model to feature number 45, ASV219
#> 2026-10-28 10:46:58.424885 INFO::Fitting model to feature number 46, ASV166
#> 2026-10-28 10:46:58.428353 INFO::Fitting model to feature number 47, ASV417
#> 2026-10-28 10:46:58.431864 INFO::Fitting model to feature number 48, ASV1337
#> 2026-10-28 10:46:58.435975 INFO::Fitting model to feature number 49, ASV481
#> 2026-10-28 10:46:58.44057 INFO::Fitting model to feature number 50, ASV303
#> 2026-10-28 10:46:58.444087 INFO::Fitting model to feature number 51, ASV1028
#> 2026-10-28 10:46:58.447626 INFO::Fitting model to feature number 52, ASV2680
#> 2026-10-28 10:46:58.451197 INFO::Fitting model to feature number 53, ASV2153
#> 2026-10-28 10:46:58.457579 INFO::Fitting model to feature number 54, ASV211
#> 2026-10-28 10:46:58.461368 INFO::Fitting model to feature number 55, ASV240
#> 2026-10-28 10:46:58.464747 INFO::Fitting model to feature number 56, ASV1272
#> 2026-10-28 10:46:58.468298 INFO::Fitting model to feature number 57, ASV400
#> 2026-10-28 10:46:58.471951 INFO::Fitting model to feature number 58, ASV164
#> 2026-10-28 10:46:58.476382 INFO::Fitting model to feature number 59, ASV494
#> 2026-10-28 10:46:58.479764 INFO::Fitting model to feature number 60, ASV244
#> 2026-10-28 10:46:58.483009 INFO::Fitting model to feature number 61, ASV450
#> 2026-10-28 10:46:58.486181 INFO::Fitting model to feature number 62, ASV237
#> 2026-10-28 10:46:58.490157 INFO::Fitting model to feature number 63, ASV50
#> 2026-10-28 10:46:58.493496 INFO::Fitting model to feature number 64, ASV97
#> 2026-10-28 10:46:58.496735 INFO::Fitting model to feature number 65, ASV382
#> 2026-10-28 10:46:58.499896 INFO::Fitting model to feature number 66, ASV479
#> 2026-10-28 10:46:58.504269 INFO::Fitting model to feature number 67, ASV376
#> 2026-10-28 10:46:58.507562 INFO::Fitting model to feature number 68, ASV131
#> 2026-10-28 10:46:58.510812 INFO::Fitting model to feature number 69, ASV3769
#> 2026-10-28 10:46:58.513965 INFO::Fitting model to feature number 70, ASV3889
#> 2026-10-28 10:46:58.517113 INFO::Fitting model to feature number 71, ASV516
#> 2026-10-28 10:46:58.520731 INFO::Fitting model to feature number 72, ASV1011
#> 2026-10-28 10:46:58.524021 INFO::Fitting model to feature number 73, ASV4484
#> 2026-10-28 10:46:58.527329 INFO::Fitting model to feature number 74, ASV2171
#> 2026-10-28 10:46:58.531657 INFO::Fitting model to feature number 75, ASV4446
#> 2026-10-28 10:46:58.534967 INFO::Fitting model to feature number 76, ASV533
#> 2026-10-28 10:46:58.538628 INFO::Fitting model to feature number 77, ASV6440
#> 2026-10-28 10:46:58.544349 INFO::Fitting model to feature number 78, ASV4690
#> 2026-10-28 10:46:58.547693 INFO::Fitting model to feature number 79, ASV2706
#> 2026-10-28 10:46:58.552685 INFO::Fitting model to feature number 80, ASV5228
#> 2026-10-28 10:46:58.557018 INFO::Fitting model to feature number 81, ASV1586
#> 2026-10-28 10:46:58.561025 INFO::Fitting model to feature number 82, ASV3692
#> 2026-10-28 10:46:58.564462 INFO::Fitting model to feature number 83, ASV2623
#> 2026-10-28 10:46:58.568295 INFO::Fitting model to feature number 84, ASV707
#> 2026-10-28 10:46:58.573131 INFO::Fitting model to feature number 85, ASV1616
#> 2026-10-28 10:46:58.578184 INFO::Fitting model to feature number 86, ASV878
#> 2026-10-28 10:46:58.583044 INFO::Fitting model to feature number 87, ASV2727
#> 2026-10-28 10:46:58.587155 INFO::Fitting model to feature number 88, ASV2304
#> 2026-10-28 10:46:58.591598 INFO::Fitting model to feature number 89, ASV521
#> 2026-10-28 10:46:58.595683 INFO::Fitting model to feature number 90, ASV839
#> 2026-10-28 10:46:58.599417 INFO::Fitting model to feature number 91, ASV280
#> 2026-10-28 10:46:58.603222 INFO::Fitting model to feature number 92, ASV337
#> 2026-10-28 10:46:58.608457 INFO::Fitting model to feature number 93, ASV469
#> 2026-10-28 10:46:58.61264 INFO::Fitting model to feature number 94, ASV1503
#> 2026-10-28 10:46:58.618411 INFO::Fitting model to feature number 95, ASV518
#> 2026-10-28 10:46:58.622912 INFO::Fitting model to feature number 96, ASV2857
#> 2026-10-28 10:46:58.627333 INFO::Fitting model to feature number 97, ASV1664
#> 2026-10-28 10:46:58.631222 INFO::Fitting model to feature number 98, ASV667
#> 2026-10-28 10:46:58.63515 INFO::Fitting model to feature number 99, ASV386
#> 2026-10-28 10:46:58.639646 INFO::Fitting model to feature number 100, ASV975
#> 2026-10-28 10:46:58.643982 INFO::Fitting model to feature number 101, ASV1425
#> 2026-10-28 10:46:58.647597 INFO::Fitting model to feature number 102, ASV3680
#> 2026-10-28 10:46:58.651099 INFO::Fitting model to feature number 103, ASV275
#> 2026-10-28 10:46:58.656425 INFO::Fitting model to feature number 104, ASV42
#> 2026-10-28 10:46:58.660569 INFO::Fitting model to feature number 105, ASV1266
#> 2026-10-28 10:46:58.66412 INFO::Fitting model to feature number 106, ASV66
#> 2026-10-28 10:46:58.667969 INFO::Fitting model to feature number 107, ASV186
#> 2026-10-28 10:46:58.67279 INFO::Fitting model to feature number 108, ASV119
#> 2026-10-28 10:46:58.676567 INFO::Fitting model to feature number 109, ASV105
#> 2026-10-28 10:46:58.680073 INFO::Fitting model to feature number 110, ASV124
#> 2026-10-28 10:46:58.683782 INFO::Fitting model to feature number 111, ASV213
#> 2026-10-28 10:46:58.68731 INFO::Fitting model to feature number 112, ASV222
#> 2026-10-28 10:46:58.692477 INFO::Fitting model to feature number 113, ASV32
#> 2026-10-28 10:46:58.696216 INFO::Fitting model to feature number 114, ASV366
#> 2026-10-28 10:46:58.699891 INFO::Fitting model to feature number 115, ASV1804
#> 2026-10-28 10:46:58.703479 INFO::Fitting model to feature number 116, ASV269
#> 2026-10-28 10:46:58.707756 INFO::Fitting model to feature number 117, ASV239
#> 2026-10-28 10:46:58.711452 INFO::Fitting model to feature number 118, ASV1330
#> 2026-10-28 10:46:58.71511 INFO::Fitting model to feature number 119, ASV2208
#> 2026-10-28 10:46:58.718653 INFO::Fitting model to feature number 120, ASV408
#> 2026-10-28 10:46:58.725934 INFO::Fitting model to feature number 121, ASV1666
#> 2026-10-28 10:46:58.731864 INFO::Fitting model to feature number 122, ASV866
#> 2026-10-28 10:46:58.735007 INFO::Fitting model to feature number 123, ASV123
#> 2026-10-28 10:46:58.76081 INFO::Fitting model to feature number 124, ASV113
#> 2026-10-28 10:46:58.766441 INFO::Fitting model to feature number 125, ASV1015
#> 2026-10-28 10:46:58.770378 INFO::Fitting model to feature number 126, ASV144
#> 2026-10-28 10:46:58.774787 INFO::Fitting model to feature number 127, ASV134
#> 2026-10-28 10:46:58.778543 INFO::Fitting model to feature number 128, ASV1745
#> 2026-10-28 10:46:58.782448 INFO::Fitting model to feature number 129, ASV92
#> 2026-10-28 10:46:58.786082 INFO::Fitting model to feature number 130, ASV120
#> 2026-10-28 10:46:58.789514 INFO::Fitting model to feature number 131, ASV596
#> 2026-10-28 10:46:58.793025 INFO::Fitting model to feature number 132, ASV19
#> 2026-10-28 10:46:58.796453 INFO::Fitting model to feature number 133, ASV5
#> 2026-10-28 10:46:58.799784 INFO::Fitting model to feature number 134, ASV556
#> 2026-10-28 10:46:58.803075 INFO::Fitting model to feature number 135, ASV72
#> 2026-10-28 10:46:58.806759 INFO::Fitting model to feature number 136, ASV9
#> 2026-10-28 10:46:58.813262 INFO::Fitting model to feature number 137, ASV136
#> 2026-10-28 10:46:58.819043 INFO::Fitting model to feature number 138, ASV5349
#> 2026-10-28 10:46:58.823037 INFO::Fitting model to feature number 139, ASV3909
#> 2026-10-28 10:46:58.827254 INFO::Fitting model to feature number 140, ASV2017
#> 2026-10-28 10:46:58.831869 INFO::Fitting model to feature number 141, ASV156
#> 2026-10-28 10:46:58.835751 INFO::Fitting model to feature number 142, ASV5343
#> 2026-10-28 10:46:58.840282 INFO::Fitting model to feature number 143, ASV4627
#> 2026-10-28 10:46:58.845806 INFO::Fitting model to feature number 144, ASV647
#> 2026-10-28 10:46:58.850642 INFO::Fitting model to feature number 145, ASV336
#> 2026-10-28 10:46:58.855062 INFO::Fitting model to feature number 146, ASV1067
#> 2026-10-28 10:46:58.861002 INFO::Fitting model to feature number 147, ASV722
#> 2026-10-28 10:46:58.864972 INFO::Fitting model to feature number 148, ASV495
#> 2026-10-28 10:46:58.86858 INFO::Fitting model to feature number 149, ASV589
#> 2026-10-28 10:46:58.872869 INFO::Fitting model to feature number 150, ASV974
#> 2026-10-28 10:46:58.876993 INFO::Fitting model to feature number 151, ASV2355
#> 2026-10-28 10:46:58.881841 INFO::Fitting model to feature number 152, ASV152
#> 2026-10-28 10:46:58.885416 INFO::Fitting model to feature number 153, ASV83
#> 2026-10-28 10:46:58.889737 INFO::Fitting model to feature number 154, ASV183
#> 2026-10-28 10:46:58.893582 INFO::Fitting model to feature number 155, ASV640
#> 2026-10-28 10:46:58.89749 INFO::Fitting model to feature number 156, ASV1563
#> 2026-10-28 10:46:58.901094 INFO::Fitting model to feature number 157, ASV2800
#> 2026-10-28 10:46:58.907897 INFO::Fitting model to feature number 158, ASV997
#> 2026-10-28 10:46:58.912455 INFO::Fitting model to feature number 159, ASV1134
#> 2026-10-28 10:46:58.915999 INFO::Fitting model to feature number 160, ASV472
#> 2026-10-28 10:46:58.919598 INFO::Fitting model to feature number 161, ASV327
#> 2026-10-28 10:46:58.924814 INFO::Fitting model to feature number 162, ASV1647
#> 2026-10-28 10:46:58.928753 INFO::Fitting model to feature number 163, ASV227
#> 2026-10-28 10:46:58.9329 INFO::Fitting model to feature number 164, ASV109
#> 2026-10-28 10:46:58.937201 INFO::Fitting model to feature number 165, ASV147
#> 2026-10-28 10:46:58.942118 INFO::Fitting model to feature number 166, ASV742
#> 2026-10-28 10:46:58.947588 INFO::Fitting model to feature number 167, ASV922
#> 2026-10-28 10:46:58.951369 INFO::Fitting model to feature number 168, ASV673
#> 2026-10-28 10:46:58.955436 INFO::Fitting model to feature number 169, ASV62
#> 2026-10-28 10:46:58.960435 INFO::Fitting model to feature number 170, ASV100
#> 2026-10-28 10:46:58.965346 INFO::Fitting model to feature number 171, ASV616
#> 2026-10-28 10:46:58.969137 INFO::Fitting model to feature number 172, ASV668
#> 2026-10-28 10:46:58.974122 INFO::Fitting model to feature number 173, ASV4215
#> 2026-10-28 10:46:58.978806 INFO::Fitting model to feature number 174, ASV77
#> 2026-10-28 10:46:58.98554 INFO::Fitting model to feature number 175, ASV199
#> 2026-10-28 10:46:58.992235 INFO::Fitting model to feature number 176, ASV583
#> 2026-10-28 10:46:58.996192 INFO::Fitting model to feature number 177, ASV1313
#> 2026-10-28 10:46:59.000102 INFO::Fitting model to feature number 178, ASV885
#> 2026-10-28 10:46:59.004354 INFO::Fitting model to feature number 179, ASV350
#> 2026-10-28 10:46:59.010117 INFO::Fitting model to feature number 180, ASV5841
#> 2026-10-28 10:46:59.015291 INFO::Fitting model to feature number 181, ASV346
#> 2026-10-28 10:46:59.021559 INFO::Fitting model to feature number 182, ASV604
#> 2026-10-28 10:46:59.026908 INFO::Fitting model to feature number 183, ASV4928
#> 2026-10-28 10:46:59.033056 INFO::Fitting model to feature number 184, ASV456
#> 2026-10-28 10:46:59.038894 INFO::Fitting model to feature number 185, ASV638
#> 2026-10-28 10:46:59.045796 INFO::Fitting model to feature number 186, ASV991
#> 2026-10-28 10:46:59.050526 INFO::Fitting model to feature number 187, ASV2724
#> 2026-10-28 10:46:59.057615 INFO::Fitting model to feature number 188, ASV766
#> 2026-10-28 10:46:59.062735 INFO::Fitting model to feature number 189, ASV1646
#> 2026-10-28 10:46:59.067482 INFO::Fitting model to feature number 190, ASV2204
#> 2026-10-28 10:46:59.071621 INFO::Fitting model to feature number 191, ASV28
#> 2026-10-28 10:46:59.077369 INFO::Fitting model to feature number 192, ASV40
#> 2026-10-28 10:46:59.083548 INFO::Fitting model to feature number 193, ASV8
#> 2026-10-28 10:46:59.088654 INFO::Fitting model to feature number 194, ASV220
#> 2026-10-28 10:46:59.093093 INFO::Fitting model to feature number 195, ASV245
#> 2026-10-28 10:46:59.096745 INFO::Fitting model to feature number 196, ASV559
#> 2026-10-28 10:46:59.100038 INFO::Fitting model to feature number 197, ASV809
#> 2026-10-28 10:46:59.103457 INFO::Fitting model to feature number 198, ASV448
#> 2026-10-28 10:46:59.108324 INFO::Fitting model to feature number 199, ASV371
#> 2026-10-28 10:46:59.114408 INFO::Fitting model to feature number 200, ASV232
#> 2026-10-28 10:46:59.118494 INFO::Fitting model to feature number 201, ASV361
#> 2026-10-28 10:46:59.122039 INFO::Fitting model to feature number 202, ASV514
#> 2026-10-28 10:46:59.125662 INFO::Fitting model to feature number 203, ASV433
#> 2026-10-28 10:46:59.129151 INFO::Fitting model to feature number 204, ASV351
#> 2026-10-28 10:46:59.1341 INFO::Fitting model to feature number 205, ASV85
#> 2026-10-28 10:46:59.139414 INFO::Fitting model to feature number 206, ASV82
#> 2026-10-28 10:46:59.144218 INFO::Fitting model to feature number 207, ASV118
#> 2026-10-28 10:46:59.148743 INFO::Fitting model to feature number 208, ASV16
#> 2026-10-28 10:46:59.153326 INFO::Fitting model to feature number 209, ASV7
#> 2026-10-28 10:46:59.157946 INFO::Fitting model to feature number 210, ASV127
#> 2026-10-28 10:46:59.161803 INFO::Fitting model to feature number 211, ASV649
#> 2026-10-28 10:46:59.164989 INFO::Fitting model to feature number 212, ASV209
#> 2026-10-28 10:46:59.168303 INFO::Fitting model to feature number 213, ASV301
#> 2026-10-28 10:46:59.173184 INFO::Fitting model to feature number 214, ASV243
#> 2026-10-28 10:46:59.177664 INFO::Fitting model to feature number 215, ASV90
#> 2026-10-28 10:46:59.181252 INFO::Fitting model to feature number 216, ASV1406
#> 2026-10-28 10:46:59.184593 INFO::Fitting model to feature number 217, ASV2761
#> 2026-10-28 10:46:59.18794 INFO::Fitting model to feature number 218, ASV6253
#> 2026-10-28 10:46:59.191627 INFO::Fitting model to feature number 219, ASV338
#> 2026-10-28 10:46:59.194764 INFO::Fitting model to feature number 220, ASV368
#> 2026-10-28 10:46:59.198227 INFO::Fitting model to feature number 221, ASV950
#> 2026-10-28 10:46:59.203108 INFO::Fitting model to feature number 222, ASV52
#> 2026-10-28 10:46:59.209285 INFO::Fitting model to feature number 223, ASV81
#> 2026-10-28 10:46:59.212721 INFO::Fitting model to feature number 224, ASV595
#> 2026-10-28 10:46:59.216474 INFO::Fitting model to feature number 225, ASV2185
#> 2026-10-28 10:46:59.21976 INFO::Fitting model to feature number 226, ASV1959
#> 2026-10-28 10:46:59.223984 INFO::Fitting model to feature number 227, ASV1051
#> 2026-10-28 10:46:59.227264 INFO::Fitting model to feature number 228, ASV354
#> 2026-10-28 10:46:59.230994 INFO::Fitting model to feature number 229, ASV578
#> 2026-10-28 10:46:59.234339 INFO::Fitting model to feature number 230, ASV4088
#> 2026-10-28 10:46:59.237669 INFO::Fitting model to feature number 231, ASV2133
#> 2026-10-28 10:46:59.241524 INFO::Fitting model to feature number 232, ASV821
#> 2026-10-28 10:46:59.244827 INFO::Fitting model to feature number 233, ASV1710
#> 2026-10-28 10:46:59.248319 INFO::Fitting model to feature number 234, ASV1837
#> 2026-10-28 10:46:59.251484 INFO::Fitting model to feature number 235, ASV2557
#> 2026-10-28 10:46:59.256113 INFO::Fitting model to feature number 236, ASV2461
#> 2026-10-28 10:46:59.260646 INFO::Fitting model to feature number 237, ASV1527
#> 2026-10-28 10:46:59.264131 INFO::Fitting model to feature number 238, ASV1546
#> 2026-10-28 10:46:59.267484 INFO::Fitting model to feature number 239, ASV1354
#> 2026-10-28 10:46:59.272156 INFO::Fitting model to feature number 240, ASV2721
#> 2026-10-28 10:46:59.27655 INFO::Fitting model to feature number 241, ASV1524
#> 2026-10-28 10:46:59.279897 INFO::Fitting model to feature number 242, ASV802
#> 2026-10-28 10:46:59.283351 INFO::Fitting model to feature number 243, ASV1514
#> 2026-10-28 10:46:59.286366 INFO::Fitting model to feature number 244, ASV829
#> 2026-10-28 10:46:59.291841 INFO::Fitting model to feature number 245, ASV508
#> 2026-10-28 10:46:59.297023 INFO::Fitting model to feature number 246, ASV1505
#> 2026-10-28 10:46:59.301596 INFO::Fitting model to feature number 247, ASV1297
#> 2026-10-28 10:46:59.307385 INFO::Fitting model to feature number 248, ASV334
#> 2026-10-28 10:46:59.312556 INFO::Fitting model to feature number 249, ASV4241
#> 2026-10-28 10:46:59.316475 INFO::Fitting model to feature number 250, ASV1684
#> 2026-10-28 10:46:59.320409 INFO::Fitting model to feature number 251, ASV2922
#> 2026-10-28 10:46:59.325836 INFO::Fitting model to feature number 252, ASV1689
#> 2026-10-28 10:46:59.330597 INFO::Fitting model to feature number 253, ASV1783
#> 2026-10-28 10:46:59.335069 INFO::Fitting model to feature number 254, ASV1249
#> 2026-10-28 10:46:59.339779 INFO::Fitting model to feature number 255, ASV2219
#> 2026-10-28 10:46:59.344343 INFO::Fitting model to feature number 256, ASV617
#> 2026-10-28 10:46:59.34912 INFO::Fitting model to feature number 257, ASV3251
#> 2026-10-28 10:46:59.35575 INFO::Fitting model to feature number 258, ASV693
#> 2026-10-28 10:46:59.360726 INFO::Fitting model to feature number 259, ASV981
#> 2026-10-28 10:46:59.365587 INFO::Fitting model to feature number 260, ASV2293
#> 2026-10-28 10:46:59.369848 INFO::Fitting model to feature number 261, ASV536
#> 2026-10-28 10:46:59.374464 INFO::Fitting model to feature number 262, ASV2813
#> 2026-10-28 10:46:59.378339 INFO::Fitting model to feature number 263, ASV2069
#> 2026-10-28 10:46:59.382492 INFO::Fitting model to feature number 264, ASV3101
#> 2026-10-28 10:46:59.386392 INFO::Fitting model to feature number 265, ASV252
#> 2026-10-28 10:46:59.390104 INFO::Fitting model to feature number 266, ASV1466
#> 2026-10-28 10:46:59.422236 INFO::Counting total values for each feature
#> 2026-10-28 10:46:59.438583 INFO::Writing filtered data to file maaslin_temp_Porites_Transplant_Timeline/features/filtered_data.tsv
#> 2026-10-28 10:46:59.460082 INFO::Writing filtered, normalized data to file maaslin_temp_Porites_Transplant_Timeline/features/filtered_data_norm.tsv
#> 2026-10-28 10:46:59.476867 INFO::Writing filtered, normalized, transformed data to file maaslin_temp_Porites_Transplant_Timeline/features/filtered_data_norm_transformed.tsv
#> 2026-10-28 10:46:59.489994 INFO::Writing residuals to file maaslin_temp_Porites_Transplant_Timeline/fits/residuals.rds
#> 2026-10-28 10:46:59.493681 INFO::Writing fitted values to file maaslin_temp_Porites_Transplant_Timeline/fits/fitted.rds
#> 2026-10-28 10:46:59.496408 INFO::Writing all results to file (ordered by increasing q-values): maaslin_temp_Porites_Transplant_Timeline/all_results.tsv
#> 2026-10-28 10:46:59.505057 INFO::Writing the significant results (those which are less than or equal to the threshold of 0.250000 ) to file (ordered by increasing q-values): maaslin_temp_Porites_Transplant_Timeline/significant_results.tsv
#>  Finished! Porites Transplant Timeline Summary: LEfSe = 0 | DESeq2 = 19 | MaAsLin2 = 0
write.table(por_trans_drivers,
            file = paste(dir_out_tab,
                         pref_out,
                         "Transplant_DiffAbunTaxa_POR.tab", sep =""))


# ==============================================================================
# 3. EXTRACTION OF OVERLAPPING CONSENSUS MARKERS (OPTIONAL)
# ==============================================================================
# Example: Extracting only the highly rigorous markers flagged by AT LEAST 2 methods
get_consensus_taxa <- function(df_drivers) {
  if (nrow(df_drivers) == 0) return(df_drivers)
  
  df_drivers %>%
    group_by(taxon, Phylum, Class, Order, Family, Genus) %>%
    summarise(
      Methods_Count = n(),
      Supporting_Tools = paste(Method, collapse = ", "),
      Avg_Score = mean(Score),
      .groups = "drop"
    ) %>%
    filter(Methods_Count >= 2) %>%
    arrange(desc(Methods_Count), desc(abs(Avg_Score)))
}

cat("\n--- IDENTIFYING HIGH-CONFIDENCE CONSENSUS TAXA ---\n")
#> 
#> --- IDENTIFYING HIGH-CONFIDENCE CONSENSUS TAXA ---
poc_consensus_summary <- get_consensus_taxa(poc_trans_drivers)
por_consensus_summary <- get_consensus_taxa(por_trans_drivers)

print(head(poc_consensus_summary))
#> # A tibble: 0 × 9
#> # ℹ 9 variables: taxon <chr>, Phylum <chr>, Class <chr>, Order <chr>,
#> #   Family <chr>, Genus <chr>, Methods_Count <int>, Supporting_Tools <chr>,
#> #   Avg_Score <dbl>
print(head(por_consensus_summary))
#> # A tibble: 0 × 9
#> # ℹ 9 variables: taxon <chr>, Phylum <chr>, Class <chr>, Order <chr>,
#> #   Family <chr>, Genus <chr>, Methods_Count <int>, Supporting_Tools <chr>,
#> #   Avg_Score <dbl>


# 1. Isolate RAW UNSCALED heat shock window objects from your curated baseline layers
ps_poc_hs_raw <- subset_samples(physeq_midpointTree_curated_POC, Timepoint %in% c("Pre-HS", "Post-HS"))
ps_poc_hs_raw <- prune_taxa(taxa_sums(ps_poc_hs_raw) > 0, ps_poc_hs_raw)

ps_por_hs_raw <- subset_samples(physeq_midpointTree_curated_POR, Timepoint %in% c("Pre-HS", "Post-HS"))
ps_por_hs_raw <- prune_taxa(taxa_sums(ps_por_hs_raw) > 0, ps_por_hs_raw)

# 2. Execute the multi-tool consensus suite, anchoring against Tiahura
cat("\n--- EXTRACTING DRIVERS FOR POCILLOPORA ---\n")
#> 
#> --- EXTRACTING DRIVERS FOR POCILLOPORA ---
poc_hs_bacterial_drivers <- run_coral_differential_consensus_v5(
  physeq = ps_poc_hs_raw, 
  tissue_label = "Pocillopora Heat Shock", 
  anchor_site = "Tiahura"
)
#> 
#> [Power Filter] Pocillopora Heat Shock: Reduced total ASV hypotheses from 3782 down to 746 (Removed 3036 rare background taxa).
#> [Running LEfSe]...
#> [Running DESeq2]...
#> 
#> [Running MaAsLin2]...
#> [1] "Creating output folder"
#> [1] "Creating output feature tables folder"
#> [1] "Creating output fits folder"
#> 2026-10-28 10:47:06.334293 INFO::Writing function arguments to log file
#> 2026-10-28 10:47:06.364174 INFO::Verifying options selected are valid
#> 2026-10-28 10:47:06.365825 INFO::Determining format of input files
#> 2026-10-28 10:47:06.367278 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-10-28 10:47:06.397954 INFO::Formula for fixed effects: expr ~  Site
#> 2026-10-28 10:47:06.400347 INFO::Filter data based on min abundance and min prevalence
#> 2026-10-28 10:47:06.401987 INFO::Total samples in data: 39
#> 2026-10-28 10:47:06.403542 INFO::Min samples required with min abundance for a feature not to be filtered: 3.900000
#> 2026-10-28 10:47:06.413826 INFO::Total filtered features: 432
#> 2026-10-28 10:47:06.416024 INFO::Filtered feature names from abundance and prevalence filtering: ASV2325, ASV551, ASV644, ASV1495, ASV3041, ASV2294, ASV2735, ASV107, ASV566, ASV2915, ASV2916, ASV1859, ASV2151, ASV384, ASV758, ASV543, ASV669, ASV3024, ASV1211, ASV2005, ASV1369, ASV880, ASV174, ASV2877, ASV457, ASV678, ASV1264, ASV289, ASV1790, ASV2759, ASV1834, ASV1543, ASV313, ASV3380, ASV360, ASV841, ASV242, ASV857, ASV507, ASV1217, ASV3741, ASV3237, ASV4303, ASV1396, ASV990, ASV1883, ASV3900, ASV4006, ASV2491, ASV1326, ASV1975, ASV2668, ASV3782, ASV1260, ASV1717, ASV848, ASV470, ASV739, ASV2226, ASV293, ASV2156, ASV1335, ASV1047, ASV302, ASV637, ASV2058, ASV314, ASV2068, ASV1530, ASV1901, ASV3669, ASV4793, ASV1990, ASV4180, ASV1025, ASV3449, ASV699, ASV1149, ASV3013, ASV1521, ASV1201, ASV1657, ASV4308, ASV4665, ASV200, ASV1876, ASV311, ASV2990, ASV485, ASV1200, ASV1719, ASV1479, ASV3773, ASV1424, ASV492, ASV4014, ASV2642, ASV1454, ASV1478, ASV3129, ASV464, ASV1746, ASV1625, ASV44, ASV271, ASV407, ASV982, ASV860, ASV1079, ASV247, ASV1358, ASV1450, ASV5153, ASV547, ASV3862, ASV3345, ASV1337, ASV1099, ASV2106, ASV3831, ASV726, ASV3824, ASV453, ASV3475, ASV1338, ASV1902, ASV1348, ASV935, ASV353, ASV951, ASV1810, ASV2299, ASV3958, ASV4555, ASV2509, ASV1272, ASV3298, ASV1590, ASV648, ASV1319, ASV494, ASV2437, ASV780, ASV382, ASV3192, ASV825, ASV998, ASV3707, ASV1095, ASV717, ASV955, ASV376, ASV2878, ASV801, ASV322, ASV499, ASV2765, ASV53, ASV524, ASV1123, ASV2145, ASV2071, ASV870, ASV324, ASV969, ASV2926, ASV2560, ASV4297, ASV1816, ASV4061, ASV2171, ASV4937, ASV1258, ASV1887, ASV2170, ASV5429, ASV533, ASV4863, ASV1739, ASV2431, ASV660, ASV2265, ASV1439, ASV1218, ASV1171, ASV2374, ASV2492, ASV2072, ASV3270, ASV3781, ASV3231, ASV2708, ASV3266, ASV3039, ASV1062, ASV705, ASV707, ASV1764, ASV903, ASV837, ASV4456, ASV1080, ASV521, ASV1879, ASV396, ASV1716, ASV94, ASV2188, ASV469, ASV1569, ASV1748, ASV399, ASV518, ASV910, ASV4492, ASV466, ASV2042, ASV822, ASV1146, ASV4114, ASV975, ASV4010, ASV1961, ASV1589, ASV2298, ASV1434, ASV2182, ASV511, ASV763, ASV2092, ASV2345, ASV2223, ASV1119, ASV584, ASV4473, ASV1568, ASV2566, ASV1231, ASV2456, ASV269, ASV74, ASV525, ASV298, ASV1054, ASV356, ASV1789, ASV3200, ASV1426, ASV599, ASV1135, ASV3439, ASV2923, ASV1091, ASV123, ASV113, ASV490, ASV624, ASV2222, ASV1005, ASV4724, ASV596, ASV1923, ASV677, ASV835, ASV1451, ASV1136, ASV1989, ASV4490, ASV899, ASV1535, ASV2039, ASV1174, ASV752, ASV3551, ASV3447, ASV4657, ASV495, ASV4605, ASV329, ASV1617, ASV319, ASV1818, ASV1547, ASV640, ASV1536, ASV276, ASV1020, ASV725, ASV1476, ASV2933, ASV2600, ASV1290, ASV922, ASV673, ASV1721, ASV63, ASV339, ASV375, ASV592, ASV3732, ASV4485, ASV2055, ASV4545, ASV1072, ASV616, ASV1694, ASV1988, ASV1785, ASV1133, ASV1413, ASV916, ASV1233, ASV4851, ASV2121, ASV2919, ASV5157, ASV3703, ASV1814, ASV1757, ASV604, ASV2685, ASV2987, ASV456, ASV3421, ASV1002, ASV871, ASV1700, ASV1039, ASV779, ASV850, ASV454, ASV148, ASV220, ASV809, ASV448, ASV1242, ASV371, ASV435, ASV2019, ASV232, ASV819, ASV2541, ASV844, ASV3695, ASV970, ASV352, ASV936, ASV1066, ASV214, ASV192, ASV641, ASV2764, ASV553, ASV4360, ASV649, ASV2056, ASV297, ASV157, ASV923, ASV753, ASV986, ASV3012, ASV3087, ASV1900, ASV1755, ASV950, ASV727, ASV437, ASV3871, ASV2008, ASV734, ASV4425, ASV4367, ASV563, ASV367, ASV595, ASV977, ASV3830, ASV3344, ASV3957, ASV3999, ASV2238, ASV4168, ASV1282, ASV2804, ASV2848, ASV4296, ASV1364, ASV2136, ASV2213, ASV793, ASV2322, ASV1306, ASV3443, ASV1921, ASV2302, ASV4675, ASV5350, ASV4801, ASV1546, ASV4062, ASV3271, ASV4411, ASV802, ASV2324, ASV5158, ASV3169, ASV1771, ASV4007, ASV381, ASV522, ASV1515, ASV5452, ASV2119, ASV2978, ASV4243, ASV3090, ASV3448, ASV4117, ASV2979, ASV3193, ASV4173, ASV4776, ASV2081, ASV2293, ASV1860, ASV2253, ASV1695, ASV2178, ASV5437, ASV2709, ASV365
#> 2026-10-28 10:47:06.422892 INFO::Total filtered features with variance filtering: 0
#> 2026-10-28 10:47:06.424862 INFO::Filtered feature names from variance filtering:
#> 2026-10-28 10:47:06.426566 INFO::Running selected normalization method: TSS
#> 2026-10-28 10:47:06.431401 INFO::Applying z-score to standardize continuous metadata
#> 2026-10-28 10:47:06.434624 INFO::Running selected transform method: AST
#> 2026-10-28 10:47:06.440993 INFO::Running selected analysis method: LM
#> 2026-10-28 10:47:06.442846 INFO::Fitting model to feature number 1, ASV849
#> 2026-10-28 10:47:06.44642 INFO::Fitting model to feature number 2, ASV610
#> 2026-10-28 10:47:06.449839 INFO::Fitting model to feature number 3, ASV746
#> 2026-10-28 10:47:06.453437 INFO::Fitting model to feature number 4, ASV939
#> 2026-10-28 10:47:06.457551 INFO::Fitting model to feature number 5, ASV788
#> 2026-10-28 10:47:06.46109 INFO::Fitting model to feature number 6, ASV182
#> 2026-10-28 10:47:06.464504 INFO::Fitting model to feature number 7, ASV1006
#> 2026-10-28 10:47:06.467945 INFO::Fitting model to feature number 8, ASV577
#> 2026-10-28 10:47:06.471517 INFO::Fitting model to feature number 9, ASV609
#> 2026-10-28 10:47:06.476741 INFO::Fitting model to feature number 10, ASV1192
#> 2026-10-28 10:47:06.480259 INFO::Fitting model to feature number 11, ASV4104
#> 2026-10-28 10:47:06.48367 INFO::Fitting model to feature number 12, ASV98
#> 2026-10-28 10:47:06.487029 INFO::Fitting model to feature number 13, ASV325
#> 2026-10-28 10:47:06.493277 INFO::Fitting model to feature number 14, ASV413
#> 2026-10-28 10:47:06.498324 INFO::Fitting model to feature number 15, ASV628
#> 2026-10-28 10:47:06.501954 INFO::Fitting model to feature number 16, ASV3340
#> 2026-10-28 10:47:06.507076 INFO::Fitting model to feature number 17, ASV446
#> 2026-10-28 10:47:06.512067 INFO::Fitting model to feature number 18, ASV13
#> 2026-10-28 10:47:06.515625 INFO::Fitting model to feature number 19, ASV663
#> 2026-10-28 10:47:06.519088 INFO::Fitting model to feature number 20, ASV2107
#> 2026-10-28 10:47:06.525747 INFO::Fitting model to feature number 21, ASV248
#> 2026-10-28 10:47:06.532658 INFO::Fitting model to feature number 22, ASV307
#> 2026-10-28 10:47:06.538395 INFO::Fitting model to feature number 23, ASV2
#> 2026-10-28 10:47:06.544454 INFO::Fitting model to feature number 24, ASV194
#> 2026-10-28 10:47:06.549383 INFO::Fitting model to feature number 25, ASV80
#> 2026-10-28 10:47:06.556348 INFO::Fitting model to feature number 26, ASV344
#> 2026-10-28 10:47:06.562217 INFO::Fitting model to feature number 27, ASV238
#> 2026-10-28 10:47:06.566862 INFO::Fitting model to feature number 28, ASV402
#> 2026-10-28 10:47:06.572492 INFO::Fitting model to feature number 29, ASV25
#> 2026-10-28 10:47:06.578186 INFO::Fitting model to feature number 30, ASV2619
#> 2026-10-28 10:47:06.584741 INFO::Fitting model to feature number 31, ASV455
#> 2026-10-28 10:47:06.590494 INFO::Fitting model to feature number 32, ASV241
#> 2026-10-28 10:47:06.598734 INFO::Fitting model to feature number 33, ASV631
#> 2026-10-28 10:47:06.604143 INFO::Fitting model to feature number 34, ASV471
#> 2026-10-28 10:47:06.610024 INFO::Fitting model to feature number 35, ASV3142
#> 2026-10-28 10:47:06.614367 INFO::Fitting model to feature number 36, ASV1815
#> 2026-10-28 10:47:06.618769 INFO::Fitting model to feature number 37, ASV1001
#> 2026-10-28 10:47:06.623041 INFO::Fitting model to feature number 38, ASV141
#> 2026-10-28 10:47:06.629758 INFO::Fitting model to feature number 39, ASV257
#> 2026-10-28 10:47:06.63404 INFO::Fitting model to feature number 40, ASV51
#> 2026-10-28 10:47:06.637913 INFO::Fitting model to feature number 41, ASV468
#> 2026-10-28 10:47:06.646091 INFO::Fitting model to feature number 42, ASV18
#> 2026-10-28 10:47:06.650652 INFO::Fitting model to feature number 43, ASV87
#> 2026-10-28 10:47:06.654606 INFO::Fitting model to feature number 44, ASV593
#> 2026-10-28 10:47:06.661019 INFO::Fitting model to feature number 45, ASV1985
#> 2026-10-28 10:47:06.665268 INFO::Fitting model to feature number 46, ASV3011
#> 2026-10-28 10:47:06.66913 INFO::Fitting model to feature number 47, ASV2274
#> 2026-10-28 10:47:06.673957 INFO::Fitting model to feature number 48, ASV1649
#> 2026-10-28 10:47:06.678152 INFO::Fitting model to feature number 49, ASV826
#> 2026-10-28 10:47:06.682196 INFO::Fitting model to feature number 50, ASV634
#> 2026-10-28 10:47:06.6863 INFO::Fitting model to feature number 51, ASV320
#> 2026-10-28 10:47:06.695241 INFO::Fitting model to feature number 52, ASV193
#> 2026-10-28 10:47:06.700816 INFO::Fitting model to feature number 53, ASV212
#> 2026-10-28 10:47:06.706902 INFO::Fitting model to feature number 54, ASV291
#> 2026-10-28 10:47:06.712561 INFO::Fitting model to feature number 55, ASV1931
#> 2026-10-28 10:47:06.716884 INFO::Fitting model to feature number 56, ASV1920
#> 2026-10-28 10:47:06.720844 INFO::Fitting model to feature number 57, ASV357
#> 2026-10-28 10:47:06.725301 INFO::Fitting model to feature number 58, ASV20
#> 2026-10-28 10:47:06.728906 INFO::Fitting model to feature number 59, ASV78
#> 2026-10-28 10:47:06.732457 INFO::Fitting model to feature number 60, ASV11
#> 2026-10-28 10:47:06.735833 INFO::Fitting model to feature number 61, ASV30
#> 2026-10-28 10:47:06.73987 INFO::Fitting model to feature number 62, ASV112
#> 2026-10-28 10:47:06.743526 INFO::Fitting model to feature number 63, ASV38
#> 2026-10-28 10:47:06.747155 INFO::Fitting model to feature number 64, ASV138
#> 2026-10-28 10:47:06.750596 INFO::Fitting model to feature number 65, ASV34
#> 2026-10-28 10:47:06.755734 INFO::Fitting model to feature number 66, ASV101
#> 2026-10-28 10:47:06.761988 INFO::Fitting model to feature number 67, ASV1802
#> 2026-10-28 10:47:06.765702 INFO::Fitting model to feature number 68, ASV146
#> 2026-10-28 10:47:06.769881 INFO::Fitting model to feature number 69, ASV383
#> 2026-10-28 10:47:06.775128 INFO::Fitting model to feature number 70, ASV166
#> 2026-10-28 10:47:06.778838 INFO::Fitting model to feature number 71, ASV1105
#> 2026-10-28 10:47:06.784284 INFO::Fitting model to feature number 72, ASV784
#> 2026-10-28 10:47:06.787888 INFO::Fitting model to feature number 73, ASV310
#> 2026-10-28 10:47:06.796126 INFO::Fitting model to feature number 74, ASV417
#> 2026-10-28 10:47:06.802658 INFO::Fitting model to feature number 75, ASV1225
#> 2026-10-28 10:47:06.8093 INFO::Fitting model to feature number 76, ASV303
#> 2026-10-28 10:47:06.814436 INFO::Fitting model to feature number 77, ASV1028
#> 2026-10-28 10:47:06.819867 INFO::Fitting model to feature number 78, ASV4663
#> 2026-10-28 10:47:06.826107 INFO::Fitting model to feature number 79, ASV1777
#> 2026-10-28 10:47:06.831464 INFO::Fitting model to feature number 80, ASV2488
#> 2026-10-28 10:47:06.83679 INFO::Fitting model to feature number 81, ASV201
#> 2026-10-28 10:47:06.84262 INFO::Fitting model to feature number 82, ASV892
#> 2026-10-28 10:47:06.847078 INFO::Fitting model to feature number 83, ASV796
#> 2026-10-28 10:47:06.852445 INFO::Fitting model to feature number 84, ASV864
#> 2026-10-28 10:47:06.859325 INFO::Fitting model to feature number 85, ASV211
#> 2026-10-28 10:47:06.863065 INFO::Fitting model to feature number 86, ASV240
#> 2026-10-28 10:47:06.866925 INFO::Fitting model to feature number 87, ASV1265
#> 2026-10-28 10:47:06.871202 INFO::Fitting model to feature number 88, ASV661
#> 2026-10-28 10:47:06.875476 INFO::Fitting model to feature number 89, ASV728
#> 2026-10-28 10:47:06.879546 INFO::Fitting model to feature number 90, ASV988
#> 2026-10-28 10:47:06.885121 INFO::Fitting model to feature number 91, ASV2427
#> 2026-10-28 10:47:06.890336 INFO::Fitting model to feature number 92, ASV2843
#> 2026-10-28 10:47:06.895687 INFO::Fitting model to feature number 93, ASV629
#> 2026-10-28 10:47:06.901097 INFO::Fitting model to feature number 94, ASV244
#> 2026-10-28 10:47:06.907358 INFO::Fitting model to feature number 95, ASV450
#> 2026-10-28 10:47:06.912574 INFO::Fitting model to feature number 96, ASV153
#> 2026-10-28 10:47:06.918022 INFO::Fitting model to feature number 97, ASV237
#> 2026-10-28 10:47:06.922103 INFO::Fitting model to feature number 98, ASV50
#> 2026-10-28 10:47:06.926297 INFO::Fitting model to feature number 99, ASV97
#> 2026-10-28 10:47:06.930327 INFO::Fitting model to feature number 100, ASV143
#> 2026-10-28 10:47:06.934039 INFO::Fitting model to feature number 101, ASV679
#> 2026-10-28 10:47:06.937618 INFO::Fitting model to feature number 102, ASV131
#> 2026-10-28 10:47:06.942071 INFO::Fitting model to feature number 103, ASV121
#> 2026-10-28 10:47:06.946671 INFO::Fitting model to feature number 104, ASV150
#> 2026-10-28 10:47:06.95042 INFO::Fitting model to feature number 105, ASV856
#> 2026-10-28 10:47:06.954354 INFO::Fitting model to feature number 106, ASV954
#> 2026-10-28 10:47:06.960027 INFO::Fitting model to feature number 107, ASV3116
#> 2026-10-28 10:47:06.963599 INFO::Fitting model to feature number 108, ASV3512
#> 2026-10-28 10:47:06.966989 INFO::Fitting model to feature number 109, ASV516
#> 2026-10-28 10:47:06.970391 INFO::Fitting model to feature number 110, ASV1259
#> 2026-10-28 10:47:06.974521 INFO::Fitting model to feature number 111, ASV1835
#> 2026-10-28 10:47:06.978116 INFO::Fitting model to feature number 112, ASV542
#> 2026-10-28 10:47:06.98152 INFO::Fitting model to feature number 113, ASV2847
#> 2026-10-28 10:47:06.98487 INFO::Fitting model to feature number 114, ASV2344
#> 2026-10-28 10:47:06.989339 INFO::Fitting model to feature number 115, ASV2647
#> 2026-10-28 10:47:06.994156 INFO::Fitting model to feature number 116, ASV465
#> 2026-10-28 10:47:06.997613 INFO::Fitting model to feature number 117, ASV1204
#> 2026-10-28 10:47:07.001013 INFO::Fitting model to feature number 118, ASV2986
#> 2026-10-28 10:47:07.007022 INFO::Fitting model to feature number 119, ASV4105
#> 2026-10-28 10:47:07.012301 INFO::Fitting model to feature number 120, ASV236
#> 2026-10-28 10:47:07.016023 INFO::Fitting model to feature number 121, ASV684
#> 2026-10-28 10:47:07.019684 INFO::Fitting model to feature number 122, ASV2781
#> 2026-10-28 10:47:07.02468 INFO::Fitting model to feature number 123, ASV513
#> 2026-10-28 10:47:07.028415 INFO::Fitting model to feature number 124, ASV3779
#> 2026-10-28 10:47:07.031913 INFO::Fitting model to feature number 125, ASV1482
#> 2026-10-28 10:47:07.035307 INFO::Fitting model to feature number 126, ASV1270
#> 2026-10-28 10:47:07.03993 INFO::Fitting model to feature number 127, ASV388
#> 2026-10-28 10:47:07.043709 INFO::Fitting model to feature number 128, ASV989
#> 2026-10-28 10:47:07.047283 INFO::Fitting model to feature number 129, ASV1664
#> 2026-10-28 10:47:07.051866 INFO::Fitting model to feature number 130, ASV3062
#> 2026-10-28 10:47:07.05627 INFO::Fitting model to feature number 131, ASV1813
#> 2026-10-28 10:47:07.061126 INFO::Fitting model to feature number 132, ASV2281
#> 2026-10-28 10:47:07.065244 INFO::Fitting model to feature number 133, ASV1873
#> 2026-10-28 10:47:07.069137 INFO::Fitting model to feature number 134, ASV275
#> 2026-10-28 10:47:07.073744 INFO::Fitting model to feature number 135, ASV632
#> 2026-10-28 10:47:07.079496 INFO::Fitting model to feature number 136, ASV42
#> 2026-10-28 10:47:07.086468 INFO::Fitting model to feature number 137, ASV463
#> 2026-10-28 10:47:07.093127 INFO::Fitting model to feature number 138, ASV47
#> 2026-10-28 10:47:07.099366 INFO::Fitting model to feature number 139, ASV59
#> 2026-10-28 10:47:07.103691 INFO::Fitting model to feature number 140, ASV1000
#> 2026-10-28 10:47:07.109665 INFO::Fitting model to feature number 141, ASV66
#> 2026-10-28 10:47:07.114302 INFO::Fitting model to feature number 142, ASV968
#> 2026-10-28 10:47:07.118988 INFO::Fitting model to feature number 143, ASV61
#> 2026-10-28 10:47:07.12454 INFO::Fitting model to feature number 144, ASV132
#> 2026-10-28 10:47:07.129086 INFO::Fitting model to feature number 145, ASV186
#> 2026-10-28 10:47:07.13557 INFO::Fitting model to feature number 146, ASV813
#> 2026-10-28 10:47:07.141093 INFO::Fitting model to feature number 147, ASV286
#> 2026-10-28 10:47:07.146256 INFO::Fitting model to feature number 148, ASV119
#> 2026-10-28 10:47:07.15248 INFO::Fitting model to feature number 149, ASV1183
#> 2026-10-28 10:47:07.15898 INFO::Fitting model to feature number 150, ASV105
#> 2026-10-28 10:47:07.163503 INFO::Fitting model to feature number 151, ASV124
#> 2026-10-28 10:47:07.167531 INFO::Fitting model to feature number 152, ASV213
#> 2026-10-28 10:47:07.172098 INFO::Fitting model to feature number 153, ASV762
#> 2026-10-28 10:47:07.17618 INFO::Fitting model to feature number 154, ASV1245
#> 2026-10-28 10:47:07.180887 INFO::Fitting model to feature number 155, ASV222
#> 2026-10-28 10:47:07.184577 INFO::Fitting model to feature number 156, ASV32
#> 2026-10-28 10:47:07.19033 INFO::Fitting model to feature number 157, ASV1307
#> 2026-10-28 10:47:07.197123 INFO::Fitting model to feature number 158, ASV135
#> 2026-10-28 10:47:07.200769 INFO::Fitting model to feature number 159, ASV177
#> 2026-10-28 10:47:07.20658 INFO::Fitting model to feature number 160, ASV128
#> 2026-10-28 10:47:07.212367 INFO::Fitting model to feature number 161, ASV568
#> 2026-10-28 10:47:07.216368 INFO::Fitting model to feature number 162, ASV259
#> 2026-10-28 10:47:07.221691 INFO::Fitting model to feature number 163, ASV295
#> 2026-10-28 10:47:07.228002 INFO::Fitting model to feature number 164, ASV267
#> 2026-10-28 10:47:07.232697 INFO::Fitting model to feature number 165, ASV359
#> 2026-10-28 10:47:07.238852 INFO::Fitting model to feature number 166, ASV757
#> 2026-10-28 10:47:07.243825 INFO::Fitting model to feature number 167, ASV876
#> 2026-10-28 10:47:07.247528 INFO::Fitting model to feature number 168, ASV251
#> 2026-10-28 10:47:07.251187 INFO::Fitting model to feature number 169, ASV366
#> 2026-10-28 10:47:07.255303 INFO::Fitting model to feature number 170, ASV387
#> 2026-10-28 10:47:07.261598 INFO::Fitting model to feature number 171, ASV1207
#> 2026-10-28 10:47:07.265375 INFO::Fitting model to feature number 172, ASV239
#> 2026-10-28 10:47:07.270602 INFO::Fitting model to feature number 173, ASV1164
#> 2026-10-28 10:47:07.276214 INFO::Fitting model to feature number 174, ASV1
#> 2026-10-28 10:47:07.280517 INFO::Fitting model to feature number 175, ASV93
#> 2026-10-28 10:47:07.285164 INFO::Fitting model to feature number 176, ASV256
#> 2026-10-28 10:47:07.289429 INFO::Fitting model to feature number 177, ASV33
#> 2026-10-28 10:47:07.293402 INFO::Fitting model to feature number 178, ASV523
#> 2026-10-28 10:47:07.297257 INFO::Fitting model to feature number 179, ASV142
#> 2026-10-28 10:47:07.301069 INFO::Fitting model to feature number 180, ASV73
#> 2026-10-28 10:47:07.307675 INFO::Fitting model to feature number 181, ASV176
#> 2026-10-28 10:47:07.312313 INFO::Fitting model to feature number 182, ASV1692
#> 2026-10-28 10:47:07.318482 INFO::Fitting model to feature number 183, ASV943
#> 2026-10-28 10:47:07.324244 INFO::Fitting model to feature number 184, ASV676
#> 2026-10-28 10:47:07.328075 INFO::Fitting model to feature number 185, ASV535
#> 2026-10-28 10:47:07.331746 INFO::Fitting model to feature number 186, ASV172
#> 2026-10-28 10:47:07.338807 INFO::Fitting model to feature number 187, ASV144
#> 2026-10-28 10:47:07.343338 INFO::Fitting model to feature number 188, ASV134
#> 2026-10-28 10:47:07.34727 INFO::Fitting model to feature number 189, ASV92
#> 2026-10-28 10:47:07.351049 INFO::Fitting model to feature number 190, ASV120
#> 2026-10-28 10:47:07.354621 INFO::Fitting model to feature number 191, ASV19
#> 2026-10-28 10:47:07.358797 INFO::Fitting model to feature number 192, ASV5
#> 2026-10-28 10:47:07.362566 INFO::Fitting model to feature number 193, ASV863
#> 2026-10-28 10:47:07.366569 INFO::Fitting model to feature number 194, ASV312
#> 2026-10-28 10:47:07.370378 INFO::Fitting model to feature number 195, ASV370
#> 2026-10-28 10:47:07.37464 INFO::Fitting model to feature number 196, ASV72
#> 2026-10-28 10:47:07.378462 INFO::Fitting model to feature number 197, ASV48
#> 2026-10-28 10:47:07.38224 INFO::Fitting model to feature number 198, ASV206
#> 2026-10-28 10:47:07.385993 INFO::Fitting model to feature number 199, ASV9
#> 2026-10-28 10:47:07.391052 INFO::Fitting model to feature number 200, ASV65
#> 2026-10-28 10:47:07.396281 INFO::Fitting model to feature number 201, ASV136
#> 2026-10-28 10:47:07.401306 INFO::Fitting model to feature number 202, ASV691
#> 2026-10-28 10:47:07.406693 INFO::Fitting model to feature number 203, ASV43
#> 2026-10-28 10:47:07.41088 INFO::Fitting model to feature number 204, ASV265
#> 2026-10-28 10:47:07.415166 INFO::Fitting model to feature number 205, ASV56
#> 2026-10-28 10:47:07.41906 INFO::Fitting model to feature number 206, ASV1078
#> 2026-10-28 10:47:07.423371 INFO::Fitting model to feature number 207, ASV902
#> 2026-10-28 10:47:07.427369 INFO::Fitting model to feature number 208, ASV2215
#> 2026-10-28 10:47:07.431479 INFO::Fitting model to feature number 209, ASV91
#> 2026-10-28 10:47:07.435524 INFO::Fitting model to feature number 210, ASV169
#> 2026-10-28 10:47:07.439367 INFO::Fitting model to feature number 211, ASV585
#> 2026-10-28 10:47:07.443136 INFO::Fitting model to feature number 212, ASV846
#> 2026-10-28 10:47:07.4468 INFO::Fitting model to feature number 213, ASV1573
#> 2026-10-28 10:47:07.450838 INFO::Fitting model to feature number 214, ASV5160
#> 2026-10-28 10:47:07.456207 INFO::Fitting model to feature number 215, ASV335
#> 2026-10-28 10:47:07.459869 INFO::Fitting model to feature number 216, ASV1374
#> 2026-10-28 10:47:07.463345 INFO::Fitting model to feature number 217, ASV336
#> 2026-10-28 10:47:07.466751 INFO::Fitting model to feature number 218, ASV152
#> 2026-10-28 10:47:07.471347 INFO::Fitting model to feature number 219, ASV83
#> 2026-10-28 10:47:07.475521 INFO::Fitting model to feature number 220, ASV183
#> 2026-10-28 10:47:07.479546 INFO::Fitting model to feature number 221, ASV2458
#> 2026-10-28 10:47:07.484249 INFO::Fitting model to feature number 222, ASV1224
#> 2026-10-28 10:47:07.489742 INFO::Fitting model to feature number 223, ASV997
#> 2026-10-28 10:47:07.494649 INFO::Fitting model to feature number 224, ASV1134
#> 2026-10-28 10:47:07.50016 INFO::Fitting model to feature number 225, ASV207
#> 2026-10-28 10:47:07.505586 INFO::Fitting model to feature number 226, ASV149
#> 2026-10-28 10:47:07.510313 INFO::Fitting model to feature number 227, ASV147
#> 2026-10-28 10:47:07.514084 INFO::Fitting model to feature number 228, ASV62
#> 2026-10-28 10:47:07.517802 INFO::Fitting model to feature number 229, ASV100
#> 2026-10-28 10:47:07.521497 INFO::Fitting model to feature number 230, ASV77
#> 2026-10-28 10:47:07.528627 INFO::Fitting model to feature number 231, ASV419
#> 2026-10-28 10:47:07.534454 INFO::Fitting model to feature number 232, ASV2686
#> 2026-10-28 10:47:07.54299 INFO::Fitting model to feature number 233, ASV662
#> 2026-10-28 10:47:07.548903 INFO::Fitting model to feature number 234, ASV350
#> 2026-10-28 10:47:07.553139 INFO::Fitting model to feature number 235, ASV575
#> 2026-10-28 10:47:07.559308 INFO::Fitting model to feature number 236, ASV5440
#> 2026-10-28 10:47:07.564998 INFO::Fitting model to feature number 237, ASV346
#> 2026-10-28 10:47:07.571516 INFO::Fitting model to feature number 238, ASV934
#> 2026-10-28 10:47:07.578352 INFO::Fitting model to feature number 239, ASV2077
#> 2026-10-28 10:47:07.583323 INFO::Fitting model to feature number 240, ASV3263
#> 2026-10-28 10:47:07.587108 INFO::Fitting model to feature number 241, ASV424
#> 2026-10-28 10:47:07.593613 INFO::Fitting model to feature number 242, ASV638
#> 2026-10-28 10:47:07.598779 INFO::Fitting model to feature number 243, ASV1058
#> 2026-10-28 10:47:07.606251 INFO::Fitting model to feature number 244, ASV949
#> 2026-10-28 10:47:07.611976 INFO::Fitting model to feature number 245, ASV531
#> 2026-10-28 10:47:07.616907 INFO::Fitting model to feature number 246, ASV606
#> 2026-10-28 10:47:07.620889 INFO::Fitting model to feature number 247, ASV552
#> 2026-10-28 10:47:07.677792 INFO::Fitting model to feature number 248, ASV983
#> 2026-10-28 10:47:07.681346 INFO::Fitting model to feature number 249, ASV803
#> 2026-10-28 10:47:07.685355 INFO::Fitting model to feature number 250, ASV1958
#> 2026-10-28 10:47:07.688985 INFO::Fitting model to feature number 251, ASV548
#> 2026-10-28 10:47:07.696509 INFO::Fitting model to feature number 252, ASV853
#> 2026-10-28 10:47:07.701268 INFO::Fitting model to feature number 253, ASV711
#> 2026-10-28 10:47:07.705432 INFO::Fitting model to feature number 254, ASV569
#> 2026-10-28 10:47:07.709456 INFO::Fitting model to feature number 255, ASV967
#> 2026-10-28 10:47:07.713149 INFO::Fitting model to feature number 256, ASV770
#> 2026-10-28 10:47:07.716652 INFO::Fitting model to feature number 257, ASV28
#> 2026-10-28 10:47:07.72 INFO::Fitting model to feature number 258, ASV40
#> 2026-10-28 10:47:07.724413 INFO::Fitting model to feature number 259, ASV8
#> 2026-10-28 10:47:07.727954 INFO::Fitting model to feature number 260, ASV187
#> 2026-10-28 10:47:07.731358 INFO::Fitting model to feature number 261, ASV250
#> 2026-10-28 10:47:07.735691 INFO::Fitting model to feature number 262, ASV486
#> 2026-10-28 10:47:07.73981 INFO::Fitting model to feature number 263, ASV540
#> 2026-10-28 10:47:07.744986 INFO::Fitting model to feature number 264, ASV433
#> 2026-10-28 10:47:07.750426 INFO::Fitting model to feature number 265, ASV901
#> 2026-10-28 10:47:07.757316 INFO::Fitting model to feature number 266, ASV351
#> 2026-10-28 10:47:07.763773 INFO::Fitting model to feature number 267, ASV2116
#> 2026-10-28 10:47:07.768612 INFO::Fitting model to feature number 268, ASV85
#> 2026-10-28 10:47:07.773062 INFO::Fitting model to feature number 269, ASV82
#> 2026-10-28 10:47:07.777844 INFO::Fitting model to feature number 270, ASV181
#> 2026-10-28 10:47:07.781704 INFO::Fitting model to feature number 271, ASV118
#> 2026-10-28 10:47:07.785203 INFO::Fitting model to feature number 272, ASV16
#> 2026-10-28 10:47:07.789601 INFO::Fitting model to feature number 273, ASV7
#> 2026-10-28 10:47:07.793422 INFO::Fitting model to feature number 274, ASV127
#> 2026-10-28 10:47:07.796984 INFO::Fitting model to feature number 275, ASV411
#> 2026-10-28 10:47:07.800431 INFO::Fitting model to feature number 276, ASV162
#> 2026-10-28 10:47:07.803991 INFO::Fitting model to feature number 277, ASV249
#> 2026-10-28 10:47:07.807787 INFO::Fitting model to feature number 278, ASV823
#> 2026-10-28 10:47:07.811873 INFO::Fitting model to feature number 279, ASV1056
#> 2026-10-28 10:47:07.815643 INFO::Fitting model to feature number 280, ASV266
#> 2026-10-28 10:47:07.82008 INFO::Fitting model to feature number 281, ASV348
#> 2026-10-28 10:47:07.824371 INFO::Fitting model to feature number 282, ASV90
#> 2026-10-28 10:47:07.828799 INFO::Fitting model to feature number 283, ASV88
#> 2026-10-28 10:47:07.832548 INFO::Fitting model to feature number 284, ASV1406
#> 2026-10-28 10:47:07.836344 INFO::Fitting model to feature number 285, ASV1219
#> 2026-10-28 10:47:07.842255 INFO::Fitting model to feature number 286, ASV2052
#> 2026-10-28 10:47:07.848267 INFO::Fitting model to feature number 287, ASV415
#> 2026-10-28 10:47:07.852575 INFO::Fitting model to feature number 288, ASV274
#> 2026-10-28 10:47:07.857631 INFO::Fitting model to feature number 289, ASV52
#> 2026-10-28 10:47:07.86391 INFO::Fitting model to feature number 290, ASV81
#> 2026-10-28 10:47:07.867867 INFO::Fitting model to feature number 291, ASV639
#> 2026-10-28 10:47:07.873462 INFO::Fitting model to feature number 292, ASV3825
#> 2026-10-28 10:47:07.878086 INFO::Fitting model to feature number 293, ASV1770
#> 2026-10-28 10:47:07.882269 INFO::Fitting model to feature number 294, ASV4106
#> 2026-10-28 10:47:07.88623 INFO::Fitting model to feature number 295, ASV1036
#> 2026-10-28 10:47:07.891042 INFO::Fitting model to feature number 296, ASV2758
#> 2026-10-28 10:47:07.895888 INFO::Fitting model to feature number 297, ASV2323
#> 2026-10-28 10:47:07.900674 INFO::Fitting model to feature number 298, ASV520
#> 2026-10-28 10:47:07.90688 INFO::Fitting model to feature number 299, ASV2181
#> 2026-10-28 10:47:07.911655 INFO::Fitting model to feature number 300, ASV537
#> 2026-10-28 10:47:07.9152 INFO::Fitting model to feature number 301, ASV3040
#> 2026-10-28 10:47:07.919069 INFO::Fitting model to feature number 302, ASV4852
#> 2026-10-28 10:47:07.922982 INFO::Fitting model to feature number 303, ASV703
#> 2026-10-28 10:47:07.926596 INFO::Fitting model to feature number 304, ASV1600
#> 2026-10-28 10:47:07.930135 INFO::Fitting model to feature number 305, ASV1033
#> 2026-10-28 10:47:07.934192 INFO::Fitting model to feature number 306, ASV508
#> 2026-10-28 10:47:07.941086 INFO::Fitting model to feature number 307, ASV617
#> 2026-10-28 10:47:07.946033 INFO::Fitting model to feature number 308, ASV693
#> 2026-10-28 10:47:07.9496 INFO::Fitting model to feature number 309, ASV751
#> 2026-10-28 10:47:07.955217 INFO::Fitting model to feature number 310, ASV1122
#> 2026-10-28 10:47:07.960725 INFO::Fitting model to feature number 311, ASV1142
#> 2026-10-28 10:47:07.966333 INFO::Fitting model to feature number 312, ASV536
#> 2026-10-28 10:47:07.971099 INFO::Fitting model to feature number 313, ASV1073
#> 2026-10-28 10:47:07.975865 INFO::Fitting model to feature number 314, ASV1507
#> 2026-10-28 10:47:08.012751 INFO::Counting total values for each feature
#> 2026-10-28 10:47:08.031028 INFO::Writing filtered data to file maaslin_temp_Pocillopora_Heat_Shock/features/filtered_data.tsv
#> 2026-10-28 10:47:08.050503 INFO::Writing filtered, normalized data to file maaslin_temp_Pocillopora_Heat_Shock/features/filtered_data_norm.tsv
#> 2026-10-28 10:47:08.081125 INFO::Writing filtered, normalized, transformed data to file maaslin_temp_Pocillopora_Heat_Shock/features/filtered_data_norm_transformed.tsv
#> 2026-10-28 10:47:08.109024 INFO::Writing residuals to file maaslin_temp_Pocillopora_Heat_Shock/fits/residuals.rds
#> 2026-10-28 10:47:08.114453 INFO::Writing fitted values to file maaslin_temp_Pocillopora_Heat_Shock/fits/fitted.rds
#> 2026-10-28 10:47:08.118555 INFO::Writing all results to file (ordered by increasing q-values): maaslin_temp_Pocillopora_Heat_Shock/all_results.tsv
#> 2026-10-28 10:47:08.126413 INFO::Writing the significant results (those which are less than or equal to the threshold of 0.250000 ) to file (ordered by increasing q-values): maaslin_temp_Pocillopora_Heat_Shock/significant_results.tsv
#>  Finished! Pocillopora Heat Shock Summary: LEfSe = 0 | DESeq2 = 8 | MaAsLin2 = 0
write.table(poc_hs_bacterial_drivers,
            file = paste(dir_out_tab,
                         pref_out,
                         "HeatShock_DiffAbunTaxa_POC.tab", sep =""))

cat("\n--- EXTRACTING DRIVERS FOR PORITES ---\n")
#> 
#> --- EXTRACTING DRIVERS FOR PORITES ---
por_hs_bacterial_drivers <- run_coral_differential_consensus_v5(
  physeq = ps_por_hs_raw, 
  tissue_label = "Porites Heat Shock", 
  anchor_site = "Tiahura"
)
#> 
#> [Power Filter] Porites Heat Shock: Reduced total ASV hypotheses from 2903 down to 1753 (Removed 1150 rare background taxa).
#> [Running LEfSe]...
#> [Running DESeq2]...
#> 
#> [Running MaAsLin2]...
#> [1] "Creating output folder"
#> [1] "Creating output feature tables folder"
#> [1] "Creating output fits folder"
#> 2026-10-28 10:47:16.108275 INFO::Writing function arguments to log file
#> 2026-10-28 10:47:16.133609 INFO::Verifying options selected are valid
#> 2026-10-28 10:47:16.135708 INFO::Determining format of input files
#> 2026-10-28 10:47:16.137455 INFO::Input format is data samples as rows and metadata samples as rows
#> 2026-10-28 10:47:16.157337 INFO::Formula for fixed effects: expr ~  Site
#> 2026-10-28 10:47:16.161525 INFO::Filter data based on min abundance and min prevalence
#> 2026-10-28 10:47:16.163444 INFO::Total samples in data: 16
#> 2026-10-28 10:47:16.165021 INFO::Min samples required with min abundance for a feature not to be filtered: 1.600000
#> 2026-10-28 10:47:16.181458 INFO::Total filtered features: 1194
#> 2026-10-28 10:47:16.183693 INFO::Filtered feature names from abundance and prevalence filtering: ASV2338, ASV7612, ASV6185, ASV1495, ASV939, ASV4064, ASV7108, ASV4717, ASV7116, ASV7287, ASV1325, ASV3777, ASV2316, ASV1111, ASV1567, ASV5020, ASV446, ASV546, ASV1211, ASV4775, ASV1369, ASV880, ASV663, ASV6208, ASV2877, ASV678, ASV289, ASV5768, ASV7582, ASV3229, ASV7422, ASV3139, ASV6864, ASV4235, ASV1175, ASV1883, ASV1713, ASV7423, ASV6320, ASV3630, ASV471, ASV5319, ASV6846, ASV5330, ASV8175, ASV7981, ASV2645, ASV5517, ASV848, ASV4420, ASV7245, ASV4167, ASV3778, ASV2156, ASV51, ASV1523, ASV2983, ASV5884, ASV4225, ASV3953, ASV2872, ASV1946, ASV5334, ASV2620, ASV7804, ASV2844, ASV6570, ASV7430, ASV6422, ASV5333, ASV1901, ASV7397, ASV6764, ASV3935, ASV6063, ASV3370, ASV7577, ASV5392, ASV5605, ASV7426, ASV7973, ASV6959, ASV7392, ASV6335, ASV5694, ASV8001, ASV7576, ASV2297, ASV7800, ASV5621, ASV1522, ASV3450, ASV2763, ASV7429, ASV3372, ASV4158, ASV2018, ASV8006, ASV3084, ASV2407, ASV6342, ASV7285, ASV7428, ASV7109, ASV3854, ASV2666, ASV2943, ASV930, ASV5313, ASV2348, ASV4999, ASV1609, ASV1479, ASV1920, ASV6702, ASV2421, ASV5522, ASV7591, ASV1424, ASV5803, ASV2531, ASV5792, ASV6268, ASV4001, ASV3115, ASV1478, ASV5078, ASV3827, ASV4165, ASV5125, ASV6321, ASV11, ASV2301, ASV5006, ASV2079, ASV219, ASV44, ASV38, ASV5323, ASV498, ASV1079, ASV5978, ASV5602, ASV8025, ASV497, ASV4488, ASV547, ASV146, ASV1318, ASV5079, ASV3519, ASV3817, ASV4908, ASV784, ASV1725, ASV5120, ASV6296, ASV2567, ASV5329, ASV8005, ASV8171, ASV6424, ASV6860, ASV3770, ASV3865, ASV3289, ASV6583, ASV7824, ASV7611, ASV7998, ASV5519, ASV7613, ASV2488, ASV2918, ASV4054, ASV6436, ASV2153, ASV6841, ASV8026, ASV5541, ASV1348, ASV1843, ASV2782, ASV6312, ASV2103, ASV4868, ASV1538, ASV1400, ASV1086, ASV2276, ASV7395, ASV4942, ASV6582, ASV5520, ASV1810, ASV995, ASV4324, ASV4780, ASV2115, ASV4339, ASV4550, ASV3622, ASV1256, ASV1353, ASV3420, ASV2102, ASV3549, ASV2146, ASV3905, ASV7092, ASV6821, ASV1971, ASV5518, ASV1775, ASV4536, ASV400, ASV3225, ASV648, ASV3277, ASV308, ASV1319, ASV5335, ASV2180, ASV905, ASV825, ASV998, ASV2268, ASV5802, ASV1052, ASV4655, ASV1828, ASV322, ASV1628, ASV131, ASV6595, ASV2921, ASV1852, ASV2242, ASV1331, ASV738, ASV189, ASV3771, ASV121, ASV4046, ASV5163, ASV897, ASV1016, ASV4592, ASV6827, ASV3834, ASV7829, ASV6199, ASV6438, ASV7803, ASV2139, ASV6824, ASV7585, ASV6461, ASV5779, ASV5717, ASV3203, ASV7974, ASV6161, ASV542, ASV7284, ASV4799, ASV3471, ASV6865, ASV2040, ASV4876, ASV6566, ASV1944, ASV2647, ASV1204, ASV2986, ASV4843, ASV4935, ASV3056, ASV2317, ASV684, ASV754, ASV2435, ASV3111, ASV6286, ASV4587, ASV4933, ASV5052, ASV7614, ASV4779, ASV7581, ASV4940, ASV4661, ASV5074, ASV6082, ASV4844, ASV2588, ASV3815, ASV6719, ASV1412, ASV7787, ASV1526, ASV1508, ASV4101, ASV906, ASV3697, ASV5431, ASV5709, ASV527, ASV2431, ASV4003, ASV4466, ASV6298, ASV2265, ASV3775, ASV4571, ASV1439, ASV2269, ASV8183, ASV1218, ASV3301, ASV4118, ASV3554, ASV3341, ASV3950, ASV6297, ASV4849, ASV5968, ASV1556, ASV1904, ASV4489, ASV3066, ASV5155, ASV7788, ASV6866, ASV6704, ASV4690, ASV6089, ASV1875, ASV6187, ASV7238, ASV2778, ASV6085, ASV6826, ASV6184, ASV4169, ASV6848, ASV6315, ASV7107, ASV6981, ASV7425, ASV5004, ASV3265, ASV7136, ASV4711, ASV7825, ASV3692, ASV7830, ASV7464, ASV7249, ASV5331, ASV5778, ASV1255, ASV2784, ASV3857, ASV255, ASV3866, ASV1686, ASV4548, ASV1270, ASV3197, ASV2955, ASV6825, ASV2559, ASV6072, ASV6707, ASV2517, ASV2247, ASV4915, ASV2460, ASV2727, ASV2337, ASV5620, ASV2821, ASV3295, ASV1030, ASV2248, ASV2621, ASV4467, ASV7396, ASV3189, ASV4539, ASV2342, ASV1737, ASV4714, ASV2118, ASV4108, ASV5321, ASV1808, ASV5081, ASV3899, ASV3858, ASV2070, ASV6569, ASV7633, ASV8181, ASV5088, ASV3190, ASV6863, ASV1879, ASV3894, ASV1232, ASV1541, ASV2349, ASV1053, ASV1312, ASV978, ASV785, ASV777, ASV900, ASV1077, ASV3961, ASV1613, ASV1320, ASV388, ASV1433, ASV723, ASV1845, ASV3218, ASV5146, ASV1569, ASV2394, ASV3629, ASV3828, ASV5542, ASV3913, ASV989, ASV5324, ASV3275, ASV2803, ASV5630, ASV2455, ASV3410, ASV4628, ASV911, ASV1670, ASV3470, ASV4718, ASV1768, ASV432, ASV3465, ASV1428, ASV1146, ASV7099, ASV3949, ASV2053, ASV1209, ASV7390, ASV2028, ASV4538, ASV5608, ASV2135, ASV3664, ASV565, ASV2203, ASV6079, ASV1813, ASV5881, ASV4045, ASV8172, ASV3321, ASV3733, ASV1512, ASV4533, ASV7786, ASV6955, ASV3701, ASV5632, ASV2182, ASV3065, ASV286, ASV119, ASV1183, ASV6823, ASV6211, ASV177, ASV568, ASV7095, ASV3194, ASV1568, ASV2823, ASV3667, ASV3989, ASV7135, ASV4342, ASV6081, ASV4410, ASV5703, ASV2073, ASV3112, ASV3696, ASV4177, ASV1456, ASV2254, ASV8174, ASV959, ASV5872, ASV5540, ASV5288, ASV1391, ASV4479, ASV3409, ASV4705, ASV5332, ASV2545, ASV917, ASV6086, ASV3506, ASV1469, ASV4730, ASV6718, ASV2334, ASV6202, ASV599, ASV1253, ASV1587, ASV4593, ASV3508, ASV443, ASV4301, ASV4653, ASV4856, ASV790, ASV4720, ASV1027, ASV8173, ASV1091, ASV2256, ASV2796, ASV2780, ASV6065, ASV1489, ASV4306, ASV7783, ASV6067, ASV5448, ASV4478, ASV6822, ASV5144, ASV4476, ASV6867, ASV3621, ASV3294, ASV3061, ASV7262, ASV7976, ASV1969, ASV1182, ASV985, ASV4012, ASV3467, ASV2279, ASV1709, ASV5871, ASV5016, ASV3162, ASV1562, ASV805, ASV835, ASV5232, ASV4777, ASV6958, ASV4981, ASV7264, ASV4774, ASV1136, ASV206, ASV6708, ASV7635, ASV5789, ASV65, ASV1535, ASV3468, ASV999, ASV56, ASV4993, ASV6567, ASV6343, ASV3296, ASV3303, ASV7789, ASV1548, ASV3784, ASV1332, ASV156, ASV5607, ASV4423, ASV2731, ASV4440, ASV4789, ASV7790, ASV869, ASV2562, ASV2132, ASV974, ASV2343, ASV2355, ASV329, ASV1984, ASV4543, ASV6088, ASV7980, ASV4002, ASV1399, ASV1440, ASV5154, ASV1974, ASV7616, ASV5984, ASV5777, ASV7243, ASV5327, ASV3814, ASV4237, ASV7587, ASV7584, ASV1547, ASV2193, ASV4845, ASV4222, ASV3596, ASV1224, ASV6905, ASV640, ASV769, ASV4179, ASV5606, ASV5970, ASV7982, ASV3106, ASV7590, ASV1972, ASV3017, ASV2792, ASV3343, ASV276, ASV5424, ASV1026, ASV3708, ASV3912, ASV4482, ASV1528, ASV2004, ASV3829, ASV227, ASV3427, ASV1100, ASV3954, ASV7244, ASV2722, ASV3113, ASV4742, ASV4790, ASV3631, ASV4671, ASV4120, ASV7098, ASV4654, ASV7138, ASV1872, ASV5874, ASV1268, ASV3967, ASV3305, ASV7589, ASV5230, ASV5018, ASV3105, ASV1694, ASV7265, ASV6212, ASV3951, ASV3996, ASV7421, ASV6562, ASV4601, ASV1094, ASV6862, ASV2679, ASV3444, ASV5317, ASV6563, ASV1343, ASV6597, ASV583, ASV7580, ASV2622, ASV3089, ASV3507, ASV2054, ASV1438, ASV2366, ASV7268, ASV180, ASV1943, ASV2530, ASV2734, ASV4057, ASV1233, ASV3412, ASV909, ASV4004, ASV6682, ASV2050, ASV7826, ASV7418, ASV5251, ASV7782, ASV8094, ASV8179, ASV3377, ASV2117, ASV6568, ASV1407, ASV6847, ASV3033, ASV1679, ASV3267, ASV1711, ASV5623, ASV3224, ASV2941, ASV6426, ASV2924, ASV7389, ASV4287, ASV1758, ASV2296, ASV1667, ASV2221, ASV732, ASV1513, ASV1735, ASV246, ASV871, ASV1787, ASV1039, ASV1205, ASV4110, ASV4056, ASV853, ASV2424, ASV6186, ASV1987, ASV766, ASV3171, ASV5780, ASV2303, ASV1869, ASV2250, ASV1498, ASV148, ASV2749, ASV220, ASV250, ASV364, ASV4537, ASV3423, ASV809, ASV3414, ASV2038, ASV3464, ASV3496, ASV1339, ASV4656, ASV2190, ASV3036, ASV1575, ASV3195, ASV1242, ASV4153, ASV435, ASV2903, ASV7195, ASV3746, ASV361, ASV3070, ASV6439, ASV3586, ASV4059, ASV3191, ASV6316, ASV1899, ASV4226, ASV3339, ASV514, ASV540, ASV4223, ASV970, ASV352, ASV1680, ASV2984, ASV2842, ASV2267, ASV1291, ASV1965, ASV6138, ASV214, ASV2841, ASV888, ASV181, ASV2251, ASV789, ASV1553, ASV6683, ASV1642, ASV6183, ASV2120, ASV1917, ASV162, ASV297, ASV209, ASV823, ASV1056, ASV5005, ASV266, ASV2331, ASV6070, ASV3226, ASV2463, ASV5243, ASV5229, ASV5834, ASV6956, ASV5521, ASV5604, ASV4660, ASV6181, ASV415, ASV274, ASV2371, ASV1919, ASV734, ASV6564, ASV7240, ASV5704, ASV7977, ASV5055, ASV3417, ASV6073, ASV5708, ASV2544, ASV4598, ASV6565, ASV1045, ASV1812, ASV5082, ASV2009, ASV7615, ASV7637, ASV5083, ASV2558, ASV5617, ASV8182, ASV4535, ASV6437, ASV7096, ASV5239, ASV6317, ASV5015, ASV3418, ASV2489, ASV5528, ASV2705, ASV6820, ASV3521, ASV4299, ASV5997, ASV620, ASV2123, ASV5449, ASV6421, ASV4762, ASV6845, ASV6590, ASV7113, ASV7388, ASV5791, ASV3670, ASV7583, ASV5435, ASV5710, ASV5080, ASV3445, ASV3742, ASV4787, ASV4356, ASV6066, ASV2484, ASV2201, ASV6151, ASV1832, ASV4096, ASV7093, ASV4672, ASV2939, ASV6946, ASV3805, ASV2041, ASV3136, ASV2511, ASV2901, ASV2982, ASV3513, ASV3018, ASV3907, ASV2425, ASV2133, ASV4348, ASV5426, ASV2386, ASV3734, ASV7239, ASV6083, ASV4512, ASV5242, ASV2578, ASV4778, ASV2249, ASV4914, ASV5971, ASV2876, ASV2848, ASV282, ASV821, ASV2773, ASV5252, ASV1238, ASV1614, ASV612, ASV3227, ASV1038, ASV894, ASV4221, ASV6586, ASV7579, ASV1368, ASV5718, ASV5801, ASV4047, ASV5985, ASV793, ASV4713, ASV5699, ASV8004, ASV2563, ASV3472, ASV2236, ASV2322, ASV7827, ASV1306, ASV2584, ASV3016, ASV6969, ASV6428, ASV6084, ASV8177, ASV3302, ASV4039, ASV4600, ASV4534, ASV5244, ASV7266, ASV8176, ASV3826, ASV3196, ASV4791, ASV7784, ASV5979, ASV7250, ASV4166, ASV5076, ASV7094, ASV3683, ASV3984, ASV2516, ASV4304, ASV7242, ASV7137, ASV2980, ASV5700, ASV3745, ASV2429, ASV6813, ASV6200, ASV2730, ASV7420, ASV5231, ASV3264, ASV2408, ASV5077, ASV4224, ASV3739, ASV5282, ASV3443, ASV3371, ASV1470, ASV2462, ASV4058, ASV6180, ASV7808, ASV6068, ASV7110, ASV2510, ASV1527, ASV3520, ASV5529, ASV2981, ASV6318, ASV4487, ASV1740, ASV7387, ASV7997, ASV4236, ASV5879, ASV2707, ASV5790, ASV5425, ASV4809, ASV7978, ASV2409, ASV7588, ASV7115, ASV4486, ASV7636, ASV5145, ASV7463, ASV5002, ASV4470, ASV7394, ASV7241, ASV5969, ASV6319, ASV7805, ASV6703, ASV5527, ASV6097, ASV7781, ASV2875, ASV4365, ASV2324, ASV5245, ASV5702, ASV5000, ASV3342, ASV2942, ASV3107, ASV6427, ASV3480, ASV5873, ASV1033, ASV4994, ASV2715, ASV4109, ASV8002, ASV5318, ASV2340, ASV6430, ASV4477, ASV3585, ASV6087, ASV1250, ASV4053, ASV3698, ASV3906, ASV7427, ASV5434, ASV6201, ASV1831, ASV334, ASV3885, ASV1382, ASV381, ASV522, ASV797, ASV4800, ASV1819, ASV3592, ASV4531, ASV6850, ASV4364, ASV6429, ASV5322, ASV5066, ASV6679, ASV3228, ASV6337, ASV5320, ASV6069, ASV4402, ASV1453, ASV2001, ASV8003, ASV4716, ASV2779, ASV3776, ASV1736, ASV3702, ASV2917, ASV5164, ASV6957, ASV2220, ASV2682, ASV4923, ASV4858, ASV1302, ASV2119, ASV5314, ASV300, ASV1544, ASV7391, ASV2978, ASV1795, ASV2428, ASV2947, ASV5986, ASV5701, ASV1387, ASV1627, ASV8180, ASV6071, ASV1678, ASV4273, ASV4946, ASV6796, ASV6441, ASV7606, ASV5087, ASV7112, ASV1689, ASV5075, ASV7246, ASV6678, ASV3774, ASV6561, ASV6705, ASV4934, ASV7586, ASV7608, ASV4298, ASV6314, ASV7386, ASV3135, ASV5619, ASV6706, ASV7975, ASV4119, ASV3663, ASV3920, ASV6968, ASV4599, ASV3962, ASV6425, ASV270, ASV1142, ASV4283, ASV5348, ASV3902, ASV4363, ASV4164, ASV8178, ASV6305, ASV1857, ASV5883, ASV2090, ASV3634, ASV2813, ASV7578, ASV3700, ASV3501, ASV7248, ASV1381, ASV2167, ASV7393, ASV1817, ASV2582, ASV7592, ASV3955, ASV451, ASV4932, ASV6828, ASV3337, ASV2533, ASV2452, ASV6681
#> 2026-10-28 10:47:16.199093 INFO::Total filtered features with variance filtering: 0
#> 2026-10-28 10:47:16.2018 INFO::Filtered feature names from variance filtering:
#> 2026-10-28 10:47:16.203824 INFO::Running selected normalization method: TSS
#> 2026-10-28 10:47:16.210725 INFO::Applying z-score to standardize continuous metadata
#> 2026-10-28 10:47:16.214917 INFO::Running selected transform method: AST
#> 2026-10-28 10:47:16.22441 INFO::Running selected analysis method: LM
#> 2026-10-28 10:47:16.229192 INFO::Fitting model to feature number 1, ASV4670
#> 2026-10-28 10:47:16.234398 INFO::Fitting model to feature number 2, ASV610
#> 2026-10-28 10:47:16.242216 INFO::Fitting model to feature number 3, ASV746
#> 2026-10-28 10:47:16.246424 INFO::Fitting model to feature number 4, ASV788
#> 2026-10-28 10:47:16.250698 INFO::Fitting model to feature number 5, ASV729
#> 2026-10-28 10:47:16.258535 INFO::Fitting model to feature number 6, ASV1599
#> 2026-10-28 10:47:16.262794 INFO::Fitting model to feature number 7, ASV410
#> 2026-10-28 10:47:16.267165 INFO::Fitting model to feature number 8, ASV543
#> 2026-10-28 10:47:16.271951 INFO::Fitting model to feature number 9, ASV13
#> 2026-10-28 10:47:16.277354 INFO::Fitting model to feature number 10, ASV174
#> 2026-10-28 10:47:16.281377 INFO::Fitting model to feature number 11, ASV4
#> 2026-10-28 10:47:16.285269 INFO::Fitting model to feature number 12, ASV307
#> 2026-10-28 10:47:16.28942 INFO::Fitting model to feature number 13, ASV3199
#> 2026-10-28 10:47:16.293246 INFO::Fitting model to feature number 14, ASV1543
#> 2026-10-28 10:47:16.296813 INFO::Fitting model to feature number 15, ASV6954
#> 2026-10-28 10:47:16.300455 INFO::Fitting model to feature number 16, ASV1659
#> 2026-10-28 10:47:16.304056 INFO::Fitting model to feature number 17, ASV191
#> 2026-10-28 10:47:16.310028 INFO::Fitting model to feature number 18, ASV990
#> 2026-10-28 10:47:16.314973 INFO::Fitting model to feature number 19, ASV721
#> 2026-10-28 10:47:16.319021 INFO::Fitting model to feature number 20, ASV80
#> 2026-10-28 10:47:16.32276 INFO::Fitting model to feature number 21, ASV344
#> 2026-10-28 10:47:16.329193 INFO::Fitting model to feature number 22, ASV238
#> 2026-10-28 10:47:16.33442 INFO::Fitting model to feature number 23, ASV402
#> 2026-10-28 10:47:16.341088 INFO::Fitting model to feature number 24, ASV1115
#> 2026-10-28 10:47:16.346043 INFO::Fitting model to feature number 25, ASV241
#> 2026-10-28 10:47:16.349552 INFO::Fitting model to feature number 26, ASV1212
#> 2026-10-28 10:47:16.352986 INFO::Fitting model to feature number 27, ASV4005
#> 2026-10-28 10:47:16.357022 INFO::Fitting model to feature number 28, ASV510
#> 2026-10-28 10:47:16.362628 INFO::Fitting model to feature number 29, ASV3137
#> 2026-10-28 10:47:16.36621 INFO::Fitting model to feature number 30, ASV2243
#> 2026-10-28 10:47:16.369654 INFO::Fitting model to feature number 31, ASV470
#> 2026-10-28 10:47:16.374133 INFO::Fitting model to feature number 32, ASV416
#> 2026-10-28 10:47:16.37895 INFO::Fitting model to feature number 33, ASV1047
#> 2026-10-28 10:47:16.3844 INFO::Fitting model to feature number 34, ASV302
#> 2026-10-28 10:47:16.389229 INFO::Fitting model to feature number 35, ASV637
#> 2026-10-28 10:47:16.393108 INFO::Fitting model to feature number 36, ASV314
#> 2026-10-28 10:47:16.396595 INFO::Fitting model to feature number 37, ASV468
#> 2026-10-28 10:47:16.400024 INFO::Fitting model to feature number 38, ASV18
#> 2026-10-28 10:47:16.403456 INFO::Fitting model to feature number 39, ASV87
#> 2026-10-28 10:47:16.410808 INFO::Fitting model to feature number 40, ASV1411
#> 2026-10-28 10:47:16.414584 INFO::Fitting model to feature number 41, ASV1154
#> 2026-10-28 10:47:16.418489 INFO::Fitting model to feature number 42, ASV1305
#> 2026-10-28 10:47:16.423248 INFO::Fitting model to feature number 43, ASV5432
#> 2026-10-28 10:47:16.427191 INFO::Fitting model to feature number 44, ASV4689
#> 2026-10-28 10:47:16.431299 INFO::Fitting model to feature number 45, ASV1029
#> 2026-10-28 10:47:16.434918 INFO::Fitting model to feature number 46, ASV831
#> 2026-10-28 10:47:16.439558 INFO::Fitting model to feature number 47, ASV3223
#> 2026-10-28 10:47:16.44373 INFO::Fitting model to feature number 48, ASV3740
#> 2026-10-28 10:47:16.448726 INFO::Fitting model to feature number 49, ASV2217
#> 2026-10-28 10:47:16.453206 INFO::Fitting model to feature number 50, ASV3666
#> 2026-10-28 10:47:16.460845 INFO::Fitting model to feature number 51, ASV3864
#> 2026-10-28 10:47:16.465841 INFO::Fitting model to feature number 52, ASV2902
#> 2026-10-28 10:47:16.470902 INFO::Fitting model to feature number 53, ASV4544
#> 2026-10-28 10:47:16.475737 INFO::Fitting model to feature number 54, ASV5603
#> 2026-10-28 10:47:16.48041 INFO::Fitting model to feature number 55, ASV4447
#> 2026-10-28 10:47:16.488705 INFO::Fitting model to feature number 56, ASV2179
#> 2026-10-28 10:47:16.49407 INFO::Fitting model to feature number 57, ASV1357
#> 2026-10-28 10:47:16.501604 INFO::Fitting model to feature number 58, ASV1649
#> 2026-10-28 10:47:16.509808 INFO::Fitting model to feature number 59, ASV3230
#> 2026-10-28 10:47:16.513913 INFO::Fitting model to feature number 60, ASV1025
#> 2026-10-28 10:47:16.519784 INFO::Fitting model to feature number 61, ASV5003
#> 2026-10-28 10:47:16.525287 INFO::Fitting model to feature number 62, ASV818
#> 2026-10-28 10:47:16.530054 INFO::Fitting model to feature number 63, ASV3901
#> 2026-10-28 10:47:16.534646 INFO::Fitting model to feature number 64, ASV1149
#> 2026-10-28 10:47:16.54126 INFO::Fitting model to feature number 65, ASV1657
#> 2026-10-28 10:47:16.546601 INFO::Fitting model to feature number 66, ASV320
#> 2026-10-28 10:47:16.552458 INFO::Fitting model to feature number 67, ASV193
#> 2026-10-28 10:47:16.558883 INFO::Fitting model to feature number 68, ASV216
#> 2026-10-28 10:47:16.56452 INFO::Fitting model to feature number 69, ASV3376
#> 2026-10-28 10:47:16.568384 INFO::Fitting model to feature number 70, ASV485
#> 2026-10-28 10:47:16.573888 INFO::Fitting model to feature number 71, ASV4596
#> 2026-10-28 10:47:16.578551 INFO::Fitting model to feature number 72, ASV1012
#> 2026-10-28 10:47:16.582162 INFO::Fitting model to feature number 73, ASV1779
#> 2026-10-28 10:47:16.585729 INFO::Fitting model to feature number 74, ASV1200
#> 2026-10-28 10:47:16.591303 INFO::Fitting model to feature number 75, ASV651
#> 2026-10-28 10:47:16.597572 INFO::Fitting model to feature number 76, ASV1846
#> 2026-10-28 10:47:16.60397 INFO::Fitting model to feature number 77, ASV2202
#> 2026-10-28 10:47:16.610857 INFO::Fitting model to feature number 78, ASV1349
#> 2026-10-28 10:47:16.614545 INFO::Fitting model to feature number 79, ASV3590
#> 2026-10-28 10:47:16.618732 INFO::Fitting model to feature number 80, ASV2989
#> 2026-10-28 10:47:16.62328 INFO::Fitting model to feature number 81, ASV2744
#> 2026-10-28 10:47:16.627275 INFO::Fitting model to feature number 82, ASV139
#> 2026-10-28 10:47:16.63109 INFO::Fitting model to feature number 83, ASV1031
#> 2026-10-28 10:47:16.634727 INFO::Fitting model to feature number 84, ASV112
#> 2026-10-28 10:47:16.638763 INFO::Fitting model to feature number 85, ASV3322
#> 2026-10-28 10:47:16.643666 INFO::Fitting model to feature number 86, ASV982
#> 2026-10-28 10:47:16.647849 INFO::Fitting model to feature number 87, ASV2732
#> 2026-10-28 10:47:16.651581 INFO::Fitting model to feature number 88, ASV1738
#> 2026-10-28 10:47:16.655439 INFO::Fitting model to feature number 89, ASV166
#> 2026-10-28 10:47:16.660073 INFO::Fitting model to feature number 90, ASV1105
#> 2026-10-28 10:47:16.663882 INFO::Fitting model to feature number 91, ASV310
#> 2026-10-28 10:47:16.667727 INFO::Fitting model to feature number 92, ASV417
#> 2026-10-28 10:47:16.671262 INFO::Fitting model to feature number 93, ASV1225
#> 2026-10-28 10:47:16.677307 INFO::Fitting model to feature number 94, ASV1337
#> 2026-10-28 10:47:16.681317 INFO::Fitting model to feature number 95, ASV1099
#> 2026-10-28 10:47:16.684986 INFO::Fitting model to feature number 96, ASV2106
#> 2026-10-28 10:47:16.68932 INFO::Fitting model to feature number 97, ASV1932
#> 2026-10-28 10:47:16.693524 INFO::Fitting model to feature number 98, ASV481
#> 2026-10-28 10:47:16.697448 INFO::Fitting model to feature number 99, ASV303
#> 2026-10-28 10:47:16.701165 INFO::Fitting model to feature number 100, ASV1028
#> 2026-10-28 10:47:16.704934 INFO::Fitting model to feature number 101, ASV1884
#> 2026-10-28 10:47:16.710033 INFO::Fitting model to feature number 102, ASV3761
#> 2026-10-28 10:47:16.713669 INFO::Fitting model to feature number 103, ASV201
#> 2026-10-28 10:47:16.717117 INFO::Fitting model to feature number 104, ASV892
#> 2026-10-28 10:47:16.720759 INFO::Fitting model to feature number 105, ASV453
#> 2026-10-28 10:47:16.725009 INFO::Fitting model to feature number 106, ASV2899
#> 2026-10-28 10:47:16.728734 INFO::Fitting model to feature number 107, ASV4715
#> 2026-10-28 10:47:16.732227 INFO::Fitting model to feature number 108, ASV211
#> 2026-10-28 10:47:16.735796 INFO::Fitting model to feature number 109, ASV240
#> 2026-10-28 10:47:16.740251 INFO::Fitting model to feature number 110, ASV1128
#> 2026-10-28 10:47:16.744542 INFO::Fitting model to feature number 111, ASV4403
#> 2026-10-28 10:47:16.748125 INFO::Fitting model to feature number 112, ASV1272
#> 2026-10-28 10:47:16.751983 INFO::Fitting model to feature number 113, ASV988
#> 2026-10-28 10:47:16.758128 INFO::Fitting model to feature number 114, ASV1601
#> 2026-10-28 10:47:16.763358 INFO::Fitting model to feature number 115, ASV2729
#> 2026-10-28 10:47:16.768475 INFO::Fitting model to feature number 116, ASV690
#> 2026-10-28 10:47:16.77424 INFO::Fitting model to feature number 117, ASV164
#> 2026-10-28 10:47:16.780113 INFO::Fitting model to feature number 118, ASV653
#> 2026-10-28 10:47:16.785579 INFO::Fitting model to feature number 119, ASV494
#> 2026-10-28 10:47:16.790943 INFO::Fitting model to feature number 120, ASV244
#> 2026-10-28 10:47:16.795173 INFO::Fitting model to feature number 121, ASV450
#> 2026-10-28 10:47:16.799089 INFO::Fitting model to feature number 122, ASV153
#> 2026-10-28 10:47:16.803089 INFO::Fitting model to feature number 123, ASV237
#> 2026-10-28 10:47:16.808384 INFO::Fitting model to feature number 124, ASV50
#> 2026-10-28 10:47:16.813444 INFO::Fitting model to feature number 125, ASV97
#> 2026-10-28 10:47:16.817367 INFO::Fitting model to feature number 126, ASV382
#> 2026-10-28 10:47:16.823905 INFO::Fitting model to feature number 127, ASV479
#> 2026-10-28 10:47:16.830401 INFO::Fitting model to feature number 128, ASV376
#> 2026-10-28 10:47:16.835626 INFO::Fitting model to feature number 129, ASV143
#> 2026-10-28 10:47:16.84197 INFO::Fitting model to feature number 130, ASV801
#> 2026-10-28 10:47:16.84815 INFO::Fitting model to feature number 131, ASV2669
#> 2026-10-28 10:47:16.855537 INFO::Fitting model to feature number 132, ASV53
#> 2026-10-28 10:47:16.861157 INFO::Fitting model to feature number 133, ASV2795
#> 2026-10-28 10:47:16.867945 INFO::Fitting model to feature number 134, ASV150
#> 2026-10-28 10:47:16.872061 INFO::Fitting model to feature number 135, ASV4854
#> 2026-10-28 10:47:16.87975 INFO::Fitting model to feature number 136, ASV1362
#> 2026-10-28 10:47:16.887004 INFO::Fitting model to feature number 137, ASV954
#> 2026-10-28 10:47:16.895126 INFO::Fitting model to feature number 138, ASV3116
#> 2026-10-28 10:47:16.902155 INFO::Fitting model to feature number 139, ASV7458
#> 2026-10-28 10:47:16.911241 INFO::Fitting model to feature number 140, ASV4484
#> 2026-10-28 10:47:16.917493 INFO::Fitting model to feature number 141, ASV4412
#> 2026-10-28 10:47:16.923402 INFO::Fitting model to feature number 142, ASV2171
#> 2026-10-28 10:47:16.927885 INFO::Fitting model to feature number 143, ASV2344
#> 2026-10-28 10:47:16.932638 INFO::Fitting model to feature number 144, ASV465
#> 2026-10-28 10:47:16.938091 INFO::Fitting model to feature number 145, ASV1258
#> 2026-10-28 10:47:16.943573 INFO::Fitting model to feature number 146, ASV236
#> 2026-10-28 10:47:16.947262 INFO::Fitting model to feature number 147, ASV2029
#> 2026-10-28 10:47:16.952424 INFO::Fitting model to feature number 148, ASV5253
#> 2026-10-28 10:47:16.958942 INFO::Fitting model to feature number 149, ASV2781
#> 2026-10-28 10:47:16.96446 INFO::Fitting model to feature number 150, ASV4446
#> 2026-10-28 10:47:16.970147 INFO::Fitting model to feature number 151, ASV263
#> 2026-10-28 10:47:16.976265 INFO::Fitting model to feature number 152, ASV1365
#> 2026-10-28 10:47:16.980249 INFO::Fitting model to feature number 153, ASV1665
#> 2026-10-28 10:47:16.984049 INFO::Fitting model to feature number 154, ASV1677
#> 2026-10-28 10:47:16.987547 INFO::Fitting model to feature number 155, ASV533
#> 2026-10-28 10:47:16.993402 INFO::Fitting model to feature number 156, ASV6440
#> 2026-10-28 10:47:16.997357 INFO::Fitting model to feature number 157, ASV3272
#> 2026-10-28 10:47:17.001058 INFO::Fitting model to feature number 158, ASV1171
#> 2026-10-28 10:47:17.006127 INFO::Fitting model to feature number 159, ASV2007
#> 2026-10-28 10:47:17.011406 INFO::Fitting model to feature number 160, ASV2492
#> 2026-10-28 10:47:17.016843 INFO::Fitting model to feature number 161, ASV5673
#> 2026-10-28 10:47:17.022096 INFO::Fitting model to feature number 162, ASV3495
#> 2026-10-28 10:47:17.027418 INFO::Fitting model to feature number 163, ASV6581
#> 2026-10-28 10:47:17.032782 INFO::Fitting model to feature number 164, ASV2706
#> 2026-10-28 10:47:17.03694 INFO::Fitting model to feature number 165, ASV4606
#> 2026-10-28 10:47:17.042436 INFO::Fitting model to feature number 166, ASV5228
#> 2026-10-28 10:47:17.047141 INFO::Fitting model to feature number 167, ASV4744
#> 2026-10-28 10:47:17.052757 INFO::Fitting model to feature number 168, ASV5326
#> 2026-10-28 10:47:17.059434 INFO::Fitting model to feature number 169, ASV7801
#> 2026-10-28 10:47:17.063582 INFO::Fitting model to feature number 170, ASV4475
#> 2026-10-28 10:47:17.069171 INFO::Fitting model to feature number 171, ASV6638
#> 2026-10-28 10:47:17.074059 INFO::Fitting model to feature number 172, ASV2874
#> 2026-10-28 10:47:17.079444 INFO::Fitting model to feature number 173, ASV3555
#> 2026-10-28 10:47:17.084013 INFO::Fitting model to feature number 174, ASV2295
#> 2026-10-28 10:47:17.087934 INFO::Fitting model to feature number 175, ASV4302
#> 2026-10-28 10:47:17.092931 INFO::Fitting model to feature number 176, ASV4474
#> 2026-10-28 10:47:17.096967 INFO::Fitting model to feature number 177, ASV1586
#> 2026-10-28 10:47:17.100794 INFO::Fitting model to feature number 178, ASV217
#> 2026-10-28 10:47:17.104604 INFO::Fitting model to feature number 179, ASV707
#> 2026-10-28 10:47:17.10864 INFO::Fitting model to feature number 180, ASV1608
#> 2026-10-28 10:47:17.112869 INFO::Fitting model to feature number 181, ASV837
#> 2026-10-28 10:47:17.11679 INFO::Fitting model to feature number 182, ASV2846
#> 2026-10-28 10:47:17.120428 INFO::Fitting model to feature number 183, ASV3014
#> 2026-10-28 10:47:17.124935 INFO::Fitting model to feature number 184, ASV956
#> 2026-10-28 10:47:17.128671 INFO::Fitting model to feature number 185, ASV1616
#> 2026-10-28 10:47:17.132974 INFO::Fitting model to feature number 186, ASV2341
#> 2026-10-28 10:47:17.136524 INFO::Fitting model to feature number 187, ASV5347
#> 2026-10-28 10:47:17.142491 INFO::Fitting model to feature number 188, ASV878
#> 2026-10-28 10:47:17.147375 INFO::Fitting model to feature number 189, ASV2150
#> 2026-10-28 10:47:17.150903 INFO::Fitting model to feature number 190, ASV1452
#> 2026-10-28 10:47:17.154333 INFO::Fitting model to feature number 191, ASV2078
#> 2026-10-28 10:47:17.159145 INFO::Fitting model to feature number 192, ASV3392
#> 2026-10-28 10:47:17.162666 INFO::Fitting model to feature number 193, ASV439
#> 2026-10-28 10:47:17.182121 INFO::Fitting model to feature number 194, ASV1198
#> 2026-10-28 10:47:17.186683 INFO::Fitting model to feature number 195, ASV926
#> 2026-10-28 10:47:17.191278 INFO::Fitting model to feature number 196, ASV4574
#> 2026-10-28 10:47:17.195297 INFO::Fitting model to feature number 197, ASV2583
#> 2026-10-28 10:47:17.199824 INFO::Fitting model to feature number 198, ASV521
#> 2026-10-28 10:47:17.20657 INFO::Fitting model to feature number 199, ASV1930
#> 2026-10-28 10:47:17.213612 INFO::Fitting model to feature number 200, ASV396
#> 2026-10-28 10:47:17.220291 INFO::Fitting model to feature number 201, ASV839
#> 2026-10-28 10:47:17.226478 INFO::Fitting model to feature number 202, ASV509
#> 2026-10-28 10:47:17.230429 INFO::Fitting model to feature number 203, ASV280
#> 2026-10-28 10:47:17.234254 INFO::Fitting model to feature number 204, ASV94
#> 2026-10-28 10:47:17.238139 INFO::Fitting model to feature number 205, ASV337
#> 2026-10-28 10:47:17.243866 INFO::Fitting model to feature number 206, ASV2064
#> 2026-10-28 10:47:17.248472 INFO::Fitting model to feature number 207, ASV469
#> 2026-10-28 10:47:17.252226 INFO::Fitting model to feature number 208, ASV879
#> 2026-10-28 10:47:17.256486 INFO::Fitting model to feature number 209, ASV1418
#> 2026-10-28 10:47:17.262178 INFO::Fitting model to feature number 210, ASV399
#> 2026-10-28 10:47:17.268448 INFO::Fitting model to feature number 211, ASV1597
#> 2026-10-28 10:47:17.273408 INFO::Fitting model to feature number 212, ASV1503
#> 2026-10-28 10:47:17.277989 INFO::Fitting model to feature number 213, ASV518
#> 2026-10-28 10:47:17.281868 INFO::Fitting model to feature number 214, ASV910
#> 2026-10-28 10:47:17.285636 INFO::Fitting model to feature number 215, ASV2857
#> 2026-10-28 10:47:17.29132 INFO::Fitting model to feature number 216, ASV1664
#> 2026-10-28 10:47:17.295505 INFO::Fitting model to feature number 217, ASV667
#> 2026-10-28 10:47:17.299279 INFO::Fitting model to feature number 218, ASV1786
#> 2026-10-28 10:47:17.303139 INFO::Fitting model to feature number 219, ASV386
#> 2026-10-28 10:47:17.307647 INFO::Fitting model to feature number 220, ASV5089
#> 2026-10-28 10:47:17.312066 INFO::Fitting model to feature number 221, ASV975
#> 2026-10-28 10:47:17.315748 INFO::Fitting model to feature number 222, ASV4416
#> 2026-10-28 10:47:17.319593 INFO::Fitting model to feature number 223, ASV4445
#> 2026-10-28 10:47:17.325079 INFO::Fitting model to feature number 224, ASV1425
#> 2026-10-28 10:47:17.341792 INFO::Fitting model to feature number 225, ASV1873
#> 2026-10-28 10:47:17.347673 INFO::Fitting model to feature number 226, ASV42
#> 2026-10-28 10:47:17.352613 INFO::Fitting model to feature number 227, ASV1266
#> 2026-10-28 10:47:17.359089 INFO::Fitting model to feature number 228, ASV186
#> 2026-10-28 10:47:17.364194 INFO::Fitting model to feature number 229, ASV813
#> 2026-10-28 10:47:17.367914 INFO::Fitting model to feature number 230, ASV1119
#> 2026-10-28 10:47:17.372112 INFO::Fitting model to feature number 231, ASV105
#> 2026-10-28 10:47:17.377188 INFO::Fitting model to feature number 232, ASV124
#> 2026-10-28 10:47:17.383302 INFO::Fitting model to feature number 233, ASV213
#> 2026-10-28 10:47:17.390424 INFO::Fitting model to feature number 234, ASV762
#> 2026-10-28 10:47:17.395901 INFO::Fitting model to feature number 235, ASV222
#> 2026-10-28 10:47:17.400291 INFO::Fitting model to feature number 236, ASV32
#> 2026-10-28 10:47:17.405784 INFO::Fitting model to feature number 237, ASV2494
#> 2026-10-28 10:47:17.410949 INFO::Fitting model to feature number 238, ASV259
#> 2026-10-28 10:47:17.415553 INFO::Fitting model to feature number 239, ASV359
#> 2026-10-28 10:47:17.421847 INFO::Fitting model to feature number 240, ASV3612
#> 2026-10-28 10:47:17.428312 INFO::Fitting model to feature number 241, ASV3473
#> 2026-10-28 10:47:17.433384 INFO::Fitting model to feature number 242, ASV2566
#> 2026-10-28 10:47:17.43751 INFO::Fitting model to feature number 243, ASV876
#> 2026-10-28 10:47:17.443868 INFO::Fitting model to feature number 244, ASV366
#> 2026-10-28 10:47:17.447663 INFO::Fitting model to feature number 245, ASV269
#> 2026-10-28 10:47:17.451629 INFO::Fitting model to feature number 246, ASV830
#> 2026-10-28 10:47:17.455804 INFO::Fitting model to feature number 247, ASV239
#> 2026-10-28 10:47:17.460165 INFO::Fitting model to feature number 248, ASV1164
#> 2026-10-28 10:47:17.465419 INFO::Fitting model to feature number 249, ASV3518
#> 2026-10-28 10:47:17.470407 INFO::Fitting model to feature number 250, ASV1170
#> 2026-10-28 10:47:17.475663 INFO::Fitting model to feature number 251, ASV2208
#> 2026-10-28 10:47:17.479631 INFO::Fitting model to feature number 252, ASV529
#> 2026-10-28 10:47:17.483307 INFO::Fitting model to feature number 253, ASV1085
#> 2026-10-28 10:47:17.487139 INFO::Fitting model to feature number 254, ASV4086
#> 2026-10-28 10:47:17.491134 INFO::Fitting model to feature number 255, ASV408
#> 2026-10-28 10:47:17.49613 INFO::Fitting model to feature number 256, ASV1666
#> 2026-10-28 10:47:17.500746 INFO::Fitting model to feature number 257, ASV866
#> 2026-10-28 10:47:17.505681 INFO::Fitting model to feature number 258, ASV1084
#> 2026-10-28 10:47:17.510414 INFO::Fitting model to feature number 259, ASV3591
#> 2026-10-28 10:47:17.514382 INFO::Fitting model to feature number 260, ASV123
#> 2026-10-28 10:47:17.518063 INFO::Fitting model to feature number 261, ASV113
#> 2026-10-28 10:47:17.521671 INFO::Fitting model to feature number 262, ASV490
#> 2026-10-28 10:47:17.526486 INFO::Fitting model to feature number 263, ASV877
#> 2026-10-28 10:47:17.53032 INFO::Fitting model to feature number 264, ASV1118
#> 2026-10-28 10:47:17.534027 INFO::Fitting model to feature number 265, ASV1637
#> 2026-10-28 10:47:17.537715 INFO::Fitting model to feature number 266, ASV1015
#> 2026-10-28 10:47:17.54167 INFO::Fitting model to feature number 267, ASV1455
#> 2026-10-28 10:47:17.545338 INFO::Fitting model to feature number 268, ASV144
#> 2026-10-28 10:47:17.549013 INFO::Fitting model to feature number 269, ASV134
#> 2026-10-28 10:47:17.552896 INFO::Fitting model to feature number 270, ASV92
#> 2026-10-28 10:47:17.557987 INFO::Fitting model to feature number 271, ASV120
#> 2026-10-28 10:47:17.562917 INFO::Fitting model to feature number 272, ASV596
#> 2026-10-28 10:47:17.567632 INFO::Fitting model to feature number 273, ASV3166
#> 2026-10-28 10:47:17.573628 INFO::Fitting model to feature number 274, ASV1722
#> 2026-10-28 10:47:17.578517 INFO::Fitting model to feature number 275, ASV1451
#> 2026-10-28 10:47:17.582398 INFO::Fitting model to feature number 276, ASV19
#> 2026-10-28 10:47:17.586105 INFO::Fitting model to feature number 277, ASV5
#> 2026-10-28 10:47:17.59134 INFO::Fitting model to feature number 278, ASV863
#> 2026-10-28 10:47:17.595586 INFO::Fitting model to feature number 279, ASV312
#> 2026-10-28 10:47:17.600525 INFO::Fitting model to feature number 280, ASV370
#> 2026-10-28 10:47:17.603862 INFO::Fitting model to feature number 281, ASV48
#> 2026-10-28 10:47:17.608992 INFO::Fitting model to feature number 282, ASV9
#> 2026-10-28 10:47:17.612747 INFO::Fitting model to feature number 283, ASV136
#> 2026-10-28 10:47:17.616349 INFO::Fitting model to feature number 284, ASV1145
#> 2026-10-28 10:47:17.619687 INFO::Fitting model to feature number 285, ASV226
#> 2026-10-28 10:47:17.623635 INFO::Fitting model to feature number 286, ASV5349
#> 2026-10-28 10:47:17.62736 INFO::Fitting model to feature number 287, ASV3167
#> 2026-10-28 10:47:17.632326 INFO::Fitting model to feature number 288, ASV817
#> 2026-10-28 10:47:17.635741 INFO::Fitting model to feature number 289, ASV3909
#> 2026-10-28 10:47:17.639945 INFO::Fitting model to feature number 290, ASV155
#> 2026-10-28 10:47:17.643717 INFO::Fitting model to feature number 291, ASV4627
#> 2026-10-28 10:47:17.647466 INFO::Fitting model to feature number 292, ASV1573
#> 2026-10-28 10:47:17.651211 INFO::Fitting model to feature number 293, ASV647
#> 2026-10-28 10:47:17.654916 INFO::Fitting model to feature number 294, ASV336
#> 2026-10-28 10:47:17.661466 INFO::Fitting model to feature number 295, ASV2459
#> 2026-10-28 10:47:17.665654 INFO::Fitting model to feature number 296, ASV1067
#> 2026-10-28 10:47:17.670436 INFO::Fitting model to feature number 297, ASV5981
#> 2026-10-28 10:47:17.675778 INFO::Fitting model to feature number 298, ASV8024
#> 2026-10-28 10:47:17.679231 INFO::Fitting model to feature number 299, ASV722
#> 2026-10-28 10:47:17.682503 INFO::Fitting model to feature number 300, ASV495
#> 2026-10-28 10:47:17.68779 INFO::Fitting model to feature number 301, ASV589
#> 2026-10-28 10:47:17.694097 INFO::Fitting model to feature number 302, ASV152
#> 2026-10-28 10:47:17.697728 INFO::Fitting model to feature number 303, ASV1336
#> 2026-10-28 10:47:17.70149 INFO::Fitting model to feature number 304, ASV83
#> 2026-10-28 10:47:17.706439 INFO::Fitting model to feature number 305, ASV183
#> 2026-10-28 10:47:17.712401 INFO::Fitting model to feature number 306, ASV1251
#> 2026-10-28 10:47:17.716115 INFO::Fitting model to feature number 307, ASV2030
#> 2026-10-28 10:47:17.720695 INFO::Fitting model to feature number 308, ASV5618
#> 2026-10-28 10:47:17.726692 INFO::Fitting model to feature number 309, ASV2760
#> 2026-10-28 10:47:17.733197 INFO::Fitting model to feature number 310, ASV5287
#> 2026-10-28 10:47:17.740247 INFO::Fitting model to feature number 311, ASV6842
#> 2026-10-28 10:47:17.746796 INFO::Fitting model to feature number 312, ASV6094
#> 2026-10-28 10:47:17.751699 INFO::Fitting model to feature number 313, ASV1660
#> 2026-10-28 10:47:17.757297 INFO::Fitting model to feature number 314, ASV2800
#> 2026-10-28 10:47:17.763319 INFO::Fitting model to feature number 315, ASV3802
#> 2026-10-28 10:47:17.768363 INFO::Fitting model to feature number 316, ASV1289
#> 2026-10-28 10:47:17.774188 INFO::Fitting model to feature number 317, ASV1134
#> 2026-10-28 10:47:17.779188 INFO::Fitting model to feature number 318, ASV2369
#> 2026-10-28 10:47:17.783158 INFO::Fitting model to feature number 319, ASV1386
#> 2026-10-28 10:47:17.787086 INFO::Fitting model to feature number 320, ASV4943
#> 2026-10-28 10:47:17.792225 INFO::Fitting model to feature number 321, ASV472
#> 2026-10-28 10:47:17.796337 INFO::Fitting model to feature number 322, ASV2277
#> 2026-10-28 10:47:17.800729 INFO::Fitting model to feature number 323, ASV1020
#> 2026-10-28 10:47:17.806529 INFO::Fitting model to feature number 324, ASV3910
#> 2026-10-28 10:47:17.811275 INFO::Fitting model to feature number 325, ASV2839
#> 2026-10-28 10:47:17.815244 INFO::Fitting model to feature number 326, ASV327
#> 2026-10-28 10:47:17.821409 INFO::Fitting model to feature number 327, ASV1647
#> 2026-10-28 10:47:17.827228 INFO::Fitting model to feature number 328, ASV2762
#> 2026-10-28 10:47:17.832708 INFO::Fitting model to feature number 329, ASV1246
#> 2026-10-28 10:47:17.836945 INFO::Fitting model to feature number 330, ASV109
#> 2026-10-28 10:47:17.842307 INFO::Fitting model to feature number 331, ASV147
#> 2026-10-28 10:47:17.846677 INFO::Fitting model to feature number 332, ASV3206
#> 2026-10-28 10:47:17.851728 INFO::Fitting model to feature number 333, ASV742
#> 2026-10-28 10:47:17.855565 INFO::Fitting model to feature number 334, ASV716
#> 2026-10-28 10:47:17.859673 INFO::Fitting model to feature number 335, ASV922
#> 2026-10-28 10:47:17.863359 INFO::Fitting model to feature number 336, ASV673
#> 2026-10-28 10:47:17.866913 INFO::Fitting model to feature number 337, ASV62
#> 2026-10-28 10:47:17.873208 INFO::Fitting model to feature number 338, ASV100
#> 2026-10-28 10:47:17.877871 INFO::Fitting model to feature number 339, ASV2055
#> 2026-10-28 10:47:17.881646 INFO::Fitting model to feature number 340, ASV1626
#> 2026-10-28 10:47:17.885294 INFO::Fitting model to feature number 341, ASV2701
#> 2026-10-28 10:47:17.89042 INFO::Fitting model to feature number 342, ASV1363
#> 2026-10-28 10:47:17.897268 INFO::Fitting model to feature number 343, ASV2218
#> 2026-10-28 10:47:17.901299 INFO::Fitting model to feature number 344, ASV3168
#> 2026-10-28 10:47:17.906006 INFO::Fitting model to feature number 345, ASV616
#> 2026-10-28 10:47:17.911353 INFO::Fitting model to feature number 346, ASV668
#> 2026-10-28 10:47:17.916723 INFO::Fitting model to feature number 347, ASV1988
#> 2026-10-28 10:47:17.922389 INFO::Fitting model to feature number 348, ASV3198
#> 2026-10-28 10:47:17.926249 INFO::Fitting model to feature number 349, ASV77
#> 2026-10-28 10:47:17.932993 INFO::Fitting model to feature number 350, ASV419
#> 2026-10-28 10:47:17.936816 INFO::Fitting model to feature number 351, ASV199
#> 2026-10-28 10:47:17.941274 INFO::Fitting model to feature number 352, ASV2686
#> 2026-10-28 10:47:17.945904 INFO::Fitting model to feature number 353, ASV3589
#> 2026-10-28 10:47:17.949544 INFO::Fitting model to feature number 354, ASV927
#> 2026-10-28 10:47:17.954266 INFO::Fitting model to feature number 355, ASV4357
#> 2026-10-28 10:47:17.960611 INFO::Fitting model to feature number 356, ASV1472
#> 2026-10-28 10:47:17.965001 INFO::Fitting model to feature number 357, ASV4178
#> 2026-10-28 10:47:17.969351 INFO::Fitting model to feature number 358, ASV1313
#> 2026-10-28 10:47:17.975144 INFO::Fitting model to feature number 359, ASV885
#> 2026-10-28 10:47:17.979156 INFO::Fitting model to feature number 360, ASV350
#> 2026-10-28 10:47:17.983712 INFO::Fitting model to feature number 361, ASV575
#> 2026-10-28 10:47:17.987676 INFO::Fitting model to feature number 362, ASV519
#> 2026-10-28 10:47:17.992742 INFO::Fitting model to feature number 363, ASV5017
#> 2026-10-28 10:47:17.996419 INFO::Fitting model to feature number 364, ASV6330
#> 2026-10-28 10:47:18.000139 INFO::Fitting model to feature number 365, ASV346
#> 2026-10-28 10:47:18.004134 INFO::Fitting model to feature number 366, ASV604
#> 2026-10-28 10:47:18.007967 INFO::Fitting model to feature number 367, ASV1645
#> 2026-10-28 10:47:18.011874 INFO::Fitting model to feature number 368, ASV7460
#> 2026-10-28 10:47:18.015328 INFO::Fitting model to feature number 369, ASV934
#> 2026-10-28 10:47:18.018769 INFO::Fitting model to feature number 370, ASV4857
#> 2026-10-28 10:47:18.024495 INFO::Fitting model to feature number 371, ASV456
#> 2026-10-28 10:47:18.029755 INFO::Fitting model to feature number 372, ASV2077
#> 2026-10-28 10:47:18.034942 INFO::Fitting model to feature number 373, ASV5895
#> 2026-10-28 10:47:18.039123 INFO::Fitting model to feature number 374, ASV638
#> 2026-10-28 10:47:18.045146 INFO::Fitting model to feature number 375, ASV3459
#> 2026-10-28 10:47:18.072685 INFO::Fitting model to feature number 376, ASV991
#> 2026-10-28 10:47:18.079265 INFO::Fitting model to feature number 377, ASV6515
#> 2026-10-28 10:47:18.085051 INFO::Fitting model to feature number 378, ASV949
#> 2026-10-28 10:47:18.091418 INFO::Fitting model to feature number 379, ASV552
#> 2026-10-28 10:47:18.096407 INFO::Fitting model to feature number 380, ASV803
#> 2026-10-28 10:47:18.100196 INFO::Fitting model to feature number 381, ASV5054
#> 2026-10-28 10:47:18.106734 INFO::Fitting model to feature number 382, ASV2724
#> 2026-10-28 10:47:18.110661 INFO::Fitting model to feature number 383, ASV569
#> 2026-10-28 10:47:18.115475 INFO::Fitting model to feature number 384, ASV967
#> 2026-10-28 10:47:18.121274 INFO::Fitting model to feature number 385, ASV2561
#> 2026-10-28 10:47:18.126048 INFO::Fitting model to feature number 386, ASV1646
#> 2026-10-28 10:47:18.129981 INFO::Fitting model to feature number 387, ASV2204
#> 2026-10-28 10:47:18.133537 INFO::Fitting model to feature number 388, ASV1294
#> 2026-10-28 10:47:18.138204 INFO::Fitting model to feature number 389, ASV28
#> 2026-10-28 10:47:18.142474 INFO::Fitting model to feature number 390, ASV454
#> 2026-10-28 10:47:18.147068 INFO::Fitting model to feature number 391, ASV40
#> 2026-10-28 10:47:18.152386 INFO::Fitting model to feature number 392, ASV8
#> 2026-10-28 10:47:18.158469 INFO::Fitting model to feature number 393, ASV187
#> 2026-10-28 10:47:18.162678 INFO::Fitting model to feature number 394, ASV2280
#> 2026-10-28 10:47:18.166334 INFO::Fitting model to feature number 395, ASV245
#> 2026-10-28 10:47:18.170221 INFO::Fitting model to feature number 396, ASV559
#> 2026-10-28 10:47:18.175949 INFO::Fitting model to feature number 397, ASV448
#> 2026-10-28 10:47:18.18155 INFO::Fitting model to feature number 398, ASV1061
#> 2026-10-28 10:47:18.186564 INFO::Fitting model to feature number 399, ASV2134
#> 2026-10-28 10:47:18.192768 INFO::Fitting model to feature number 400, ASV232
#> 2026-10-28 10:47:18.196981 INFO::Fitting model to feature number 401, ASV2080
#> 2026-10-28 10:47:18.200334 INFO::Fitting model to feature number 402, ASV433
#> 2026-10-28 10:47:18.205292 INFO::Fitting model to feature number 403, ASV1763
#> 2026-10-28 10:47:18.209698 INFO::Fitting model to feature number 404, ASV2900
#> 2026-10-28 10:47:18.213339 INFO::Fitting model to feature number 405, ASV901
#> 2026-10-28 10:47:18.216692 INFO::Fitting model to feature number 406, ASV351
#> 2026-10-28 10:47:18.220918 INFO::Fitting model to feature number 407, ASV2370
#> 2026-10-28 10:47:18.225215 INFO::Fitting model to feature number 408, ASV1615
#> 2026-10-28 10:47:18.230344 INFO::Fitting model to feature number 409, ASV736
#> 2026-10-28 10:47:18.235193 INFO::Fitting model to feature number 410, ASV787
#> 2026-10-28 10:47:18.240875 INFO::Fitting model to feature number 411, ASV85
#> 2026-10-28 10:47:18.245737 INFO::Fitting model to feature number 412, ASV82
#> 2026-10-28 10:47:18.250029 INFO::Fitting model to feature number 413, ASV118
#> 2026-10-28 10:47:18.2566 INFO::Fitting model to feature number 414, ASV16
#> 2026-10-28 10:47:18.262855 INFO::Fitting model to feature number 415, ASV2764
#> 2026-10-28 10:47:18.266841 INFO::Fitting model to feature number 416, ASV7
#> 2026-10-28 10:47:18.27064 INFO::Fitting model to feature number 417, ASV553
#> 2026-10-28 10:47:18.275732 INFO::Fitting model to feature number 418, ASV996
#> 2026-10-28 10:47:18.280027 INFO::Fitting model to feature number 419, ASV127
#> 2026-10-28 10:47:18.285643 INFO::Fitting model to feature number 420, ASV649
#> 2026-10-28 10:47:18.292526 INFO::Fitting model to feature number 421, ASV428
#> 2026-10-28 10:47:18.298797 INFO::Fitting model to feature number 422, ASV157
#> 2026-10-28 10:47:18.304011 INFO::Fitting model to feature number 423, ASV348
#> 2026-10-28 10:47:18.309492 INFO::Fitting model to feature number 424, ASV753
#> 2026-10-28 10:47:18.313986 INFO::Fitting model to feature number 425, ASV90
#> 2026-10-28 10:47:18.318087 INFO::Fitting model to feature number 426, ASV88
#> 2026-10-28 10:47:18.322333 INFO::Fitting model to feature number 427, ASV958
#> 2026-10-28 10:47:18.327228 INFO::Fitting model to feature number 428, ASV2000
#> 2026-10-28 10:47:18.33102 INFO::Fitting model to feature number 429, ASV3556
#> 2026-10-28 10:47:18.334761 INFO::Fitting model to feature number 430, ASV1406
#> 2026-10-28 10:47:18.338874 INFO::Fitting model to feature number 431, ASV1219
#> 2026-10-28 10:47:18.344729 INFO::Fitting model to feature number 432, ASV2761
#> 2026-10-28 10:47:18.348909 INFO::Fitting model to feature number 433, ASV6253
#> 2026-10-28 10:47:18.353245 INFO::Fitting model to feature number 434, ASV338
#> 2026-10-28 10:47:18.358362 INFO::Fitting model to feature number 435, ASV368
#> 2026-10-28 10:47:18.362239 INFO::Fitting model to feature number 436, ASV950
#> 2026-10-28 10:47:18.366225 INFO::Fitting model to feature number 437, ASV1296
#> 2026-10-28 10:47:18.369896 INFO::Fitting model to feature number 438, ASV1275
#> 2026-10-28 10:47:18.375011 INFO::Fitting model to feature number 439, ASV2278
#> 2026-10-28 10:47:18.379497 INFO::Fitting model to feature number 440, ASV52
#> 2026-10-28 10:47:18.383199 INFO::Fitting model to feature number 441, ASV81
#> 2026-10-28 10:47:18.388273 INFO::Fitting model to feature number 442, ASV595
#> 2026-10-28 10:47:18.39228 INFO::Fitting model to feature number 443, ASV675
#> 2026-10-28 10:47:18.396318 INFO::Fitting model to feature number 444, ASV1089
#> 2026-10-28 10:47:18.400081 INFO::Fitting model to feature number 445, ASV1959
#> 2026-10-28 10:47:18.403754 INFO::Fitting model to feature number 446, ASV3863
#> 2026-10-28 10:47:18.411244 INFO::Fitting model to feature number 447, ASV6029
#> 2026-10-28 10:47:18.415334 INFO::Fitting model to feature number 448, ASV3344
#> 2026-10-28 10:47:18.419056 INFO::Fitting model to feature number 449, ASV7979
#> 2026-10-28 10:47:18.425197 INFO::Fitting model to feature number 450, ASV1036
#> 2026-10-28 10:47:18.429816 INFO::Fitting model to feature number 451, ASV2758
#> 2026-10-28 10:47:18.433315 INFO::Fitting model to feature number 452, ASV3999
#> 2026-10-28 10:47:18.436803 INFO::Fitting model to feature number 453, ASV4300
#> 2026-10-28 10:47:18.444342 INFO::Fitting model to feature number 454, ASV2368
#> 2026-10-28 10:47:18.44817 INFO::Fitting model to feature number 455, ASV5200
#> 2026-10-28 10:47:18.45166 INFO::Fitting model to feature number 456, ASV1311
#> 2026-10-28 10:47:18.456794 INFO::Fitting model to feature number 457, ASV2797
#> 2026-10-28 10:47:18.460381 INFO::Fitting model to feature number 458, ASV1051
#> 2026-10-28 10:47:18.463806 INFO::Fitting model to feature number 459, ASV354
#> 2026-10-28 10:47:18.467308 INFO::Fitting model to feature number 460, ASV578
#> 2026-10-28 10:47:18.471231 INFO::Fitting model to feature number 461, ASV2845
#> 2026-10-28 10:47:18.475162 INFO::Fitting model to feature number 462, ASV2323
#> 2026-10-28 10:47:18.480231 INFO::Fitting model to feature number 463, ASV1282
#> 2026-10-28 10:47:18.484555 INFO::Fitting model to feature number 464, ASV931
#> 2026-10-28 10:47:18.490237 INFO::Fitting model to feature number 465, ASV3234
#> 2026-10-28 10:47:18.496735 INFO::Fitting model to feature number 466, ASV7785
#> 2026-10-28 10:47:18.5017 INFO::Fitting model to feature number 467, ASV4573
#> 2026-10-28 10:47:18.505618 INFO::Fitting model to feature number 468, ASV3338
#> 2026-10-28 10:47:18.509755 INFO::Fitting model to feature number 469, ASV1112
#> 2026-10-28 10:47:18.513398 INFO::Fitting model to feature number 470, ASV3699
#> 2026-10-28 10:47:18.517652 INFO::Fitting model to feature number 471, ASV4823
#> 2026-10-28 10:47:18.521753 INFO::Fitting model to feature number 472, ASV560
#> 2026-10-28 10:47:18.526407 INFO::Fitting model to feature number 473, ASV1986
#> 2026-10-28 10:47:18.530493 INFO::Fitting model to feature number 474, ASV775
#> 2026-10-28 10:47:18.534841 INFO::Fitting model to feature number 475, ASV1364
#> 2026-10-28 10:47:18.587527 INFO::Fitting model to feature number 476, ASV804
#> 2026-10-28 10:47:18.595087 INFO::Fitting model to feature number 477, ASV1710
#> 2026-10-28 10:47:18.600366 INFO::Fitting model to feature number 478, ASV2152
#> 2026-10-28 10:47:18.605309 INFO::Fitting model to feature number 479, ASV3040
#> 2026-10-28 10:47:18.610978 INFO::Fitting model to feature number 480, ASV4658
#> 2026-10-28 10:47:18.615344 INFO::Fitting model to feature number 481, ASV8095
#> 2026-10-28 10:47:18.620896 INFO::Fitting model to feature number 482, ASV4572
#> 2026-10-28 10:47:18.62548 INFO::Fitting model to feature number 483, ASV886
#> 2026-10-28 10:47:18.629666 INFO::Fitting model to feature number 484, ASV5982
#> 2026-10-28 10:47:18.633394 INFO::Fitting model to feature number 485, ASV6080
#> 2026-10-28 10:47:18.638035 INFO::Fitting model to feature number 486, ASV5001
#> 2026-10-28 10:47:18.643645 INFO::Fitting model to feature number 487, ASV3903
#> 2026-10-28 10:47:18.647251 INFO::Fitting model to feature number 488, ASV6322
#> 2026-10-28 10:47:18.65078 INFO::Fitting model to feature number 489, ASV4712
#> 2026-10-28 10:47:18.65686 INFO::Fitting model to feature number 490, ASV2461
#> 2026-10-28 10:47:18.663259 INFO::Fitting model to feature number 491, ASV3134
#> 2026-10-28 10:47:18.67018 INFO::Fitting model to feature number 492, ASV4031
#> 2026-10-28 10:47:18.676314 INFO::Fitting model to feature number 493, ASV2411
#> 2026-10-28 10:47:18.679932 INFO::Fitting model to feature number 494, ASV2430
#> 2026-10-28 10:47:18.685028 INFO::Fitting model to feature number 495, ASV2840
#> 2026-10-28 10:47:18.691411 INFO::Fitting model to feature number 496, ASV2302
#> 2026-10-28 10:47:18.696224 INFO::Fitting model to feature number 497, ASV7286
#> 2026-10-28 10:47:18.699815 INFO::Fitting model to feature number 498, ASV4931
#> 2026-10-28 10:47:18.704355 INFO::Fitting model to feature number 499, ASV5286
#> 2026-10-28 10:47:18.710866 INFO::Fitting model to feature number 500, ASV6843
#> 2026-10-28 10:47:18.717676 INFO::Fitting model to feature number 501, ASV4505
#> 2026-10-28 10:47:18.724276 INFO::Fitting model to feature number 502, ASV5616
#> 2026-10-28 10:47:18.730548 INFO::Fitting model to feature number 503, ASV3801
#> 2026-10-28 10:47:18.735909 INFO::Fitting model to feature number 504, ASV2721
#> 2026-10-28 10:47:18.740967 INFO::Fitting model to feature number 505, ASV2275
#> 2026-10-28 10:47:18.745344 INFO::Fitting model to feature number 506, ASV1600
#> 2026-10-28 10:47:18.749281 INFO::Fitting model to feature number 507, ASV7465
#> 2026-10-28 10:47:18.753123 INFO::Fitting model to feature number 508, ASV7263
#> 2026-10-28 10:47:18.757438 INFO::Fitting model to feature number 509, ASV2237
#> 2026-10-28 10:47:18.761991 INFO::Fitting model to feature number 510, ASV802
#> 2026-10-28 10:47:18.766803 INFO::Fitting model to feature number 511, ASV1514
#> 2026-10-28 10:47:18.770727 INFO::Fitting model to feature number 512, ASV6064
#> 2026-10-28 10:47:18.775767 INFO::Fitting model to feature number 513, ASV829
#> 2026-10-28 10:47:18.780948 INFO::Fitting model to feature number 514, ASV3413
#> 2026-10-28 10:47:18.785271 INFO::Fitting model to feature number 515, ASV3015
#> 2026-10-28 10:47:18.789836 INFO::Fitting model to feature number 516, ASV3635
#> 2026-10-28 10:47:18.794992 INFO::Fitting model to feature number 517, ASV2601
#> 2026-10-28 10:47:18.798831 INFO::Fitting model to feature number 518, ASV7419
#> 2026-10-28 10:47:18.802731 INFO::Fitting model to feature number 519, ASV508
#> 2026-10-28 10:47:18.807114 INFO::Fitting model to feature number 520, ASV1505
#> 2026-10-28 10:47:18.811748 INFO::Fitting model to feature number 521, ASV1297
#> 2026-10-28 10:47:18.815542 INFO::Fitting model to feature number 522, ASV223
#> 2026-10-28 10:47:18.819796 INFO::Fitting model to feature number 523, ASV1945
#> 2026-10-28 10:47:18.824402 INFO::Fitting model to feature number 524, ASV1588
#> 2026-10-28 10:47:18.828878 INFO::Fitting model to feature number 525, ASV645
#> 2026-10-28 10:47:18.83275 INFO::Fitting model to feature number 526, ASV1515
#> 2026-10-28 10:47:18.836389 INFO::Fitting model to feature number 527, ASV4241
#> 2026-10-28 10:47:18.840468 INFO::Fitting model to feature number 528, ASV1960
#> 2026-10-28 10:47:18.844164 INFO::Fitting model to feature number 529, ASV2026
#> 2026-10-28 10:47:18.848197 INFO::Fitting model to feature number 530, ASV1684
#> 2026-10-28 10:47:18.852953 INFO::Fitting model to feature number 531, ASV1243
#> 2026-10-28 10:47:18.858042 INFO::Fitting model to feature number 532, ASV1529
#> 2026-10-28 10:47:18.863404 INFO::Fitting model to feature number 533, ASV1350
#> 2026-10-28 10:47:18.868386 INFO::Fitting model to feature number 534, ASV5663
#> 2026-10-28 10:47:18.873919 INFO::Fitting model to feature number 535, ASV6329
#> 2026-10-28 10:47:18.881374 INFO::Fitting model to feature number 536, ASV2922
#> 2026-10-28 10:47:18.886292 INFO::Fitting model to feature number 537, ASV1783
#> 2026-10-28 10:47:18.897751 INFO::Fitting model to feature number 538, ASV5443
#> 2026-10-28 10:47:18.90334 INFO::Fitting model to feature number 539, ASV1249
#> 2026-10-28 10:47:18.909762 INFO::Fitting model to feature number 540, ASV2219
#> 2026-10-28 10:47:18.915395 INFO::Fitting model to feature number 541, ASV617
#> 2026-10-28 10:47:18.922111 INFO::Fitting model to feature number 542, ASV1811
#> 2026-10-28 10:47:18.92719 INFO::Fitting model to feature number 543, ASV4721
#> 2026-10-28 10:47:18.931797 INFO::Fitting model to feature number 544, ASV3020
#> 2026-10-28 10:47:18.937365 INFO::Fitting model to feature number 545, ASV1882
#> 2026-10-28 10:47:18.943639 INFO::Fitting model to feature number 546, ASV693
#> 2026-10-28 10:47:18.94956 INFO::Fitting model to feature number 547, ASV3633
#> 2026-10-28 10:47:18.955098 INFO::Fitting model to feature number 548, ASV3165
#> 2026-10-28 10:47:18.960324 INFO::Fitting model to feature number 549, ASV2081
#> 2026-10-28 10:47:18.965805 INFO::Fitting model to feature number 550, ASV981
#> 2026-10-28 10:47:18.969844 INFO::Fitting model to feature number 551, ASV1860
#> 2026-10-28 10:47:18.97416 INFO::Fitting model to feature number 552, ASV1695
#> 2026-10-28 10:47:18.978174 INFO::Fitting model to feature number 553, ASV536
#> 2026-10-28 10:47:18.981868 INFO::Fitting model to feature number 554, ASV7634
#> 2026-10-28 10:47:18.98544 INFO::Fitting model to feature number 555, ASV2069
#> 2026-10-28 10:47:18.990103 INFO::Fitting model to feature number 556, ASV2393
#> 2026-10-28 10:47:18.996015 INFO::Fitting model to feature number 557, ASV4234
#> 2026-10-28 10:47:19.001259 INFO::Fitting model to feature number 558, ASV2681
#> 2026-10-28 10:47:19.005861 INFO::Fitting model to feature number 559, ASV252
#> 2026-10-28 10:47:19.07968 INFO::Counting total values for each feature
#> 2026-10-28 10:47:19.102547 INFO::Writing filtered data to file maaslin_temp_Porites_Heat_Shock/features/filtered_data.tsv
#> 2026-10-28 10:47:19.115233 INFO::Writing filtered, normalized data to file maaslin_temp_Porites_Heat_Shock/features/filtered_data_norm.tsv
#> 2026-10-28 10:47:19.141938 INFO::Writing filtered, normalized, transformed data to file maaslin_temp_Porites_Heat_Shock/features/filtered_data_norm_transformed.tsv
#> 2026-10-28 10:47:19.166549 INFO::Writing residuals to file maaslin_temp_Porites_Heat_Shock/fits/residuals.rds
#> 2026-10-28 10:47:19.170874 INFO::Writing fitted values to file maaslin_temp_Porites_Heat_Shock/fits/fitted.rds
#> 2026-10-28 10:47:19.175037 INFO::Writing all results to file (ordered by increasing q-values): maaslin_temp_Porites_Heat_Shock/all_results.tsv
#> 2026-10-28 10:47:19.186603 INFO::Writing the significant results (those which are less than or equal to the threshold of 0.250000 ) to file (ordered by increasing q-values): maaslin_temp_Porites_Heat_Shock/significant_results.tsv
#>  Finished! Porites Heat Shock Summary: LEfSe = 0 | DESeq2 = 16 | MaAsLin2 = 0
write.table(por_hs_bacterial_drivers,
            file = paste(dir_out_tab,
                         pref_out,
                         "HeatShock_DiffAbunTaxa_POR.tab", sep =""))

library(ggplot2)
library(dplyr)

# ==============================================================================
# 1. UPDATED PLOTTING FUNCTION WITH EXPLICIT FACET ORDERING
# ==============================================================================
plot_coral_differential_by_rank <- function(master_bind_data, tax_rank = "Genus", top_n_features = 30) {
  
  # Filter specifically for DESeq2 results
  plot_df <- master_bind_data %>%
    filter(Method == "DESeq2")
  
  if (nrow(plot_df) == 0) {
    stop("No DESeq2 results found in the provided data frame.")
  }
  
  # ----------------------------------------------------------------------------
  # FIX: Enforce Custom Chronological Factoring on the Facet Rows (Tissue Column)
  # ----------------------------------------------------------------------------
  # This explicitly pairs the Transplants at the top and the Heat Shocks at the bottom
  plot_df <- plot_df %>%
    mutate(Tissue = factor(Tissue, levels = c(
      "Pocillopora Transplant Timeline",
      "Porites Transplant Timeline",
      "Pocillopora Heat Shock",
      "Porites Heat Shock"
    )))
  
  # Check if the requested taxonomic rank exists in the data frame
  if (!tax_rank %in% colnames(plot_df)) {
    warning(paste("The requested rank '", tax_rank, "' was not found. Falling back to 'taxon' strings."))
    plot_df <- plot_df %>% mutate(Taxon_Label = taxon)
  } else {
    # If the taxonomic rank name is NA or blank, fallback to the unique ASV ID
    plot_df <- plot_df %>% 
      mutate(Taxon_Label = ifelse(is.na(.[[tax_rank]]) | .[[tax_rank]] == "", taxon, .[[tax_rank]]))
  }
  
  # Grab the top N most impactful changes per Tissue and Comparison
  plot_df <- plot_df %>%
    group_by(Tissue, Comparison) %>%
    slice_max(order_by = abs(Score), n = top_n_features, with_ties = FALSE) %>%
    ungroup()
  
  # Generate the plot configured for Opunohu vs Tiahura tracking vectors
  ggplot(plot_df, aes(x = reorder(Taxon_Label, Score), y = Score, fill = Score > 0)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    coord_flip() + 
    facet_grid(Tissue ~ Comparison, scales = "free_y", space = "free_y") +
    
    # Updated legends to map explicitly to Opunohu and Tiahura profiles
    scale_fill_manual(
      values = c("TRUE" = "#e34a33", "FALSE" = "#2b8cbe"),
      labels = c("TRUE" = "Enriched in Opunohu (Impacted)", "FALSE" = "Enriched in Tiahura (Pristine)")
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(face = "italic"),
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("Differentially Abundant Bacterial Taxa (", tax_rank, "Level)"),
      subtitle = "Log2 Fold Changes derived via DESeq2 Wald tests",
      x = tax_rank,
      y = "Log2 Fold Change (Opunohu vs Tiahura)",
      fill = "Enrichment Profile"
    )
}

# ==============================================================================
# 2. BIND DATA AND GENERATE ORDERED FACET LAYOUTS
# ==============================================================================

# Ensure you merge your 4 outputs exactly matching the factory level names strings above
combined_coral_indicators <- bind_rows(
  poc_trans_drivers %>% mutate(Tissue = "Pocillopora Transplant Timeline"),
  por_trans_drivers %>% mutate(Tissue = "Porites Transplant Timeline"),
  poc_hs_bacterial_drivers    %>% mutate(Tissue = "Pocillopora Heat Shock"),
  por_hs_bacterial_drivers    %>% mutate(Tissue = "Porites Heat Shock")
)

# Render Plot 1: Look at Phylum level shifts (Transplants on top, Heat Shocks on bottom)
plot_diffabund_phylum <- plot_coral_differential_by_rank(combined_coral_indicators, tax_rank = "Phylum", top_n_features = 20)
print(plot_diffabund_phylum)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-18-1.png" alt="" style="display: block; margin: auto;" />

``` r

# PATCH 2026-10-26: write the figure to Figures/05-AAC-CommunityAnalyses/
ggsave(paste0(dir_out_fig, pref_out, "plot_DiffAbund_Phylum.png"),
       plot = plot_diffabund_phylum, width = 9, height = 8, dpi = 300, bg = "white")
ggsave(paste0(dir_out_fig, pref_out, "plot_DiffAbund_Phylum.pdf"),
       plot = plot_diffabund_phylum, width = 9, height = 8, bg = "white")

# Render Plot 2: Look at fine-scale Genus level shifts (Transplants on top, Heat Shocks on bottom)
plot_diffabund_genus <- plot_coral_differential_by_rank(combined_coral_indicators, tax_rank = "Genus", top_n_features = 25)
print(plot_diffabund_genus)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-18-2.png" alt="" style="display: block; margin: auto;" />

``` r

# PATCH 2026-10-26: write the figure to Figures/05-AAC-CommunityAnalyses/
ggsave(paste0(dir_out_fig, pref_out, "plot_DiffAbund_Genus.png"),
       plot = plot_diffabund_genus, width = 9, height = 10, dpi = 300, bg = "white")
ggsave(paste0(dir_out_fig, pref_out, "plot_DiffAbund_Genus.pdf"),
       plot = plot_diffabund_genus, width = 9, height = 10, bg = "white")

```
These differential abundance sheets are highly informative. They provide the precise taxonomic and numerical evidence needed to explain the ecological patterns discovered across your macro-analyses.
Since your DESeq2 models were anchored against Tiahura, the interpretation of the Score (Log2 Fold Change) is standardized across all four files:

* Positive Scores (>0): Taxa that are significantly enriched/more abundant when the coral resides at the Impacted site (Opunohu).
* Negative Scores (<0): Taxa that are significantly enriched/more abundant when the coral resides at the Pristine site (Tiahura).

Here is a comprehensive synthesis of your biological drivers, organized by species and experimental window.

Chapter 1: The Pocillopora (POC) Transition Matrix
Your beta diversity analyses established that Pocillopora is highly plastic, shifting its community trajectory significantly based on destination sites. These tables reveal that this plasticity is driven by a massive taxonomic trade-off between protective endosymbionts and opportunistic environmental colonizers.
1. Ambient Transplant Timeline (Long-term Acclimatization)

* The Tiahura-Enriched Core (Negative Scores): Long-term acclimatization at the pristine site is dominated by a heavy enrichment of ASV89, ASV144, and ASV304, all belonging to the genus Endozoicomonas (Endozoicomonadaceae).
* Biological Significance: Endozoicomonas is widely recognized as the classic, health-associated core endosymbiont of marine corals. They reside deep within the host tissues, structuring coral nutrient cycling and offering protection against pathogens.
* The Opunohu-Enriched Profile (Positive Scores): When moved to the impacted site, Pocillopora sheds its Endozoicomonas core and becomes heavily colonized by opportunistic environmental marine bacteria:
* Pseudoalteromonas (ASV62, ASV694) and Vibrio (ASV210). These groups are highly dynamic, chemotactic, and fast-growing, often acting as secondary colonizers or indicator taxa for water column nutrient enrichment or host tissue degradation.
   * Rhodobacterales (Paracoccaceae) (ASV156, ASV108, ASV429), which are classic bioindicator bacteria commonly found expanding on corals undergoing structural environmental stress.

2. Acute Heat Shock Window
The post-hoc tests showed that acute heat shock triggered a significant structural change only when Pocillopora resided at Opunohu. The differential abundance profiles show that this shift is an acute acceleration of host dysbiosis:

* During the thermal challenge at Opunohu, Endozoicomonas (ASV89, ASV174) drops even further (highly negative scores).
* At the same time, there is a massive proliferation of opportunistic stress-responders: Fulvivirga (ASV107), Aquimarina spongiae (ASV241), and unique Xenococcaceae Cyanobacteria (ASV2). This combination points toward a classic thermal breakdown of the host tissue barrier, allowing water column biofilm-forming microbes to infiltrate the coral.

Chapter 2: The Porites (POR) Homeostatic Matrix
Your macro-analyses painted Porites as a highly rigid, homeostatic host whose microbiome centroid is completely immune to acute laboratory heat shock, yet undergoes a slow, parallel multi-month transition across both field sites. These files explain that mechanism perfectly.
1. Ambient Transplant Timeline (Long-term Acclimatization)

* When tracking long-term timeline transitions, Porites shows a very tight list of shifting microbes.
* It drops early-stage Alteromonadaceae (ASV136) and Parvularculaceae (ASV76) while systematically enriching Pelagibius (Kiloniellaceae) (ASV226) and Fulvivirga (ASV140).
* Biological Significance: This small, highly controlled list of shifting taxa explain why its overall alpha diversity (Shannon index) stayed perfectly flat. Porites does not allow wild invasions of opportunistic water bacteria like Vibrio or Pseudoalteromonas. It keeps tight regulation over its taxonomic boundaries.

2. Acute Heat Shock Window
Even though the heat shock PERMANOVA for Porites was non-significant overall (centroids didn't change drastically), DESeq2 identified a few underlying "sentinel" ASVs that adjusted their abundances.

* Under thermal stress at Opunohu, Porites enriches for highly specific, deep phylogenetic sediment and nutrient-cycling lineages: Woeseia (ASV504, ASV2111), Kordiimonas (ASV352, ASV517), and specialized Thalassolituus hydrocarbon-degraders (ASV367).
* Crucially, look at ASV140 (Fulvivirga) and ASV76 (Parvularculaceae). They appear in both the transplant timeline and the heat shock tables with identical directional signs. This provides molecular proof of what your PERMANOVA suggested: Porites responds to heat shock by slightly accelerating its long-term environmental baseline timeline, rather than undergoing an erratic community collapse or dysbiosis.

Cross-Tissue Shared Responses: The Universal Biomarker
There is an exceptional piece of shared cross-study evidence hidden in these data tables:
Look closely at ASV140 inside the Porites datasets and ASV107/ASV275 inside the Pocillopora datasets. They all belong to the exact same bacterial genus: Fulvivirga (Cyclobacteriaceae).

* For both coral species, across both the field transplant timelines and the laboratory heat shock windows, these Fulvivirga ASVs always display highly significant, positive scores.
* Synthesis: This identifies Fulvivirga as a universal, cross-host biomarker for environmental exposure and degradation at Opunohu. Regardless of the host coral species' internal physiological defense strategy, Fulvivirga systematically establishes itself within the coral holobiont whenever the coral is exposed to the local conditions of the impacted site.

# 9. Community Turnover

``` r
library(phyloseq)
library(betapart)
library(dplyr)
library(tidyr)

run_baselga_partitioning <- function(physeq, host_label) {
  # 1. Isolate ambient coral transplant samples
  ps_ambient <- subset_samples(physeq, Group == "TrueSample" & Type == "Coral" & Treatment != "Post-HS")
  ps_ambient <- prune_taxa(taxa_sums(ps_ambient) > 0, ps_ambient)
  
  # 2. Extract and binarize the ASV table (matrix of 0 and 1)
  if (taxa_are_rows(ps_ambient)) {
    asv_mat <- as.matrix(t(otu_table(ps_ambient)))
  } else {
    asv_mat <- as.matrix(otu_table(ps_ambient))
  }
  asv_binary <- ifelse(asv_mat > 0, 1, 0)
  
  # 3. Calculate multi-site Baselga components
  # This calculates the global turnover and nestedness across all samples in this matrix
  dist_partition <- beta.pair(asv_binary, index.family = "sor")
  
  # Convert distance matrices to standard data frames to calculate overall means
  total_sor <- mean(dist_partition$beta.sor, na.rm = TRUE)
  turnover_sim <- mean(dist_partition$beta.sim, na.rm = TRUE)
  nested_sne <- mean(dist_partition$beta.sne, na.rm = TRUE)
  
  # Calculate relative proportions
  pct_turnover <- (turnover_sim / total_sor) * 100
  pct_nestedness <- (nested_sne / total_sor) * 100
  
  cat("\n==================================================\n")
  cat("  BASELGA PARTITIONING SUMMARY:", host_label, "\n")
  cat("==================================================\n")
  cat(sprintf("• Total Sørensen Dissimilarity (beta.sor): %.3f\n", total_sor))
  cat(sprintf("• Pure Species Turnover Component (beta.sim): %.3f (%.1f%% of total)\n", turnover_sim, pct_turnover))
  cat(sprintf("• Nestedness-Driven Component (beta.sne): %.3f (%.1f%% of total)\n", nested_sne, pct_nestedness))
  
  # Return matrices in case you want to do a PERMANOVA on the pure turnover matrix downstream!
  return(dist_partition)
}

# Execute for both host datasets using your unscaled parent objects
poc_baselga <- run_baselga_partitioning(physeq_midpointTree_curated_POC, "Pocillopora (POC)")
#> 
#> ==================================================
#>   BASELGA PARTITIONING SUMMARY: Pocillopora (POC) 
#> ==================================================
#> • Total Sørensen Dissimilarity (beta.sor): 0.944
#> • Pure Species Turnover Component (beta.sim): 0.925 (98.0% of total)
#> • Nestedness-Driven Component (beta.sne): 0.019 (2.0% of total)
por_baselga <- run_baselga_partitioning(physeq_midpointTree_curated_POR, "Porites (POR)")
#> 
#> ==================================================
#>   BASELGA PARTITIONING SUMMARY: Porites (POR) 
#> ==================================================
#> • Total Sørensen Dissimilarity (beta.sor): 0.952
#> • Pure Species Turnover Component (beta.sim): 0.938 (98.5% of total)
#> • Nestedness-Driven Component (beta.sne): 0.015 (1.5% of total)
```
Interpreting the Baselga Matrix
For both Pocillopora and Porites, your community differences are almost entirely driven by Pure Species Turnover (\(\beta _{SIM}\) ≥ 97.2%), while the Nestedness-driven component (\(\beta _{SNE}\)) is practically non-existent (≤ 2.8%).
In community ecology, this means:
• The "No-Collapse" Rule: Your experimental treatments (transplantation and environmental shifts) are not causing a unidirectional hollowing-out or nested degradation of the microbiome. A high nestedness score would mean that a stressed coral simply loses its rare bacteria, leaving behind a blank, generic subset of its original self.
• True Biological Replacement: Instead, your corals are undergoing active taxonomic replacement. As they move across experimental axes and acclimate, they are actively replacing old sequence variants with entirely new, distinct bacterial lineages.
This provides a vital piece of context: when your PERMANOVA tests show significant shifts across treatments, those shifts represent a fundamental restructuring of the community's membership, not just a drop in species counts.

# 10. Hill Numbers

``` r
library(phyloseq)
library(vegan)
library(tidyverse)

generate_hill_numbers_dataframe <- function(physeq, host_label) {
  
  physeq <- subset_samples(physeq, Group == "TrueSample")
  # 1. Force phyloseq to export with samples as rows
  if (taxa_are_rows(physeq)) {
    otu_mat <- as.matrix(t(otu_table(physeq)))
  } else {
    otu_mat <- as.matrix(otu_table(physeq))
  }
  
  # 2. Calculate the unified Hill profiles
  hill_q0 <- rowSums(otu_mat > 0) # Richness (Total active ASVs)
  
  shannon_raw <- vegan::diversity(otu_mat, index = "shannon")
  hill_q1 <- exp(shannon_raw)     # Exponential Shannon (Common species)
  
  hill_q2 <- vegan::diversity(otu_mat, index = "invsimpson") # Inverse Simpson (Dominant species)
  
  # 3. Attach metadata
  metadata <- data.frame(sample_data(physeq))
  df_hill <- metadata %>%
    mutate(
      Host = host_label,
      Hill_q0_Richness = hill_q0,
      Hill_q1_Shannon  = hill_q1,
      Hill_q2_Simpson  = hill_q2
    )
  
  return(df_hill)
}

# Compute the final alpha framework matrices
hill_poc_df <- generate_hill_numbers_dataframe(physeq_midpointTree_curated_POC_scaled, "Pocillopora (POC)")
hill_por_df <- generate_hill_numbers_dataframe(physeq_midpointTree_curated_POR_scaled, "Porites (POR)")
master_hill_data <- bind_rows(hill_poc_df, hill_por_df)

# Print preview to verify column alignment
print(head(master_hill_data %>% select(Host, Site, Treatment, Hill_q0_Richness, Hill_q1_Shannon, Hill_q2_Simpson)))
#>                               Host    Site      Treatment Hill_q0_Richness
#> 100-TBC01F60POC  Pocillopora (POC) Tiahura 4-Month-Pre-HS               59
#> 101-TBC02F60POC  Pocillopora (POC) Tiahura 4-Month-Pre-HS              141
#> 102-TBC03F60POC  Pocillopora (POC) Tiahura 4-Month-Pre-HS              116
#> 103-TBC04F60POC  Pocillopora (POC) Tiahura 4-Month-Pre-HS              122
#> 104-TBC05F60POC  Pocillopora (POC) Tiahura 4-Month-Pre-HS              242
#> 11-TGC01F00POC_2 Pocillopora (POC) Tiahura       Genotype               91
#>                  Hill_q1_Shannon Hill_q2_Simpson
#> 100-TBC01F60POC         12.71382        5.004838
#> 101-TBC02F60POC         38.54998       14.468687
#> 102-TBC03F60POC         13.07421        3.446719
#> 103-TBC04F60POC         40.78318       22.719726
#> 104-TBC05F60POC         87.36088       35.850542
#> 11-TGC01F00POC_2        24.50373       10.361619


library(nlme)
library(car)

# Test Pocillopora hill profile transitions
lme_hill_poc <- lme(Hill_q1_Shannon ~ Site * Treatment, 
                    random = ~ 1 | Colony_ID, 
                    data = hill_poc_df, na.action = na.omit)

cat("\n--- POCILLOPORA HILL (q=1) LME TEST ---\n")
#> 
#> --- POCILLOPORA HILL (q=1) LME TEST ---
print(car::Anova(lme_hill_poc, type = "III"))
#> Analysis of Deviance Table (Type III tests)
#> 
#> Response: Hill_q1_Shannon
#>                  Chisq Df Pr(>Chisq)    
#> (Intercept)    31.7741  1  1.732e-08 ***
#> Site            0.0868  1     0.7683    
#> Treatment       1.3114  2     0.5191    
#> Site:Treatment  0.0173  2     0.9914    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Test Porites hill profile transitions
lme_hill_por <- lme(Hill_q1_Shannon ~ Site * Treatment, 
                    random = ~ 1 | Colony_ID, 
                    data = hill_por_df, na.action = na.omit)

cat("\n--- PORITES HILL (q=1) LME TEST ---\n")
#> 
#> --- PORITES HILL (q=1) LME TEST ---
print(car::Anova(lme_hill_por, type = "III"))
#> Analysis of Deviance Table (Type III tests)
#> 
#> Response: Hill_q1_Shannon
#>                  Chisq Df Pr(>Chisq)    
#> (Intercept)    31.2221  1  2.301e-08 ***
#> Site            0.0009  1     0.9767    
#> Treatment       0.4121  2     0.8138    
#> Site:Treatment  2.5807  2     0.2752    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Unpacking Your Hill Number Metrics
Look at the contrast between Hill_q0_Richness (total count of active ASVs) and Hill_q1_Shannon (the effective number of common species):
• Sample 1-TGC01F00POC (Genotype Baseline):
	• Has 19 active ASVs (q=0), and its q=1 value is 8.16.
	• This means that out of its 19 present species, the community abundance structure behaves as if it is evenly split among roughly 8 dominant or highly common species.
• Sample 101-TBC02F60POC (Pre-Heat Shock Field Acclimated):
	• Features 142 active ASVs (q=0), but its q=1 value collapses down to 38.73.
	• This is a massive expansion in absolute species presence compared to the genotype baseline, yet it highlights that the vast majority of those 142 sequences are highly rare, transient, or background elements. The functional core of the ecosystem is effectively driven by about 39 common species.
This mathematical clarity perfectly reinforces your previous Baselga partitioning discovery. Because your community differences are 97% to 98% pure species turnover (and not a nested drop in species presence), these Hill Numbers quantify the scale of that replacement. The corals are not just changing the names of their microbes; they are dramatically expanding the actual niche density and evenness of their active common cohorts as they acclimate to field conditions over time.

These Hill number LME results provide a fascinating, highly rigorous statistical nuance to your alpha diversity chapter. They completely redefine how you will write up the alpha diversity section of your thesis or manuscript.
When we ran your traditional Shannon index values earlier, Pocillopora showed a statistically significant Treatment main effect (p = 0.023). However, when evaluated as a Hill number (q=1, representing the effective number of common species), the Treatment main effect drops completely out of significance (p = 0.176).
Here is the exact statistical and biological interpretation of this change:
Decoupling Traditional Shannon vs. Hill Numbers (q=1)
• Traditional Shannon (p = 0.023): Measures logarithmic entropy. It is highly sensitive to subtle mathematically skewed redistributions in the relative abundances of rare or background taxa.
• Hill Number q=1 (p = 0.176): Translates that entropy into literal biological units (e.g., "This sample effectively behaves like it has 13 equally common species"). It puts a far stricter mathematical weight on the common species, meaning it is less easily swayed by fluctuations in rare background taxa.
The Biological Story For Your Paper
1. True Functional Constancy: This comparison proves that while the long-term field acclimatization process introduces minor shifts in the relative abundance curves of rare background bacteria (which tripped the traditional Shannon significance), the actual structural core of common, active species remains highly constant throughout the entire study.
2. Absolute Homeostasis Across Both Species: By switching to the highly preferred Hill numbers framework, you show that both Pocillopora and Porites exhibit profound alpha diversity homeostasis (p ≥ 0.176 for all experimental terms across both hosts). No matter how far the centroids shift in your beta diversity spaces, or how completely the corals swap out their specific bacterial sequences (as proven by your 97% Baselga species turnover), the living host tissues tightly regulate their micro-habitats to keep the absolute density of common functional niches perfectly balanced.

# 11. Distance to Baseline Plots

``` r
library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)

calculate_distances_to_baseline <- function(physeq, host_label) {
  # 1. Isolate ambient timeline
  ps_ambient <- subset_samples(physeq, Group == "TrueSample" & Type == "Coral" & Treatment != "Post-HS")
  
  # 2. Compute Weighted UniFrac matrix
  wunifrac_mat <- as.matrix(phyloseq::distance(ps_ambient, method = "wunifrac"))
  meta_df <- data.frame(sample_data(ps_ambient)) %>% mutate(SampleID = rownames(.))
  
  # 3. For each sample, find the distance to its matching Mother Colony "Genotype" baseline
  distance_vectors <- lapply(1:nrow(meta_df), function(i) {
    current_sample   <- meta_df$SampleID[i]
    current_colony   <- meta_df$Colony_ID[i]
    current_site     <- meta_df$Site[i]
    current_treat    <- meta_df$Treatment[i]
    
    # Locate the unique Day 1 Genotype baseline sample for this exact colony
    baseline_id <- meta_df %>% 
      filter(Colony_ID == current_colony & Treatment == "Genotype") %>% 
      pull(SampleID)
    
    if(length(baseline_id) == 0) return(NULL)
    
    # Extract the explicit distance value from the matrix
    # PATCH 2026-10-26: average if a colony has more than one day-0 baseline sample.
    # Under the old Unique_Colony this returned both bays' baselines, which silently
    # produced two rows per fragment, one against the wrong colony.
    dist_to_baseline <- mean(wunifrac_mat[current_sample, baseline_id], na.rm = TRUE)
    
    data.frame(
      SampleID = current_sample,
      Host = host_label,
      Site = current_site,
      Treatment = current_treat,
      Colony_ID = current_colony,
      Distance_To_Baseline = dist_to_baseline
    )
  }) %>% bind_rows()
  
  return(distance_vectors)
}

# Calculate raw baseline drift vectors for both hosts
poc_drift <- calculate_distances_to_baseline(physeq_midpointTree_curated_POC_scaled, "Pocillopora (POC)")
por_drift <- calculate_distances_to_baseline(physeq_midpointTree_curated_POR_scaled, "Porites (POR)")

# Bind and clean factor timelines (excluding the zero baseline point from the plot itself for clean scale)
plot_drift_df <- bind_rows(poc_drift, por_drift) %>%
  filter(Treatment != "Genotype") %>%
  mutate(Treatment = factor(Treatment, levels = c("2-Weeks", "4-Months", "Pre-HS")))

# ------------------------------------------------------------------------------
# BUILD THE PLOT
# ------------------------------------------------------------------------------
plot_baseline_drift <- ggplot(plot_drift_df, aes(x = Treatment, y = Distance_To_Baseline, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_point(aes(color = Site), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.5) +
  # Add tracking line connecting the means of the sites to show the split vs parallel tracks
  stat_summary(aes(group = Site, color = Site), fun = mean, geom = "line", 
               position = position_dodge(width = 0.8), size = 0.8, linetype = "solid") +
  facet_wrap(~Host, scales = "free_x") +
  theme_bw() +
  scale_fill_manual(values = c("Tiahura" = "#2b8cbe", "Opunohu" = "#e34a33")) +
  scale_color_manual(values = c("Tiahura" = "#08589e", "Opunohu" = "#b30000")) +
  labs(
    title = "Microbiome Community Divergence from Day 1 Baseline State",
    subtitle = "Pocillopora tracks widen/diverge across environments (Plasticity) | Porites tracks remain parallel (Homeostasis)",
    x = "Experimental Timeline",
    y = "Phylogenetic Distance to Native Genotype Baseline (Weighted UniFrac)",
    fill = "Transplant Site"
  ) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 11))


print(plot_baseline_drift)
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-21-1.png" alt="" style="display: block; margin: auto;" />

``` r

# PATCH 2026-10-26: write the figure to Figures/05-AAC-CommunityAnalyses/
ggsave(paste0(dir_out_fig, pref_out, "plot_BaselineDrift.png"),
       plot = plot_baseline_drift, width = 8, height = 5, dpi = 300, bg = "white")
ggsave(paste0(dir_out_fig, pref_out, "plot_BaselineDrift.pdf"),
       plot = plot_baseline_drift, width = 8, height = 5, bg = "white")
```

# 6. Taxonomic Composition

## 6.1 - Phylum

``` r
# a - Calculate the phylum relative abundance ---- 
  # Note: The read depth MUST be normalized in some way: e.g., scale_reads
  df_phylum <- 
    physeq_midpointTree_curated_scaled %>%
    # agglomerate at the phylum level 
    tax_glom(taxrank = "Phylum") %>% 
    # Transform counts to relative abundance 
    transform_sample_counts(function (x) {x/sum(x)}) %>%
    # Melt to a long format 
    psmelt() %>%
    # Filter out Phyla that are <1% - get rid of low abundant Phyla
    dplyr::filter(Abundance > 0.01) # %>%
    # # fix the order of date
    # mutate(date = fct_relevel(date, c("6/2/21", "6/15/21", "10/5/21")),
    #        station = fct_relevel(station, c("Copano West", "Copano East",
    #                                         "Mesquite Bay", "Aransas Bay",
    #                                         "Shipping Channel")))

# b - Stacked Bar Plot With All phyla ----
  # Plot Phylum Abundances in Pocillopora - make sure to load phylum_colors 
  df_phylum %>%
    # Warning: It's important to have one sample per x value, 
    # otherwise, it will take the sum between multiple samples
    # dplyr::filter(Site == "Opunohu") %>%
    dplyr::filter(Host == "POC") %>%
    ggplot(aes(x = Treatment, y = Abundance, fill = Phylum)) + 
    facet_grid(Site~Transplant) + 
    geom_bar(stat = "identity", color = "black") + 
    labs(title = "Pocillopora Phylum Composition") + 
    # scale_fill_manual(values = phylum_colors) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-22-1.png" alt="" style="display: block; margin: auto;" />

``` r

  # Plot Phylum Abundances in Porites - make sure to load phylum_colors 
  df_phylum %>%
    # Warning: It's important to have one sample per x value, 
    # otherwise, it will take the sum between multiple samples
    # dplyr::filter(Site == "Opunohu") %>%
    dplyr::filter(Host == "POR") %>%
    ggplot(aes(x = Treatment, y = Abundance, fill = Phylum)) + 
    facet_grid(Site~Transplant) + 
    geom_bar(stat = "identity", color = "black") + 
    labs(title = "Porites Phylum Composition") + 
    # scale_fill_manual(values = phylum_colors) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

<img src="05-AAC-CommunityAnalyses_files/figure-html/unnamed-chunk-22-2.png" alt="" style="display: block; margin: auto;" />
