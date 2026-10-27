# 07 - AAC - Pocillopora meandrina heat-shock RNAseq

**Status: drafted, not run.** The count data do not exist yet.

## What runs where

| Stage | Where | Script |
|---|---|---|
| Concatenate the two runs per library, trim, build the salmon index, quantify | workstation or BioHPC - **not** in this repo, and not in a sandbox | `Code/07-AAC-RNAseq_quantification.sh` |
| tximport, QC, DESeq2, figures | RStudio, from this folder | `07-AAC-RNAseq.Rmd` |

The read set is **41 GB** (68 single-end FASTQs) and the reference is 347 Mb, so neither
can be held here; `.gitignore` blocks both. Only the per-sample `quant.sf` files,
`tx2gene.tsv` and `mapping_rates.tsv` are copied back into `Data/07-AAC-RNAseq/`.

## Library layout - read this before touching the design

34 libraries, all *Pocillopora meandrina*, **single-end (R1 only)**, each sequenced twice:

```
<SampleID>_<S>_L001_R1_001.fastq.gz        run 1
<SampleID>-run2_<S>_L001_R1_001.fastq.gz   run 2
```

The two runs are **technical replicates of one library** and are concatenated before
quantification. Treating them as separate samples would inflate n from 34 to 68 and
invalidate every p-value.

## Design

|  | held at Opunohu | held at Tiahura | total |
|---|---|---|---|
| Pre-HS (F60) | 10 | 10 | 20 |
| Post-HS (F80) | 5 | 9 | 14 |

Pre-HS is complete across all 10 colonies at both bays. Post-HS is missing six fragments,
and the loss is **bay-biased** - 5 Opunohu against 9 Tiahura. Destination and interaction
terms therefore carry survivorship confounding that no model removes.

## Reference

*Pocillopora meandrina* **PMEA_v1**, accession **GCA_942486045.1**, submitted by Genoscope
CEA under BioProject PRJEB51619 (Tara Pacific). Contig-level: 347,233,126 bp in 252
contigs, contig N50 4.75 Mb, with a Genoscope annotation of 32,095 protein-coding genes.
This is the only *P. meandrina* assembly in NCBI. The annotation carries **no GO or KEGG
mapping**, so functional enrichment needs an orthology step (eggNOG-mapper or
InterProScan on the provided protein fasta) before the exported ranked tables can be used.

## Why a decoy-aware salmon index

These are holobiont libraries: reads not from the coral are largely *Cladocopium*. Indexing
the host transcriptome with the whole host genome as decoy keeps those reads from being
force-assigned to host transcripts. The host mapping rate is then retained as a crude
**symbiont-fraction proxy** and tested against timepoint in the Rmd rather than discarded
as a QC number - a fall in host mapping rate after heat shock would itself be evidence of
a shift in symbiont load.

## Link to script 06

Script 06 found that acute heat shock increased among-fragment variability in the
Symbiodiniaceae community ~1.75-fold in this same host without moving the mean community,
while the bacterial community did not destabilise. Two questions carry over, and the Rmd
addresses both: whether the canonical heat-shock machinery was actually induced (the ITS2
interpretation currently rests on treatment rather than measured stress), and whether host
expression shows the same rise in among-fragment dispersion.
