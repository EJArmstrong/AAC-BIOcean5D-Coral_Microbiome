#!/usr/bin/env bash
# =============================================================================
# 07-AAC-RNAseq_quantification.sh
# Pocillopora meandrina heat-shock RNAseq: raw reads -> gene-level counts
#
# Runs OUTSIDE this repository (BioHPC / workstation shell). Nothing it writes
# belongs in git except the small per-sample quant.sf files and the QC summary;
# FASTQs, the reference and the salmon index are gitignored - see .gitignore.
#
# Library layout: 34 libraries, SINGLE-END (R1 only), each sequenced TWICE.
#   <SampleID>_<S>_L001_R1_001.fastq.gz        (run 1)
#   <SampleID>-run2_<S>_L001_R1_001.fastq.gz   (run 2)
# The two runs are technical replicates of one library and are CONCATENATED
# before quantification - do not treat them as biological replicates.
# =============================================================================
set -euo pipefail

# ------------------------------- configuration -------------------------------
RAW=/path/to/BIOcean5D/Data_Sequences/True_Sequences_RNAseq   # 68 .fastq.gz, ~41 GB
WORK=/workdir/$USER/POC_RNAseq                                 # scratch, not in repo
THREADS=12

CAT=$WORK/01-concatenated
TRIM=$WORK/02-trimmed
REF=$WORK/ref
IDX=$WORK/salmon_index
QUANT=$WORK/03-salmon
mkdir -p "$CAT" "$TRIM" "$REF" "$QUANT"

# ------------------------- 1. concatenate the two runs -----------------------
# Sample ID is the filename up to the first underscore, minus any -run2 suffix.
for f in "$RAW"/*_L001_R1_001.fastq.gz; do
  b=$(basename "$f")
  s=${b%%_*}; s=${s%-run2}
  echo "$s"
done | sort -u > "$WORK/sample_ids.txt"
wc -l < "$WORK/sample_ids.txt"          # expect 34

while read -r s; do
  if [[ ! -s "$CAT/${s}_R1.fastq.gz" ]]; then
    # gzip streams concatenate directly; order is irrelevant for quantification
    cat "$RAW/${s}_"*_L001_R1_001.fastq.gz "$RAW/${s}-run2_"*_L001_R1_001.fastq.gz \
        > "$CAT/${s}_R1.fastq.gz"
  fi
done < "$WORK/sample_ids.txt"

# ----------------------------- 2. QC and trimming ----------------------------
# fastp on single-end reads: adapter auto-detection, polyG trim (NovaSeq/NextSeq),
# quality trimming, and a per-sample JSON that feeds the QC table in the Rmd.
while read -r s; do
  [[ -s "$TRIM/${s}_R1.fastq.gz" ]] && continue
  fastp -i "$CAT/${s}_R1.fastq.gz" -o "$TRIM/${s}_R1.fastq.gz" \
        --trim_poly_g --cut_tail --cut_tail_mean_quality 20 --length_required 36 \
        --thread "$THREADS" \
        --json "$TRIM/${s}.fastp.json" --html "$TRIM/${s}.fastp.html"
done < "$WORK/sample_ids.txt"

# --------------------------- 3. reference transcriptome ----------------------
# Pocillopora meandrina PMEA_v1, Genoscope CEA (Tara Pacific), BioProject
# PRJEB51619: 347.2 Mb in 252 contigs, contig N50 4.75 Mb, 32,095 annotated
# protein-coding genes. NOTE: this is the ONLY P. meandrina assembly in NCBI.
if [[ ! -s "$REF/pmea_transcripts.fa" ]]; then
  cd "$REF"
  datasets download genome accession GCA_942486045.1 \
      --include genome,gff3,cds,protein --filename pmea.zip
  unzip -o pmea.zip
  GENOME=$(find . -name "*_genomic.fna" | head -1)
  GFF=$(find . -name "genomic.gff" | head -1)
  # Transcript fasta keyed by transcript ID (gffread writes the mRNA ID).
  gffread "$GFF" -g "$GENOME" -w pmea_transcripts.fa
  # tx2gene: mRNA ID -> parent gene ID, consumed by tximport in the Rmd.
  awk -F'\t' '$3=="mRNA"{
        split($9,a,";"); id=""; par="";
        for(i in a){ if(a[i]~/^ID=/) id=substr(a[i],4); if(a[i]~/^Parent=/) par=substr(a[i],8) }
        print id"\t"par }' "$GFF" | sort -u > tx2gene.tsv
  wc -l tx2gene.tsv
  # Decoys: the full genome, so reads from unannotated/intronic regions and from
  # the Cladocopium symbiont are far less likely to be force-assigned to a host
  # transcript. This matters here - the holobiont carries a large symbiont
  # transcriptome that has no host locus.
  grep "^>" "$GENOME" | cut -d" " -f1 | sed 's/^>//' > decoys.txt
  cat pmea_transcripts.fa "$GENOME" > gentrome.fa
fi

# ------------------------------ 4. salmon index ------------------------------
if [[ ! -d "$IDX" ]]; then
  salmon index -t "$REF/gentrome.fa" -d "$REF/decoys.txt" -i "$IDX" \
               -k 31 -p "$THREADS"
fi

# ------------------------------ 5. quantify ----------------------------------
# -l A        : let salmon infer library type
# -r          : single-end reads
# --seqBias   : sequence-specific bias (do NOT use --gcBias, which is paired-end only)
# --validateMappings is the default from salmon 1.0 onward.
while read -r s; do
  [[ -s "$QUANT/$s/quant.sf" ]] && continue
  salmon quant -i "$IDX" -l A -r "$TRIM/${s}_R1.fastq.gz" \
               -p "$THREADS" --seqBias --validateMappings \
               -o "$QUANT/$s"
done < "$WORK/sample_ids.txt"

# ------------------------- 6. mapping-rate QC summary ------------------------
# Mapping rate is a QC metric AND a biological one here: host mapping rate falls
# as the symbiont fraction rises, so it is a crude proxy for symbiont load and
# should be checked against Pre-HS / Post-HS before it is treated as noise.
{
  echo -e "Sample_ID\tnum_processed\tnum_mapped\tpercent_mapped"
  while read -r s; do
    m=$QUANT/$s/aux_info/meta_info.json
    python3 - "$s" "$m" <<'PY'
import json, sys
s, m = sys.argv[1], sys.argv[2]
d = json.load(open(m))
print(f"{s}\t{d['num_processed']}\t{d['num_mapped']}\t{100*d['num_mapped']/d['num_processed']:.2f}")
PY
  done < "$WORK/sample_ids.txt"
} > "$QUANT/mapping_rates.tsv"

# ---------------------------- 7. copy into the repo --------------------------
# Only the small files. quant.sf is ~1-2 MB per sample (34 samples, ~50 MB total):
# if that is more than you want in git history, keep --no-quant-sf and instead copy
# only the gene-level matrix written by the Rmd.
REPO=/path/to/AAC-BIOcean5D-Coral_Microbiome
mkdir -p "$REPO/Data/07-AAC-RNAseq/salmon"
while read -r s; do
  mkdir -p "$REPO/Data/07-AAC-RNAseq/salmon/$s"
  cp "$QUANT/$s/quant.sf" "$REPO/Data/07-AAC-RNAseq/salmon/$s/"
done < "$WORK/sample_ids.txt"
cp "$QUANT/mapping_rates.tsv" "$REPO/Data/07-AAC-RNAseq/"
cp "$REF/tx2gene.tsv"         "$REPO/Data/07-AAC-RNAseq/"

echo "Done. Next: knit Analysis/07-AAC-RNAseq/07-AAC-RNAseq.Rmd"
