#!/bin/bash
# Unknown Bacterial Isolate Analysis Pipeline
# Oxford Nanopore Sequencing Data
# Author: Sahre Hilal Sezer
# Date: 2026-07-27

set -e

# ── VARIABLES ────────────────────────────────
READS="data/unknown_isolate.fastq.gz"
ASSEMBLY="results/assembly/assembly.fasta"
THREADS=4

# ── STEP 1: QC ───────────────────────────────
echo "[1/6] Running NanoPlot QC..."
NanoPlot --fastq $READS \
         --outdir results/qc \
         --threads $THREADS

# ── STEP 2: ASSEMBLY ─────────────────────────
echo "[2/6] Running Flye assembly..."
flye --nano-hq $READS \
     --out-dir results/assembly \
     --threads $THREADS

# ── STEP 3: ASSEMBLY QC ──────────────────────
echo "[3/6] Running QUAST assembly QC..."
python -m quast $ASSEMBLY \
       -o results/qc/quast

# ── STEP 4: TAXONOMIC IDENTIFICATION ─────────
# BLASTn was performed via NCBI web interface:
# https://blast.ncbi.nlm.nih.gov/
# Query: contig_4 (chromosome) against NCBI nt database
# Result: Klebsiella pneumoniae strain 18471 (OZ185751.1), 98.14% identity
echo "[4/6] Taxonomic ID: performed via NCBI BLASTn web (see results/taxonomy/blast_result.txt)"

# ── STEP 5: MLST ─────────────────────────────
echo "[5/6] Running MLST..."
mkdir -p results/typing
mlst --scheme klebsiella $ASSEMBLY > results/typing/mlst_result.txt

# ── STEP 6: AMR ANALYSIS ─────────────────────
# ResFinder was run via web interface:
# https://cge.food.dtu.dk/services/ResFinder/
# PlasmidFinder was run via web interface:
# https://cge.food.dtu.dk/services/PlasmidFinder/
# Results saved to: results/amr/resfinder/
echo "[6/6] AMR + Plasmid typing: performed via web tools (see results/amr/resfinder/)"

# ── EXTENDED TYPING (Kleborate) ───────────────
echo "Running Kleborate extended typing..."
mkdir -p results/typing
kleborate -a $ASSEMBLY -o results/typing/ \
  -m klebsiella_pneumo_complex__mlst,klebsiella_pneumo_complex__amr,\
klebsiella_pneumo_complex__kaptive,klebsiella_pneumo_complex__virulence_score,\
klebsiella_pneumo_complex__resistance_score,klebsiella__ybst,klebsiella__cbst,\
klebsiella__abst,klebsiella_pneumo_complex__wzi

# ── VISUALIZATIONS & REPORT ───────────────────
echo "Generating visualizations..."
PYTHONPATH=code python code/visualize.py

echo "Generating HTML clinical report..."
PYTHONPATH=code python code/report_generator.py

echo ""
echo "Pipeline complete."
echo "Key outputs:"
echo "  QC report      : results/qc/NanoPlot-report.html"
echo "  Assembly QC    : results/qc/quast/report.html"
echo "  MLST result    : results/typing/mlst_result.txt"
echo "  Kleborate      : results/typing/klebsiella_pneumo_complex__mlst_output.txt"
echo "  AMR summary    : results/amr/amr_summary.html"
echo "  Clinical report: results/report/final_report.html"
echo "  Findings       : results/amr/findings.md"
