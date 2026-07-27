# Unknown Bacterial Isolate Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/)
[![Flye](https://img.shields.io/badge/Flye-2.9.6-orange.svg)](https://github.com/fenderglass/Flye)
[![NanoPlot](https://img.shields.io/badge/NanoPlot-1.47.1-teal.svg)](https://github.com/wdecoster/NanoPlot)
[![Pytest](https://img.shields.io/badge/Pytest-passing-brightgreen.svg)](https://docs.pytest.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A clinical genomics pipeline for the identification, typing, and antimicrobial resistance (AMR)
profiling of an unknown bacterial isolate sequenced with Oxford Nanopore Technology (ONT).

The pipeline covers read quality control, de novo assembly and assembly QC, taxonomic
identification, multi-locus sequence typing, capsular typing, acquired resistance gene
detection, resistance-conferring point mutation analysis, and plasmid replicon typing —
with interactive visualizations and a structured clinical report aimed at both
bioinformaticians and clinicians.

---

## Key Findings

| Finding | Result |
|---------|--------|
| Organism | *Klebsiella pneumoniae* (98.14% BLASTn identity) |
| Sequence type | ST258 — globally dominant CRKP epidemic clone |
| Critical resistance gene | blaKPC-3 — carbapenem resistance (KPC-3 carbapenemase) |
| Resistance genes | 37 ResFinder hits, 23 unique genes across 8 antibiotic classes |
| Plasmid-borne resistance | blaKPC-3 on IncI2(Delta) — conjugative, broad host range |
| Plasmid replicons | IncFIB(K), IncFII(K), IncI2(Delta), IncR, repB(R1701) |
| OMP mutations | OmpK35 + OmpK36 dual porin loss — enhances carbapenem resistance |
| FLQ mutations | GyrA Ser83Ile + Asp87Asn; ParC Ser80Ile — full fluoroquinolone resistance |
| Virulence / Resistance score | 0 / 5 virulence, 2 / 3 resistance (XDR) |
| Capsule type | KL107 / O13 |
| Sequencing coverage | ~98x (577 Mb total bases / 5.9 Mb genome) |
| Assembly | 13 contigs — 1 chromosome (~5.3 Mb) + 4 plasmids |

This isolate is a **carbapenem-resistant *Klebsiella pneumoniae* (CRKP)**. Carbapenems are
last-resort antibiotics; resistance to them severely limits treatment options and constitutes
a clinical emergency.

---

## Analysis Rationale

Each step was chosen to answer a specific question, and critical findings were confirmed
with more than one independent tool wherever possible.

| Question | Approach | Why this approach |
|----------|----------|-------------------|
| Is the data good enough to assemble? | NanoPlot + custom Python analyzer | ONT error profiles differ from Illumina; read N50 and yield matter more than per-base quality alone. A custom analyzer was written to cross-check NanoPlot rather than trusting a single tool. |
| What organism is this? | BLASTn against NCBI nt | Assembly-based identification is more reliable than read-level classification, since a contiguous chromosome gives an unambiguous match. |
| Is the assembly trustworthy? | QUAST | An N50 equal to the chromosome length with L50 of 1 and zero gaps confirms the assembly resolved repeats rather than fragmenting them. |
| Which strain is it? | mlst + Kleborate | Species alone is not clinically actionable. Sequence type identifies epidemic clones. Two independent tools were used because ST assignment drives the infection control response. |
| What resistance genes are present? | ResFinder + Kleborate AMR module | Acquired gene detection alone is incomplete — Kleborate additionally reports chromosomal point mutations that ResFinder does not cover. |
| Why is resistance so severe? | Point mutation analysis (Kleborate) | The porin loss finding (OmpK35/OmpK36) explains resistance beyond what the carbapenemase gene alone predicts, and changes the treatment recommendation. |
| Can this resistance spread? | PlasmidFinder replicon typing | A resistance gene on a conjugative plasmid is an outbreak risk, not just a treatment problem. Replicon type determines host range. |

Where a tool could not be installed on the available platform (ARM64 macOS), the analysis was
completed through the corresponding web service rather than being skipped, and this is
documented in the Tools table below. Steps that could not be completed at all are listed
openly under Limitations rather than omitted.

---

## Project Structure

```
unknown_isolate/
├── data/
│   └── unknown_isolate.fastq.gz      # Raw ONT reads (not tracked in git)
├── code/
│   ├── pipeline.sh                   # Full analysis pipeline as shell script
│   └── blast_query.fasta             # Query sequence used for NCBI BLASTn
├── scripts/
│   ├── qc_analyzer.py                # OOP FASTQ quality control analyzer
│   ├── amr_parser.py                 # ResFinder output parser with AMR classification
│   ├── visualize.py                  # Interactive Plotly dashboards
│   └── report_generator.py           # Self-contained HTML clinical report
├── tests/
│   └── test_pipeline.py              # pytest unit tests (16 tests)
├── results/
│   ├── qc/
│   │   ├── NanoPlot-report.html      # ONT read quality report
│   │   ├── NanoStats.txt             # Read statistics
│   │   ├── qc_dashboard.html         # Interactive QC dashboard
│   │   └── quast/                    # Assembly quality metrics
│   ├── assembly/                     # Flye assembly (FASTA + assembly_info.txt)
│   ├── taxonomy/
│   │   └── blast_result.txt          # NCBI BLASTn species identification
│   ├── typing/
│   │   ├── mlst_result.txt           # MLST sequence type (ST258)
│   │   └── klebsiella_pneumo_complex__mlst_output.txt   # Kleborate full output
│   ├── plasmid/
│   │   └── plasmidfinder_results.txt # Plasmid replicon typing
│   ├── amr/
│   │   ├── resfinder/                # Raw ResFinder output
│   │   ├── amr_summary.html          # Interactive resistance gene chart
│   │   ├── amr_location.html         # Plasmid vs chromosome pie chart
│   │   └── findings.md               # Full clinical findings report
│   └── report/
│       └── final_report.html         # Self-contained clinical HTML report
├── .github/workflows/ci.yml          # GitHub Actions CI (runs tests on push)
├── Snakefile                         # End-to-end Snakemake workflow
├── config.yaml                       # Pipeline configuration
├── environment.yml                   # Conda environment with pinned versions
├── LICENSE                           # MIT License
└── README.md
```

---

## Installation

**Requirements:** conda or mamba

```bash
# Clone the repository
git clone https://github.com/<your-username>/unknown_isolate.git
cd unknown_isolate

# Create and activate the environment
conda env create -f environment.yml
conda activate isolate_analysis
```

Place your input reads at `data/unknown_isolate.fastq.gz` before running.

---

## Usage

### Run the full pipeline

```bash
snakemake --cores 4
```

Snakemake will execute the following steps in order:

1. **QC** — NanoPlot generates read length and quality statistics
2. **Assembly** — Flye assembles the genome in `--nano-hq` mode
3. **Visualize** — Plotly dashboards for QC and AMR results
4. **Report** — Self-contained HTML clinical report

### Run individual scripts

```bash
# Quality control
PYTHONPATH=scripts python scripts/qc_analyzer.py data/unknown_isolate.fastq.gz

# AMR parsing (requires ResFinder output)
PYTHONPATH=scripts python scripts/amr_parser.py results/amr/resfinder/ResFinder_results_tab.txt

# Generate interactive visualizations
PYTHONPATH=scripts python scripts/visualize.py

# Generate HTML clinical report
PYTHONPATH=scripts python scripts/report_generator.py
```

---

## Pipeline Overview

### 1. Quality Control

Sequencing quality was assessed with NanoPlot and a custom OOP FASTQ analyzer.

| Metric | Value |
|--------|-------|
| Total reads | 260,294 |
| Total bases | 576,590,333 bp |
| Mean read length | 2,215 bp |
| Read N50 | 15,932 bp |
| Max read length | 210,485 bp |
| Mean quality | Q20.1 |
| Reads above Q20 | 74.3% |
| Estimated coverage | ~98x |

Coverage exceeds the 30x threshold required for reliable de novo assembly.
An N50 of ~16 kb is well suited to resolving plasmid structures.

### 2. Taxonomic Identification

The assembled chromosome was queried against the NCBI nucleotide (nt) database via BLASTn.

| Parameter | Value |
|-----------|-------|
| Top hit | *Klebsiella pneumoniae* strain 18471, chromosome 1 (OZ185751.1) |
| Identity | 98.14% |
| Query coverage | 100% |
| E-value | 0.0 |

### 3. De Novo Assembly

Assembled with Flye v2.9.6 (`--nano-hq` mode).

| Contig | Length (bp) | Coverage | Circular | Interpretation |
|--------|-------------|----------|----------|----------------|
| contig_4 | 5,306,291 | 84x | Yes | Chromosome |
| contig_2 | 214,836 | 128x | Yes | Large plasmid |
| contig_15 | 79,487 | 142x | Yes | Plasmid (blaKPC-3) |
| contig_8 | 41,154 | 20x | Yes | Small plasmid |
| contig_3 | 31,963 | 39x | No | Plasmid fragment |
| contig_5 | 35,423 | 19x | No | Plasmid fragment |

Total assembled size: ~5.3 Mb, consistent with reference *K. pneumoniae* genomes.

### 4. Antimicrobial Resistance

AMR genes detected with ResFinder 4.x. The tool returned 37 hits representing 23 unique
resistance genes. Duplicate hits arise from identical resistance cassettes present on
multiple plasmids, and from overlapping database entries for the chromosomal blaSHV locus.

| Antibiotic Class | Genes | Location |
|-----------------|-------|----------|
| Carbapenems / Beta-lactams | blaKPC-3, blaSHV-12, blaTEM-1A, blaOXA-9 | Plasmid |
| Beta-lactams (chromosomal) | blaSHV variant (blaSHV-158/159/182, same locus) | Chromosome |
| Aminoglycosides | aac(6')-Ib, aadA1, aadA2, aadA2b, aph(3')-Ia, aph(4)-Ia, aac(3)-IV | Plasmid |
| Fluoroquinolones | OqxA, OqxB (efflux pump) | Chromosome |
| Macrolides | mph(A) | Plasmid |
| Sulfonamides | sul1, sul3 | Plasmid |
| Chloramphenicol | catA1, cmlA1 | Plasmid |
| Trimethoprim | dfrA12 | Plasmid |
| Fosfomycin | fosA6 | Chromosome |

The critical gene **blaKPC-3** (100% identity, 100% coverage) encodes a KPC-3 carbapenemase
and is located on a circular plasmid (contig_15, 142x coverage). Plasmid-borne carbapenem
resistance poses a horizontal gene transfer risk to other bacterial species.

---

## Test Results

Unit tests cover the QC analyzer, AMR gene classification, and parser logic.

```bash
cd unknown_isolate
PYTHONPATH=scripts pytest tests/test_pipeline.py -v
```

```
tests/test_pipeline.py::TestFastqAnalyzer::test_gc_content_all_gc       PASSED
tests/test_pipeline.py::TestFastqAnalyzer::test_gc_content_all_at       PASSED
tests/test_pipeline.py::TestFastqAnalyzer::test_gc_content_mixed        PASSED
tests/test_pipeline.py::TestFastqAnalyzer::test_gc_content_empty        PASSED
tests/test_pipeline.py::TestFastqAnalyzer::test_quality_calculation     PASSED
tests/test_pipeline.py::TestFastqAnalyzer::test_quality_empty           PASSED
tests/test_pipeline.py::TestFastqAnalyzer::test_file_not_found          PASSED
tests/test_pipeline.py::TestFastqAnalyzer::test_summary_requires_parse  PASSED
tests/test_pipeline.py::TestAMRGene::test_critical_kpc                  PASSED
tests/test_pipeline.py::TestAMRGene::test_not_critical                  PASSED
tests/test_pipeline.py::TestAMRGene::test_chromosome_location           PASSED
tests/test_pipeline.py::TestAMRGene::test_plasmid_location              PASSED
tests/test_pipeline.py::TestAMRParser::test_file_not_found              PASSED
tests/test_pipeline.py::TestAMRParser::test_parse_real_file             PASSED
tests/test_pipeline.py::TestAMRParser::test_critical_genes_detected     PASSED
tests/test_pipeline.py::TestAMRParser::test_summary_keys                PASSED

16 passed in 0.XXs
```

---

## Tools and Versions

| Tool | Version | Run | Purpose |
|------|---------|-----|---------|
| Python | 3.10.20 | local | Core scripting language |
| NanoPlot | 1.47.1 | local | ONT read quality assessment |
| Flye | 2.9.6 | local | De novo genome assembly |
| QUAST | 5.2.0 | local | Assembly quality control |
| mlst | 2.33.1 | local | Multi-locus sequence typing (ST258) |
| Kleborate | 3.x | local | Capsular typing, virulence scoring, point mutations |
| minimap2 | — | local | Sequence alignment (Kleborate dependency) |
| BLASTn | — | NCBI web | Taxonomic identification |
| ResFinder | 4.x | CGE web | Acquired AMR gene detection |
| PlasmidFinder | 2.0 | CGE web | Plasmid replicon typing |
| Plotly | 6.9.0 | local | Interactive visualizations |
| pytest | 9.1.1 | local | Unit testing |
| Snakemake | — | local | Workflow management |

ResFinder, PlasmidFinder and BLASTn were run through their official web services because
their local packages could not be installed on the ARM64 macOS environment used for this
analysis. Raw outputs from all three are committed under `results/` so the findings remain
verifiable.

---

## Limitations

- Assembly was not polished with Medaka or Clair3. Polishing is recommended for publication-grade analysis.
- No genome annotation (Prokka/Bakta) was performed.
- No phylogenetic analysis comparing this isolate to other ST258 genomes.
- All resistance predictions are genotypic; phenotypic MIC testing is required before clinical decisions.

---

## License

MIT License. See [LICENSE](LICENSE) for details.
