# Findings Report: Unknown Bacterial Isolate

## Note to Professor Kılıç

Dear Professor Kılıç,

The bacterial isolate from your patient has been identified as **Klebsiella pneumoniae ST258**,
a gram-negative bacterium and one of the most clinically dangerous hospital-acquired pathogens worldwide.

**The critical findings are:**

1. This strain carries *blaKPC-3*, a gene that destroys carbapenems — the antibiotics used as a
   last resort when everything else fails.
2. Beyond the enzyme itself, the bacterium has additionally lost two proteins (OmpK35 and OmpK36)
   that normally allow antibiotics to enter the cell. This double barrier means the bacteria is
   resistant even to some newer drugs designed to overcome KPC-type resistance.
3. The resistance genes are on a plasmid — a small, mobile piece of DNA that can spread from
   this bacterium to other bacterial species in the same patient or ward.
4. ST258 is a globally recognized outbreak clone responsible for large hospital epidemics.

**What this means for the patient:** Standard antibiotics and most last-resort options will not
work. The clinical team should urgently consult an infectious disease specialist. Given the dual
porin loss combined with KPC-3, ceftazidime-avibactam alone may have reduced efficacy.
**Cefiderocol** or **aztreonam-avibactam** may be considered, but phenotypic susceptibility
testing (MIC determination) is essential before any treatment decision.

**Infection control is equally urgent:** This strain can spread its resistance to other bacteria.
Strict isolation measures should be implemented immediately.

---

## 1. Data Quality (QC)

| Metric | Value |
|--------|-------|
| Total reads | 260,294 |
| Total bases | 576,590,333 bp |
| Mean read length | 2,215 bp |
| Read length N50 | 15,932 bp |
| Max read length | 210,485 bp |
| Mean read quality | Q20.1 |
| Reads above Q20 | 74.3% |
| Estimated coverage | ~98x |

Coverage was estimated as total sequenced bases divided by the final assembled
genome size (576,590,333 bp / 5,913,563 bp = 97.5x).

**Assessment:** High-quality ONT data. Coverage is well above the 30x threshold 
required for reliable assembly. N50 of ~16 kb is excellent for complete genome assembly.

---

## 2. Organism Identification

**Species:** *Klebsiella pneumoniae*  
**Method:** BLASTn against NCBI nt database  
**Identity:** 98.14%  
**E-value:** 0.0  
**Query coverage:** 100%  
**Top hit:** Klebsiella pneumoniae strain 18471, chromosome 1 (OZ185751.1)

### MLST (Multi-Locus Sequence Typing)

**Sequence Type: ST258**

| Locus | Allele |
|-------|--------|
| gapA | 3 |
| infB | 3 |
| mdh | 1 |
| pgi | 1 |
| phoE | 1 |
| rpoB | 1 |
| tonB | 79 |

ST258 is the globally dominant carbapenem-resistant *Klebsiella pneumoniae* (CRKP) clone.
It is responsible for large-scale hospital outbreaks worldwide and is strongly associated
with KPC-type carbapenemases — consistent with the blaKPC-3 gene detected in this isolate.
This sequence type is a recognized high-risk epidemic clone and should be reported to
infection control and, where applicable, public health authorities.

---

## 3. Genome Assembly

**Assembler:** Flye v2.9.6 (--nano-hq mode)  
**Assembly QC:** QUAST v5.2.0

| Contig | Length (bp) | Coverage | Circular | Interpretation |
|--------|-------------|----------|----------|----------------|
| contig_4 | 5,306,291 | 84x | Yes | Chromosome |
| contig_2 | 214,836 | 128x | Yes | Large plasmid |
| contig_15 | 79,487 | 142x | Yes | Plasmid (carries blaKPC-3) |
| contig_8 | 41,154 | 20x | Yes | Small plasmid |
| contig_3 | 31,963 | 39x | No | Plasmid fragment |
| contig_5 | 35,423 | 19x | No | Plasmid fragment |

Total assembled size: 5,913,563 bp (~5.9 Mb, consistent with *K. pneumoniae* reference genomes)

### QUAST Assembly Quality Metrics

| Metric | Value | Assessment |
|--------|-------|------------|
| Total contigs | 13 | — |
| Largest contig | 5,306,291 bp | Chromosome assembled as single contig |
| Total length | 5,913,563 bp | — |
| N50 | 5,306,291 bp | Exceptional — equals chromosome length |
| L50 | 1 | Only 1 contig needed to cover 50% of assembly |
| GC content | 56.79% | Consistent with *K. pneumoniae* reference (~57%) |
| N's per 100 kbp | 0.00 | No gaps — gapless assembly |

An N50 equal to the chromosome length and L50 of 1 indicates a near-complete, gapless
chromosome assembly. This is a direct result of the long-read ONT sequencing and
high coverage (~98x), which enabled the assembler to resolve repetitive regions.

---

## 4. Antimicrobial Resistance (AMR)

**Tool:** ResFinder 4.x (CARD/ResFinder database)

ResFinder returned **37 hits** across all contigs, representing **23 unique resistance genes**.
The difference arises from two sources: (1) identical resistance cassettes present on multiple
plasmids (contig_3, contig_5, and contig_8 carry the same aminoglycoside/sulfonamide/chloramphenicol
cassette, likely from related plasmid replicons); (2) three overlapping blaSHV database entries
(blaSHV-158, -159, -182) mapping to the same chromosomal locus.

### Critical Finding
| Gene | Resistance | Location | Identity |
|------|-----------|----------|----------|
| **blaKPC-3** | Carbapenems (imipenem, meropenem, ertapenem) | Plasmid (contig_15) | 100% |

This is a **carbapenem-resistant Klebsiella pneumoniae (CRKP)**. Carbapenems are 
last-line antibiotics. Resistance to these agents severely limits treatment options.

### Full Resistance Profile (23 unique genes)

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

### Kleborate Extended Typing Results

**Tool:** Kleborate v3.x | **Modules:** MLST, AMR, Kaptive, Virulence, Resistance scoring

| Category | Finding |
|----------|---------|
| Sequence Type | ST258 (confirmed) |
| K antigen locus | KL107 (99.97% identity) |
| O antigen locus | OL13 / O13 (98.26% identity) |
| WZI type | wzi154 |
| Virulence score | 0 / 5 — non-hypervirulent (no yersiniabactin, colibactin, aerobactin) |
| Resistance score | 2 / 3 — XDR (extensively drug-resistant) |

#### Outer Membrane Protein (OMP) Mutations — Critical

| Gene | Mutation | Effect |
|------|----------|--------|
| OmpK35 | p.Glu42fs (frameshift) | Complete loss of function — porin absent |
| OmpK36 | p.Val59fs (frameshift) | Complete loss of function — porin absent |

The simultaneous loss of both OmpK35 and OmpK36 porins is clinically critical.
In *K. pneumoniae*, carbapenems enter the cell through these outer membrane channels.
When combined with KPC-3 carbapenemase, the **dual porin loss** results in extremely
high-level carbapenem resistance that may not be overcome even by novel
beta-lactam/beta-lactamase inhibitor combinations such as ceftazidime-avibactam.

#### Fluoroquinolone Resistance Mutations

| Gene | Mutation | Effect |
|------|----------|--------|
| GyrA | p.Ser83Ile | Primary fluoroquinolone resistance mutation (QRDR) |
| GyrA | p.Asp87Asn | Secondary fluoroquinolone resistance mutation (QRDR) |
| ParC | p.Ser80Ile | Additional resistance — reduces inhibitor binding |

Combined with the OqxA/OqxB efflux pump detected by ResFinder, this strain carries
**both chromosomal target mutations and efflux-mediated resistance** to fluoroquinolones,
rendering this entire antibiotic class ineffective.

#### SHV Beta-Lactamase Mutations

| Mutation | Effect |
|----------|--------|
| SHV p.Gly238Ser | Classic ESBL-conferring mutation |
| SHV p.Glu240Lys | ESBL-conferring mutation |

These mutations in the chromosomal SHV enzyme confer extended-spectrum
beta-lactamase (ESBL) activity in addition to the plasmid-borne SHV-12 and KPC-3.

### Plasmid Replicon Typing (PlasmidFinder 2.0)

| Replicon | Identity | Contig | Clinical Significance |
|----------|----------|--------|-----------------------|
| IncFIB(K) | 100% | contig_2 | Common K. pneumoniae plasmid; associated with ST258 |
| IncFII(K) | 100% | contig_2 | K. pneumoniae-specific IncF replicon; ST258 hallmark |
| IncFII(Yp) | 95.63% | contig_1 | IncF family; broad host range |
| IncI2(Delta) | 98.42% | contig_15 | **Carries blaKPC-3**; conjugative, broad host range |
| IncR | 100% | contig_5, contig_6, contig_8 | Associated with resistance gene accumulation |
| repB(R1701) | 100% | contig_7 | Small mobilizable replicon |

### Plasmid Risk Assessment

The **blaKPC-3** carbapenemase gene resides on an **IncI2(Delta)** plasmid (contig_15, 142x coverage,
circular). IncI2 plasmids are conjugative with a broad host range across Enterobacteriaceae,
meaning this resistance element can be transferred horizontally to *E. coli*, *Salmonella*,
and other clinically relevant species.

The co-occurrence of **IncFIB(K) + IncFII(K)** on contig_2 is a well-documented signature
of *K. pneumoniae* ST258, further confirming the epidemic clone identification.
The IncR replicons on multiple small contigs are frequently associated with integron-mediated
resistance gene accumulation, consistent with the dense AMR gene cassette observed in this isolate.

This isolate represents a high-priority infection control risk: carbapenem resistance on a
mobile, conjugative plasmid in a globally recognized epidemic clone (ST258).

---

## 5. Tools and Database Versions

| Tool | Version | Purpose | Database |
|------|---------|---------|----------|
| NanoPlot | 1.47.1 | Read QC | — |
| Flye | 2.9.6 | De novo assembly | — |
| QUAST | 5.2.0 | Assembly quality control | — |
| BLASTn | Web (NCBI) | Taxonomic identification | nt (accessed 2026-07-27) |
| mlst | 2.33.1 | Multi-locus sequence typing | PubMLST klebsiella scheme |
| ResFinder | 4.x (web) | AMR gene detection | ResFinder DB (accessed 2026-07-27) |
| PlasmidFinder | 2.0 (web) | Plasmid replicon typing | Enterobacteriales DB (accessed 2026-07-27) |
| Kleborate | 3.x | Extended K. pneumoniae typing | — |
| minimap2 | — | Sequence alignment (Kleborate dep.) | — |

---

## 6. Limitations

- **Assembly not polished.** Flye output was not polished with Medaka or Clair3. Long-read assemblies may carry residual indel errors; polishing is recommended for publication-grade analysis.
- **No genome annotation.** Functional gene annotation (Prokka/Bakta) was not performed. Annotation would enable identification of genomic islands, integrons, and mobile genetic element boundaries.
- **No phylogenetic analysis.** Comparison with publicly available ST258 genomes (e.g., via Parsnp or Mashtree) would place this isolate in a broader epidemiological context.
- **Phenotypic susceptibility testing not available.** All resistance predictions are genotypic. MIC determination is required to confirm clinical breakpoints before treatment decisions.

---

## 7. Recommendations

1. **Urgent — clinical:** Notify the infectious disease team and infection control immediately. This is a carbapenem-resistant *K. pneumoniae* ST258 — a notifiable high-risk pathogen in most jurisdictions.
2. **Treatment:** Do not rely on carbapenems, fluoroquinolones, or standard beta-lactams. Given dual OmpK35/OmpK36 porin loss combined with KPC-3, ceftazidime-avibactam efficacy may be compromised. Consider **cefiderocol** or **aztreonam-avibactam** pending MIC results.
3. **Susceptibility testing:** Perform phenotypic MIC determination for all relevant antibiotic classes before finalizing treatment.
4. **Infection control:** Implement contact precautions immediately. The blaKPC-3 gene is on a conjugative IncI2(Delta) plasmid — horizontal transfer to other ward pathogens is a real risk.
5. **Public health reporting:** ST258 CRKP is a reportable pathogen in many countries. Notify the relevant public health authority if required by local regulations.
