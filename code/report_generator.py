#!/usr/bin/env python3
"""
Clinical Report Generator
Produces a self-contained HTML report for the isolate analysis.
Author: Sahre Hilal Sezer
"""

import sys
from pathlib import Path
from datetime import date

sys.path.insert(0, "code")
from qc_analyzer import FastqAnalyzer
from amr_parser import AMRParser


def parse_nanostats(path: str) -> dict:
    """
    Parses NanoPlot's NanoStats.txt into a dictionary.

    NanoStats covers the complete read set, whereas FastqAnalyzer is run on a
    subset for speed. Whole-dataset figures are used in the report so that the
    numbers match the NanoPlot output and findings.md.
    """
    stats = {}
    with open(path) as f:
        for line in f:
            if ":" in line:
                key, _, value = line.partition(":")
                stats[key.strip()] = value.strip()
    return stats


def generate_report(output_path: str):
    """Generates a self-contained HTML clinical report."""

    # Whole-dataset QC figures from NanoPlot
    ns = parse_nanostats("results/qc/NanoStats.txt")

    def ns_num(key, default=0.0):
        try:
            return float(ns.get(key, default).replace(",", ""))
        except (AttributeError, ValueError):
            return default

    total_reads = int(ns_num("Number of reads"))
    total_bases = int(ns_num("Total bases"))
    read_n50 = int(ns_num("Read length N50"))
    mean_qual = ns_num("Mean read quality")

    # Percentage of reads above Q20, reported by NanoPlot as ">Q20"
    q20_line = ns.get(">Q20", "")
    q20_pct = q20_line.split("(")[1].split(")")[0] if "(" in q20_line else "n/a"

    # Coverage estimate: total sequenced bases / assembled genome size
    genome_size = 5_913_563
    coverage = round(total_bases / genome_size) if total_bases else 0

    # Independent QC pass over a read subset (validates the custom analyzer)
    analyzer = FastqAnalyzer("data/unknown_isolate.fastq.gz")
    analyzer.parse(max_reads=50000)
    qc = analyzer.summary()

    # Collect AMR data - chromosome/plasmid assignment read from assembly info
    parser = AMRParser("results/amr/resfinder/ResFinder_results_tab.txt")
    parser.load_assembly_info("results/assembly/assembly_info.txt")
    parser.parse()
    amr = parser.summary()

    unique_genes = len({g.gene for g in parser.genes})

    # Build resistance table rows
    gene_rows = ""
    for gene in sorted(parser.genes, key=lambda g: (not g.is_critical, g.gene)):
        badge = (
            '<span style="background:#c00;color:#fff;padding:2px 8px;'
            'border-radius:4px;font-size:11px;margin-left:6px">CRITICAL</span>'
            if gene.is_critical else ""
        )
        phenotype = gene.phenotype[:70] + ("..." if len(gene.phenotype) > 70 else "")
        gene_rows += f"""
        <tr>
            <td><b>{gene.gene}</b>{badge}</td>
            <td>{gene.identity:.2f}%</td>
            <td>{gene.coverage:.1f}%</td>
            <td>{gene.contig}</td>
            <td>{gene.location_type}</td>
            <td style="font-size:12px">{phenotype}</td>
        </tr>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Clinical Isolate Analysis Report</title>
<style>
  body {{ font-family: Arial, Helvetica, sans-serif; max-width: 1050px;
          margin: 40px auto; color: #222; line-height: 1.6; padding: 0 20px; }}
  h1 {{ color: #1F4E79; border-bottom: 3px solid #2E75B6; padding-bottom: 8px; }}
  h2 {{ color: #2E75B6; margin-top: 40px; }}
  h3 {{ color: #444; margin-top: 24px; font-size: 16px; }}
  .alert {{ background: #fff3cd; border-left: 5px solid #c00;
            padding: 16px; margin: 20px 0; border-radius: 4px; }}
  .critical {{ color: #c00; font-weight: bold; }}
  .card {{ background: #f8f9fa; border-radius: 8px;
           padding: 20px; margin: 16px 0; border: 1px solid #e4e7eb; }}
  table {{ width: 100%; border-collapse: collapse; margin-top: 16px;
           font-size: 14px; }}
  th {{ background: #1F4E79; color: white; padding: 10px; text-align: left; }}
  td {{ padding: 8px 10px; border-bottom: 1px solid #ddd; }}
  tr:hover {{ background: #f0f4ff; }}
  .stat {{ display: inline-block; background: #2E75B6; color: white;
           border-radius: 8px; padding: 12px 24px; margin: 8px 8px 8px 0;
           text-align: center; min-width: 90px; }}
  .stat-num {{ font-size: 26px; font-weight: bold; display: block; }}
  .stat-lbl {{ font-size: 12px; opacity: 0.85; }}
  footer {{ margin-top: 60px; color: #888; font-size: 12px;
            border-top: 1px solid #ddd; padding-top: 16px; }}
</style>
</head>
<body>

<h1>Clinical Isolate Analysis Report</h1>
<p><b>Date:</b> {date.today()} &nbsp;|&nbsp;
   <b>Analyst:</b> Sahre Hilal Sezer &nbsp;|&nbsp;
   <b>Sample:</b> Unknown clinical isolate (Oxford Nanopore sequencing)</p>

<div class="alert">
  <b>CRITICAL FINDING:</b> This isolate is
  <span class="critical">Klebsiella pneumoniae ST258</span> carrying
  <span class="critical">blaKPC-3</span> on a conjugative plasmid, combined with
  <span class="critical">dual OmpK35/OmpK36 porin loss</span>.
  This confers extremely high-level carbapenem resistance that may not be
  overcome by ceftazidime-avibactam. ST258 is the globally dominant
  carbapenem-resistant K. pneumoniae outbreak clone.
  Urgent clinical and infection control action is required.
</div>

<h2>1. Organism Identification and Typing</h2>
<div class="card">
  <b>Species:</b> <i>Klebsiella pneumoniae</i><br>
  <b>Method:</b> BLASTn vs NCBI nt database<br>
  <b>Identity:</b> 98.14% &nbsp;|&nbsp; <b>E-value:</b> 0.0 &nbsp;|&nbsp;
  <b>Coverage:</b> 100%<br>
  <b>Top hit:</b> K. pneumoniae strain 18471, chromosome 1 (OZ185751.1)
</div>

<h3>Multi-Locus Sequence Typing (MLST)</h3>
<div class="card">
  <b>Sequence Type: <span class="critical">ST258</span></b>
  &nbsp;(confirmed independently by mlst v2.33.1 and Kleborate)<br>
  <b>Allelic profile:</b> gapA(3) infB(3) mdh(1) pgi(1) phoE(1) rpoB(1) tonB(79)<br><br>
  ST258 is the globally dominant carbapenem-resistant <i>K. pneumoniae</i> (CRKP)
  clone, responsible for large-scale hospital outbreaks worldwide and strongly
  associated with KPC-type carbapenemases.
</div>

<h3>Capsular and Antigen Typing (Kleborate / Kaptive)</h3>
<div class="card">
  <b>K locus:</b> KL107 (99.97% identity) &nbsp;|&nbsp;
  <b>O locus:</b> OL13 / O13 (98.26% identity) &nbsp;|&nbsp;
  <b>WZI type:</b> wzi154<br>
  <b>Virulence score:</b> 0 / 5 (non-hypervirulent) &nbsp;|&nbsp;
  <b>Resistance score:</b> 2 / 3 (extensively drug-resistant, XDR)
</div>

<h2>2. Sequencing Quality (QC)</h2>
<p style="font-size:13px;color:#666">
  Figures below are whole-dataset values reported by NanoPlot v1.47.1.
</p>
<div>
  <div class="stat">
    <span class="stat-num">{total_reads:,}</span>
    <span class="stat-lbl">Total Reads</span>
  </div>
  <div class="stat">
    <span class="stat-num">{total_bases/1e6:,.0f}</span>
    <span class="stat-lbl">Total Bases (Mb)</span>
  </div>
  <div class="stat">
    <span class="stat-num">{read_n50:,}</span>
    <span class="stat-lbl">Read N50 (bp)</span>
  </div>
  <div class="stat">
    <span class="stat-num">Q{mean_qual:.1f}</span>
    <span class="stat-lbl">Mean Quality</span>
  </div>
  <div class="stat">
    <span class="stat-num">{q20_pct}</span>
    <span class="stat-lbl">Reads &ge; Q20</span>
  </div>
  <div class="stat">
    <span class="stat-num">~{coverage}x</span>
    <span class="stat-lbl">Coverage</span>
  </div>
</div>
<div class="card" style="font-size:13px">
  <b>Independent validation:</b> a custom Python QC analyzer
  (<code>code/qc_analyzer.py</code>) was run over a {qc['total_reads']:,}-read
  subset as a cross-check of the NanoPlot output, giving a read N50 of
  {qc['n50']:,} bp, mean quality Q{qc['mean_quality']}, and mean GC content
  {qc['mean_gc']}%. Mean quality is computed by averaging per-base error
  probabilities rather than raw Phred scores, matching NanoPlot's method.<br><br>
  <b>Assessment:</b> Coverage of ~{coverage}x is well above the 30x threshold
  required for reliable bacterial assembly, and the read N50 of
  {read_n50/1000:.1f} kb is sufficient to resolve plasmid structures.
</div>

<h2>3. Genome Assembly</h2>
<div class="card">
  <b>Assembler:</b> Flye v2.9.6 (--nano-hq) &nbsp;|&nbsp;
  <b>Assembly QC:</b> QUAST v5.2.0<br>
  <b>Total contigs:</b> 13 &nbsp;|&nbsp;
  <b>Total length:</b> 5,913,563 bp &nbsp;|&nbsp;
  <b>GC content:</b> 56.79%<br>
  <b>Chromosome:</b> 5,306,291 bp (contig_4, circular) &nbsp;|&nbsp;
  <b>N50:</b> 5,306,291 bp &nbsp;|&nbsp; <b>L50:</b> 1 &nbsp;|&nbsp;
  <b>Gaps:</b> 0 N per 100 kbp<br><br>
  An N50 equal to the chromosome length with L50 of 1 and zero gaps indicates a
  near-complete, gapless chromosome assembly.<br><br>
  <b>Key plasmid:</b> contig_15 (79,487 bp, circular, 142x) &mdash;
  carries <span class="critical">blaKPC-3</span>
</div>

<h2>4. Antimicrobial Resistance</h2>
<div class="card">
  <b>ResFinder hits:</b> {amr['total_genes']} &nbsp;|&nbsp;
  <b>Unique genes:</b> {unique_genes} &nbsp;|&nbsp;
  <b>On plasmids:</b> {amr['plasmid_genes']} &nbsp;|&nbsp;
  <b>On chromosome:</b> {amr['chromosome_genes']} &nbsp;|&nbsp;
  <b>Resistance classes:</b> {amr['resistance_classes']}<br><br>
  <span class="critical">Critical genes: {', '.join(amr['critical_genes'])}</span>
</div>

<h3>Outer Membrane Protein Mutations (Critical)</h3>
<div class="card">
  <table>
    <tr><th>Gene</th><th>Mutation</th><th>Effect</th></tr>
    <tr><td>OmpK35</td><td>p.Glu42fs (frameshift)</td>
        <td>Complete loss of function &mdash; porin absent</td></tr>
    <tr><td>OmpK36</td><td>p.Val59fs (frameshift)</td>
        <td>Complete loss of function &mdash; porin absent</td></tr>
  </table>
  <p>Carbapenems enter <i>K. pneumoniae</i> through the OmpK35 and OmpK36 outer
  membrane channels. Simultaneous loss of both porins, combined with KPC-3
  carbapenemase, produces extremely high-level carbapenem resistance that may
  not be overcome by ceftazidime-avibactam.</p>
</div>

<h3>Fluoroquinolone Resistance Mutations</h3>
<div class="card">
  <table>
    <tr><th>Gene</th><th>Mutation</th><th>Effect</th></tr>
    <tr><td>GyrA</td><td>p.Ser83Ile</td>
        <td>Primary QRDR resistance mutation</td></tr>
    <tr><td>GyrA</td><td>p.Asp87Asn</td>
        <td>Secondary QRDR resistance mutation</td></tr>
    <tr><td>ParC</td><td>p.Ser80Ile</td>
        <td>Additional target-site resistance</td></tr>
  </table>
  <p>Combined with the OqxA/OqxB efflux pump, this strain carries both
  chromosomal target mutations and efflux-mediated resistance, rendering the
  entire fluoroquinolone class ineffective.</p>
</div>

<h3>Plasmid Replicon Typing (PlasmidFinder 2.0)</h3>
<div class="card">
  <table>
    <tr><th>Replicon</th><th>Identity</th><th>Contig</th>
        <th>Clinical Significance</th></tr>
    <tr><td>IncFIB(K)</td><td>100%</td><td>contig_2</td>
        <td>Common K. pneumoniae plasmid; ST258-associated</td></tr>
    <tr><td>IncFII(K)</td><td>100%</td><td>contig_2</td>
        <td>K. pneumoniae-specific IncF replicon; ST258 hallmark</td></tr>
    <tr><td>IncFII(Yp)</td><td>95.63%</td><td>contig_1</td>
        <td>IncF family; broad host range</td></tr>
    <tr><td><b>IncI2(Delta)</b></td><td>98.42%</td><td>contig_15</td>
        <td><span class="critical">Carries blaKPC-3</span>; conjugative,
        broad host range</td></tr>
    <tr><td>IncR</td><td>100%</td><td>contig_5, 6, 8</td>
        <td>Associated with resistance gene accumulation</td></tr>
    <tr><td>repB(R1701)</td><td>100%</td><td>contig_7</td>
        <td>Small mobilizable replicon</td></tr>
  </table>
  <p>The blaKPC-3 gene resides on an IncI2(Delta) conjugative plasmid with broad
  host range across Enterobacteriaceae, meaning resistance can transfer
  horizontally to <i>E. coli</i>, <i>Salmonella</i>, and other species.</p>
</div>

<h3>Full Resistance Gene Table</h3>
<table>
  <tr>
    <th>Gene</th><th>Identity</th><th>Coverage</th>
    <th>Contig</th><th>Location</th><th>Phenotype</th>
  </tr>
  {gene_rows}
</table>

<h2>5. Recommendations</h2>
<div class="card">
  <ol>
    <li><b>Urgent:</b> Notify the infectious disease team and infection control
        immediately. This is a carbapenem-resistant <i>K. pneumoniae</i> ST258,
        a notifiable high-risk pathogen in most jurisdictions.</li>
    <li><b>Treatment:</b> Do not rely on carbapenems, fluoroquinolones, or
        standard beta-lactams. Given dual OmpK35/OmpK36 porin loss combined with
        KPC-3, ceftazidime-avibactam efficacy may be compromised. Consider
        <b>cefiderocol</b> or <b>aztreonam-avibactam</b> pending MIC results.</li>
    <li><b>Susceptibility testing:</b> Perform phenotypic MIC determination for
        all relevant antibiotic classes before finalizing treatment.</li>
    <li><b>Infection control:</b> Implement contact precautions immediately.
        The blaKPC-3 gene is on a conjugative IncI2(Delta) plasmid &mdash;
        horizontal transfer to other ward pathogens is a real risk.</li>
    <li><b>Public health reporting:</b> ST258 CRKP is reportable in many
        countries. Notify the relevant authority if required locally.</li>
    <li>Screen close contacts and the ward environment for transmission.</li>
  </ol>
</div>

<h2>Note to Professor Kılıç</h2>
<div class="card">
  Dear Professor Kılıç,<br><br>

  The patient's bacterial isolate has been identified as
  <i>Klebsiella pneumoniae</i> <b>ST258</b>, one of the most clinically dangerous
  hospital-acquired pathogens worldwide.<br><br>

  <b>The critical findings are:</b>
  <ol>
    <li>The strain carries <b>blaKPC-3</b>, a gene that destroys carbapenems &mdash;
        the antibiotics used as a last resort when everything else fails.</li>
    <li>Beyond the enzyme itself, the bacterium has lost two proteins
        (<b>OmpK35</b> and <b>OmpK36</b>) that normally allow antibiotics to enter
        the cell. This double barrier means it resists even some newer drugs
        designed to overcome KPC-type resistance.</li>
    <li>The resistance genes sit on a <b>plasmid</b> &mdash; a small, mobile piece
        of DNA that can spread to other bacterial species in the same patient
        or ward.</li>
    <li><b>ST258</b> is a globally recognised outbreak clone responsible for large
        hospital epidemics.</li>
  </ol>

  <b>What this means for the patient:</b> Standard antibiotics and most
  last-resort options will not work. The clinical team should urgently consult an
  infectious disease specialist. Given the dual porin loss combined with KPC-3,
  ceftazidime-avibactam alone may have reduced efficacy. Cefiderocol or
  aztreonam-avibactam may be considered, but phenotypic susceptibility testing
  is essential before any treatment decision.<br><br>

  <b>Infection control is equally urgent:</b> This strain can pass its resistance
  to other bacteria. Strict isolation measures should be implemented immediately.
</div>

<h2>Limitations</h2>
<div class="card">
  <ul>
    <li>Assembly was not polished with Medaka or Clair3; residual indel errors
        may remain.</li>
    <li>No functional genome annotation (Prokka/Bakta) was performed.</li>
    <li>No phylogenetic comparison against other ST258 genomes.</li>
    <li>All resistance predictions are genotypic. Phenotypic MIC testing is
        required to confirm clinical breakpoints.</li>
  </ul>
</div>

<footer>
  Pipeline: NanoPlot &middot; Flye &middot; QUAST &middot; BLASTn &middot;
  mlst &middot; Kleborate &middot; ResFinder &middot; PlasmidFinder &middot;
  Snakemake<br>
  Databases accessed: {date.today()} &nbsp;|&nbsp;
  Analyst: Sahre Hilal Sezer
</footer>

</body>
</html>"""

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)
    print(f"Report saved: {output_path}")


if __name__ == "__main__":
    generate_report("results/report/final_report.html")
