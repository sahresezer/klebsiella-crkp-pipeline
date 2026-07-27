#!/usr/bin/env python3
"""
Visualization module for isolate analysis pipeline.
Generates interactive Plotly charts for QC and AMR results.
Author: Sahre Hilal Sezer
"""

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from pathlib import Path
from amr_parser import AMRParser
from qc_analyzer import FastqAnalyzer


def plot_qc_dashboard(stats: dict, reads: list, output_dir: str):
    """Creates a multi-panel QC dashboard."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "Read Length Distribution",
            "Quality Score Distribution",
            "GC Content Distribution",
            "Read Length vs Quality"
        )
    )

    lengths = [r.length for r in reads]
    qualities = [r.mean_quality for r in reads]
    gcs = [r.gc_content for r in reads]

    # Read length histogram
    fig.add_trace(go.Histogram(
        x=lengths, nbinsx=50, name="Read Length",
        marker_color="#2E75B6", opacity=0.8
    ), row=1, col=1)

    # Quality score histogram
    fig.add_trace(go.Histogram(
        x=qualities, nbinsx=40, name="Quality Score",
        marker_color="#375623", opacity=0.8
    ), row=1, col=2)

    # GC content histogram
    fig.add_trace(go.Histogram(
        x=gcs, nbinsx=40, name="GC Content",
        marker_color="#C55A11", opacity=0.8
    ), row=2, col=1)

    # Length vs Quality scatter
    fig.add_trace(go.Scatter(
        x=lengths[:5000], y=qualities[:5000],
        mode="markers", name="Length vs Quality",
        marker=dict(color="#7030A0", size=3, opacity=0.5)
    ), row=2, col=2)

    fig.update_layout(
        title=f"<b>QC Dashboard — Unknown Isolate</b><br>"
              f"<sup>Plotted from a {stats['total_reads']:,}-read subset | "
              f"subset N50: {stats['n50']:,} bp | "
              f"subset mean Q: {stats['mean_quality']} | "
              f"whole-dataset figures in NanoPlot-report.html</sup>",
        height=700,
        showlegend=False,
        template="plotly_white"
    )

    # plotly.js is loaded from CDN rather than embedded: embedding adds ~4.5 MB
    # per file, which would triple the size of this repository for no benefit.
    fig.write_html(str(output_dir / "qc_dashboard.html"), include_plotlyjs="cdn")
    print(f"QC dashboard saved: {output_dir}/qc_dashboard.html")


def plot_amr_summary(amr_path: str, output_dir: str,
                     assembly_info: str = "results/assembly/assembly_info.txt"):
    """Creates AMR resistance summary chart."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    parser = AMRParser(amr_path)
    parser.load_assembly_info(assembly_info)
    parser.parse()

    # Resistance classes bar chart
    classes = parser.resistance_classes()
    df = pd.DataFrame([
        {"Antibiotic": cls, "Gene Count": len(genes), "Genes": ", ".join(genes)}
        for cls, genes in classes.items()
    ]).sort_values("Gene Count", ascending=True)

    fig = px.bar(
        df, x="Gene Count", y="Antibiotic",
        orientation="h",
        color="Gene Count",
        color_continuous_scale="Reds",
        hover_data=["Genes"],
        title="<b>Antimicrobial Resistance Profile</b><br>"
              "<sup>Klebsiella pneumoniae — Unknown Isolate</sup>"
    )

    # Highlight critical finding
    fig.add_annotation(
        text="CRITICAL: blaKPC-3 detected — carbapenem resistance",
        xref="paper", yref="paper",
        x=0.5, y=1.08, showarrow=False,
        font=dict(size=13, color="red"),
        bgcolor="lightyellow",
        bordercolor="red", borderwidth=1
    )

    fig.update_layout(
        height=700, template="plotly_white",
        coloraxis_showscale=False
    )

    fig.write_html(str(output_dir / "amr_summary.html"), include_plotlyjs="cdn")
    print(f"AMR chart saved: {output_dir}/amr_summary.html")

    # Plasmid vs Chromosome pie chart
    location = parser.by_location()
    fig2 = go.Figure(go.Pie(
        labels=["Plasmid", "Chromosome"],
        values=[len(location["Plasmid"]), len(location["Chromosome"])],
        hole=0.4,
        marker_colors=["#C55A11", "#2E75B6"]
    ))
    fig2.update_layout(
        title="<b>Resistance Gene Location</b><br>"
              "<sup>Plasmid-borne genes can spread to other bacteria</sup>",
        template="plotly_white"
    )
    fig2.write_html(str(output_dir / "amr_location.html"), include_plotlyjs="cdn")
    print(f"AMR location chart saved: {output_dir}/amr_location.html")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, "code")

    print("Generating QC visualizations...")
    analyzer = FastqAnalyzer("data/unknown_isolate.fastq.gz")
    analyzer.parse(max_reads=10000)
    stats = analyzer.summary()
    plot_qc_dashboard(stats, analyzer.reads, "results/qc")

    print("Generating AMR visualizations...")
    plot_amr_summary(
        "results/amr/resfinder/ResFinder_results_tab.txt",
        "results/amr"
    )
    print("Done!")
