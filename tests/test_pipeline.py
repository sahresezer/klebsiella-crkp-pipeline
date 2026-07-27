#!/usr/bin/env python3
"""
Unit tests for the isolate analysis pipeline.
Author: Sahre Hilal Sezer
"""

import sys
import pytest
sys.path.insert(0, "code")

from qc_analyzer import FastqAnalyzer, ReadMetrics
from amr_parser import AMRParser, AMRGene


# ── QC ANALYZER TESTS ─────────────────────────────────────

class TestFastqAnalyzer:

    def test_gc_content_all_gc(self):
        analyzer = FastqAnalyzer.__new__(FastqAnalyzer)
        assert analyzer._calc_gc("GGCC") == 100.0

    def test_gc_content_all_at(self):
        analyzer = FastqAnalyzer.__new__(FastqAnalyzer)
        assert analyzer._calc_gc("AATT") == 0.0

    def test_gc_content_mixed(self):
        analyzer = FastqAnalyzer.__new__(FastqAnalyzer)
        assert analyzer._calc_gc("ATGC") == 50.0

    def test_gc_content_empty(self):
        analyzer = FastqAnalyzer.__new__(FastqAnalyzer)
        assert analyzer._calc_gc("") == 0.0

    def test_quality_calculation(self):
        analyzer = FastqAnalyzer.__new__(FastqAnalyzer)
        # ASCII 'I' = 73, 73-33 = 40
        assert analyzer._calc_quality("IIII") == 40.0

    def test_quality_empty(self):
        analyzer = FastqAnalyzer.__new__(FastqAnalyzer)
        assert analyzer._calc_quality("") == 0.0

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            FastqAnalyzer("nonexistent.fastq.gz")

    def test_summary_requires_parse(self):
        analyzer = FastqAnalyzer.__new__(FastqAnalyzer)
        analyzer.reads = []
        with pytest.raises(ValueError):
            analyzer.summary()


# ── AMR PARSER TESTS ──────────────────────────────────────

class TestAMRGene:

    def test_critical_kpc(self):
        gene = AMRGene("blaKPC-3", 100.0, 100.0, "contig_15",
                       "Imipenem, Meropenem", "HM769262")
        assert gene.is_critical is True

    def test_not_critical(self):
        gene = AMRGene("aadA1", 100.0, 100.0, "contig_15",
                       "Streptomycin", "JQ480156")
        assert gene.is_critical is False

    def test_chromosome_location(self):
        gene = AMRGene("OqxA", 100.0, 100.0, "contig_4",
                       "Ciprofloxacin", "EU370913", is_chromosome=True)
        assert gene.location_type == "Chromosome"

    def test_plasmid_location(self):
        gene = AMRGene("blaKPC-3", 100.0, 100.0, "contig_15",
                       "Imipenem", "HM769262", is_chromosome=False)
        assert gene.location_type == "Plasmid"


class TestAMRParser:

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            AMRParser("nonexistent.txt")

    def test_parse_real_file(self):
        parser = AMRParser(
            "results/amr/resfinder/ResFinder_results_tab.txt"
        )
        parser.parse()
        assert len(parser.genes) > 0

    def test_critical_genes_detected(self):
        parser = AMRParser(
            "results/amr/resfinder/ResFinder_results_tab.txt"
        )
        parser.parse()
        critical = [g.gene for g in parser.critical_genes()]
        assert "blaKPC-3" in critical

    def test_summary_keys(self):
        parser = AMRParser(
            "results/amr/resfinder/ResFinder_results_tab.txt"
        )
        parser.parse()
        summary = parser.summary()
        assert "total_genes" in summary
        assert "critical_genes" in summary
        assert summary["total_genes"] > 0
