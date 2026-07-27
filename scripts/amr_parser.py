#!/usr/bin/env python3
"""
AMR Result Parser for ResFinder Output
Parses and interprets antimicrobial resistance gene findings
Author: Sahre Hilal Sezer
"""

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass
class AMRGene:
    """Represents a single AMR gene finding."""
    gene: str
    identity: float
    coverage: float
    contig: str
    phenotype: str
    accession: str
    is_chromosome: bool = False  # set by AMRParser.load_assembly_info()

    @property
    def is_critical(self) -> bool:
        """Flags carbapenem resistance genes as critical."""
        critical = ["blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP"]
        return any(c in self.gene for c in critical)

    @property
    def location_type(self) -> str:
        """Returns Chromosome or Plasmid based on assembly info."""
        return "Chromosome" if self.is_chromosome else "Plasmid"


class AMRParser:
    """
    Parses ResFinder tab-delimited output and provides
    structured access to resistance gene findings.
    """

    def __init__(self, results_path: str):
        self.path = Path(results_path)
        if not self.path.exists():
            raise FileNotFoundError(f"File not found: {results_path}")
        self.genes: list[AMRGene] = []

    def load_assembly_info(self, assembly_info_path: str,
                           min_chromosome_size: int = 1_000_000):
        """
        Reads Flye assembly_info.txt and marks chromosome contigs.
        Contigs >= min_chromosome_size (default 1 Mb) are treated as chromosome.
        Must be called before or after parse(); updates genes in place.
        """
        path = Path(assembly_info_path)
        if not path.exists():
            raise FileNotFoundError(f"Assembly info not found: {assembly_info_path}")

        chromosome_contigs: set[str] = set()
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    name, length = parts[0], int(parts[1])
                    if length >= min_chromosome_size:
                        chromosome_contigs.add(name)

        self._chromosome_contigs = chromosome_contigs
        for gene in self.genes:
            gene.is_chromosome = gene.contig in chromosome_contigs
        return self

    def parse(self):
        """Parses ResFinder tab file into AMRGene objects."""
        chromosome_contigs = getattr(self, "_chromosome_contigs", set())
        with open(self.path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                try:
                    contig = row["Contig"].strip()
                    self.genes.append(AMRGene(
                        gene=row["Resistance gene"],
                        identity=float(row["Identity"]),
                        coverage=float(row["Coverage"]),
                        contig=contig,
                        phenotype=row["Phenotype"],
                        accession=row["Accession no."],
                        is_chromosome=contig in chromosome_contigs,
                    ))
                except (KeyError, ValueError):
                    continue
        return self

    def critical_genes(self) -> list[AMRGene]:
        """Returns only critical (last-line) resistance genes."""
        return [g for g in self.genes if g.is_critical]

    def by_location(self) -> dict:
        """Groups genes by chromosome vs plasmid."""
        result = {"Chromosome": [], "Plasmid": []}
        for gene in self.genes:
            result[gene.location_type].append(gene)
        return result

    def resistance_classes(self) -> dict:
        """Groups genes by antibiotic class."""
        classes = {}
        for gene in self.genes:
            phenotypes = [p.strip() for p in gene.phenotype.split(",")]
            for p in phenotypes:
                if p not in classes:
                    classes[p] = []
                if gene.gene not in classes[p]:
                    classes[p].append(gene.gene)
        return classes

    def summary(self) -> dict:
        """Returns summary statistics."""
        return {
            "total_genes": len(self.genes),
            "critical_genes": [g.gene for g in self.critical_genes()],
            "plasmid_genes": len(self.by_location()["Plasmid"]),
            "chromosome_genes": len(self.by_location()["Chromosome"]),
            "resistance_classes": len(self.resistance_classes()),
        }


if __name__ == "__main__":
    import sys
    import json

    parser = AMRParser(sys.argv[1])
    parser.parse()

    print("=== CRITICAL GENES ===")
    for g in parser.critical_genes():
        print(f"  {g.gene} | {g.location_type} ({g.contig}) | Identity: {g.identity}%")

    print("\n=== SUMMARY ===")
    print(json.dumps(parser.summary(), indent=2))

    print("\n=== RESISTANCE CLASSES ===")
    for cls, genes in parser.resistance_classes().items():
        print(f"  {cls}: {', '.join(genes)}")
