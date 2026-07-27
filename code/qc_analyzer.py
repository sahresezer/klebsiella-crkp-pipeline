#!/usr/bin/env python3
"""
QC Analyzer for Oxford Nanopore Sequencing Data
OOP-based FASTQ quality control analysis
Author: Sahre Hilal Sezer
"""

import gzip
import math
import statistics
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ReadMetrics:
    """Stores metrics for a single sequencing read."""
    read_id: str
    length: int
    gc_content: float
    mean_quality: float


class FastqAnalyzer:
    """
    OOP-based FASTQ analyzer for Oxford Nanopore data.
    Supports both .fastq and .fastq.gz formats.
    """

    def __init__(self, fastq_path: str):
        self.path = Path(fastq_path)
        if not self.path.exists():
            raise FileNotFoundError(f"File not found: {fastq_path}")
        self.reads: list[ReadMetrics] = []

    def _open_file(self):
        """Opens gzipped or plain FASTQ files."""
        if self.path.suffix == ".gz":
            return gzip.open(self.path, "rt")
        return open(self.path, "r")

    def _calc_gc(self, sequence: str) -> float:
        """Calculates GC content percentage."""
        if not sequence:
            return 0.0
        gc = sequence.upper().count("G") + sequence.upper().count("C")
        return round(gc / len(sequence) * 100, 2)

    def _calc_quality(self, qual_string: str) -> float:
        """
        Converts a Phred+33 ASCII quality string to a mean quality score.

        Phred scores are logarithmic, so averaging them arithmetically
        overestimates read quality. The correct approach is to convert each
        score to its error probability, average those, then convert back.
        This matches how NanoPlot and other ONT QC tools report mean quality.
        """
        if not qual_string:
            return 0.0
        error_probs = [10 ** (-(ord(c) - 33) / 10) for c in qual_string]
        mean_error = statistics.mean(error_probs)
        if mean_error == 0:
            return 0.0
        return round(-10 * math.log10(mean_error), 2)

    def parse(self, max_reads: int = None):
        """Parses FASTQ file and stores read metrics."""
        self.reads = []
        with self._open_file() as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                sequence = f.readline().strip()
                f.readline()  # '+' line
                quality = f.readline().strip()

                read_id = header[1:].split()[0]
                self.reads.append(ReadMetrics(
                    read_id=read_id,
                    length=len(sequence),
                    gc_content=self._calc_gc(sequence),
                    mean_quality=self._calc_quality(quality),
                ))

                if max_reads and len(self.reads) >= max_reads:
                    break
        return self

    def summary(self) -> dict:
        """Returns summary statistics."""
        if not self.reads:
            raise ValueError("No reads parsed. Run parse() first.")

        lengths = [r.length for r in self.reads]
        qualities = [r.mean_quality for r in self.reads]
        gcs = [r.gc_content for r in self.reads]

        lengths_sorted = sorted(lengths)
        total = sum(lengths_sorted)
        cumsum, n50 = 0, 0
        for l in reversed(lengths_sorted):
            cumsum += l
            if cumsum >= total / 2:
                n50 = l
                break

        return {
            "total_reads": len(self.reads),
            "total_bases": total,
            "mean_length": round(statistics.mean(lengths), 1),
            "median_length": statistics.median(lengths),
            "n50": n50,
            "max_length": max(lengths),
            "min_length": min(lengths),
            "mean_quality": round(statistics.mean(qualities), 1),
            "median_quality": round(statistics.median(qualities), 1),
            "mean_gc": round(statistics.mean(gcs), 1),
            "reads_q20": sum(1 for q in qualities if q >= 20),
            "reads_q20_pct": round(
                sum(1 for q in qualities if q >= 20) / len(qualities) * 100, 1
            ),
        }


if __name__ == "__main__":
    import sys
    import json

    analyzer = FastqAnalyzer(sys.argv[1])
    analyzer.parse(max_reads=50000)
    stats = analyzer.summary()
    print(json.dumps(stats, indent=2))
