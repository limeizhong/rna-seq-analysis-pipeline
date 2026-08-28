#!/usr/bin/env python3
"""Filter a raw tab-delimited gene-count matrix for DESeq2."""
import argparse
import csv
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-count", type=int, default=10)
    p.add_argument("--min-samples", type=int, required=True)
    a = p.parse_args()
    if a.min_count < 0 or a.min_samples < 1:
        p.error("--min-count must be >= 0 and --min-samples must be >= 1")

    with open(a.input, newline="") as src:
        rows = list(csv.reader(src, delimiter="\t"))
    if len(rows) < 2 or len(rows[0]) < 3:
        raise SystemExit("Count matrix requires a header and at least two sample columns")
    header, data = rows[0], rows[1:]
    if not header[0] or len(set(header[1:])) != len(header[1:]):
        raise SystemExit("First header cell must name genes; sample names must be unique")
    kept, zeros = [], 0
    for line_no, row in enumerate(data, 2):
        if len(row) != len(header):
            raise SystemExit(f"Line {line_no}: expected {len(header)} columns")
        try:
            values = [int(x) for x in row[1:]]
        except ValueError as exc:
            raise SystemExit(f"Line {line_no}: counts must be integers") from exc
        if any(x < 0 for x in values):
            raise SystemExit(f"Line {line_no}: counts cannot be negative")
        zeros += not any(values)
        if sum(x >= a.min_count for x in values) >= a.min_samples:
            kept.append(row)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with open(a.output, "w", newline="") as dst:
        csv.writer(dst, delimiter="\t", lineterminator="\n").writerows([header, *kept])
    print(f"Input genes: {len(data)}; all-zero: {zeros}; retained: {len(kept)}")
    print(f"Rule: count >= {a.min_count} in >= {a.min_samples} samples")


if __name__ == "__main__":
    main()
