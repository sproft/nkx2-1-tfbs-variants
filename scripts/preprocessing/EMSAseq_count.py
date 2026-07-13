#This is a Python script for counting the number of unique sequences in each processed sequencing fastq file
#Freyda Lenihan-Geels, 2025

import gzip
import sys
import argparse
from collections import Counter

def count_sequences(input_fastq, output_txt):
    sequence_counts = Counter()

    with gzip.open(input_fastq, "rt") as f:
        while True:
            try:
                _ = next(f)  # Header line
                seq = next(f).strip()  # Sequence line
                _ = next(f)  # '+' line
                _ = next(f)  # Quality line

                if len(seq) == 14:
                    sequence = seq[:14]
                    sequence_counts[sequence] += 1
            except StopIteration:
                break

    with open(output_txt, "w") as out:
        for seq, count in sequence_counts.items():
            out.write(f"{seq}\t{count}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count 14bp sequences from a FASTQ file.")
    parser.add_argument("input", help="Input FASTQ.gz file")
    parser.add_argument("-o", "--output", required=True, help="Output .txt file for sequence counts")

    args = parser.parse_args()
    count_sequences(args.input, args.output)

