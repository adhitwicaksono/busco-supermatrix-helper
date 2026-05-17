#!/usr/bin/env python3

from pathlib import Path
import argparse
from collections import defaultdict, Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Group renamed BUSCO FASTA sequences by BUSCO ID across samples."
    )

    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="Renamed BUSCO FASTA files with headers formatted as sample|BUSCO_ID"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for grouped BUSCO FASTA files"
    )

    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Optional sample order. Example: AGIS1.0 CI CL GJ IN BSM KM MPE PP"
    )

    parser.add_argument(
        "--shared_only",
        action="store_true",
        help="Keep only BUSCO IDs present in all samples"
    )

    parser.add_argument(
        "--summary",
        default="BUSCO_grouping_summary.txt",
        help="Output summary file"
    )

    parser.add_argument(
        "--manifest",
        default="BUSCO_grouped_manifest.tsv",
        help="Output manifest table"
    )

    return parser.parse_args()


def read_fasta(path):
    records = {}
    current = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(seq_parts)

                current = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if current is not None:
        records[current] = "".join(seq_parts)

    return records


def main():
    args = parse_args()

    input_files = [Path(f) for f in args.input]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # BUSCO_ID -> sample -> sequence
    busco_groups = defaultdict(dict)

    # Track samples
    detected_samples = []

    # Track duplicate sample|BUSCO_ID entries
    duplicate_entries = []

    for fasta_file in input_files:
        records = read_fasta(fasta_file)

        seen_in_file = Counter()

        for header, seq in records.items():
            if "|" not in header:
                raise ValueError(
                    f"Header does not contain '|': {header}\n"
                    "Expected format: sample|BUSCO_ID"
                )

            sample, busco_id = header.split("|", 1)

            detected_samples.append(sample)
            seen_in_file[(sample, busco_id)] += 1

            if sample in busco_groups[busco_id]:
                duplicate_entries.append((sample, busco_id, fasta_file.name))
                continue

            busco_groups[busco_id][sample] = seq

        for (sample, busco_id), count in seen_in_file.items():
            if count > 1:
                duplicate_entries.append((sample, busco_id, fasta_file.name))

    detected_samples = sorted(set(detected_samples))

    if args.samples:
        sample_order = args.samples
    else:
        sample_order = detected_samples

    expected_sample_count = len(sample_order)

    written = 0
    skipped_missing = 0
    manifest_lines = [
        "BUSCO_ID\tN_samples\tSamples_present\tOutput_file"
    ]

    for busco_id in sorted(busco_groups.keys()):
        sample_to_seq = busco_groups[busco_id]
        present_samples = [s for s in sample_order if s in sample_to_seq]

        if args.shared_only and len(present_samples) != expected_sample_count:
            skipped_missing += 1
            continue

        out_fasta = outdir / f"{busco_id}.fasta"

        with open(out_fasta, "w") as out:
            for sample in present_samples:
                out.write(f">{sample}\n")
                seq = sample_to_seq[sample]
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")

        written += 1
        manifest_lines.append(
            f"{busco_id}\t{len(present_samples)}\t"
            f"{','.join(present_samples)}\t{out_fasta.name}"
        )

    summary_lines = [
        f"Input FASTA files: {len(input_files)}",
        f"Detected samples: {', '.join(detected_samples)}",
        f"Sample order used: {', '.join(sample_order)}",
        f"Total BUSCO IDs detected: {len(busco_groups)}",
        f"Grouped BUSCO FASTA files written: {written}",
        f"BUSCO IDs skipped because not shared by all samples: {skipped_missing}",
        f"Duplicate sample|BUSCO_ID entries detected: {len(duplicate_entries)}",
        f"Output directory: {outdir}",
        ""
    ]

    with open(args.summary, "w") as f:
        f.write("\n".join(summary_lines))

    with open(args.manifest, "w") as f:
        f.write("\n".join(manifest_lines) + "\n")

    print("\n".join(summary_lines))

    if duplicate_entries:
        dup_file = Path("BUSCO_grouping_duplicates.tsv")
        with open(dup_file, "w") as f:
            f.write("Sample\tBUSCO_ID\tFile\n")
            for sample, busco_id, fname in duplicate_entries:
                f.write(f"{sample}\t{busco_id}\t{fname}\n")
        print(f"Duplicate details written to: {dup_file}")


if __name__ == "__main__":
    main()