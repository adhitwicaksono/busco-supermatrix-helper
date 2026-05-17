#!/usr/bin/env python3

from pathlib import Path
import argparse
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rename BUSCO FASTA headers using BUSCO full_table.tsv."
    )

    parser.add_argument(
        "--sample",
        required=True,
        help="Sample name to use in FASTA headers, e.g. CI, AGIS1.0, CL"
    )

    parser.add_argument(
        "--full_table",
        required=True,
        help="BUSCO full_table.tsv file"
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="BUSCO sequence FASTA file with coordinate-style headers"
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output renamed FASTA"
    )

    parser.add_argument(
        "--summary",
        default=None,
        help="Optional output summary text file"
    )

    parser.add_argument(
        "--unmatched",
        default=None,
        help="Optional output file listing unmatched FASTA headers"
    )

    return parser.parse_args()


def read_busco_table(full_table):
    """
    Build coordinate-to-BUSCO mapping from BUSCO full_table.tsv.

    Expected BUSCO full_table columns:
    BUSCO_ID, Status, Sequence, Gene Start, Gene End, Strand, Score, Length, ...
    """

    coord_to_busco = {}
    busco_ids = []

    with open(full_table) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 6:
                continue

            busco_id = parts[0]
            status = parts[1]
            seqid = parts[2]
            start = parts[3]
            end = parts[4]
            strand = parts[5]

            # Keep only Complete BUSCOs.
            # Exclude Duplicated, Fragmented, Missing.
            if status != "Complete":
                continue

            # FASTA headers may appear as either start-end or end-start,
            # especially for minus-strand hits.
            key1 = f"{seqid}:{start}-{end}|{strand}"
            key2 = f"{seqid}:{end}-{start}|{strand}"

            coord_to_busco[key1] = busco_id
            coord_to_busco[key2] = busco_id

            busco_ids.append(busco_id)

    return coord_to_busco, busco_ids


def rename_fasta(sample, fasta, coord_to_busco, out_fasta):
    input_entries = 0
    kept_entries = 0
    unmatched_headers = []
    renamed_buscos = []

    keep = False

    with open(fasta) as inp, open(out_fasta, "w") as out:
        for line in inp:
            line = line.rstrip("\n")

            if line.startswith(">"):
                input_entries += 1

                header = line[1:].strip().split()[0]

                if header in coord_to_busco:
                    busco_id = coord_to_busco[header]
                    out.write(f">{sample}|{busco_id}\n")
                    keep = True
                    kept_entries += 1
                    renamed_buscos.append(busco_id)
                else:
                    keep = False
                    unmatched_headers.append(header)

            else:
                if keep:
                    out.write(line + "\n")

    return input_entries, kept_entries, unmatched_headers, renamed_buscos


def main():
    args = parse_args()

    full_table = Path(args.full_table)
    fasta = Path(args.fasta)
    out_fasta = Path(args.out)

    coord_to_busco, complete_busco_ids = read_busco_table(full_table)

    input_entries, kept_entries, unmatched_headers, renamed_buscos = rename_fasta(
        sample=args.sample,
        fasta=fasta,
        coord_to_busco=coord_to_busco,
        out_fasta=out_fasta
    )

    duplicate_renamed = [
        busco for busco, count in Counter(renamed_buscos).items()
        if count > 1
    ]

    summary_lines = [
        f"Sample: {args.sample}",
        f"Input FASTA entries: {input_entries}",
        f"Complete BUSCO records used from table: {len(set(complete_busco_ids))}",
        f"Renamed/kept FASTA entries: {kept_entries}",
        f"Unique BUSCO IDs retained: {len(set(renamed_buscos))}",
        f"Unmatched/discarded FASTA entries: {len(unmatched_headers)}",
        f"Duplicated renamed BUSCO IDs: {len(duplicate_renamed)}",
        ""
    ]

    summary_text = "\n".join(summary_lines)

    print(summary_text)

    if args.summary:
        with open(args.summary, "w") as f:
            f.write(summary_text)

    if args.unmatched:
        with open(args.unmatched, "w") as f:
            for header in unmatched_headers:
                f.write(header + "\n")


if __name__ == "__main__":
    main()