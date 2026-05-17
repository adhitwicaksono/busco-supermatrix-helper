from pathlib import Path

aligned_dir = Path("aligned")
output_fasta = Path("BUSCO_supermatrix.fasta")
partition_file = Path("BUSCO_partitions.txt")

samples = ["AGIS1.0", "CI", "CL", "GJ", "IN", "BSM", "KM", "MPE", "PP"]

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

supermatrix = {s: [] for s in samples}
partitions = []
start = 1

files = sorted(aligned_dir.glob("*.aln.fasta"))

for aln_file in files:
    busco_id = aln_file.stem.replace(".aln", "")
    records = read_fasta(aln_file)

    missing = [s for s in samples if s not in records]
    if missing:
        print(f"Skipping {busco_id}, missing samples: {missing}")
        continue

    lengths = {len(records[s]) for s in samples}
    if len(lengths) != 1:
        print(f"Skipping {busco_id}, unequal alignment lengths: {lengths}")
        continue

    aln_len = lengths.pop()
    end = start + aln_len - 1
    partitions.append((busco_id, start, end))

    for s in samples:
        supermatrix[s].append(records[s])

    start = end + 1

with open(output_fasta, "w") as out:
    for s in samples:
        seq = "".join(supermatrix[s])
        out.write(f">{s}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

with open(partition_file, "w") as out:
    for busco_id, start, end in partitions:
        out.write(f"DNA, {busco_id} = {start}-{end}\n")

print(f"Concatenated BUSCO alignments: {len(partitions)}")
print(f"Final alignment length: {start - 1} bp")
print(f"Output FASTA: {output_fasta}")
print(f"Partition file: {partition_file}")