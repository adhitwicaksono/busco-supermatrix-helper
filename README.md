# BUSCO Supermatrix Helper

Lightweight Python helper scripts for preparing BUSCO-derived sequences for phylogenomic analysis.

This workflow was developed to solve a practical problem encountered during comparative genome analysis: BUSCO sequence FASTA files may contain coordinate-based headers rather than BUSCO identifiers, making it difficult to group orthologous BUSCO sequences across multiple samples.

The scripts in this repository help convert BUSCO outputs into ortholog-wise FASTA files suitable for MAFFT alignment, concatenation, and downstream phylogenetic tree inference.

---

## Overview

The workflow converts BUSCO outputs from this structure:

```text
Sample_A_BUSCO.fasta
  >chr1:1000-2500|+
  ATG...

  >chr2:5000-7200|-
  ATG...

Sample_B_BUSCO.fasta
  >chr3:1200-2600|+
  ATG...

  >chr5:8000-9500|-
  ATG...
```

into this structure:

```text
BUSCO_ID_001.fasta
  >Sample_A
  ATG...

  >Sample_B
  ATG...

BUSCO_ID_002.fasta
  >Sample_A
  ATG...

  >Sample_B
  ATG...
```

Each BUSCO ortholog can then be aligned separately before concatenation into a phylogenomic supermatrix.

---

## Workflow

```text
BUSCO full_table.tsv + BUSCO sequence FASTA
        ↓
rename_busco_headers.py
        ↓
sample|BUSCO_ID FASTA files
        ↓
group_buscos_by_id.py
        ↓
one FASTA per BUSCO ortholog
        ↓
MAFFT alignment per BUSCO
        ↓
concatenate_busco_alignments.py
        ↓
supermatrix FASTA + partition file
        ↓
IQ-TREE / FastTree / other phylogenetic tools
```

---

## Repository structure

```text
busco-supermatrix-helper/
├── README.md
├── scripts/
│   ├── rename_busco_headers.py
│   ├── group_buscos_by_id.py
│   └── concatenate_busco_alignments.py
├── examples/
│   ├── example_full_table.tsv
│   ├── example_busco_sequences.fasta
│   └── expected_output/
├── docs/
│   └── workflow_diagram.png
├── LICENSE
└── CITATION.cff
```

---

## Requirements

The Python scripts require only standard Python libraries.

Recommended:

```bash
python >= 3.8
mafft
iqtree2
```

Optional:

```bash
fasttree
```

---

## Input files

For each sample, the workflow requires:

1. BUSCO `full_table.tsv`
2. BUSCO nucleotide or protein FASTA file

Example:

```text
CI_full_table.tsv
CI_busco_sequences.fasta
```

The BUSCO FASTA headers are expected to contain genomic coordinate-style identifiers, such as:

```text
>NC_089035.1:42474564-42477069|-
```

The BUSCO `full_table.tsv` is used to map these coordinates back to BUSCO identifiers.

---

## Step 1 — Rename BUSCO FASTA headers

Use `rename_busco_headers.py` to convert coordinate-based FASTA headers into:

```text
>sample|BUSCO_ID
```

Example:

```bash
python3 scripts/rename_busco_headers.py \
  --sample CI \
  --full_table CI_full_table.tsv \
  --fasta CI_busco_sequences.fasta \
  --out CI_BUSCO_complete_renamed.fasta \
  --summary CI_BUSCO_rename_summary.txt \
  --unmatched CI_BUSCO_unmatched_headers.txt
```

Example output headers:

```text
>CI|2at3193
>CI|5at3193
>CI|10at3193
```

Only BUSCO records annotated as `Complete` are retained.

Duplicated, fragmented, missing, and unmatched sequences are excluded.

---

## Step 2 — Group BUSCO sequences by BUSCO ID

After renaming all sample FASTA files, group orthologous BUSCO sequences across samples using `group_buscos_by_id.py`.

Example:

```bash
python3 scripts/group_buscos_by_id.py \
  --input AGIS1.0_BUSCO_complete_renamed.fasta CI_BUSCO_complete_renamed.fasta CL_BUSCO_complete_renamed.fasta GJ_BUSCO_complete_renamed.fasta IN_BUSCO_complete_renamed.fasta BSM_BUSCO_complete_renamed.fasta KM_BUSCO_complete_renamed.fasta MPE_BUSCO_complete_renamed.fasta PP_BUSCO_complete_renamed.fasta \
  --outdir grouped \
  --samples AGIS1.0 CI CL GJ IN BSM KM MPE PP \
  --shared_only \
  --summary BUSCO_grouping_summary.txt \
  --manifest BUSCO_grouped_manifest.tsv
```

This produces one FASTA file per shared BUSCO ID:

```text
grouped/
├── 2at3193.fasta
├── 5at3193.fasta
├── 10at3193.fasta
└── ...
```

Each grouped FASTA contains one sequence per sample:

```text
>AGIS1.0
ATG...

>CI
ATG...

>CL
ATG...
```

---

## Step 3 — Align each BUSCO ortholog with MAFFT

Run MAFFT on each grouped BUSCO FASTA.

For nucleotide or protein BUSCO sequences, the same basic MAFFT command can be used:

```bash
mkdir -p aligned

for f in grouped/*.fasta
do
  base=$(basename "$f" .fasta)
  mafft --auto --thread 1 "$f" > aligned/${base}.aln.fasta
done
```

This produces:

```text
aligned/
├── 2at3193.aln.fasta
├── 5at3193.aln.fasta
├── 10at3193.aln.fasta
└── ...
```

---

## Step 4 — Concatenate aligned BUSCOs

After MAFFT alignment, concatenate all aligned BUSCOs into one supermatrix using `concatenate_busco_alignments.py`.

Example:

```bash
python3 scripts/concatenate_busco_alignments.py \
  --aligned_dir aligned \
  --samples AGIS1.0 CI CL GJ IN BSM KM MPE PP \
  --out BUSCO_supermatrix.fasta \
  --partition BUSCO_partitions.txt \
  --summary BUSCO_concatenation_summary.txt
```

Expected final FASTA structure:

```text
>AGIS1.0
ATG...

>CI
ATG...

>CL
ATG...
```

The resulting supermatrix can be used for phylogenetic inference.

---

## Step 5 — Build a phylogenetic tree

### IQ-TREE example for nucleotide BUSCO supermatrix

```bash
iqtree2 \
  -s BUSCO_supermatrix.fasta \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T AUTO
```

### IQ-TREE example for protein BUSCO supermatrix

```bash
iqtree2 \
  -s BUSCO_supermatrix.fasta \
  -st AA \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T AUTO
```

### FastTree example for nucleotide alignment

```bash
FastTree -nt -gtr BUSCO_supermatrix.fasta > BUSCO_tree.nwk
```

---

## Outputs

Main outputs:

```text
*_BUSCO_complete_renamed.fasta
*_BUSCO_rename_summary.txt
*_BUSCO_unmatched_headers.txt
BUSCO_grouping_summary.txt
BUSCO_grouped_manifest.tsv
grouped/*.fasta
aligned/*.aln.fasta
BUSCO_supermatrix.fasta
BUSCO_partitions.txt
BUSCO_tree.nwk
```

---

## Important notes

### Do not align all BUSCOs from one sample together

This is incorrect:

```text
CI_BUSCO_complete_renamed.fasta → MAFFT
```

because the file contains many unrelated BUSCO genes from the same sample.

Instead, BUSCOs must first be grouped by BUSCO ID across samples.

Correct:

```text
2at3193.fasta → MAFFT
5at3193.fasta → MAFFT
10at3193.fasta → MAFFT
```

Each grouped FASTA should contain the same BUSCO ortholog across all samples.

---

## Recommended use cases

This workflow is useful when:

- BUSCO sequence headers are coordinate-based rather than BUSCO-ID-based
- BUSCO outputs were generated separately for multiple genomes
- users want manual control over ortholog grouping and alignment
- users are combining Galaxy-based BUSCO analysis with local MAFFT/IQ-TREE workflows
- users want to prepare nucleotide or protein BUSCO supermatrices for phylogenomics

---

## Limitations

This workflow does not replace full phylogenomics pipelines such as BuscoPhylo, BUSCO_Phylogenomics, OrthoFinder, or Phylociraptor.

It is a lightweight helper workflow for reformatting and preparing BUSCO-derived sequences.

Potential limitations:

- It assumes BUSCO coordinate headers can be matched to `full_table.tsv`
- It retains only BUSCO records marked as `Complete`
- It does not automatically resolve paralogy beyond excluding non-complete or unmatched records
- It does not perform alignment trimming
- It does not perform tree inference internally
- It assumes that BUSCO IDs shared across samples represent comparable orthologous loci

For complex genomes, additional filtering may be needed, especially for:

- high paralogy
- high missingness
- contamination
- horizontal gene transfer
- long-branch attraction
- highly fragmented assemblies

---

## Example application

This workflow was first tested on comparative rice genome data using BUSCO outputs from:

```text
AGIS1.0
CI
CL
GJ
IN
BSM
KM
MPE
PP
```

It was used to prepare shared BUSCO nucleotide and protein datasets for MAFFT alignment, supermatrix concatenation, and IQ-TREE phylogenetic reconstruction.

---

## AI-assisted development statement

Parts of this workflow were developed with assistance from ChatGPT (OpenAI), which was used for workflow planning, coding support, and debugging of Python scripts.

All scripts, outputs, and interpretations were reviewed, tested, and verified by the authors. ChatGPT was not used as an author and did not independently generate scientific conclusions.

---

## Citation and acknowledgment

Please cite the relevant tools used in this workflow:

- BUSCO
- MAFFT
- IQ-TREE
- FastTree, if used

Related tools and workflows include:

- BuscoPhylo
- BUSCO_Phylogenomics
- OrthoFinder
- Phylociraptor

This repository is intended as a lightweight helper workflow, not a replacement for these established phylogenomics tools.

---

## License

MIT License is recommended for this repository.

---

## Author

Developed by Adhityo Wicaksono. Project led by Yekti Asih Purwestri (Universitas Gadjah Mada, Indonesia).
