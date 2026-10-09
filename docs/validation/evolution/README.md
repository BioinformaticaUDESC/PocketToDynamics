# Evolutionary conservation companion analysis

This reproducible companion analysis projects a descriptive sequence-conservation profile for *Escherichia coli* K-12 adenylate kinase (UniProt P69441) onto the Ap5A-bound 1AKE structure. It is separate from the M01–M21 screening and MD workflow and does not alter pocket detection, candidate eligibility, or docking ranks.

## Data snapshot

UniProtKB was queried on 2026-10-09 with:

```text
ec:2.7.4.3 AND reviewed:true AND taxonomy_id:2 AND length:[190 TO 230]
```

The raw TSV and query metadata are included here. The deterministic selector retains one reviewed sequence per species by minimizing length difference from the 214-residue reference (accession breaks ties), then removes exact duplicate sequences. The resulting 405-sequence FASTA and metadata are archived. This selection is a pragmatic enzyme-annotation set; orthology was not independently established.

## Reproduction

Run from the repository root in the configured environment (Biopython, RDKit and Matplotlib required):

```bash
python scripts/validation/evolution/prepare_sequences.py \
  --input docs/validation/evolution/uniprot_adenylate_kinases.tsv \
  --reference P69441 \
  --output work/evolution/analysis_set.fasta

python scripts/validation/evolution/build_reference_msa.py \
  --input work/evolution/analysis_set.fasta \
  --reference P69441 \
  --output work/evolution/reference_guided_msa.fasta

python scripts/validation/evolution/map_conservation.py \
  --alignment work/evolution/reference_guided_msa.fasta \
  --structure docs/validation/adk/inputs/1AKE.pdb \
  --ligand docs/validation/adk/inputs/reference.sdf \
  --chain A --reference P69441 --contact-cutoff 5.0 \
  --out work/evolution/results

python scripts/validation/evolution/plot_conservation.py \
  --input work/evolution/results/adk_conservation.tsv \
  --output work/evolution/results/evolution_conservation.png
```

`build_reference_msa.py` uses global pairwise alignments to the reference with BLOSUM62, gap-open −10 and gap-extension −0.5. It retains sequences with at least 75% reference coverage and 25% identity over paired positions, then expands insertions into a reference-anchored matrix. This is not a phylogeny-aware or iterative multiple-sequence alignment. Conservation is one minus Shannon entropy normalized by log2(20), with gaps excluded from entropy and reported as coverage. Structural mapping uses a global sequence alignment to chain A. A residue is in the Ap5A contact set if any receptor heavy atom is within 5 Å of any ligand heavy atom.

## Recorded result

The archived output contains 405 sequences, 610 alignment columns, all 214 reference positions mapped, and 45 Ap5A-contact residues. Mean conservation was 0.7878 for the contact set and 0.5169 for the remaining 169 mapped positions. These values are descriptive; there is no phylogenetic weighting, correction for sequence dependence, inferential test, or affinity prediction.

Files `analysis_set.fasta`, `reference_guided_msa.fasta`, `adk_conservation.tsv`, `summary.json`, and the query metadata preserve the inputs and outputs. `SHA256SUMS` records their SHA-256 checksums. The matching manuscript figure is in `docs/paper/figures/evolutionary_conservation.png` and `.svg`.
