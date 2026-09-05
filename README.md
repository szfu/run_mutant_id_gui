# Mutant ID Analyzer

A small Python GUI for analyzing Sanger sequencing `.ab1` files or FASTA mutant sequences from mutant identification experiments.

The program:

- reads one or more `.ab1` chromatogram files or FASTA sequence files
- accepts a pasted PCR amplicon DNA/FASTA sequence in the GUI text box
- aligns each read to a user-supplied PCR amplicon reference sequence
- detects substitutions, insertions, and deletions
- marks the sgRNA position on the PCR reference sequence
- highlights sequence differences in the DNA alignment
- translates the coding sequence in codons and stops at the first stop codon
- supports an optional gene-structure or protein-domain schematic
- supports multiple mutant samples in one run
- exports publication-style PNG and editable PDF figures

<img width="1210" height="534" alt="image" src="https://github.com/user-attachments/assets/e3ee45f5-5ad4-4dd8-ac58-0aface8a82cf" />

## Files

```text
mutant_id_gui.py             GUI application
mutant_id_core.py            sequence loading, alignment, variant calling, translation, and plotting
requirements-mutant-id.txt   Python dependencies
run_mutant_id_gui.bat        Windows launcher
README.md                    English documentation
```

## Requirements

- Python 3.10 or newer
- Tkinter support

Python packages:

```text
biopython
matplotlib
numpy
pillow
```

Install them with:

```bash
python -m pip install -r requirements-mutant-id.txt
```

## How to Run

On Windows, double-click:

```text
run_mutant_id_gui.bat
```

Or start it from a terminal:

```bash
python mutant_id_gui.py
```

## Input Workflow

1. Paste the PCR amplicon reference sequence or FASTA text directly into the **PCR amplicon reference sequence** box. This can be the sequence between the sequencing primers; it does not need to be the full gene or full CDS.
2. Enter the sgRNA sequence.
3. Paste the full CDS sequence in the **Full CDS sequence** box. This is used only for cDNA-coordinate annotation and protein translation.
4. Enter how many nucleotides to show before and after the sgRNA first base. These two values add up to the local WT/mutant DNA display length, which wraps at 60 reference bases per row.
5. Enter how many amino acids to show before and after the amino acid containing the sgRNA first base. These two values add up to the displayed protein window, which wraps at 60 amino acids per row.
6. Set **Reference-mutant row gap** to control the vertical distance between Reference/Mutant DNA rows and WT/mutant protein rows within each block.
7. Set **60-base block gap** to control the vertical space before the next 60-reference-base DNA block. The same setting is also used between wrapped protein blocks.
8. Set **DNA-to-protein panel gap** to control the vertical space between the DNA panel and the protein panel.
9. Set **sgRNA line-to-DNA gap** to control the vertical distance between the sgRNA line and the first DNA row.
10. Set **DNA font size**, **Protein font size**, and **Structure label font size** independently. The structure setting controls region names such as `NBS`, `LRR domain`, `TM`, and `Kinase domain`.
11. Enter **Sequence labels**: the first line is the reference label, and each following line is the label for one mutant in input-file order. Mutant labels are shown in italics.
12. Choose **Gene structure** or **Protein domain** if you want the schematic panel.
13. Enter the total CDS length or protein length in **Structure total length**.
14. Enter the sgRNA start position in the same structure coordinates. For protein domains, use the amino-acid coordinate; for gene structures, use the CDS base-pair coordinate. For example, enter `249` when the sgRNA begins at amino acid 249. Do not enter the total protein/CDS length here.
15. Add one or more structure rows with a name, start, end, and color.
16. Use CDS bp positions for gene structure and amino-acid positions for protein domain.
17. Gene structure exons are shown as black boxes; protein-domain regions use the selected colors.
18. Add one or more `.ab1` sequencing files or FASTA mutant-sequence files.
19. Click **Run analysis**.

## Default Display Settings

- DNA bases before / after sgRNA: `20 / 40`
- Protein amino acids before / after sgRNA: `20 / 40`
- Reference-mutant row gap: `0.4`
- 60-base block gap: `0.5`
- DNA-to-protein panel gap: `0.25`
- sgRNA line-to-DNA gap: `0.2`
- DNA font size: `10`
- Protein font size: `10`
- Structure label font size: `10`

## What the Program Does

For each `.ab1` or FASTA file, the software:

1. Reads the called sequence from each Sanger `.ab1` file or the first sequence record from each FASTA file. Chromatogram channels are accepted as part of the AB1 file format but are not plotted in the final figure.
2. Automatically checks both forward and reverse-complement orientations and keeps the better alignment.
3. Aligns the read to the supplied PCR amplicon reference sequence.
4. Calls sequence differences as substitutions, insertions, or deletions.
5. Locates the sgRNA directly in the full CDS and independently checks that the PCR amplicon covers the same target.
6. Independently maps the PCR amplicon reference onto the supplied full CDS, in either orientation.
7. Maps called sequence differences from PCR-reference coordinates to full-CDS coordinates, then translates the mutated CDS codon by codon and stops at the first stop codon.
8. Generates a summary figure with:
   - a local PCR-reference/mutant DNA alignment panel for the user-selected nucleotide window; insertion and deletion columns are shown with `-`, changed bases are colored red, the sgRNA interval is marked above the alignment, and long sequences wrap at 60 reference bases per line
   - aligned WT and mutant protein sequences, limited to the user-selected amino-acid window and wrapped at 60 columns per line
   - an optional gene-structure or protein-domain schematic at the top when structure regions are provided, with user-defined names, spans, colors, and sgRNA position; the sgRNA is linked to the sequence panel with a connector line

## Output

Results are written to:

```text
mutant_id_output/<timestamp>/<sample_name>/
```

Each sample folder contains:

```text
<sample>_mutation_summary.png
<sample>_mutation_summary.pdf
<sample>_variants.csv
<sample>_report.txt
```

A batch summary file is also created:

```text
batch_summary.csv
publication_mutation_summary.png
publication_mutation_summary.pdf
```

The `publication_mutation_summary` figure contains a shared local PCR-reference/mutant DNA alignment panel and aligned WT/mutant protein sequences for the selected amino-acid window. It does not display the full PCR amplicon, the full translated protein, a CDS schematic, or Sanger chromatograms.

## Notes

- Use a clean PCR-product reference sequence without spaces or numbering. It may be shorter than the full CDS.
- FASTA headers beginning with `>` are ignored automatically when pasted into either sequence box.
- Paste the full CDS from the same gene/reference version. The PCR product and CDS must share a sufficient overlapping sequence, but neither input needs to contain the other completely.
- The program locates the sgRNA from the supplied full CDS and sgRNA sequence; no manual CDS coordinate is required.
- The sgRNA must map to one unambiguous position in the full CDS.
- The sgRNA line above the DNA alignment is positioned from the sgRNA sequence match within the PCR reference.
- Reports include both PCR-reference coordinates and mapped full-CDS `c.` coordinates. Differences outside the mapped coding overlap are reported as outside the mapped CDS.
- Uncovered ends of a longer reference sequence are not counted as deletions.
- The program detects exact sgRNA matches only.
- `.ab1` files and FASTA files are used as the source of called mutant bases for alignment and variant detection; chromatograms are not included in exported figures.

## Suggested Methods Text

```text
Sanger sequencing reads in .ab1 format or base-called mutant sequences in FASTA format were analyzed using a custom Python graphical tool. The called sequence from each read or FASTA record was aligned against a user-supplied PCR amplicon reference sequence, and the program automatically evaluated both forward and reverse-complement orientations. Sequence differences were classified as substitutions, insertions, or deletions. The sgRNA sequence was located directly in the user-supplied full CDS sequence and independently matched to the PCR amplicon reference. User-selected nucleotide flanks around the sgRNA defined the local PCR-reference/mutant DNA region displayed in the summary figure. Local reference and mutant sequences were aligned with gap characters used to represent insertions or deletions, and all non-matching bases were displayed in red. Only variants inside this selected PCR-reference window were projected onto the complete CDS for protein translation. The resulting mutant protein sequence was globally aligned to the WT protein sequence, and only the user-selected local amino-acid window was displayed. Translation was stopped at the first stop codon. Publication-style summary figures were exported as PNG and editable PDF files.
```

## Sharing the Program

If you share this tool with someone else, include:

- the Python files
- the requirements file
- the Windows launcher
- this README

Do not include private sequencing data unless it is intended to be shared.
