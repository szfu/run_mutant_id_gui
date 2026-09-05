# Mutant ID Analyzer

Mutant ID Analyzer is a small Python GUI for CRISPR/Cas mutant identification from Sanger sequencing `.ab1` files or FASTA mutant sequences. It aligns mutant sequencing reads to a pasted PCR amplicon reference sequence, calls local DNA variants, projects the variants onto a full CDS, translates the mutant protein, and exports publication-style DNA/protein mutation summary figures.

The figure style is designed for manuscript mutation panels: wild-type and mutant DNA sequence alignment, wild-type and mutant protein sequence alignment, changed bases and amino acids highlighted in red, optional gene-structure or protein-domain schematic, and editable PDF output.
<img width="1226" height="555" alt="image" src="https://github.com/user-attachments/assets/71cc528e-d8a5-4309-8d3e-4ab532d9bf4b" />


## Features

- Paste PCR amplicon reference sequence directly into the GUI.
- Paste full CDS sequence for coding-coordinate annotation and protein translation.
- Enter sgRNA sequence; the program locates it automatically in the PCR reference and CDS.
- Analyze one or more Sanger `.ab1` files or FASTA mutant sequence files in one run.
- Automatically test forward and reverse-complement read orientation.
- Detect substitutions, insertions, and deletions.
- Display local DNA sequence around the sgRNA first base.
- Display local WT and mutant protein sequence around the sgRNA-targeted amino acid.
- Highlight mutant DNA bases and amino acids that differ from WT.
- Show insertions or deletions with `-` gap characters.
- Draw optional gene-structure or protein-domain schematics.
- Customize sequence labels, domain labels, colors, font sizes, and vertical spacing.
- Export PNG and editable PDF figures.
- Export per-sample variant tables and a batch summary table.

## Files

```text
mutant_id_gui.py             GUI application
mutant_id_core.py            sequence loading, alignment, variant calling, translation, and plotting
requirements-mutant-id.txt   Python dependencies
run_mutant_id_gui.bat        Windows launcher
README.md                    GitHub help document
```

## Requirements

- Python 3.10 or newer
- Tkinter support, included with most standard Python installations on Windows

Python packages:

```text
biopython>=1.83
matplotlib>=3.8
numpy>=1.26
pillow>=10.0
```

Install dependencies:

```bash
python -m pip install -r requirements-mutant-id.txt
```

If another Python environment already contains incompatible scientific packages, create a clean virtual environment first:

```bash
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements-mutant-id.txt
```

## How to Run

On Windows, double-click:

```text
run_mutant_id_gui.bat
```

Or run from a terminal:

```bash
python mutant_id_gui.py
```

## Input Workflow

1. Paste the PCR amplicon reference sequence into **PCR amplicon reference sequence**. FASTA text is accepted; header lines beginning with `>` are ignored.
2. Enter the sgRNA sequence.
3. Paste the full CDS sequence into **Full CDS sequence**. This sequence is used for CDS-coordinate mapping and protein translation.
4. Set **Sequence bp before / after sgRNA**. These values are anchored on the first base of the sgRNA, not on the full sgRNA span. The default is `20 / 40`, so the displayed DNA window contains 60 reference bases.
5. Set **Protein aa before / after sgRNA**. These values are anchored on the amino acid containing the sgRNA first base. The default is `20 / 40`.
6. Adjust sequence layout parameters if needed:
   - **Reference-mutant row gap**
   - **60-base block gap**
   - **DNA-to-protein panel gap**
   - **sgRNA line-to-DNA gap**
   - **DNA font size**
   - **Protein font size**
   - **Structure label font size**
7. Enter **Sequence labels**. The first line is the reference label. Each following line is used for one mutant sample in input-file order. Mutant labels are shown in italics.
8. Optional: choose **Gene structure** or **Protein domain** and enter structure information.
9. Add one or more `.ab1` sequencing files or FASTA mutant-sequence files.
10. Click **Run analysis**.

## Structure or Domain Panel

The schematic panel is optional.

For **Gene structure**:

- Enter the total gene/CDS length in base pairs.
- Enter the sgRNA start position in the same base-pair coordinate system.
- Add exon regions with name, start, end, and color fields.
- Exons are drawn as black boxes on a horizontal line.

For **Protein domain**:

- Enter the total protein length in amino acids.
- Enter the sgRNA start position in amino-acid coordinates.
- Add domain regions with name, start, end, and color fields.
- Domains are drawn as colored boxes on a horizontal line.

The sgRNA start position is connected from the structure panel to the corresponding DNA sequence position. The program uses the start coordinate only for this connector.

## Default Display Settings

```text
DNA bases before / after sgRNA:      20 / 40
Protein amino acids before / after:  20 / 40
Reference-mutant row gap:            0.4
60-base block gap:                   0.5
DNA-to-protein panel gap:            0.25
sgRNA line-to-DNA gap:               0.2
DNA font size:                       10
Protein font size:                   10
Structure label font size:           10
Line width:                          60 characters per row
```

## Output

Results are written to:

```text
mutant_id_output/<timestamp>/
```

Each sample folder contains:

```text
<sample>_mutation_summary.png
<sample>_mutation_summary.pdf
<sample>_variants.csv
<sample>_report.txt
```

The batch output folder also contains:

```text
batch_summary.csv
publication_mutation_summary.png
publication_mutation_summary.pdf
```

The PDF files use editable text where supported by the PDF editor. They are suitable for further figure assembly in tools such as Adobe Illustrator.

## Analysis Method

For each Sanger sequencing file, the program reads the base-called sequence from the `.ab1` file using Biopython. For FASTA input, the program reads the first sequence record from the file. The read or FASTA sequence is aligned against the user-supplied PCR amplicon reference sequence, and both forward and reverse-complement orientations are tested. The orientation with the better alignment score is used for variant calling.

Sequence differences are classified as substitutions, insertions, or deletions according to the read-to-reference alignment. The sgRNA sequence is independently located in the PCR amplicon reference and in the full CDS. The displayed DNA window is defined by user-selected flanks around the first base of the sgRNA. Within this local window, WT and mutant DNA sequences are displayed with gap characters for indels and red highlighting for non-matching bases.

The PCR amplicon reference is mapped to the full CDS. Variants inside the displayed PCR-reference window are projected onto the full CDS, producing a mutant CDS sequence. The WT and mutant CDS sequences are translated codon by codon, and translation stops at the first stop codon. The WT and mutant protein sequences are aligned, and the local protein window around the sgRNA-targeted amino acid is displayed. Amino acids that differ from WT are highlighted in red.

## Suggested Methods Text

```text
Sanger sequencing reads in .ab1 format or base-called mutant sequences in FASTA format were analyzed using a custom Python graphical tool. The base-called sequence from each read or FASTA record was aligned to a user-supplied PCR amplicon reference sequence, and both forward and reverse-complement orientations were evaluated automatically. Sequence differences were classified as substitutions, insertions, or deletions from the read-to-reference alignment. The sgRNA sequence was located in both the PCR amplicon reference and the full CDS sequence. A user-defined local nucleotide window around the first base of the sgRNA was used for display of wild-type and mutant DNA alignments. Insertions and deletions were represented by gap characters, and non-matching bases were highlighted in red. Variants within the selected PCR-reference window were projected onto the full CDS, and the resulting mutant CDS was translated codon by codon until the first stop codon. Wild-type and mutant protein sequences were aligned, and a user-defined local amino-acid window around the sgRNA-targeted amino acid was displayed. Publication-style figures were exported as PNG and editable PDF files.
```

## Notes and Limitations

- The PCR amplicon reference can be shorter than the full CDS, but it must cover the sgRNA target region.
- The full CDS should come from the same gene/reference version as the PCR amplicon.
- The sgRNA is detected by exact sequence matching.
- If the sgRNA matches multiple positions, inspect the input sequence carefully before interpreting the result.
- The tool uses the base-called sequence from `.ab1` files or the sequence from FASTA files. Sanger chromatogram traces are read from the AB1 file format but are not plotted in the current publication summary figure.
- Very noisy or mixed sequencing reads may require manual review in dedicated Sanger trace software.
- This tool is intended for research visualization and screening support, not clinical or diagnostic use.

## Troubleshooting

### The BAT file closes immediately

Open a terminal in the program folder and run:

```bash
run_mutant_id_gui.bat
```

or:

```bash
python mutant_id_gui.py
```

This keeps the error message visible.

### Missing package error

Install dependencies:

```bash
python -m pip install -r requirements-mutant-id.txt
```

### NumPy or binary compatibility error

Use a clean virtual environment:

```bash
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements-mutant-id.txt
python mutant_id_gui.py
```

### sgRNA not found

Check that:

- the sgRNA sequence is in the PCR amplicon reference
- the sgRNA sequence is in the full CDS
- the sequence was pasted in the correct orientation
- there are no extra non-DNA characters

### Protein sequence looks incorrect

Check that the CDS starts at the correct ATG and is in frame. The CDS input should be the coding sequence, not the genomic sequence with introns.

## Sharing

To share the program, include:

```text
mutant_id_gui.py
mutant_id_core.py
requirements-mutant-id.txt
run_mutant_id_gui.bat
README.md
```

Do not include private sequencing data, unpublished experimental data, or personal output folders unless they are intended to be shared.

## License

No license file is included by default. Add a license before public release if other users are allowed to reuse, modify, or redistribute the code.
