from __future__ import annotations

import csv
import os
import threading
from datetime import datetime
from pathlib import Path
from tkinter import (
    BOTH,
    BOTTOM,
    END,
    LEFT,
    RIGHT,
    TOP,
    VERTICAL,
    Button,
    Canvas,
    Entry,
    Frame,
    Label,
    LabelFrame,
    Listbox,
    StringVar,
    Text,
    Tk,
    W,
    X,
    Y,
    colorchooser,
    filedialog,
    messagebox,
    ttk,
)

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib-cache"))

try:
    from PIL import Image, ImageTk
except ImportError as exc:  # pragma: no cover - user-facing error
    raise SystemExit("Pillow is required. Install dependencies with: python -m pip install -r requirements-mutant-id.txt") from exc

from mutant_id_core import (
    SampleAnalysis,
    analyze_sample,
    clean_dna,
    render_batch_publication_figure,
)
from matplotlib.colors import is_color_like, to_hex


class ScrollableFrame(Frame):
    def __init__(self, master: Tk | Frame, **kwargs):
        super().__init__(master, **kwargs)
        canvas = Canvas(self, highlightthickness=0, width=560)
        scrollbar = ttk.Scrollbar(self, orient=VERTICAL, command=canvas.yview)
        self.inner = Frame(canvas)
        self.inner.bind(
            "<Configure>",
            lambda _event: canvas.configure(scrollregion=canvas.bbox("all")),
        )
        inner_window = canvas.create_window((0, 0), window=self.inner, anchor="nw")
        canvas.bind(
            "<Configure>",
            lambda event: canvas.itemconfigure(inner_window, width=event.width),
        )
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=LEFT, fill=BOTH, expand=True)
        scrollbar.pack(side=RIGHT, fill=Y)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)


class MutantIDApp:
    def __init__(self) -> None:
        self.root = Tk()
        self.root.title("Mutant ID Analyzer")
        self.root.geometry("1660x920")

        self.reference_text: Text | None = None
        self.labels_text: Text | None = None
        self.cds_text: Text | None = None
        self.guide_var = StringVar()
        self.sequence_before_var = StringVar(value="20")
        self.sequence_after_var = StringVar(value="40")
        self.protein_before_var = StringVar(value="20")
        self.protein_after_var = StringVar(value="40")
        self.sequence_row_spacing_var = StringVar(value="0.4")
        self.dna_block_gap_var = StringVar(value="0.5")
        self.dna_protein_gap_var = StringVar(value="0.25")
        self.sgrna_line_gap_var = StringVar(value="0.2")
        self.dna_font_size_var = StringVar(value="10")
        self.protein_font_size_var = StringVar(value="10")
        self.structure_label_font_size_var = StringVar(value="10")
        self.structure_mode_var = StringVar(value="Gene structure")
        self.structure_total_var = StringVar()
        self.structure_sgrna_start_var = StringVar()
        self.output_var = StringVar(value=str(Path(__file__).resolve().parent / "mutant_id_output"))
        self.status_var = StringVar(value="Ready.")

        self.file_paths: list[Path] = []
        self.results: list[SampleAnalysis] = []
        self.structure_rows: list[dict[str, object]] = []
        self.preview_image = None
        self.last_output_root: Path | None = None
        self.batch_publication_png: Path | None = None
        self.structure_rows_frame: Frame | None = None

        self._build_ui()

    def _build_ui(self) -> None:
        outer = Frame(self.root)
        outer.pack(fill=BOTH, expand=True)

        left = ScrollableFrame(outer)
        left.pack(side=LEFT, fill=Y, padx=8, pady=8)

        right = Frame(outer)
        right.pack(side=RIGHT, fill=BOTH, expand=True, padx=8, pady=8)

        self._build_left(left.inner)
        self._build_right(right)
        self._build_bottom()

    def _build_left(self, parent: Frame) -> None:
        ref_box = LabelFrame(parent, text="PCR amplicon reference sequence")
        ref_box.pack(fill=X, pady=(0, 8))
        ref_row = Frame(ref_box)
        ref_row.pack(fill=BOTH, expand=True)

        self.reference_text = Text(ref_row, height=14, width=45, wrap="word")
        ref_scroll = ttk.Scrollbar(ref_row, orient=VERTICAL, command=self.reference_text.yview)
        self.reference_text.configure(yscrollcommand=ref_scroll.set)
        self.reference_text.pack(side=LEFT, fill=BOTH, expand=True)
        ref_scroll.pack(side=RIGHT, fill=Y)

        ref_buttons = Frame(ref_box)
        ref_buttons.pack(fill=X, pady=(6, 0))
        Button(ref_buttons, text="Paste from clipboard", command=self.paste_reference_from_clipboard).pack(side=LEFT)
        Button(ref_buttons, text="Clear reference", command=self.clear_reference).pack(side=LEFT, padx=6)

        guide_box = LabelFrame(parent, text="Analysis inputs")
        guide_box.pack(fill=X, pady=6)

        row1 = Frame(guide_box)
        row1.pack(fill=X, pady=4)
        Label(row1, text="sgRNA sequence").pack(side=LEFT)
        Entry(row1, textvariable=self.guide_var, width=32).pack(side=RIGHT, fill=X, expand=True)

        row2 = Frame(guide_box)
        row2.pack(fill=BOTH, expand=True, pady=4)
        Label(row2, text="Full CDS sequence").pack(side=LEFT)
        self.cds_text = Text(row2, height=5, width=32, wrap="word")
        self.cds_text.pack(side=RIGHT, fill=BOTH, expand=True)

        row3 = Frame(guide_box)
        row3.pack(fill=X, pady=4)
        Label(row3, text="Sequence bp before / after sgRNA").pack(side=LEFT)
        sequence_window = Frame(row3)
        sequence_window.pack(side=RIGHT)
        Entry(sequence_window, textvariable=self.sequence_before_var, width=6).pack(side=LEFT)
        Label(sequence_window, text=" / ").pack(side=LEFT)
        Entry(sequence_window, textvariable=self.sequence_after_var, width=6).pack(side=LEFT)

        row4 = Frame(guide_box)
        row4.pack(fill=X, pady=4)
        Label(row4, text="Protein aa before / after sgRNA").pack(side=LEFT)
        protein_window = Frame(row4)
        protein_window.pack(side=RIGHT)
        Entry(protein_window, textvariable=self.protein_before_var, width=6).pack(side=LEFT)
        Label(protein_window, text=" / ").pack(side=LEFT)
        Entry(protein_window, textvariable=self.protein_after_var, width=6).pack(side=LEFT)

        row5 = Frame(guide_box)
        row5.pack(fill=X, pady=4)
        Label(row5, text="Reference-mutant row gap (0.4-2.0)").pack(side=LEFT)
        Entry(row5, textvariable=self.sequence_row_spacing_var, width=10).pack(side=RIGHT)

        row6 = Frame(guide_box)
        row6.pack(fill=X, pady=4)
        Label(row6, text="60-base block gap (0.05-2.0)").pack(side=LEFT)
        Entry(row6, textvariable=self.dna_block_gap_var, width=10).pack(side=RIGHT)

        row7 = Frame(guide_box)
        row7.pack(fill=X, pady=4)
        Label(row7, text="DNA-to-protein panel gap (0.0-2.0)").pack(side=LEFT)
        Entry(row7, textvariable=self.dna_protein_gap_var, width=10).pack(side=RIGHT)

        row8 = Frame(guide_box)
        row8.pack(fill=X, pady=4)
        Label(row8, text="sgRNA line-to-DNA gap (0.05-2.0)").pack(side=LEFT)
        Entry(row8, textvariable=self.sgrna_line_gap_var, width=10).pack(side=RIGHT)

        row9 = Frame(guide_box)
        row9.pack(fill=X, pady=4)
        Label(row9, text="DNA font size (4-14)").pack(side=LEFT)
        Entry(row9, textvariable=self.dna_font_size_var, width=10).pack(side=RIGHT)

        row10 = Frame(guide_box)
        row10.pack(fill=X, pady=4)
        Label(row10, text="Protein font size (4-14)").pack(side=LEFT)
        Entry(row10, textvariable=self.protein_font_size_var, width=10).pack(side=RIGHT)

        row11 = Frame(guide_box)
        row11.pack(fill=X, pady=4)
        Label(row11, text="Structure label font size (4-18)").pack(side=LEFT)
        Entry(row11, textvariable=self.structure_label_font_size_var, width=10).pack(side=RIGHT)

        labels_box = LabelFrame(guide_box, text="Sequence labels")
        labels_box.pack(fill=X, pady=4)
        Label(
            labels_box,
            text="First line: reference. Following lines: mutants in input file order.",
            foreground="#555555",
            wraplength=390,
            justify=LEFT,
        ).pack(anchor=W, padx=4, pady=(3, 0))
        self.labels_text = Text(labels_box, height=4, width=32, wrap="word")
        self.labels_text.pack(fill=X, padx=4, pady=(2, 4))
        self.labels_text.insert("1.0", "Reference\nMutant")

        row12 = Frame(guide_box)
        row12.pack(fill=X, pady=4)
        Label(row12, text="Output folder").pack(side=LEFT)
        Entry(row12, textvariable=self.output_var).pack(side=RIGHT, fill=X, expand=True)

        row13 = Frame(guide_box)
        row13.pack(fill=X, pady=4)
        Label(row13, text="Structure type").pack(side=LEFT)
        structure_choice = ttk.Combobox(
            row13,
            textvariable=self.structure_mode_var,
            values=("Gene structure", "Protein domain"),
            state="readonly",
            width=20,
        )
        structure_choice.pack(side=RIGHT)

        row14 = Frame(guide_box)
        row14.pack(fill=X, pady=4)
        Label(row14, text="Structure total length").pack(side=LEFT)
        Entry(row14, textvariable=self.structure_total_var, width=10).pack(side=RIGHT)

        row15 = Frame(guide_box)
        row15.pack(fill=X, pady=4)
        Label(row15, text="sgRNA start position").pack(side=LEFT)
        Entry(row15, textvariable=self.structure_sgrna_start_var, width=10).pack(side=RIGHT)

        Label(
            guide_box,
            text=(
                "Structure total length: full protein aa or full CDS bp. "
                "sgRNA start position: the first sgRNA base in that same coordinate system "
                "(protein domain = aa; gene structure = CDS bp)."
            ),
            foreground="#555555",
            wraplength=390,
            justify=LEFT,
        ).pack(anchor=W, pady=(0, 4))

        region_box = LabelFrame(guide_box, text="Structure regions")
        region_box.pack(fill=X, pady=(4, 4))
        header = Frame(region_box)
        header.pack(fill=X, padx=2, pady=(2, 1))
        Label(header, text="Name", width=16, anchor=W).pack(side=LEFT)
        Label(header, text="Start", width=8, anchor=W).pack(side=LEFT)
        Label(header, text="End", width=8, anchor=W).pack(side=LEFT)
        Label(header, text="Color", width=10, anchor=W).pack(side=LEFT)
        self.structure_rows_frame = Frame(region_box)
        self.structure_rows_frame.pack(fill=X, padx=2, pady=(0, 2))
        region_buttons = Frame(region_box)
        region_buttons.pack(fill=X, padx=2, pady=(0, 4))
        Button(region_buttons, text="Add region", command=self.add_structure_row).pack(side=LEFT)
        Button(region_buttons, text="Clear regions", command=self.clear_structure_rows).pack(side=LEFT, padx=6)
        self.add_structure_row("Exon1", "1", "", "#1f77b4")

        Label(
            guide_box,
            text="The sgRNA first base is the display anchor. The before and after values add up to the displayed DNA length; sequences wrap at 60 bases per row.",
            foreground="#555555",
            wraplength=390,
            justify=LEFT,
        ).pack(anchor=W, pady=(0, 4))

        buttons = Frame(parent)
        buttons.pack(fill=X, pady=6)
        Button(buttons, text="Add AB1/FASTA files", command=self.add_files).pack(fill=X, pady=2)
        Button(buttons, text="Remove selected", command=self.remove_selected).pack(fill=X, pady=2)
        Button(buttons, text="Clear file list", command=self.clear_files).pack(fill=X, pady=2)
        Button(buttons, text="Run analysis", command=self.run_analysis).pack(fill=X, pady=6)

        file_box = LabelFrame(parent, text="AB1 or FASTA input files")
        file_box.pack(fill=BOTH, expand=True, pady=8)
        self.file_listbox = Listbox(file_box, height=10, selectmode="extended")
        self.file_listbox.pack(side=LEFT, fill=BOTH, expand=True)
        file_scroll = ttk.Scrollbar(file_box, orient=VERTICAL, command=self.file_listbox.yview)
        file_scroll.pack(side=RIGHT, fill=Y)
        self.file_listbox.configure(yscrollcommand=file_scroll.set)

        tip = Label(
            parent,
            text="Tip: paste the PCR product sequence (between your sequencing primers) into the reference box. Ctrl+V also works.",
            wraplength=390,
            justify=LEFT,
            foreground="#555555",
        )
        tip.pack(fill=X, pady=(6, 0))

    def _build_right(self, parent: Frame) -> None:
        top = Frame(parent)
        top.pack(fill=BOTH, expand=True)

        table_frame = LabelFrame(top, text="Analysis results")
        table_frame.pack(side=TOP, fill=BOTH, expand=True)

        columns = ("file", "orientation", "variants", "guide", "effect")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=11)
        headings = {
            "file": "Sample",
            "orientation": "Orientation",
            "variants": "Variants",
            "guide": "Guide",
            "effect": "Protein effect",
        }
        widths = {
            "file": 280,
            "orientation": 130,
            "variants": 200,
            "guide": 220,
            "effect": 420,
        }
        for col in columns:
            self.tree.heading(col, text=headings[col])
            self.tree.column(col, width=widths[col], anchor=W)
        self.tree.pack(fill=BOTH, expand=True)
        self.tree.bind("<<TreeviewSelect>>", self.on_select_result)

        preview_frame = LabelFrame(top, text="Figure preview")
        preview_frame.pack(side=BOTTOM, fill=BOTH, expand=True, pady=(8, 0))
        self.preview_label = Label(preview_frame, anchor="center")
        self.preview_label.pack(fill=BOTH, expand=True)

    def _build_bottom(self) -> None:
        status = Frame(self.root)
        status.pack(fill=X, side=BOTTOM)
        Label(status, textvariable=self.status_var, anchor=W).pack(fill=X)

    def set_status(self, text: str) -> None:
        self.status_var.set(text)
        self.root.update_idletasks()

    def paste_reference_from_clipboard(self) -> None:
        try:
            content = self.root.clipboard_get()
        except Exception as exc:
            messagebox.showerror("Clipboard error", f"Could not read clipboard content:\n{exc}")
            return
        if not content.strip():
            messagebox.showwarning("Clipboard empty", "Clipboard does not contain any text.")
            return
        self.reference_text.delete("1.0", END)
        self.reference_text.insert("1.0", content)
        self.set_status("Pasted PCR amplicon reference sequence from clipboard.")

    def clear_reference(self) -> None:
        self.reference_text.delete("1.0", END)
        self.set_status("Cleared PCR amplicon reference sequence.")

    def add_files(self) -> None:
        paths = filedialog.askopenfilenames(
            title="Select one or more AB1 or FASTA files",
            filetypes=[
                ("AB1 or FASTA files", "*.ab1 *.fa *.fasta *.fas *.fna *.ffn *.txt"),
                ("AB1 files", "*.ab1"),
                ("FASTA files", "*.fa *.fasta *.fas *.fna *.ffn *.txt"),
                ("All files", "*.*"),
            ],
        )
        if not paths:
            return
        for item in paths:
            path = Path(item)
            if path not in self.file_paths:
                self.file_paths.append(path)
                self.file_listbox.insert(END, str(path))
        self.set_status(f"Added {len(paths)} file(s).")

    def remove_selected(self) -> None:
        selected = list(self.file_listbox.curselection())
        for idx in reversed(selected):
            self.file_paths.pop(idx)
            self.file_listbox.delete(idx)
        self.set_status("Removed selected files.")

    def clear_files(self) -> None:
        self.file_paths.clear()
        self.file_listbox.delete(0, END)
        self.set_status("Cleared file list.")

    @staticmethod
    def _parse_nonnegative(raw: str, label: str) -> int:
        value = int(raw.strip())
        if value < 0:
            raise ValueError(f"{label} must be zero or greater.")
        return value

    @staticmethod
    def _parse_positive(raw: str, label: str) -> int:
        value = int(raw.strip())
        if value <= 0:
            raise ValueError(f"{label} must be greater than zero.")
        return value

    @staticmethod
    def _parse_decimal_range(raw: str, label: str, minimum: float, maximum: float) -> float:
        value = float(raw.strip())
        if not minimum <= value <= maximum:
            raise ValueError(f"{label} must be between {minimum:g} and {maximum:g}.")
        return value

    @staticmethod
    def _default_structure_color(index: int) -> str:
        palette = ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]
        return palette[index % len(palette)]

    def add_structure_row(
        self,
        name: str = "",
        start: str = "",
        end: str = "",
        color: str | None = None,
    ) -> None:
        if self.structure_rows_frame is None:
            return
        row = Frame(self.structure_rows_frame)
        row.pack(fill=X, pady=2)
        name_var = StringVar(value=name)
        start_var = StringVar(value=start)
        end_var = StringVar(value=end)
        color_var = StringVar(value=color or self._default_structure_color(len(self.structure_rows)))

        Entry(row, textvariable=name_var, width=16).pack(side=LEFT, padx=(0, 4))
        Entry(row, textvariable=start_var, width=8).pack(side=LEFT, padx=(0, 4))
        Entry(row, textvariable=end_var, width=8).pack(side=LEFT, padx=(0, 4))
        Entry(row, textvariable=color_var, width=10).pack(side=LEFT, padx=(0, 4))
        color_preview = Label(row, width=2, relief="sunken")
        color_preview.pack(side=LEFT, padx=(0, 4))

        row_info: dict[str, object] = {
            "frame": row,
            "name_var": name_var,
            "start_var": start_var,
            "end_var": end_var,
            "color_var": color_var,
        }

        def choose_color() -> None:
            current = color_var.get().strip() or "#1f77b4"
            _rgb, hex_color = colorchooser.askcolor(color=current, title="Choose domain color")
            if hex_color:
                color_var.set(hex_color)

        def refresh_color_preview(*_args) -> None:
            raw_color = color_var.get().strip()
            if is_color_like(raw_color):
                color_preview.configure(background=to_hex(raw_color))
            else:
                color_preview.configure(background="#ffffff")

        color_var.trace_add("write", refresh_color_preview)
        refresh_color_preview()

        def remove_row() -> None:
            row.destroy()
            if row_info in self.structure_rows:
                self.structure_rows.remove(row_info)

        Button(row, text="Pick", command=choose_color).pack(side=LEFT, padx=(0, 4))
        Button(row, text="Remove", command=remove_row).pack(side=RIGHT)
        self.structure_rows.append(row_info)

    def clear_structure_rows(self) -> None:
        for row_info in self.structure_rows:
            frame = row_info.get("frame")
            if isinstance(frame, Frame):
                frame.destroy()
        self.structure_rows.clear()

    def _read_sequence_labels(self) -> tuple[str, list[str]]:
        raw = self.labels_text.get("1.0", END)
        labels = [line.strip() for line in raw.splitlines() if line.strip()]
        reference_label = labels[0] if labels else "Reference"
        return reference_label, labels[1:]

    def _read_cds_sequence(self) -> str:
        raw = self.cds_text.get("1.0", END)
        cds = clean_dna(raw)
        if not cds:
            raise ValueError("CDS sequence is empty.")
        return cds

    def _read_structure_config(self) -> tuple[str | None, int | None, list[tuple[str, int, int, str]], tuple[int, int] | None]:
        mode = self.structure_mode_var.get().strip().lower()
        structure_mode: str | None = None
        if mode.startswith("gene"):
            structure_mode = "gene"
        elif mode.startswith("protein"):
            structure_mode = "protein"

        total_raw = self.structure_total_var.get().strip()
        total_length: int | None = None
        if total_raw:
            total_length = self._parse_positive(total_raw, "Structure total length")

        sgrna_start_raw = self.structure_sgrna_start_var.get().strip()
        sgrna_span: tuple[int, int] | None = None
        if sgrna_start_raw:
            sgrna_start = self._parse_positive(sgrna_start_raw, "sgRNA start")
            sgrna_span = (sgrna_start, sgrna_start)

        regions: list[tuple[str, int, int, str]] = []
        for idx, row_info in enumerate(self.structure_rows, start=1):
            name_var = row_info["name_var"]
            start_var = row_info["start_var"]
            end_var = row_info["end_var"]
            color_var = row_info["color_var"]
            if not isinstance(name_var, StringVar) or not isinstance(start_var, StringVar) or not isinstance(end_var, StringVar) or not isinstance(color_var, StringVar):
                continue
            name = name_var.get().strip() or f"Region {idx}"
            start = self._parse_positive(start_var.get(), f"{name} start")
            end = self._parse_positive(end_var.get(), f"{name} end")
            if end < start:
                start, end = end, start
            color = color_var.get().strip() or self._default_structure_color(idx - 1)
            if not is_color_like(color):
                raise ValueError(
                    f"{name} color must be a named matplotlib color or a hex value such as #4f81bd."
                )
            color = to_hex(color)
            regions.append((name, start, end, color))

        if structure_mode and regions and total_length is None:
            raise ValueError("Structure total length is required when structure regions are entered.")
        if structure_mode and regions and sgrna_span is None:
            raise ValueError("sgRNA start position in the structure is required when structure regions are entered.")
        return structure_mode, total_length, regions, sgrna_span

    def _read_reference(self) -> str:
        raw = self.reference_text.get("1.0", END)
        reference = clean_dna(raw)
        if not reference:
            raise ValueError("PCR amplicon reference sequence is empty.")
        return reference

    def _clear_results(self) -> None:
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.results.clear()
        self.preview_label.configure(image="", text="")
        self.preview_image = None
        self.batch_publication_png = None

    def _load_preview(self, path: Path) -> None:
        image = Image.open(path)
        max_w, max_h = 960, 520
        ratio = min(max_w / image.width, max_h / image.height, 1.0)
        resized = image.resize((max(1, int(image.width * ratio)), max(1, int(image.height * ratio))))
        self.preview_image = ImageTk.PhotoImage(resized)
        self.preview_label.configure(image=self.preview_image)

    def _variant_summary(self, analysis: SampleAnalysis) -> str:
        if not analysis.variants:
            return "none"
        head = []
        for variant in analysis.variants[:3]:
            if variant.kind == "substitution":
                head.append(f"{variant.kind}@{variant.ref_start}")
            elif variant.kind == "insertion":
                head.append(f"ins@{variant.ref_start}")
            else:
                head.append(f"del@{variant.ref_start}-{variant.ref_end}")
        if len(analysis.variants) > 3:
            head.append(f"+{len(analysis.variants) - 3} more")
        return ", ".join(head)

    def _guide_summary(self, analysis: SampleAnalysis) -> str:
        if not analysis.guide_hits:
            return "not found"
        return "; ".join(f"{hit.start}-{hit.end}({hit.strand})" for hit in analysis.guide_hits[:3])

    def _save_batch_summary(self, output_root: Path) -> None:
        summary_path = output_root / "batch_summary.csv"
        with summary_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample_label", "file", "orientation", "variants", "guide_hits", "effect", "figure_png", "figure_pdf", "report_txt"])
            for result in self.results:
                writer.writerow(
                    [
                        result.sample_label,
                        result.input_path.name,
                        result.orientation,
                        self._variant_summary(result),
                        self._guide_summary(result),
                        result.effect_summary,
                        result.figure_png,
                        result.figure_pdf,
                        result.report_txt,
                    ]
                )

    def run_analysis(self) -> None:
        try:
            reference = self._read_reference()
            guide = self.guide_var.get().strip()
            if not guide:
                raise ValueError("sgRNA sequence is empty.")
            cds_sequence = self._read_cds_sequence()
            protein_flank_aa = (
                self._parse_nonnegative(self.protein_before_var.get(), "Protein aa before sgRNA"),
                self._parse_nonnegative(self.protein_after_var.get(), "Protein aa after sgRNA"),
            )
            sequence_flank_bp = (
                self._parse_nonnegative(self.sequence_before_var.get(), "Sequence bp before sgRNA"),
                self._parse_nonnegative(self.sequence_after_var.get(), "Sequence bp after sgRNA"),
            )
            sequence_row_spacing = self._parse_decimal_range(
                self.sequence_row_spacing_var.get(),
                "Sequence row spacing",
                0.4,
                2.0,
            )
            dna_block_gap = self._parse_decimal_range(
                self.dna_block_gap_var.get(),
                "60-base block gap",
                0.05,
                2.0,
            )
            dna_protein_gap = self._parse_decimal_range(
                self.dna_protein_gap_var.get(),
                "DNA-to-protein gap",
                0.0,
                2.0,
            )
            sgrna_line_gap = self._parse_decimal_range(
                self.sgrna_line_gap_var.get(),
                "sgRNA line-to-DNA gap",
                0.05,
                2.0,
            )
            dna_font_size = self._parse_decimal_range(
                self.dna_font_size_var.get(),
                "DNA font size",
                4.0,
                14.0,
            )
            protein_font_size = self._parse_decimal_range(
                self.protein_font_size_var.get(),
                "Protein font size",
                4.0,
                14.0,
            )
            structure_label_font_size = self._parse_decimal_range(
                self.structure_label_font_size_var.get(),
                "Structure label font size",
                4.0,
                18.0,
            )
            structure_mode, structure_total_length, structure_regions, structure_sgrna_span = self._read_structure_config()
            reference_label, mutant_labels = self._read_sequence_labels()
        except Exception as exc:
            messagebox.showerror("Input error", str(exc))
            return

        if not self.file_paths:
            messagebox.showerror("Input error", "Please add one or more AB1 or FASTA files first.")
            return

        output_root = Path(self.output_var.get().strip() or (Path(__file__).resolve().parent / "mutant_id_output"))
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_root = output_root / timestamp
        output_root.mkdir(parents=True, exist_ok=True)
        self.last_output_root = output_root
        self._clear_results()
        self.set_status("Running analysis...")

        def worker() -> None:
            failures: list[str] = []
            for idx, path in enumerate(self.file_paths, start=1):
                try:
                    sample_label = (
                        mutant_labels[idx - 1]
                        if idx <= len(mutant_labels)
                        else f"Mutant {idx}"
                    )
                    result = analyze_sample(
                        reference,
                        guide,
                        cds_sequence,
                        path,
                        output_root,
                        sample_label=sample_label,
                        protein_flank_aa=protein_flank_aa,
                        sequence_flank_bp=sequence_flank_bp,
                        structure_mode=structure_mode,
                        structure_total_length=structure_total_length,
                        structure_regions=structure_regions,
                        structure_sgrna_span=structure_sgrna_span,
                        sequence_row_spacing=sequence_row_spacing,
                        dna_block_gap=dna_block_gap,
                        dna_protein_gap=dna_protein_gap,
                        sgrna_line_gap=sgrna_line_gap,
                        dna_font_size=dna_font_size,
                        protein_font_size=protein_font_size,
                        structure_label_font_size=structure_label_font_size,
                        reference_label=reference_label,
                        mutant_label=sample_label,
                    )
                    self.results.append(result)
                    self.root.after(0, self._append_tree_row, result)
                    self.root.after(0, self.set_status, f"Finished {idx}/{len(self.file_paths)}: {path.name}")
                except Exception as exc:  # pragma: no cover - user-facing
                    failures.append(f"{path.name}: {exc}")
                    self.root.after(0, self.set_status, f"Failed {path.name}")
            if self.results:
                try:
                    png = output_root / "publication_mutation_summary.png"
                    pdf = output_root / "publication_mutation_summary.pdf"
                    render_batch_publication_figure(
                        png,
                        pdf,
                        reference,
                        clean_dna(guide),
                        self.results,
                        protein_flank_aa=protein_flank_aa,
                        structure_mode=structure_mode,
                        structure_total_length=structure_total_length,
                        structure_regions=structure_regions,
                        structure_sgrna_span=structure_sgrna_span,
                        sequence_row_spacing=sequence_row_spacing,
                        dna_block_gap=dna_block_gap,
                        dna_protein_gap=dna_protein_gap,
                        sgrna_line_gap=sgrna_line_gap,
                        dna_font_size=dna_font_size,
                        protein_font_size=protein_font_size,
                        structure_label_font_size=structure_label_font_size,
                        reference_label=reference_label,
                    )
                    self.batch_publication_png = png
                except Exception as exc:
                    failures.append(f"publication summary figure: {exc}")
            self.root.after(0, self._analysis_finished, output_root, failures)

        threading.Thread(target=worker, daemon=True).start()

    def _append_tree_row(self, result: SampleAnalysis) -> None:
        item_id = self.tree.insert(
            "",
            END,
            values=(
                result.sample_label,
                result.orientation,
                self._variant_summary(result),
                self._guide_summary(result),
                result.effect_summary,
            ),
        )
        if not self.tree.selection():
            self.tree.selection_set(item_id)
            self.tree.focus(item_id)
            self._load_preview(result.figure_png)
            self.preview_label.configure(text="")

    def _analysis_finished(self, output_root: Path, failures: list[str]) -> None:
        self._save_batch_summary(output_root)
        if self.batch_publication_png and self.batch_publication_png.exists():
            self._load_preview(self.batch_publication_png)
        elif self.results:
            self._load_preview(self.results[0].figure_png)
        if failures:
            messagebox.showwarning(
                "Analysis finished with warnings",
                "Some files failed:\n\n" + "\n".join(failures),
            )
        else:
            messagebox.showinfo(
                "Analysis finished",
                f"Results saved to:\n{output_root}",
            )
        self.set_status(f"Finished. Output saved to {output_root}")

    def on_select_result(self, _event=None) -> None:
        selection = self.tree.selection()
        if not selection:
            return
        index = self.tree.index(selection[0])
        if 0 <= index < len(self.results):
            result = self.results[index]
            self._load_preview(result.figure_png)
            self.set_status(result.report_txt.as_posix())

    def run(self) -> None:
        self.root.mainloop()


def main() -> None:
    app = MutantIDApp()
    app.run()


if __name__ == "__main__":
    main()
