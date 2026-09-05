from __future__ import annotations

import csv
import os
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".matplotlib-cache"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np

from Bio import BiopythonDeprecationWarning, BiopythonWarning

warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
warnings.filterwarnings("ignore", category=BiopythonWarning)

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


BASE_COLORS = {
    "A": "#2ca02c",
    "C": "#1f77b4",
    "G": "#111111",
    "T": "#d62728",
    "N": "#7f7f7f",
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
SEQUENCE_CHAR_STEP = 0.14
SEQUENCE_CHUNK_WIDTH = 60
SEQUENCE_LEFT_MARGIN = -1.8


@dataclass
class AB1Read:
    path: Path
    name: str
    sequence: str
    traces: dict[str, np.ndarray]
    peak_positions: np.ndarray
    channel_order: str


@dataclass
class AlignmentResult:
    reference: str
    query: str
    aligned_reference: str
    aligned_query: str
    score: float
    reverse_complement: bool
    reference_start: int
    query_start: int


@dataclass
class GuideHit:
    start: int
    end: int
    strand: str
    sequence: str


@dataclass
class VariantCall:
    kind: str
    ref_start: int
    ref_end: int
    ref_bases: str
    query_bases: str
    query_indices: list[int]


@dataclass
class SampleAnalysis:
    input_path: Path
    sample_label: str
    output_dir: Path
    orientation: str
    aligned_read: str
    aligned_reference: str
    reference: str
    query: str
    cds_sequence: str
    cds_start: int
    cds_end: int
    reverse_cds: bool
    guide_hits: list[GuideHit]
    variants: list[VariantCall]
    ref_protein: str
    mut_protein: str
    effect_summary: str
    figure_png: Path
    figure_pdf: Path
    report_txt: Path
    variants_csv: Path
    cds_map_start: int
    cds_map_end: int
    cds_map_reverse: bool
    unmapped_variants: int
    guide_cds_start: int | None
    guide_cds_end: int | None
    display_reference_start: int
    display_reference_end: int


def clean_dna(raw: str) -> str:
    lines = [line for line in raw.splitlines() if not line.lstrip().startswith(">")]
    sequence = re.sub(r"[^A-Za-z]", "", "\n".join(lines)).upper().replace("U", "T")
    return "".join(base for base in sequence if base in "ACGTN")


def reverse_complement(seq: str) -> str:
    return clean_dna(seq).translate(COMPLEMENT)[::-1]


def safe_stem(path: Path) -> str:
    stem = path.stem.strip()
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", stem)
    return stem or "sample"


def load_ab1(path: Path) -> AB1Read:
    record = SeqIO.read(str(path), "abi")
    sequence = clean_dna(str(record.seq))
    abif = record.annotations.get("abif_raw", {})

    channel_order = abif.get("FWO_1", b"ACGT")
    if isinstance(channel_order, bytes):
        channel_order = channel_order.decode("ascii", errors="ignore")
    channel_order = "".join(base for base in channel_order.upper() if base in "ACGT")
    if len(channel_order) != 4:
        channel_order = "ACGT"

    traces: dict[str, np.ndarray] = {}
    for idx, base in enumerate(channel_order):
        key = f"DATA{9 + idx}"
        traces[base] = np.asarray(abif.get(key, []), dtype=np.float64)

    if not traces:
        raise ValueError(f"Could not extract chromatogram traces from {path.name}.")

    peak_positions = abif.get("PLOC2")
    if peak_positions is None:
        peak_positions = abif.get("PLOC1")
    if peak_positions is None:
        longest = max((len(arr) for arr in traces.values()), default=len(sequence))
        if len(sequence) > 0:
            peak_positions = np.linspace(0, longest - 1, len(sequence), dtype=int)
        else:
            peak_positions = np.asarray([], dtype=int)
    else:
        peak_positions = np.asarray(peak_positions, dtype=int)

    return AB1Read(
        path=path,
        name=record.name or path.stem,
        sequence=sequence,
        traces=traces,
        peak_positions=peak_positions,
        channel_order=channel_order,
    )


def load_fasta(path: Path) -> AB1Read:
    records = list(SeqIO.parse(str(path), "fasta"))
    if records:
        record = records[0]
        sequence = clean_dna(str(record.seq))
        name = record.id or path.stem
    else:
        sequence = clean_dna(path.read_text(encoding="utf-8", errors="ignore"))
        name = path.stem
    if not sequence:
        raise ValueError(f"Could not read a DNA sequence from {path.name}.")
    return AB1Read(
        path=path,
        name=name,
        sequence=sequence,
        traces={},
        peak_positions=np.arange(1, len(sequence) + 1, dtype=int),
        channel_order="",
    )


def load_sequence_read(path: Path) -> AB1Read:
    suffix = path.suffix.lower()
    if suffix == ".ab1":
        return load_ab1(path)
    if suffix in {".fa", ".fasta", ".fas", ".fna", ".ffn", ".txt"}:
        return load_fasta(path)
    try:
        return load_ab1(path)
    except Exception as ab1_error:
        try:
            return load_fasta(path)
        except Exception as fasta_error:
            raise ValueError(
                f"Could not read {path.name} as AB1 or FASTA. "
                f"AB1 error: {ab1_error}; FASTA error: {fasta_error}"
            ) from fasta_error


def orient_read(read: AB1Read, reverse: bool) -> AB1Read:
    if not reverse:
        return read

    seq = reverse_complement(read.sequence)
    traces = {}
    if read.traces:
        for base in "ACGT":
            comp = {"A": "T", "T": "A", "C": "G", "G": "C"}[base]
            traces[base] = np.asarray(read.traces.get(comp, np.asarray([]))[::-1], dtype=np.float64)

    trace_length = max((len(arr) for arr in traces.values()), default=0)
    if trace_length:
        peaks = np.asarray([trace_length - 1 - int(pos) for pos in read.peak_positions], dtype=int)[::-1]
    elif len(read.peak_positions) == len(read.sequence):
        peaks = np.arange(1, len(seq) + 1, dtype=int)
    else:
        peaks = read.peak_positions[::-1]

    return AB1Read(
        path=read.path,
        name=read.name,
        sequence=seq,
        traces=traces,
        peak_positions=peaks,
        channel_order="ACGT" if traces else "",
    )


def align_sequences(reference: str, query: str) -> tuple[str, str, float, int, int]:
    alignment = pairwise2.align.localms(reference, query, 2, -1, -5, -1, one_alignment_only=True)
    if not alignment:
        raise ValueError("Unable to align the sequencing read to the reference sequence.")
    best = alignment[0]
    aligned_reference = best.seqA[best.start : best.end]
    aligned_query = best.seqB[best.start : best.end]
    reference_start = sum(base != "-" for base in best.seqA[: best.start]) + 1
    query_start = sum(base != "-" for base in best.seqB[: best.start]) + 1
    return aligned_reference, aligned_query, float(best.score), reference_start, query_start


def choose_orientation(reference: str, read: AB1Read) -> tuple[AB1Read, AlignmentResult]:
    forward_query = read.sequence
    reverse_read = orient_read(read, True)

    forward_ref, forward_qry, forward_score, forward_ref_start, forward_query_start = align_sequences(reference, forward_query)
    reverse_ref, reverse_qry, reverse_score, reverse_ref_start, reverse_query_start = align_sequences(reference, reverse_read.sequence)

    if reverse_score > forward_score:
        oriented = reverse_read
        result = AlignmentResult(
            reference=reference,
            query=oriented.sequence,
            aligned_reference=reverse_ref,
            aligned_query=reverse_qry,
            score=reverse_score,
            reverse_complement=True,
            reference_start=reverse_ref_start,
            query_start=reverse_query_start,
        )
    else:
        oriented = read
        result = AlignmentResult(
            reference=reference,
            query=forward_query,
            aligned_reference=forward_ref,
            aligned_query=forward_qry,
            score=forward_score,
            reverse_complement=False,
            reference_start=forward_ref_start,
            query_start=forward_query_start,
        )
    return oriented, result


def find_guide_hits(reference: str, guide: str) -> list[GuideHit]:
    guide = clean_dna(guide)
    if not guide:
        return []

    hits: list[GuideHit] = []
    candidates = [(guide, "+"), (reverse_complement(guide), "-")]
    for pattern, strand in candidates:
        start = reference.find(pattern)
        while start != -1:
            hits.append(GuideHit(start=start + 1, end=start + len(pattern), strand=strand, sequence=pattern))
            start = reference.find(pattern, start + 1)

    hits.sort(key=lambda item: (item.start, item.end, item.strand))
    return hits


def locate_cds(reference: str, cds: str) -> tuple[int, int, bool]:
    """Return 1-based inclusive CDS coordinates and whether the CDS is reverse-oriented."""
    reference = clean_dna(reference)
    cds = clean_dna(cds)
    if not cds:
        raise ValueError("CDS sequence is empty.")

    forward_start = reference.find(cds)
    if forward_start >= 0:
        return forward_start + 1, forward_start + len(cds), False

    reverse = reverse_complement(cds)
    reverse_start = reference.find(reverse)
    if reverse_start >= 0:
        return reverse_start + 1, reverse_start + len(reverse), True

    raise ValueError(
        "The CDS sequence could not be found in the PCR reference sequence. "
        "This helper is kept for legacy workflows; the GUI maps PCR reference and full CDS by local overlap."
    )


@dataclass
class CDSMapping:
    reference_length: int
    cds_length: int
    score: float
    # Coordinates below are in the PCR-product reference, not in the full CDS.
    reference_start: int
    reference_end: int
    reverse_reference: bool
    reference_to_cds: dict[int, int]

    @property
    def cds_start(self) -> int:
        return min(self.reference_to_cds.values())

    @property
    def cds_end(self) -> int:
        return max(self.reference_to_cds.values())


def guide_cds_span(cds: str, guide: str) -> tuple[int, int]:
    """Locate the sgRNA directly in the full coding sequence."""
    hits = find_guide_hits(cds, guide)
    if not hits:
        raise ValueError(
            "The sgRNA sequence was not found in the full CDS, in either orientation. "
            "Confirm that the full CDS and sgRNA use the same sequence version."
        )
    unique_spans = {(hit.start, hit.end) for hit in hits}
    if len(unique_spans) != 1:
        raise ValueError(
            "The sgRNA sequence matches more than one position in the full CDS. "
            "Use an unambiguous sgRNA sequence or verify the supplied CDS."
        )
    return next(iter(unique_spans))


def protein_window_to_cds_span(
    cds_length: int,
    guide_cds_start: int,
    guide_cds_end: int,
    before_aa: int,
    after_aa: int,
) -> tuple[int, int]:
    guide_aa_start = max(1, ((guide_cds_start - 1) // 3) + 1)
    aa_start, aa_end = _window_from_anchor(
        (cds_length + 2) // 3,
        guide_aa_start,
        before_aa,
        after_aa,
    )
    return (aa_start - 1) * 3 + 1, min(cds_length, aa_end * 3)


def _window_from_anchor(
    sequence_length: int,
    anchor: int,
    before: int,
    after: int,
) -> tuple[int, int]:
    """Build a fixed-width display window with the anchor as the after-side start."""
    if sequence_length <= 0:
        return 1, 0
    width = min(sequence_length, max(1, max(0, before) + max(0, after)))
    anchor = max(1, min(sequence_length, anchor))
    start = max(1, anchor - max(0, before))
    end = min(sequence_length, start + width - 1)
    if end - start + 1 < width:
        start = max(1, end - width + 1)
    return start, end


def reference_window_from_guide(
    reference_length: int,
    guide_hit: GuideHit,
    before_bp: int,
    after_bp: int,
) -> tuple[int, int]:
    return _window_from_anchor(
        reference_length,
        guide_hit.start,
        before_bp,
        after_bp,
    )


def variants_in_reference_window(
    variants: list[VariantCall],
    window_start: int,
    window_end: int,
) -> list[VariantCall]:
    selected: list[VariantCall] = []
    for variant in variants:
        if variant.kind == "insertion":
            anchor = variant.ref_start - 1
            if window_start - 1 <= anchor <= window_end:
                selected.append(variant)
        elif variant.ref_start <= window_end and variant.ref_end >= window_start:
            selected.append(variant)
    return selected


def reference_positions_for_cds_span(
    mapping: CDSMapping,
    cds_start: int,
    cds_end: int,
) -> list[int]:
    return sorted(
        ref_pos
        for ref_pos, cds_pos in mapping.reference_to_cds.items()
        if cds_start <= cds_pos <= cds_end
    )


def _local_reference_cds_alignment(reference: str, cds: str) -> tuple[str, str, float, dict[int, int]]:
    alignment = pairwise2.align.localms(reference, cds, 2, -1, -5, -1, one_alignment_only=True)
    if not alignment:
        raise ValueError("Unable to map the PCR-product reference sequence to the CDS.")
    best = alignment[0]
    ref_to_cds: dict[int, int] = {}
    ref_pos = 0
    cds_pos = 0
    exact_matches = 0
    for column, (ref_base, cds_base) in enumerate(zip(best.seqA, best.seqB)):
        if ref_base != "-":
            ref_pos += 1
        if cds_base != "-":
            cds_pos += 1
        if best.start <= column < best.end and ref_base != "-" and cds_base != "-":
            ref_to_cds[ref_pos] = cds_pos
            if ref_base == cds_base:
                exact_matches += 1
    mapped = sorted(ref_to_cds)
    minimum_overlap = min(20, max(8, min(len(reference), len(cds)) // 3))
    if len(mapped) < minimum_overlap or exact_matches < minimum_overlap:
        raise ValueError(
            "The PCR-product reference and full CDS do not share enough matching sequence "
            "for reliable coordinate mapping. Confirm that both inputs come from the same transcript."
        )
    return best.seqA, best.seqB, float(best.score), ref_to_cds


def map_reference_to_cds(reference: str, cds: str) -> CDSMapping:
    candidates = []
    for reverse in (False, True):
        oriented_reference = reverse_complement(reference) if reverse else reference
        seq_a, seq_b, score, oriented_map = _local_reference_cds_alignment(oriented_reference, cds)
        mapped_positions = sorted(oriented_map)
        if reverse:
            original_map = {
                len(reference) - oriented_ref_pos + 1: cds_pos
                for oriented_ref_pos, cds_pos in oriented_map.items()
            }
        else:
            original_map = dict(oriented_map)
        candidates.append(
            CDSMapping(
                reference_length=len(reference),
                cds_length=len(cds),
                score=score,
                reference_start=min(original_map),
                reference_end=max(original_map),
                reverse_reference=reverse,
                reference_to_cds=original_map,
            )
        )
    return max(candidates, key=lambda item: (item.score, len(item.reference_to_cds)))


def apply_variants_to_cds(
    cds: str,
    reference: str,
    variants: list[VariantCall],
    mapping: CDSMapping,
) -> tuple[str, list[str], int | None]:
    """Apply amplicon-reference variants to the full CDS in CDS coordinates."""
    sequence = cds
    operations: list[tuple[int, int, str, str]] = []
    unmapped: list[str] = []
    affected_positions: list[int] = []

    def mapped_position(ref_pos: int) -> int | None:
        return mapping.reference_to_cds.get(ref_pos)

    for variant in variants:
        if variant.kind == "substitution":
            cds_pos = mapped_position(variant.ref_start)
            if cds_pos is None:
                unmapped.append(f"substitution@{variant.ref_start}")
                continue
            base = reverse_complement(variant.query_bases) if mapping.reverse_reference else variant.query_bases
            operations.append((cds_pos, cds_pos, "substitution", base))
            affected_positions.append(cds_pos)
            continue

        if variant.kind == "deletion":
            cds_positions = [
                mapped_position(ref_pos)
                for ref_pos in range(variant.ref_start, variant.ref_end + 1)
            ]
            if any(pos is None for pos in cds_positions):
                unmapped.append(f"deletion@{variant.ref_start}-{variant.ref_end}")
                continue
            first = min(pos for pos in cds_positions if pos is not None)
            last = max(pos for pos in cds_positions if pos is not None)
            operations.append((first, last, "deletion", ""))
            affected_positions.extend([first, last])
            continue

        if variant.kind == "insertion":
            anchor = variant.ref_start - 1
            if mapping.reverse_reference:
                next_cds = mapped_position(variant.ref_start)
                anchor_cds = mapped_position(anchor) if anchor > 0 else None
                if anchor_cds is not None:
                    # PCR reference runs opposite to the CDS. An insertion after
                    # this PCR base belongs immediately before its CDS counterpart.
                    insert_index = anchor_cds - 1
                    affected_position = anchor_cds
                elif next_cds is not None:
                    insert_index = next_cds
                    affected_position = next_cds
                else:
                    unmapped.append(f"insertion@{variant.ref_start}")
                    continue
                inserted = reverse_complement(variant.query_bases)
            else:
                anchor_cds = mapped_position(anchor) if anchor > 0 else None
                next_cds = mapped_position(variant.ref_start)
                if anchor_cds is not None:
                    insert_index = anchor_cds
                    affected_position = anchor_cds
                elif next_cds is not None:
                    insert_index = max(0, next_cds - 1)
                    affected_position = next_cds
                else:
                    unmapped.append(f"insertion@{variant.ref_start}")
                    continue
                inserted = variant.query_bases
            operations.append((insert_index, insert_index, "insertion", inserted))
            affected_positions.append(affected_position)
            continue

        unmapped.append(f"{variant.kind}@{variant.ref_start}")

    for start, end, kind, payload in sorted(operations, key=lambda item: (item[0], item[1]), reverse=True):
        if kind == "substitution":
            index = start - 1
            if 0 <= index < len(sequence):
                sequence = sequence[:index] + payload + sequence[index + 1 :]
        elif kind == "deletion":
            sequence = sequence[: start - 1] + sequence[end:]
        elif kind == "insertion":
            index = max(0, min(len(sequence), start))
            sequence = sequence[:index] + payload + sequence[index:]

    return sequence, unmapped, min(affected_positions) if affected_positions else None


def call_variants(
    aligned_reference: str,
    aligned_query: str,
    reference_start: int = 1,
    query_start: int = 1,
) -> list[VariantCall]:
    variants: list[VariantCall] = []
    ref_pos = reference_start - 1
    query_pos = query_start - 1
    i = 0
    n = len(aligned_reference)
    covered_columns = [idx for idx, base in enumerate(aligned_query) if base != "-"]
    first_query_col = covered_columns[0] if covered_columns else 0
    last_query_col = covered_columns[-1] if covered_columns else n - 1

    while i < n:
        ref_base = aligned_reference[i]
        query_base = aligned_query[i]

        if ref_base == query_base:
            if ref_base != "-":
                ref_pos += 1
            if query_base != "-":
                query_pos += 1
            i += 1
            continue

        if ref_base == "-" and query_base != "-":
            anchor = ref_pos
            inserted: list[str] = []
            indices: list[int] = []
            while i < n and aligned_reference[i] == "-" and aligned_query[i] != "-":
                query_pos += 1
                inserted.append(aligned_query[i])
                indices.append(query_pos)
                i += 1
            variants.append(
                VariantCall(
                    kind="insertion",
                    ref_start=anchor + 1,
                    ref_end=anchor,
                    ref_bases="",
                    query_bases="".join(inserted),
                    query_indices=indices,
                )
            )
            continue

        if ref_base != "-" and query_base == "-":
            if i < first_query_col or i > last_query_col:
                while i < n and aligned_reference[i] != "-" and aligned_query[i] == "-":
                    ref_pos += 1
                    i += 1
                continue
            start = ref_pos + 1
            deleted: list[str] = []
            ref_end = start
            while i < n and aligned_reference[i] != "-" and aligned_query[i] == "-":
                ref_pos += 1
                deleted.append(aligned_reference[i])
                ref_end = ref_pos
                i += 1
            variants.append(
                VariantCall(
                    kind="deletion",
                    ref_start=start,
                    ref_end=ref_end,
                    ref_bases="".join(deleted),
                    query_bases="",
                    query_indices=[],
                )
            )
            continue

        if ref_base != "-" and query_base != "-":
            ref_pos += 1
            query_pos += 1
            variants.append(
                VariantCall(
                    kind="substitution",
                    ref_start=ref_pos,
                    ref_end=ref_pos,
                    ref_bases=ref_base,
                    query_bases=query_base,
                    query_indices=[query_pos],
                )
            )
            i += 1
            continue

        i += 1

    return variants


def build_ref_to_query_map(
    aligned_reference: str,
    aligned_query: str,
    reference_start: int = 1,
    query_start: int = 1,
) -> tuple[dict[int, list[int]], dict[int, int]]:
    ref_to_query: dict[int, list[int]] = {}
    query_to_ref: dict[int, int] = {}
    ref_pos = reference_start - 1
    query_pos = query_start - 1
    for ref_base, query_base in zip(aligned_reference, aligned_query):
        if ref_base != "-":
            ref_pos += 1
        if query_base != "-":
            query_pos += 1
        if query_base != "-":
            query_to_ref[query_pos] = ref_pos
        if ref_base != "-" and query_base != "-":
            ref_to_query.setdefault(ref_pos, []).append(query_pos)
    return ref_to_query, query_to_ref


def extract_coding_sequence_from_alignment(
    aligned_reference: str,
    aligned_query: str,
    cds_start: int,
    cds_end: int,
    reverse_cds: bool = False,
) -> tuple[str, str]:
    ref_pos = 0
    ref_coding: list[str] = []
    query_coding: list[str] = []

    for ref_base, query_base in zip(aligned_reference, aligned_query):
        if ref_base == "-" and query_base != "-":
            if cds_start <= ref_pos < cds_end:
                query_coding.append(query_base)
            continue
        if ref_base != "-":
            ref_pos += 1
        if not (cds_start <= ref_pos <= cds_end):
            continue
        ref_coding.append(ref_base)
        if query_base != "-":
            query_coding.append(query_base)

    ref_sequence = "".join(ref_coding)
    query_sequence = "".join(query_coding)
    if reverse_cds:
        return reverse_complement(ref_sequence), reverse_complement(query_sequence)
    return ref_sequence, query_sequence


def translate_dna(sequence: str) -> str:
    sequence = clean_dna(sequence)
    if not sequence:
        return ""
    remainder = len(sequence) % 3
    if remainder:
        sequence = sequence[: len(sequence) - remainder]
    translated = str(Seq(sequence).translate(table=1, to_stop=False))
    stop_index = translated.find("*")
    if stop_index != -1:
        return translated[: stop_index + 1]
    return translated


def first_difference_position(ref_aa: str, mut_aa: str) -> int | None:
    limit = min(len(ref_aa), len(mut_aa))
    for idx in range(limit):
        if ref_aa[idx] != mut_aa[idx]:
            return idx + 1
    if len(ref_aa) != len(mut_aa):
        return limit + 1
    return None


def summarize_effect(
    variants: list[VariantCall],
    ref_aa: str,
    mut_aa: str,
    first_cds_position: int | None,
) -> str:
    if not variants:
        return "No sequence difference detected relative to the supplied reference."

    indels = [item for item in variants if item.kind in {"insertion", "deletion"}]
    frameshift = any((len(item.query_bases) - len(item.ref_bases)) % 3 != 0 for item in indels if item.kind != "substitution")

    if first_cds_position is None:
        return "Sequence difference detected outside the mapped CDS region."
    aa_pos = max(1, ((first_cds_position - 1) // 3) + 1)

    if frameshift:
        stop_pos = mut_aa.find("*")
        if stop_pos != -1:
            return f"Frameshift likely begins near aa {aa_pos}; mutant protein stops at aa {stop_pos + 1}."
        return f"Frameshift likely begins near aa {aa_pos}; no stop codon detected in the translated window."

    if indels:
        return f"In-frame indel detected near aa {aa_pos}."

    diff = first_difference_position(ref_aa, mut_aa)
    if diff is not None:
        return f"Protein sequence differs first at aa {diff}."
    return "Sequence difference detected, but protein-level change could not be resolved from the current window."


def format_pcr_variant_notation(variant: VariantCall) -> str:
    prefix = (
        f"PCR.{variant.ref_start}"
        if variant.ref_start == variant.ref_end
        else f"PCR.{variant.ref_start}_{variant.ref_end}"
    )
    if variant.kind == "substitution":
        return f"{prefix}{variant.ref_bases}>{variant.query_bases}"
    if variant.kind == "deletion":
        return f"{prefix}del{variant.ref_bases}"
    if variant.kind == "insertion":
        return f"{prefix}ins{variant.query_bases}"
    return f"{prefix}{variant.kind}"


def format_variant_notation(variant: VariantCall, mapping: CDSMapping) -> str:
    """Format a PCR-reference variant in full-CDS cDNA coordinates."""
    def cds_position(ref_pos: int) -> int | None:
        return mapping.reference_to_cds.get(ref_pos)

    def orient_bases(bases: str) -> str:
        return reverse_complement(bases) if mapping.reverse_reference else bases

    if variant.kind == "substitution":
        position = cds_position(variant.ref_start)
        if position is None:
            return f"{format_pcr_variant_notation(variant)} (outside mapped CDS)"
        return f"c.{position}{orient_bases(variant.ref_bases)}>{orient_bases(variant.query_bases)}"

    if variant.kind == "deletion":
        positions = [
            cds_position(ref_pos)
            for ref_pos in range(variant.ref_start, variant.ref_end + 1)
        ]
        if any(position is None for position in positions):
            return f"{format_pcr_variant_notation(variant)} (outside mapped CDS)"
        start, end = min(positions), max(positions)
        prefix = f"c.{start}" if start == end else f"c.{start}_{end}"
        return f"{prefix}del{orient_bases(variant.ref_bases)}"

    if variant.kind == "insertion":
        anchor = variant.ref_start - 1
        anchor_cds = cds_position(anchor) if anchor > 0 else None
        next_cds = cds_position(variant.ref_start)
        if mapping.reverse_reference:
            left, right = next_cds, anchor_cds
        else:
            left, right = anchor_cds, next_cds

        if left is not None and right is not None:
            start, end = sorted((left, right))
        elif left is not None:
            start, end = left, left + 1
        elif right is not None:
            start, end = max(1, right - 1), right
        else:
            return f"{format_pcr_variant_notation(variant)} (outside mapped CDS)"
        return f"c.{start}_{end}ins{orient_bases(variant.query_bases)}"

    return format_pcr_variant_notation(variant)


def _peak_positions_for_query(read: AB1Read) -> np.ndarray:
    peaks = np.asarray(read.peak_positions, dtype=int)
    if len(peaks) == len(read.sequence):
        return peaks
    if len(read.sequence) == 0:
        return np.asarray([], dtype=int)
    trace_length = max((len(arr) for arr in read.traces.values()), default=len(read.sequence))
    return np.linspace(0, max(0, trace_length - 1), len(read.sequence), dtype=int)


def build_variant_windows(
    variants: list[VariantCall],
    ref_to_query: dict[int, list[int]],
    query_to_ref: dict[int, int],
    window_start: int,
    window_end: int,
) -> list[int]:
    selected: set[int] = set()
    for ref_pos in range(window_start, window_end + 1):
        for query_index in ref_to_query.get(ref_pos, []):
            selected.add(query_index)
    for query_index, ref_pos in query_to_ref.items():
        if window_start <= ref_pos <= window_end:
            selected.add(query_index)
    for variant in variants:
        if variant.kind == "insertion":
            selected.update(variant.query_indices)
        elif variant.kind == "substitution":
            selected.update(variant.query_indices)
        elif variant.kind == "deletion":
            for ref_pos in range(max(1, variant.ref_start - 1), variant.ref_end + 2):
                for query_index in ref_to_query.get(ref_pos, []):
                    selected.add(query_index)
    return sorted(idx for idx in selected if idx > 0)


def select_trace_bounds(
    read: AB1Read,
    query_indices: list[int],
    margin: int = 40,
) -> tuple[int, int]:
    peaks = _peak_positions_for_query(read)
    if len(query_indices) == 0 or len(peaks) == 0:
        return 0, int(peaks.max() if len(peaks) else 0)

    positions = [int(peaks[idx - 1]) for idx in query_indices if 1 <= idx <= len(peaks)]
    if not positions:
        return 0, int(peaks.max())
    left = max(0, min(positions) - margin)
    right = min(int(peaks.max()) + 1, max(positions) + margin)
    if right <= left:
        right = left + 100
    return left, right


def _trace_event_query_indices(
    events: list[tuple[int | None, int, str, str]],
    ref_to_query: dict[int, list[int]],
) -> list[int]:
    """Return peak indices for local alignment columns, including deletion flanks."""
    indices: set[int] = set()
    for query_index, ref_pos, _wt_base, _mut_base in events:
        if query_index is not None:
            indices.add(query_index)
            continue
        for nearby_ref_pos in (ref_pos - 1, ref_pos, ref_pos + 1):
            indices.update(ref_to_query.get(nearby_ref_pos, []))
    return sorted(index for index in indices if index > 0)


def render_analysis_figure(
    output_png: Path,
    output_pdf: Path,
    sample_name: str,
    reference: str,
    read: AB1Read,
    alignment: AlignmentResult,
    guide_hits: list[GuideHit],
    variants: list[VariantCall],
    ref_to_query: dict[int, list[int]],
    query_to_ref: dict[int, int],
    ref_protein: str,
    mut_protein: str,
    effect_summary: str,
    cds_mapping: CDSMapping,
    guide_cds_start: int | None = None,
    guide_cds_end: int | None = None,
    protein_flank_aa: tuple[int, int] = (8, 8),
    sequence_flank_bp: tuple[int, int] = (20, 20),
    display_reference_span: tuple[int, int] | None = None,
    structure_mode: str | None = None,
    structure_total_length: int | None = None,
    structure_regions: list[tuple[str, int, int, str]] | None = None,
    structure_sgrna_span: tuple[int, int] | None = None,
    sequence_row_spacing: float = 0.50,
    dna_block_gap: float = 0.28,
    dna_protein_gap: float = 0.25,
    sgrna_line_gap: float = 0.16,
    dna_font_size: float = 6.0,
    protein_font_size: float = 6.0,
    structure_label_font_size: float = 10.0,
    reference_label: str = "Reference",
    mutant_label: str = "Mutant",
) -> None:
    if display_reference_span is None:
        if not guide_hits:
            display_reference_span = (1, min(len(reference), SEQUENCE_CHUNK_WIDTH))
        else:
            display_reference_span = reference_window_from_guide(
                len(reference),
                guide_hits[0],
                sequence_flank_bp[0],
                sequence_flank_bp[1],
            )
    display_reference_positions = list(
        range(display_reference_span[0], display_reference_span[1] + 1)
    )
    dna_columns, dna_rows = _build_dna_alignment_columns_from_alignments(
        reference,
        [alignment],
        display_reference_span[0],
        display_reference_span[1],
    )
    protein_columns, protein_rows, protein_aa_start, protein_aa_end = _protein_window_alignment(
        [(reference_label, ref_protein), (mutant_label, mut_protein)],
        guide_cds_start,
        guide_cds_end,
        protein_flank_aa[0],
        protein_flank_aa[1],
        width=SEQUENCE_CHUNK_WIDTH,
    )
    _render_sequence_only_figure(
        output_png,
        output_pdf,
        sample_name,
        dna_columns,
        dna_rows,
        [reference_label, mutant_label],
        guide_hits,
        display_reference_span,
        protein_columns,
        protein_rows,
        protein_aa_start,
        protein_aa_end,
        structure_title=(
            "Gene structure"
            if structure_mode == "gene" and structure_regions
            else "Protein domain"
            if structure_mode == "protein" and structure_regions
            else None
        ),
        structure_total_length=structure_total_length,
        structure_regions=structure_regions,
        structure_unit_label="bp" if structure_mode == "gene" else "aa",
        structure_mode=structure_mode,
        structure_sgrna_span=structure_sgrna_span,
        sequence_row_spacing=sequence_row_spacing,
        dna_block_gap=dna_block_gap,
        dna_protein_gap=dna_protein_gap,
        sgrna_line_gap=sgrna_line_gap,
        dna_font_size=dna_font_size,
        protein_font_size=protein_font_size,
        structure_label_font_size=structure_label_font_size,
    )


def _short_label(value: Path | str, max_len: int = 18) -> str:
    label = value.stem if isinstance(value, Path) else str(value)
    if len(label) <= max_len:
        return label
    return label[: max_len - 1] + "."


def _variant_label(variant: VariantCall) -> str:
    if variant.kind == "insertion":
        return f"+{variant.query_bases}"
    if variant.kind == "deletion":
        return f"-{variant.ref_bases}"
    if variant.kind == "substitution":
        return f"{variant.ref_bases}>{variant.query_bases}"
    return variant.kind


def _variant_trace_positions(
    read: AB1Read,
    variant: VariantCall,
    ref_to_query: dict[int, list[int]],
    query_to_ref: dict[int, int],
) -> list[int]:
    peaks = _peak_positions_for_query(read)
    query_indices = list(variant.query_indices)
    if not query_indices and variant.kind == "deletion":
        query_indices = [
            idx
            for idx, ref_pos in query_to_ref.items()
            if max(1, variant.ref_start - 1) <= ref_pos <= variant.ref_end + 1
        ][:2]
    if not query_indices:
        for ref_pos in range(variant.ref_start, variant.ref_end + 1):
            query_indices.extend(ref_to_query.get(ref_pos, []))
    return [int(peaks[idx - 1]) for idx in query_indices if 1 <= idx <= len(peaks)]


def _protein_window_alignment(
    sample_rows: list[tuple[str, str]],
    guide_cds_start: int | None,
    guide_cds_end: int | None,
    before_aa: int,
    after_aa: int,
    width: int = SEQUENCE_CHUNK_WIDTH,
) -> tuple[list[tuple[int, int]], list[tuple[str, str]], int, int]:
    if not sample_rows:
        return [], [], 0, 0
    wt_label, wt_protein = sample_rows[0]
    mutant_rows = sample_rows[1:]
    if not wt_protein:
        return [], [], 0, 0

    if guide_cds_start is None:
        aa_start, aa_end = _window_from_anchor(
            len(wt_protein),
            1,
            before_aa,
            after_aa,
        )
    else:
        guide_aa_start = max(1, ((guide_cds_start - 1) // 3) + 1)
        aa_start, aa_end = _window_from_anchor(
            len(wt_protein),
            guide_aa_start,
            before_aa,
            after_aa,
        )
        if aa_start > len(wt_protein):
            return [], [], aa_start, aa_end

    per_sample: list[tuple[str, dict[int, str], dict[int, list[str]], int | None]] = []
    max_insertions: dict[int, int] = {}
    for label, mutant_protein in mutant_rows:
        mutant_stop = mutant_protein.find("*")
        wt_stop = wt_protein.find("*")
        premature_stop = mutant_stop >= 0 and (
            wt_stop < 0 or mutant_stop < wt_stop
        )

        # A premature stop is a biological endpoint, not an ordinary short
        # alignment. Keep the shared N-terminal sequence in register and put
        # the stop at the first translated mutant position instead of letting
        # a global aligner move it to the end of the WT sequence.
        if premature_stop:
            prefix_length = 0
            while (
                prefix_length < len(wt_protein)
                and prefix_length < mutant_stop
                and wt_protein[prefix_length] == mutant_protein[prefix_length]
            ):
                prefix_length += 1

            aa_by_ref: dict[int, str] = {}
            insertions: dict[int, list[str]] = {}
            for ref_position in range(1, prefix_length + 1):
                aa_by_ref[ref_position] = mutant_protein[ref_position - 1]
            for offset, mut_aa in enumerate(mutant_protein[prefix_length:]):
                ref_position = prefix_length + offset + 1
                if ref_position <= len(wt_protein):
                    aa_by_ref[ref_position] = mut_aa
                else:
                    insertions.setdefault(len(wt_protein), []).append(mut_aa)
                    max_insertions[len(wt_protein)] = max(
                        max_insertions.get(len(wt_protein), 0),
                        len(insertions[len(wt_protein)]),
                    )
            per_sample.append((label, aa_by_ref, insertions, None))
            continue

        alignment = pairwise2.align.globalms(
            wt_protein,
            mutant_protein,
            2,
            -1,
            -5,
            -1,
            one_alignment_only=True,
        )
        if not alignment:
            continue
        best = alignment[0]
        ref_pos = 0
        aa_by_ref: dict[int, str] = {}
        insertions: dict[int, list[str]] = {}
        for ref_aa, mut_aa in zip(best.seqA, best.seqB):
            if ref_aa != "-":
                ref_pos += 1
                if aa_start <= ref_pos <= aa_end:
                    aa_by_ref[ref_pos] = mut_aa if mut_aa != "-" else "-"
                continue
            if mut_aa != "-" and aa_start - 1 <= ref_pos <= aa_end:
                insertions.setdefault(ref_pos, []).append(mut_aa)
                max_insertions[ref_pos] = max(max_insertions.get(ref_pos, 0), len(insertions[ref_pos]))
        per_sample.append((label, aa_by_ref, insertions, None))

    columns: list[tuple[int, int]] = []
    for aa_pos in range(aa_start, aa_end + 1):
        for insertion_index in range(max_insertions.get(aa_pos - 1, 0)):
            columns.append((aa_pos - 1, insertion_index))
        columns.append((aa_pos, -1))
    for insertion_index in range(max_insertions.get(aa_end, 0)):
        columns.append((aa_end, insertion_index))

    sequences: list[tuple[str, str]] = []
    wt_chars = []
    for aa_pos, insertion_index in columns:
        wt_chars.append("-" if insertion_index >= 0 else wt_protein[aa_pos - 1])
    sequences.append((wt_label, "".join(wt_chars)))

    for label, aa_by_ref, insertions, _ in per_sample:
        chars = []
        for aa_pos, insertion_index in columns:
            if insertion_index >= 0:
                inserted = insertions.get(aa_pos, [])
                chars.append(inserted[insertion_index] if insertion_index < len(inserted) else "-")
            else:
                chars.append(aa_by_ref.get(aa_pos, "-"))
        sequences.append((label, "".join(chars)))

    return columns, sequences, aa_start, aa_end


def _draw_cds_schematic(
    ax,
    cds_length: int,
    guide_start: int | None,
    guide_end: int | None,
    guide_label: str,
) -> None:
    ax.set_xlim(1, max(2, cds_length))
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(
        plt.Rectangle(
            (1, 0.42),
            max(1, cds_length - 1),
            0.18,
            facecolor="white",
            edgecolor="#111111",
            linewidth=1.2,
        )
    )
    ax.text(1, 0.70, "ATG", ha="left", va="bottom", fontsize=8)
    ax.text(cds_length, 0.70, "TGA", ha="right", va="bottom", fontsize=8)
    if guide_start is None:
        return
    start = max(1, min(cds_length, guide_start))
    end = max(start, min(cds_length, guide_end or guide_start))
    ax.add_patch(
        plt.Rectangle(
            (start, 0.38),
            max(1, end - start + 1),
            0.26,
            facecolor="#111111",
            edgecolor="#111111",
            linewidth=0.8,
        )
    )
    ax.text(
        (start + end) / 2,
        0.25,
        guide_label or "sgRNA",
        ha="center",
        va="top",
        fontsize=8,
        family="monospace",
    )


def _guide_trace_span(
    read: AB1Read,
    guide_hits: list[GuideHit],
    ref_to_query: dict[int, list[int]],
) -> tuple[int, int] | None:
    if not guide_hits:
        return None
    peaks = _peak_positions_for_query(read)
    query_indices: list[int] = []
    hit = guide_hits[0]
    for ref_pos in range(hit.start, hit.end + 1):
        query_indices.extend(ref_to_query.get(ref_pos, []))
    positions = [int(peaks[idx - 1]) for idx in query_indices if 1 <= idx <= len(peaks)]
    if not positions:
        return None
    return min(positions), max(positions)


def _trace_sequence_events(
    read: AB1Read,
    alignment: AlignmentResult,
    reference_positions: list[int],
) -> list[tuple[int | None, int, str, str]]:
    if not reference_positions:
        return []

    selected = set(reference_positions)
    first_ref = min(reference_positions)
    last_ref = max(reference_positions)
    events: list[tuple[int | None, int, str, str]] = []
    ref_pos = alignment.reference_start - 1
    query_pos = alignment.query_start - 1

    for ref_base, query_base in zip(alignment.aligned_reference, alignment.aligned_query):
        if ref_base == "-":
            if query_base != "-":
                query_pos += 1
                if first_ref - 1 <= ref_pos <= last_ref:
                    events.append((query_pos, ref_pos, "-", query_base))
            continue

        ref_pos += 1
        if query_base != "-":
            query_pos += 1
        if ref_pos not in selected:
            continue
        if query_base != "-":
            events.append((query_pos, ref_pos, ref_base, query_base))
        else:
            events.append((None, ref_pos, ref_base, "-"))

    return events


def _draw_mutant_sequence_on_trace(
    ax,
    read: AB1Read,
    alignment: AlignmentResult,
    reference_positions: list[int],
    ref_to_query: dict[int, list[int]],
    left: int,
    right: int,
    ymax: float,
) -> None:
    peaks = _peak_positions_for_query(read)
    events = _trace_sequence_events(read, alignment, reference_positions)
    if not events:
        return

    sequence_y = ymax * 1.22
    ax.text(
        0.005,
        sequence_y,
        "Mutant bases",
        transform=ax.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=6.5,
    )

    for query_index, ref_pos, reference_base, mutant_base in events:
        if query_index is not None and 1 <= query_index <= len(peaks):
            peak_x = int(peaks[query_index - 1])
        else:
            left_indices = ref_to_query.get(ref_pos - 1, [])
            right_indices = ref_to_query.get(ref_pos + 1, [])
            left_peaks = [int(peaks[idx - 1]) for idx in left_indices if 1 <= idx <= len(peaks)]
            right_peaks = [int(peaks[idx - 1]) for idx in right_indices if 1 <= idx <= len(peaks)]
            if left_peaks and right_peaks:
                peak_x = int((max(left_peaks) + min(right_peaks)) / 2)
            elif left_peaks:
                peak_x = max(left_peaks) + 8
            elif right_peaks:
                peak_x = min(right_peaks) - 8
            else:
                continue
        if not left <= peak_x <= right:
            continue
        changed = reference_base != mutant_base
        ax.text(
            peak_x,
            sequence_y,
            mutant_base,
            ha="center",
            va="center",
            fontsize=7,
            family="monospace",
            color="#d62728" if changed else "#111111",
        )


def _reference_based_sequence_map(
    aligned_reference: str,
    aligned_query: str,
    reference_start: int = 1,
) -> tuple[dict[int, str], dict[int, list[str]]]:
    reference_bases: dict[int, str] = {}
    insertions: dict[int, list[str]] = {}
    ref_pos = reference_start - 1
    for ref_base, query_base in zip(aligned_reference, aligned_query):
        if ref_base == "-":
            if query_base != "-":
                insertions.setdefault(ref_pos, []).append(query_base)
            continue
        ref_pos += 1
        reference_bases[ref_pos] = query_base if query_base != "-" else "-"
    return reference_bases, insertions


def _build_dna_alignment_columns(
    reference: str,
    prepared_rows: list[tuple[SampleAnalysis, AB1Read, AlignmentResult, dict[int, list[int]], dict[int, int], int, int]],
    window_start: int,
    window_end: int,
) -> tuple[list[tuple[int, int]], list[str]]:
    alignments = [item[2] for item in prepared_rows]
    return _build_dna_alignment_columns_from_alignments(
        reference,
        alignments,
        window_start,
        window_end,
    )


def _build_dna_alignment_columns_from_alignments(
    reference: str,
    alignments: list[AlignmentResult],
    window_start: int,
    window_end: int,
) -> tuple[list[tuple[int, int]], list[str]]:
    maps = [
        _reference_based_sequence_map(
            alignment.aligned_reference,
            alignment.aligned_query,
            alignment.reference_start,
        )
        for alignment in alignments
    ]
    max_insertions: dict[int, int] = {}
    for _, insertions in maps:
        for anchor, bases in insertions.items():
            if window_start - 1 <= anchor <= window_end:
                max_insertions[anchor] = max(max_insertions.get(anchor, 0), len(bases))

    columns: list[tuple[int, int]] = []
    for ref_pos in range(window_start, window_end + 1):
        for insertion_index in range(max_insertions.get(ref_pos - 1, 0)):
            columns.append((ref_pos - 1, insertion_index))
        columns.append((ref_pos, -1))
    for insertion_index in range(max_insertions.get(window_end, 0)):
        columns.append((window_end, insertion_index))

    reference_row = []
    for ref_pos, insertion_index in columns:
        if insertion_index >= 0:
            reference_row.append("-")
        else:
            reference_row.append(reference[ref_pos - 1])

    rows = ["".join(reference_row)]
    for alignment in alignments:
        base_map, insertion_map = _reference_based_sequence_map(
            alignment.aligned_reference,
            alignment.aligned_query,
            alignment.reference_start,
        )
        row: list[str] = []
        for ref_pos, insertion_index in columns:
            if insertion_index >= 0:
                inserted = insertion_map.get(ref_pos, [])
                row.append(inserted[insertion_index] if insertion_index < len(inserted) else "-")
            else:
                row.append(base_map.get(ref_pos, "-"))
        rows.append("".join(row))
    return columns, rows


def _draw_colored_sequence(
    ax,
    label: str,
    sequence: str,
    reference_row: str,
    y: float,
    x_offset: float = 0.0,
    fontsize: float = 7.0,
    char_step: float = SEQUENCE_CHAR_STEP,
    label_fontstyle: str = "normal",
) -> None:
    ax.text(
        x_offset - 1.0,
        y,
        label,
        ha="right",
        va="center",
        fontsize=fontsize,
        fontstyle=label_fontstyle,
    )
    for idx, base in enumerate(sequence):
        ref_base = reference_row[idx] if idx < len(reference_row) else "-"
        changed = base != ref_base
        if changed and base == "-":
            color = "#d62728"
        elif changed:
            color = "#d62728"
        else:
            color = "#111111"
        ax.text(
            x_offset + idx * char_step,
            y,
            base,
            ha="center",
            va="center",
            fontsize=fontsize,
            family="monospace",
            color=color,
        )


def _compact_sequence_fontsize(column_count: int) -> float:
    """Keep sequence characters readable while filling a fixed-width panel."""
    if column_count <= 0:
        return 5.8
    return max(5.0, min(6.0, 6.0 * SEQUENCE_CHUNK_WIDTH / column_count))


def _chunk_ranges(total_length: int, chunk_width: int) -> list[tuple[int, int]]:
    if total_length <= 0:
        return []
    return [
        (start, min(total_length, start + chunk_width))
        for start in range(0, total_length, chunk_width)
    ]


def _chunk_column_ranges_by_reference(
    columns: list[tuple[int, int]],
    chunk_width: int,
) -> list[tuple[int, int]]:
    if not columns:
        return []
    ranges: list[tuple[int, int]] = []
    chunk_start = 0
    ref_count = 0
    for idx, (_ref_pos, insertion_index) in enumerate(columns):
        # Count reference bases only. Inserted bases are kept in the same
        # block as their anchor and do not consume the 60-base quota.
        if insertion_index < 0 and ref_count == chunk_width:
            ranges.append((chunk_start, idx))
            chunk_start = idx
            ref_count = 0
        if insertion_index < 0:
            ref_count += 1
    ranges.append((chunk_start, len(columns)))
    return ranges


def _draw_dna_alignment_panel(
    ax,
    columns: list[tuple[int, int]],
    rows: list[str],
    labels: list[str],
    guide_hits: list[GuideHit],
    window_start: int,
    window_end: int,
    fontsize: float = 7.0,
    panel_width: int | None = None,
    chunk_width: int = SEQUENCE_CHUNK_WIDTH,
    sequence_row_spacing: float = 0.50,
    dna_block_gap: float = 0.28,
    sgrna_line_gap: float = 0.16,
) -> None:
    if not columns or not rows:
        ax.axis("off")
        return

    sequence_row_spacing = max(0.4, min(2.0, float(sequence_row_spacing)))
    dna_block_gap = max(0.05, min(2.0, float(dna_block_gap)))
    sgrna_line_gap = max(0.05, min(2.0, float(sgrna_line_gap)))
    n_rows = min(len(rows), len(labels))
    chunks = _chunk_column_ranges_by_reference(columns, chunk_width)
    row_step = 0.72 * sequence_row_spacing
    chunk_height = max(0, n_rows - 1) * row_step + dna_block_gap
    total_height = max(1.0, len(chunks) * chunk_height)
    display_width = panel_width or chunk_width
    x_limit = max(1.0, (display_width - 1) * SEQUENCE_CHAR_STEP + 2.0)
    ax.set_xlim(SEQUENCE_LEFT_MARGIN, x_limit)
    top_padding = max(0.42 * sequence_row_spacing, sgrna_line_gap + 0.22)
    ax.set_ylim(-0.18 * sequence_row_spacing, total_height + top_padding)
    ax.axis("off")
    reference_row = rows[0]
    for chunk_index, (start, end) in enumerate(chunks):
        chunk_columns = columns[start:end]
        top_y = total_height - chunk_index * chunk_height - 0.26 * sequence_row_spacing
        ref_positions = [position for position, insertion_index in chunk_columns if insertion_index < 0]

        if guide_hits and ref_positions:
            hit = guide_hits[0]
            guide_indices = [
                idx
                for idx, (ref_pos, insertion_index) in enumerate(chunk_columns)
                if insertion_index < 0 and hit.start <= ref_pos <= hit.end
            ]
            if guide_indices:
                y_guide = top_y + sgrna_line_gap
                ax.plot(
                    [
                        (min(guide_indices) - 0.25) * SEQUENCE_CHAR_STEP,
                        (max(guide_indices) + 0.25) * SEQUENCE_CHAR_STEP,
                    ],
                    [y_guide, y_guide],
                    color="#111111",
                    linewidth=1.0,
                )
                ax.text(
                    ((min(guide_indices) + max(guide_indices)) / 2) * SEQUENCE_CHAR_STEP,
                    y_guide + max(0.08, sgrna_line_gap * 0.45),
                    "sgRNA",
                    ha="center",
                    va="bottom",
                    fontsize=fontsize,
                )

        for idx in range(n_rows):
            y = top_y - idx * row_step
            _draw_colored_sequence(
                ax,
                labels[idx],
                rows[idx][start:end],
                reference_row[start:end],
                y,
                x_offset=0,
                fontsize=fontsize,
                char_step=SEQUENCE_CHAR_STEP,
                label_fontstyle="italic" if idx > 0 else "normal",
            )


def _draw_protein_alignment_panel(
    ax,
    columns: list[tuple[int, int]],
    sequences: list[tuple[str, str]],
    aa_start: int,
    aa_end: int,
    fontsize: float = 7.0,
    panel_width: int | None = None,
    chunk_width: int = SEQUENCE_CHUNK_WIDTH,
    sequence_row_spacing: float = 0.50,
    dna_block_gap: float = 0.28,
) -> None:
    if not columns or not sequences:
        ax.axis("off")
        ax.text(0.0, 1.0, "Protein alignment unavailable.", ha="left", va="top", fontsize=fontsize)
        return

    sequence_row_spacing = max(0.4, min(2.0, float(sequence_row_spacing)))
    dna_block_gap = max(0.05, min(2.0, float(dna_block_gap)))
    chunks = _chunk_column_ranges_by_reference(columns, chunk_width)
    row_count = len(sequences)
    row_step = 0.72 * sequence_row_spacing
    block_height = max(0, row_count - 1) * row_step + dna_block_gap
    total_height = max(1.0, len(chunks) * block_height)
    display_width = panel_width or chunk_width
    x_limit = max(1.0, (display_width - 1) * SEQUENCE_CHAR_STEP + 2.0)
    ax.set_xlim(SEQUENCE_LEFT_MARGIN, x_limit)
    ax.set_ylim(-0.18 * sequence_row_spacing, total_height + 0.30 * sequence_row_spacing)
    ax.axis("off")

    reference_row = sequences[0][1]
    for chunk_index, (start, end) in enumerate(chunks):
        chunk_columns = columns[start:end]
        top_y = total_height - chunk_index * block_height - 0.20 * sequence_row_spacing
        for row_index, (label, sequence) in enumerate(sequences):
            y = top_y - row_index * row_step
            _draw_colored_sequence(
                ax,
                label,
                sequence[start:end],
                reference_row[start:end],
                y,
                x_offset=0,
                fontsize=fontsize,
                char_step=SEQUENCE_CHAR_STEP,
                label_fontstyle="italic" if row_index > 0 else "normal",
            )


def _guide_alignment_point(
    columns: list[tuple[int, int]],
    guide_hits: list[GuideHit],
    n_rows: int,
    chunk_width: int = SEQUENCE_CHUNK_WIDTH,
    sequence_row_spacing: float = 0.50,
    dna_block_gap: float = 0.28,
    sgrna_line_gap: float = 0.16,
) -> tuple[float, float] | None:
    if not guide_hits or not columns:
        return None
    sequence_row_spacing = max(0.4, min(2.0, float(sequence_row_spacing)))
    dna_block_gap = max(0.05, min(2.0, float(dna_block_gap)))
    sgrna_line_gap = max(0.05, min(2.0, float(sgrna_line_gap)))
    chunks = _chunk_column_ranges_by_reference(columns, chunk_width)
    row_step = 0.72 * sequence_row_spacing
    chunk_height = max(0, n_rows - 1) * row_step + dna_block_gap
    total_height = max(1.0, len(chunks) * chunk_height)
    hit = guide_hits[0]
    for chunk_index, (start, end) in enumerate(chunks):
        chunk_columns = columns[start:end]
        guide_indices = [
            idx
            for idx, (ref_pos, insertion_index) in enumerate(chunk_columns)
            if insertion_index < 0 and hit.start <= ref_pos <= hit.end
        ]
        if guide_indices:
            top_y = total_height - chunk_index * chunk_height - 0.26 * sequence_row_spacing
            return (
                min(guide_indices) * SEQUENCE_CHAR_STEP,
                top_y + sgrna_line_gap,
            )
    return None


def parse_structure_regions(raw: str) -> list[tuple[str, int, int, str]]:
    regions: list[tuple[str, int, int, str]] = []
    for idx, line in enumerate(raw.splitlines(), start=1):
        line = line.strip()
        if not line:
            continue
        parts = re.split(r"[,\s;]+", line)
        if len(parts) == 2:
            label = f"Region {idx}"
            start_text, end_text = parts
            color = "#1f77b4"
        elif len(parts) == 3:
            label = parts[0]
            start_text, end_text = parts[1], parts[2]
            color = "#1f77b4"
        elif len(parts) >= 4:
            label = parts[0]
            start_text, end_text, color = parts[1], parts[2], parts[3]
        else:
            raise ValueError(
                "Structure regions must be entered as 'start end', 'label start end', or 'label start end color' on separate lines."
            )
        start = int(start_text)
        end = int(end_text)
        if start <= 0 or end <= 0:
            raise ValueError("Structure region positions must be positive integers.")
        if end < start:
            start, end = end, start
        regions.append((label, start, end, color))
    return regions


def _draw_structure_panel(
    ax,
    total_length: int,
    regions: list[tuple[str, int, int, str]],
    title: str,
    unit_label: str,
    structure_mode: str | None = None,
    sgrna_span: tuple[int, int] | None = None,
    display_columns: int = SEQUENCE_CHUNK_WIDTH,
    structure_label_font_size: float = 10.0,
) -> None:
    display_columns = max(1, int(display_columns))
    sequence_extent = max(0.0, (display_columns - 1) * SEQUENCE_CHAR_STEP)
    ax.set_xlim(SEQUENCE_LEFT_MARGIN, max(1.0, sequence_extent + 2.0))
    ax.set_ylim(0, 1)
    ax.axis("off")

    # The structure is deliberately rendered on the same physical width as
    # one DNA row. Structure coordinates are normalized into that 60-column
    # display grid instead of stretching the full protein/CDS length across
    # the page.
    def position_to_x(position: int) -> float:
        if total_length <= 1:
            return 0.0
        clipped = max(1, min(total_length, int(position)))
        return (clipped - 1) / (total_length - 1) * sequence_extent

    palette = ["#4f81bd", "#c0504d", "#4f81bd", "#f79646", "#9bbb59", "#8064a2"]
    box_y = 0.40
    box_h = 0.30
    line_y = box_y + box_h / 2.0
    ax.plot(
        [0.0, sequence_extent],
        [line_y, line_y],
        color="#111111",
        linewidth=0.9,
        zorder=0,
    )
    for idx, region in enumerate(regions):
        if len(region) == 4:
            label, start, end, color = region
        else:
            label, start, end = region  # type: ignore[misc]
            color = palette[idx % len(palette)]
        start = max(1, min(total_length, start))
        end = max(start, min(total_length, end))
        if structure_mode == "gene":
            facecolor = "#111111"
            edgecolor = "#111111"
        else:
            facecolor = color
            edgecolor = "#111111"
        box_start = position_to_x(start)
        box_end = position_to_x(end)
        ax.add_patch(
            plt.Rectangle(
                (box_start, box_y),
                max(SEQUENCE_CHAR_STEP, box_end - box_start + SEQUENCE_CHAR_STEP),
                box_h,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=0.8,
                alpha=1.0,
                zorder=2,
            )
        )
        ax.text(
            (box_start + box_end) / 2,
            0.82,
            label,
            ha="center",
            va="bottom",
            fontsize=structure_label_font_size,
            zorder=3,
        )
    if sgrna_span is not None:
        sgrna_start, _sgrna_end = sgrna_span
        sgrna_start = max(1, min(total_length, sgrna_start))
        sgrna_start_x = position_to_x(sgrna_start)
        ax.plot(
            [sgrna_start_x, sgrna_start_x],
            [line_y - 0.04, line_y + 0.04],
            color="#111111",
            linewidth=0.9,
            zorder=3,
        )


def _render_sequence_only_figure(
    output_png: Path,
    output_pdf: Path,
    title: str,
    dna_columns: list[tuple[int, int]],
    dna_rows: list[str],
    dna_labels: list[str],
    guide_hits: list[GuideHit],
    reference_span: tuple[int, int],
    protein_columns: list[tuple[int, int]],
    protein_rows: list[tuple[str, str]],
    protein_aa_start: int,
    protein_aa_end: int,
    structure_title: str | None = None,
    structure_total_length: int | None = None,
    structure_regions: list[tuple[str, int, int, str]] | None = None,
    structure_unit_label: str = "bp",
    structure_mode: str | None = None,
    structure_sgrna_span: tuple[int, int] | None = None,
    sequence_row_spacing: float = 0.50,
    dna_block_gap: float = 0.28,
    dna_protein_gap: float = 0.25,
    sgrna_line_gap: float = 0.16,
    dna_font_size: float = 6.0,
    protein_font_size: float = 6.0,
    structure_label_font_size: float = 10.0,
) -> None:
    dna_chunk_ranges = _chunk_column_ranges_by_reference(dna_columns, SEQUENCE_CHUNK_WIDTH)
    protein_chunk_ranges = _chunk_column_ranges_by_reference(protein_columns, SEQUENCE_CHUNK_WIDTH)
    dna_chunks = max(1, len(dna_chunk_ranges))
    protein_chunks = max(1, len(protein_chunk_ranges))
    dna_display_columns = max(
        SEQUENCE_CHUNK_WIDTH,
        max((end - start for start, end in dna_chunk_ranges), default=SEQUENCE_CHUNK_WIDTH),
    )
    protein_display_columns = max(
        SEQUENCE_CHUNK_WIDTH,
        max((end - start for start, end in protein_chunk_ranges), default=SEQUENCE_CHUNK_WIDTH),
    )
    sequence_row_spacing = max(0.4, min(2.0, float(sequence_row_spacing)))
    dna_block_gap = max(0.05, min(2.0, float(dna_block_gap)))
    dna_protein_gap = max(0.0, min(2.0, float(dna_protein_gap)))
    sgrna_line_gap = max(0.05, min(2.0, float(sgrna_line_gap)))
    dna_font_size = max(4.0, min(14.0, float(dna_font_size)))
    protein_font_size = max(4.0, min(14.0, float(protein_font_size)))
    structure_label_font_size = max(4.0, min(18.0, float(structure_label_font_size)))
    structure_enabled = bool(structure_title and structure_total_length and structure_regions)
    structure_height = 0.60 if structure_enabled else 0.0
    dna_height = max(
        0.40,
        dna_chunks * (max(0, len(dna_rows) - 1) * 0.72 * sequence_row_spacing + dna_block_gap),
    )
    protein_height = max(
        0.40,
        protein_chunks * (max(0, len(protein_rows) - 1) * 0.72 * sequence_row_spacing + dna_block_gap),
    )
    spacer_height = max(0.01, dna_protein_gap)
    fig_height = max(
        1.45,
        0.10 + structure_height + dna_height + spacer_height + protein_height,
    )
    fig = plt.figure(figsize=(10, fig_height))
    if structure_enabled:
        gs = fig.add_gridspec(
            4,
            1,
            height_ratios=[structure_height, dna_height, spacer_height, protein_height],
            hspace=0.01,
        )
    else:
        gs = fig.add_gridspec(
            3,
            1,
            height_ratios=[dna_height, spacer_height, protein_height],
            hspace=0.01,
        )
    fig.subplots_adjust(left=0.04, right=0.99, top=0.98, bottom=0.04)

    if structure_enabled:
        ax_structure = fig.add_subplot(gs[0, 0])
        _draw_structure_panel(
            ax_structure,
            int(structure_total_length),
            structure_regions or [],
            structure_title or "",
            structure_unit_label,
            structure_mode=structure_mode,
            sgrna_span=structure_sgrna_span,
            display_columns=dna_display_columns,
            structure_label_font_size=structure_label_font_size,
        )
        dna_row_index = 1
        spacer_row_index = 2
        protein_row_index = 3
    else:
        dna_row_index = 0
        spacer_row_index = 1
        protein_row_index = 2

    ax_dna = fig.add_subplot(gs[dna_row_index, 0])
    _draw_dna_alignment_panel(
        ax_dna,
        dna_columns,
        dna_rows,
        dna_labels,
        guide_hits,
        reference_span[0],
        reference_span[1],
        fontsize=dna_font_size,
        panel_width=dna_display_columns,
        chunk_width=SEQUENCE_CHUNK_WIDTH,
        sequence_row_spacing=sequence_row_spacing,
        dna_block_gap=dna_block_gap,
        sgrna_line_gap=sgrna_line_gap,
    )

    if structure_enabled and structure_sgrna_span is not None:
        guide_point = _guide_alignment_point(
            dna_columns,
            guide_hits,
            len(dna_rows),
            sequence_row_spacing=sequence_row_spacing,
            dna_block_gap=dna_block_gap,
            sgrna_line_gap=sgrna_line_gap,
        )
        if guide_point is not None:
            sgrna_start, _sgrna_end = structure_sgrna_span
            structure_length = max(1, int(structure_total_length or 1))
            structure_extent = max(
                0.0,
                (dna_display_columns - 1) * SEQUENCE_CHAR_STEP,
            )
            if structure_length <= 1:
                sgrna_x = 0.0
            else:
                structure_start = max(1, min(structure_length, int(sgrna_start)))
                sgrna_x = (
                    (structure_start - 1)
                    / (structure_length - 1)
                    * structure_extent
                )
            connector = ConnectionPatch(
                xyA=(sgrna_x, 0.55),
                xyB=guide_point,
                coordsA="data",
                coordsB="data",
                axesA=ax_structure,
                axesB=ax_dna,
                arrowstyle="-",
                linewidth=0.9,
                linestyle="-",
                color="#555555",
            )
            fig.add_artist(connector)

    ax_spacer = fig.add_subplot(gs[spacer_row_index, 0])
    ax_spacer.axis("off")

    ax_protein = fig.add_subplot(gs[protein_row_index, 0])
    _draw_protein_alignment_panel(
        ax_protein,
        protein_columns,
        protein_rows,
        protein_aa_start,
        protein_aa_end,
        fontsize=protein_font_size,
        panel_width=protein_display_columns,
        chunk_width=SEQUENCE_CHUNK_WIDTH,
        sequence_row_spacing=sequence_row_spacing,
        dna_block_gap=dna_block_gap,
    )

    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def render_batch_publication_figure(
    output_png: Path,
    output_pdf: Path,
    reference: str,
    guide: str,
    analyses: list[SampleAnalysis],
    protein_flank_aa: tuple[int, int] = (8, 8),
    structure_mode: str | None = None,
    structure_total_length: int | None = None,
    structure_regions: list[tuple[str, int, int, str]] | None = None,
    structure_sgrna_span: tuple[int, int] | None = None,
    sequence_row_spacing: float = 0.50,
    dna_block_gap: float = 0.28,
    dna_protein_gap: float = 0.25,
    sgrna_line_gap: float = 0.16,
    dna_font_size: float = 6.0,
    protein_font_size: float = 6.0,
    structure_label_font_size: float = 10.0,
    reference_label: str = "Reference",
) -> None:
    if not analyses:
        return

    guide_hits = find_guide_hits(reference, guide)
    guide_cds_start = analyses[0].guide_cds_start
    guide_cds_end = analyses[0].guide_cds_end
    alignments: list[AlignmentResult] = []
    for analysis in analyses:
        read = load_sequence_read(analysis.input_path)
        _, alignment = choose_orientation(reference, read)
        alignments.append(alignment)

    dna_columns, dna_rows = _build_dna_alignment_columns_from_alignments(
        reference,
        alignments,
        analyses[0].display_reference_start,
        analyses[0].display_reference_end,
    )
    dna_labels = [reference_label]
    dna_labels.extend(_short_label(item.sample_label) for item in analyses)
    protein_rows = [(reference_label, analyses[0].ref_protein)]
    protein_rows.extend((_short_label(item.sample_label), item.mut_protein) for item in analyses)
    protein_columns, protein_rows, protein_aa_start, protein_aa_end = _protein_window_alignment(
        protein_rows,
        guide_cds_start,
        guide_cds_end,
        protein_flank_aa[0],
        protein_flank_aa[1],
        width=SEQUENCE_CHUNK_WIDTH,
    )
    structure_title = None
    resolved_structure_total_length = structure_total_length
    structure_unit_label = "bp"
    if structure_mode == "protein" and structure_regions:
        structure_title = "Protein domain"
        resolved_structure_total_length = resolved_structure_total_length or len(analyses[0].ref_protein)
        structure_unit_label = "aa"
    elif structure_mode == "gene" and structure_regions:
        structure_title = "Gene structure"
        resolved_structure_total_length = resolved_structure_total_length or len(analyses[0].cds_sequence)
        structure_unit_label = "bp"
    _render_sequence_only_figure(
        output_png,
        output_pdf,
        "Publication mutation summary",
        dna_columns,
        dna_rows,
        dna_labels,
        guide_hits,
        (
            analyses[0].display_reference_start,
            analyses[0].display_reference_end,
        ),
        protein_columns,
        protein_rows,
        protein_aa_start,
        protein_aa_end,
        structure_title=structure_title,
        structure_total_length=resolved_structure_total_length,
        structure_regions=structure_regions,
        structure_unit_label=structure_unit_label,
        structure_mode=structure_mode,
        structure_sgrna_span=structure_sgrna_span,
        sequence_row_spacing=sequence_row_spacing,
        dna_block_gap=dna_block_gap,
        dna_protein_gap=dna_protein_gap,
        sgrna_line_gap=sgrna_line_gap,
        dna_font_size=dna_font_size,
        protein_font_size=protein_font_size,
        structure_label_font_size=structure_label_font_size,
    )


def write_variants_csv(output: Path, variants: list[VariantCall], mapping: CDSMapping) -> None:
    with output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "kind",
                "pcr_reference_start",
                "pcr_reference_end",
                "reference_bases",
                "query_bases",
                "pcr_notation",
                "cds_notation",
            ]
        )
        for variant in variants:
            writer.writerow(
                [
                    variant.kind,
                    variant.ref_start,
                    variant.ref_end,
                    variant.ref_bases,
                    variant.query_bases,
                    format_pcr_variant_notation(variant),
                    format_variant_notation(variant, mapping),
                ]
            )


def analyze_sample(
    reference: str,
    guide: str,
    cds_sequence: str,
    read_path: Path,
    output_root: Path,
    sample_label: str | None = None,
    protein_flank_aa: tuple[int, int] = (8, 8),
    sequence_flank_bp: tuple[int, int] = (20, 20),
    structure_mode: str | None = None,
    structure_total_length: int | None = None,
    structure_regions: list[tuple[str, int, int, str]] | None = None,
    structure_sgrna_span: tuple[int, int] | None = None,
    sequence_row_spacing: float = 0.50,
    dna_block_gap: float = 0.28,
    dna_protein_gap: float = 0.25,
    sgrna_line_gap: float = 0.16,
    dna_font_size: float = 6.0,
    protein_font_size: float = 6.0,
    structure_label_font_size: float = 10.0,
    reference_label: str = "Reference",
    mutant_label: str = "Mutant",
) -> SampleAnalysis:
    reference = clean_dna(reference)
    if not reference:
        raise ValueError("PCR amplicon reference sequence is empty.")
    guide = clean_dna(guide)
    if not guide:
        raise ValueError("sgRNA sequence is empty.")
    cds_sequence = clean_dna(cds_sequence)
    mapping = map_reference_to_cds(reference, cds_sequence)

    read = load_sequence_read(read_path)
    oriented_read, alignment = choose_orientation(reference, read)
    variants = call_variants(
        alignment.aligned_reference,
        alignment.aligned_query,
        alignment.reference_start,
        alignment.query_start,
    )
    guide_hits = find_guide_hits(reference, guide)
    if not guide_hits:
        raise ValueError(
            "The sgRNA sequence was not found in the PCR amplicon reference. "
            "Use a PCR reference that covers the sgRNA target region."
        )
    display_reference_start, display_reference_end = reference_window_from_guide(
        len(reference),
        guide_hits[0],
        sequence_flank_bp[0],
        sequence_flank_bp[1],
    )
    variants = variants_in_reference_window(
        variants,
        display_reference_start,
        display_reference_end,
    )
    guide_cds_start, guide_cds_end = guide_cds_span(cds_sequence, guide)
    display_reference_positions = list(
        range(display_reference_start, display_reference_end + 1)
    )
    ref_to_query, query_to_ref = build_ref_to_query_map(
        alignment.aligned_reference,
        alignment.aligned_query,
        alignment.reference_start,
        alignment.query_start,
    )
    if not _trace_sequence_events(oriented_read, alignment, display_reference_positions):
        raise ValueError(
            "This sequencing read does not cover the selected sgRNA-centered reference nucleotide window."
        )
    mutant_cds, unmapped_variants, first_cds_position = apply_variants_to_cds(
        cds_sequence,
        reference,
        variants,
        mapping,
    )
    ref_protein = translate_dna(cds_sequence)
    mut_protein = translate_dna(mutant_cds)
    effect_summary = summarize_effect(
        variants,
        ref_protein,
        mut_protein,
        first_cds_position,
    )

    sample_dir = output_root / safe_stem(read_path)
    sample_dir.mkdir(parents=True, exist_ok=True)
    label = (
        sample_label.strip()
        if sample_label and sample_label.strip()
        else mutant_label.strip() or read_path.stem
    )

    figure_png = sample_dir / f"{safe_stem(read_path)}_mutation_summary.png"
    figure_pdf = sample_dir / f"{safe_stem(read_path)}_mutation_summary.pdf"
    report_txt = sample_dir / f"{safe_stem(read_path)}_report.txt"
    variants_csv = sample_dir / f"{safe_stem(read_path)}_variants.csv"

    render_analysis_figure(
        output_png=figure_png,
        output_pdf=figure_pdf,
        sample_name=label,
        reference=reference,
        read=oriented_read,
        alignment=alignment,
        guide_hits=guide_hits,
        variants=variants,
        ref_to_query=ref_to_query,
        query_to_ref=query_to_ref,
        ref_protein=ref_protein,
        mut_protein=mut_protein,
        effect_summary=effect_summary,
        cds_mapping=mapping,
        guide_cds_start=guide_cds_start,
        guide_cds_end=guide_cds_end,
        protein_flank_aa=protein_flank_aa,
        sequence_flank_bp=sequence_flank_bp,
        display_reference_span=(display_reference_start, display_reference_end),
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
        reference_label=reference_label.strip() or "Reference",
        mutant_label=label,
    )
    write_variants_csv(variants_csv, variants, mapping)

    report_lines = [
        f"Input file: {read_path}",
        f"Sample label: {label}",
        f"Orientation: {'reverse-complement' if alignment.reverse_complement else 'forward'}",
        f"Alignment score: {alignment.score:.2f}",
        f"PCR reference length: {len(reference)}",
        f"Read length: {len(oriented_read.sequence)}",
        f"PCR reference overlap: PCR {mapping.reference_start}-{mapping.reference_end} -> "
        f"CDS {mapping.cds_start}-{mapping.cds_end} "
        f"({'reverse-complement' if mapping.reverse_reference else 'forward'})",
        f"sgRNA CDS position: {guide_cds_start or 'n/a'}-{guide_cds_end or 'n/a'}",
        f"Displayed PCR reference nucleotide window: {display_reference_start}-{display_reference_end}",
        f"CDS length: {len(cds_sequence)}",
        f"Unmapped variants: {len(unmapped_variants)}",
        "",
        f"Reference protein: {ref_protein or 'n/a'}",
        f"Mutant protein: {mut_protein or 'n/a'}",
        f"Effect summary: {effect_summary}",
        "",
        "Guide hits:",
    ]
    if guide_hits:
        for hit in guide_hits:
            report_lines.append(f"  - {hit.start}-{hit.end} ({hit.strand}) {hit.sequence}")
    else:
        report_lines.append("  - not found")
    report_lines.append("")
    report_lines.append("Variants:")
    if variants:
        for variant in variants:
            report_lines.append(
                f"  - {format_variant_notation(variant, mapping)}  "
                f"[{format_pcr_variant_notation(variant)}; {variant.kind}]"
            )
    else:
        report_lines.append("  - none")
    report_lines.append("")
    report_lines.append("Aligned reference:")
    report_lines.append(alignment.aligned_reference)
    report_lines.append("")
    report_lines.append("Aligned query:")
    report_lines.append(alignment.aligned_query)
    report_txt.write_text("\n".join(report_lines), encoding="utf-8")

    return SampleAnalysis(
        input_path=read_path,
        sample_label=label,
        output_dir=sample_dir,
        orientation="reverse-complement" if alignment.reverse_complement else "forward",
        aligned_read=alignment.aligned_query,
        aligned_reference=alignment.aligned_reference,
        reference=reference,
        query=oriented_read.sequence,
        cds_sequence=cds_sequence,
        cds_start=1,
        cds_end=len(cds_sequence),
        reverse_cds=False,
        guide_hits=guide_hits,
        variants=variants,
        ref_protein=ref_protein,
        mut_protein=mut_protein,
        effect_summary=effect_summary,
        figure_png=figure_png,
        figure_pdf=figure_pdf,
        report_txt=report_txt,
        variants_csv=variants_csv,
        cds_map_start=mapping.cds_start,
        cds_map_end=mapping.cds_end,
        cds_map_reverse=mapping.reverse_reference,
        unmapped_variants=len(unmapped_variants),
        guide_cds_start=guide_cds_start,
        guide_cds_end=guide_cds_end,
        display_reference_start=display_reference_start,
        display_reference_end=display_reference_end,
    )
