"""Rotate circular mitochondrial assemblies to the provided reference start."""

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq


ANCHOR_LENGTH = 300
MIN_ANCHOR_COVERAGE = 0.90


def _reference_anchor(reference_fasta):
    reference = str(SeqIO.read(reference_fasta, "fasta").seq).upper()
    if len(reference) < ANCHOR_LENGTH:
        raise ValueError("reference is shorter than the 300 bp rotation anchor")
    return reference[:ANCHOR_LENGTH]


def _local_alignment_candidates(anchor, sequence, strand):
    alignments = pairwise2.align.localms(anchor, sequence + sequence, 2, -1, -5, -1)
    if not alignments:
        return []

    best_score = alignments[0][2]
    candidates = {}
    for alignment in alignments:
        aligned_anchor, aligned_sequence, score, alignment_start, alignment_end = alignment
        if score != best_score:
            break

        query_index = -1
        target_index = -1
        query_indices = []
        target_base_for_anchor_start = None

        for column, (query_base, target_base) in enumerate(
            zip(aligned_anchor, aligned_sequence)
        ):
            if target_base != "-":
                target_index += 1
            if column < alignment_start or column >= alignment_end:
                if query_base != "-":
                    query_index += 1
                continue
            if query_base != "-":
                query_index += 1
                query_indices.append(query_index)
                if query_index == 0 and target_base != "-":
                    target_base_for_anchor_start = target_index

        if (
            not query_indices
            or min(query_indices) != 0
            or len(set(query_indices)) / len(anchor) < MIN_ANCHOR_COVERAGE
            or target_base_for_anchor_start is None
        ):
            continue

        start = target_base_for_anchor_start % len(sequence)
        candidates[(strand, start)] = {
            "start": start,
            "strand": strand,
            "score": score,
            "coverage": len(set(query_indices)) / len(anchor),
        }

    return list(candidates.values())


def find_reference_anchor(sequence, reference_fasta):
    """Find the reference base-1 assembly coordinate on either strand."""
    sequence = str(sequence).upper()
    anchor = _reference_anchor(reference_fasta)
    candidates = (
        _local_alignment_candidates(anchor, sequence, "+")
        + _local_alignment_candidates(anchor, str(Seq(sequence).reverse_complement()), "-")
    )
    candidates.sort(
        key=lambda candidate: (candidate["score"], candidate["coverage"]),
        reverse=True,
    )
    if not candidates:
        raise ValueError(
            "could not find a sufficiently complete reference 1-300 anchor"
        )
    best = candidates[0]
    ambiguous = [
        candidate for candidate in candidates
        if candidate["score"] == best["score"]
        and (candidate["strand"], candidate["start"]) != (best["strand"], best["start"])
    ]
    if ambiguous:
        raise ValueError("reference anchor is ambiguous")
    return best


def _rotate_sequence(sequence, start):
    return sequence[start:] + sequence[:start]


def _map_interval(start, end, rotation, length):
    """Map a 0-based half-open interval, splitting at the new origin."""
    if end - start == length:
        return [(0, length)]
    if rotation == 0:
        return [(start, end)]
    if start < rotation < end:
        return [(start - rotation + length, length), (0, end - rotation)]

    new_start = (start - rotation) % length
    return [(new_start, new_start + end - start)]


def rotate_gff(input_gff, output_gff, rotation, length, reverse=False):
    with open(input_gff) as source, open(output_gff, "w") as target:
        for line in source:
            if not line.strip() or line.startswith("#"):
                target.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            start, end = int(fields[3]) - 1, int(fields[4])
            if reverse:
                start, end = length - end, length - start
                fields[6] = "+" if fields[6] == "-" else "-"

            for interval_start, interval_end in _map_interval(start, end, rotation, length):
                fields[3] = str(interval_start + 1)
                fields[4] = str(interval_end)
                target.write("\t".join(fields) + "\n")


def rotate_genbank(input_gb, output_gb, rotation, reverse=False, id_suffix="_rotated"):
    record = SeqIO.read(input_gb, "genbank")
    if reverse:
        record = record.reverse_complement(id=record.id)

    shifted = record[rotation:] + record[:rotation]
    shifted.id = record.id + id_suffix
    shifted.name = record.name + id_suffix
    shifted.annotations["molecule_type"] = record.annotations.get("molecule_type", "DNA")
    with open(output_gb, "w") as handle:
        SeqIO.write(shifted, handle, "genbank")


def rotate_annotations(input_fasta, input_annotation, output_fasta, output_annotation,
                       reference_fasta, annotation_type="genbank"):
    record = SeqIO.read(input_fasta, "fasta")
    candidate = find_reference_anchor(record.seq, reference_fasta)
    sequence = str(record.seq)
    reverse = candidate["strand"] == "-"
    if reverse:
        sequence = str(Seq(sequence).reverse_complement())

    # Follow the existing pipeline naming convention (see rotation.py/rotation_mitos.py)
    # so downstream cd-hit cluster parsing (getReprContig.py) can recover the original
    # contig_id by stripping "_rotated" or "_rc_rotated" from the FASTA header.
    id_suffix = "_rc_rotated" if reverse else "_rotated"

    rotation = candidate["start"]
    record.seq = Seq(_rotate_sequence(sequence, rotation))
    record.id = record.id + id_suffix
    record.description = ""
    with open(output_fasta, "w") as handle:
        SeqIO.write(record, handle, "fasta")

    if annotation_type == "genbank":
        rotate_genbank(input_annotation, output_annotation, rotation, reverse, id_suffix)
    else:
        rotate_gff(input_annotation, output_annotation, rotation, len(sequence), reverse)

    return candidate
