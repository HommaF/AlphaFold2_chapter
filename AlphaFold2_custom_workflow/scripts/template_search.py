#!/usr/bin/env python

import collections
import dataclasses
import itertools
import re
import sys
import string
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Set


msa_for_templates = sys.argv[1]

mgnify_max_hits = 501
uniref_max_hits = 10000


def _keep_line(line: str, seqnames: Set[str]) -> bool:
    """Function to decide which lines to keep."""

    if not line.strip():
        return True
    if line.strip() == '//': # End tag
        return True
    if line.startswith('# STOCKHOLM'): # Start tag
        return True
    if line.startswith('#=GC RF'): # Reference Annotation Line
        return True
    if line[:4] == '#=GS': # Description lines - keep if sequence in list
        _, seqname, _ = line.split(maxsplit=2)
        return seqname in seqnames
    elif line.startswith('#'): # Other markup - filter out
        return False
    else: # Alignment data - keep if sequence in list.
        seqname = line.partition(' ')[0]
        return seqname in seqnames



def truncate_stockholm_msa(stockholm_msa: str, max_sequences: int) -> str:
    """Truncates a stockholm file to a maximum number of sequences."""
    seqnames = set()
    filtered_lines = []
    for line in stockholm_msa.splitlines():
        if line.strip() and not line.startswith(('#', '//')):
            # Ignore blank lines, markup and end symbols - remainder are alignment
            # sequence parts.
            seqname = line.partition(' ')[0]
            seqnames.add(seqname)
            if len(seqnames) >= max_sequences:
                break

    for line in stockholm_msa.splitlines():
        if _keep_line(line, seqnames):
            filtered_lines.append(line)

    return '\n'.join(filtered_lines) + '\n'




def deduplicate_stockholm_msa(stockholm_msa: str) -> str:
    """Remove duplicate sequences (ignoring insertions wrt query)."""
    sequence_dict = collections.defaultdict(str)

    # First we must extract all sequences from the MSA.
    for line in stockholm_msa.splitlines():
        # Only consider the alignments - ignore reference annotation, empty lines, descriptions or markup
        if line.strip() and not line.startswith(('#', '//')):
            line = line.strip()
            seqname, alignment = line.split()
            sequence_dict[seqname] += alignment

    seen_sequences = set()
    seqnames = set()
    # First alignment is the query
    query_align = next(iter(sequence_dict.values()))
    mask = [c != '_' for c in query_align] # Mask is False for insertions
    for seqname, alignment in sequence_dict.items():
        # Apply mask to remove all insertions from the string.
        masked_alignment = ''.join(itertools.compress(alignment, mask))
        if masked_alignment in seen_sequences:
            continue
        else:
            seen_sequences.add(masked_alignment)
            seqnames.add(seqname)

    filtered_lines = []
    for line in stockholm_msa.splitlines():
        if _keep_line(line, seqnames):
            filtered_lines.append(line)

    return '\n'.join(filtered_lines) + '\n'




def remove_empty_columns_from_stockholm_msa(stockholm_msa: str) -> str:
    """Removes empty columns (dashes-only) from a Stockholm MSA."""

    processed_lines = {}
    unprocessed_lines = {}

    for i, line in enumerate(stockholm_msa.splitlines()):
        if line.startswith('#=GC RF'):
            reference_annotation_i = i
            reference_annotation_line = line
            # Reached the end of this chunk of the alignment. Process chunk.

            _, _, first_alignment = line.rpartition(' ')
            mask = []
            for j in range(len(first_alignment)):
                for _, unprocessed_line in unprocessed_lines.items():
                    prefix, _, alignment = unprocessed_line.rpartition(' ')
                    if alignment[j] != '-':
                        mask.append(True)
                        break
                    else: # Every row contained a hyphen - empty column.
                        mask.append(False)

            # Add reference annotation for processing with mask.
            unprocessed_lines[reference_annotation_i] = reference_annotation_line

            if not any(mask): # All columns were empty. Output empty lines for chunk.
                for line_index in unprocessed_lines:
                    processed_lines[line_index] = ''
            else:
                for line_index, unprocessed_line in unprocessed_lines.items():
                    prefix, _, alignment = unprocessed_line.rpartition(' ')
                    masked_alignment = ''.join(itertools.compress(alignment, mask))
                    processed_lines[line_index] = f'{prefix} {masked_alignment}'

            # Clear raw_alignments
            unprocessed_lines = {}

        elif line.strip() and not line.startswith(('#', '//')):
            unprocessed_lines[i] = line
        else:
            processed_lines[i] = line
    return '\n'.join((processed_lines[i] for i in range(len(processed_lines))))
        


msa_for templates = truncate_stockholm_msa(msa_for_templates, max_sequences=uniref.max_hits)
msa_for_templates = deduplicate_stockholm_msa(msa_for_templates)
msa_for_templates = remove_empty_columns_from_stockholm_msa(msa_for_templates)
