"""
UMI utility functions for standalone script execution.

This module provides UMI (Unique Molecular Identifier) processing utilities,
including UMI clustering, deduplication, and analysis for single-cell
RNA sequencing data processing.
"""

import itertools
from collections import Counter


def cluster_umis(umis, threshold=1):
    """
    Deterministically map low-count UMI neighbors to higher-count parents.
    """
    if not umis:
        return {}

    counts = Counter()
    if isinstance(umis, (dict, Counter)):
        counts.update(umis)
    else:
        counts.update(umis)

    unique_umis = sorted(counts.keys(), key=lambda x: (-counts[x], x))
    parent_map = {u: u for u in unique_umis}

    if threshold == 0:
        return parent_map

    umi_set = set(unique_umis)
    visited = set()
    bases = ("A", "C", "G", "T", "N")

    for parent in unique_umis:
        if parent in visited:
            continue

        parent_len = len(parent)
        p_chars = list(parent)

        for dist in range(1, threshold + 1):
            for positions in itertools.combinations(range(parent_len), dist):
                originals = [p_chars[pos] for pos in positions]
                replacements = []
                for orig in originals:
                    replacements.append([b for b in bases if b != orig])

                for repl_tuple in itertools.product(*replacements):
                    for pos, repl in zip(positions, repl_tuple):
                        p_chars[pos] = repl

                    child = "".join(p_chars)
                    if child in umi_set and child not in visited:
                        parent_map[child] = parent
                        visited.add(child)

                    for pos, orig in zip(positions, originals):
                        p_chars[pos] = orig

    return parent_map
