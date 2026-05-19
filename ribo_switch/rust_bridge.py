"""Thin type-conversion bridge to the ribo_rs Rust extension.

ribo_rs is a required dependency — no Python fallback.
This module is used by nsga2.py and scorer.py which pass high-level
Sequence/Structure objects.  For direct ribo_rs access, import ribo_rs
directly (see brpf.py, verify.py).
"""
from __future__ import annotations
from typing import Iterable

import numpy as _np
import ribo_rs as _ribo_rs

from ribo_switch.types import Sequence, Energy, Structure
from ribo_switch.turner import TurnerParams
from ribo_switch.fold import FoldResult


def _seq_bytes(seq: Sequence) -> list[int]:
    """Convert a Sequence to the raw byte list that ribo_rs expects."""
    return [int(b) for b in seq.bases]


def eval_energy(seq: Sequence, struct: Structure, params: TurnerParams) -> Energy:
    """Evaluate the free energy of *seq* in conformation *struct*."""
    return _ribo_rs.eval_energy(_seq_bytes(seq), struct.pair_table)


def fold_mfe(seq: Sequence, params: TurnerParams) -> FoldResult:
    """Compute the minimum-free-energy fold for *seq*."""
    mfe_e, mfe_db = _ribo_rs.fold_mfe(_seq_bytes(seq))
    n = len(seq.bases)
    _dummy = _np.empty((n, n), dtype=_np.int64)
    return FoldResult(mfe_energy=mfe_e, mfe_structure=mfe_db, v=_dummy, w=_dummy)


def evaluate_candidate(
    seq: Sequence,
    s_on: Structure,
    s_off: Structure,
    params: TurnerParams,
) -> tuple[int, int, int, str]:
    """Evaluate e_on, e_off, MFE, and MFE structure for a candidate."""
    return _ribo_rs.evaluate_candidate(
        _seq_bytes(seq), s_on.pair_table, s_off.pair_table
    )


def evaluate_batch(
    seqs: Iterable[Sequence],
    s_on: Structure,
    s_off: Structure,
    params: TurnerParams,
) -> list[tuple[int, int, int, str]]:
    """Batch-evaluate candidates using Rayon parallelism."""
    seq_lists = [_seq_bytes(s) for s in seqs]
    if not seq_lists:
        return []
    return _ribo_rs.evaluate_batch(seq_lists, s_on.pair_table, s_off.pair_table)
