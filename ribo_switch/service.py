from __future__ import annotations

import base64
from typing import Literal

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, field_validator

from ribo_switch.main import _candidate_to_json, _summary_to_kcal
from ribo_switch.nsga2 import nsga2, summarize_pareto_front
from ribo_switch.turner import TurnerParams

app = FastAPI(title="SMART Riboswitch Design Service")

_params = TurnerParams.turner2004()

_DOT_BRACKET = frozenset(".(){}")


def _validate_dot_bracket(v: str) -> str:
    if not v or not set(v) <= {'.', '(', ')'}:
        raise ValueError("must be a non-empty dot-bracket string containing only '.', '(', ')'")
    return v


# ── /design ────────────────────────────────────────────────────────────────

class DesignRequest(BaseModel):
    structure_on: str
    structure_off: str
    population: int = 100
    generations: int = 200
    mutation_rate: float = 0.1
    seed: int | None = None
    bp_distance_obj: bool = False

    @field_validator("structure_on", "structure_off")
    @classmethod
    def must_be_dot_bracket(cls, v: str) -> str:
        return _validate_dot_bracket(v)


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/design")
def design(req: DesignRequest):
    if len(req.structure_on) != len(req.structure_off):
        raise HTTPException(
            status_code=400,
            detail=(
                f"structure_on and structure_off must have equal length "
                f"({len(req.structure_on)} vs {len(req.structure_off)})"
            ),
        )
    try:
        front = nsga2(
            structure_on=req.structure_on,
            structure_off=req.structure_off,
            population_size=req.population,
            n_generations=req.generations,
            mutation_rate=req.mutation_rate,
            params=_params,
            seed=req.seed,
            include_structure_objective=req.bp_distance_obj,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))

    if not front:
        raise HTTPException(status_code=500, detail="NSGA-II returned no candidates")

    summary = summarize_pareto_front(front)
    return {
        "structure_on": req.structure_on,
        "structure_off": req.structure_off,
        "units": {"energy": "kcal/mol"},
        "summary": _summary_to_kcal(summary),
        "candidates": [_candidate_to_json(c) for c in front],
    }


# ── /visualize ─────────────────────────────────────────────────────────────

class VisualizeRequest(BaseModel):
    structure_on: str
    structure_off: str
    sequence: str | None = None
    graph_mode: Literal["combined", "on", "off"] = "combined"
    algorithm: Literal["naview", "radiate", "circular", "line"] = "radiate"

    @field_validator("structure_on", "structure_off")
    @classmethod
    def must_be_dot_bracket(cls, v: str) -> str:
        return _validate_dot_bracket(v)


def _b64(data: bytes) -> str:
    return base64.b64encode(data).decode()


@app.post("/visualize")
def visualize(req: VisualizeRequest):
    if len(req.structure_on) != len(req.structure_off):
        raise HTTPException(
            status_code=400,
            detail=(
                f"structure_on and structure_off must have equal length "
                f"({len(req.structure_on)} vs {len(req.structure_off)})"
            ),
        )
    if req.sequence is not None and len(req.sequence) != len(req.structure_on):
        raise HTTPException(
            status_code=400,
            detail=(
                f"sequence length ({len(req.sequence)}) must match "
                f"structure length ({len(req.structure_on)})"
            ),
        )

    try:
        from ribo_switch.viz import render_arc, render_circos, render_varna
    except ImportError as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Visualization dependencies not installed: {exc}. "
                   "Install with: pip install '.[viz]'",
        )

    try:
        arc_png = render_arc(req.structure_on, req.structure_off, req.graph_mode)
        circos_png = render_circos(req.structure_on, req.structure_off, req.graph_mode)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Arc/circos render failed: {exc}")

    varna_on = varna_off = None
    if req.sequence is not None:
        try:
            varna_on = _b64(render_varna(
                req.sequence, req.structure_on,
                title="S_ON", algorithm=req.algorithm,
            ))
            varna_off = _b64(render_varna(
                req.sequence, req.structure_off,
                title="S_OFF", algorithm=req.algorithm,
            ))
        except Exception as exc:
            # VARNA not configured, JAR missing, or subprocess timed out —
            # return nulls gracefully rather than 500
            varna_on = varna_off = None

    return {
        "arc":      _b64(arc_png),
        "circos":   _b64(circos_png),
        "varna_on":  varna_on,
        "varna_off": varna_off,
    }
