from __future__ import annotations

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, field_validator

from ribo_switch.main import _candidate_to_json, _summary_to_kcal
from ribo_switch.nsga2 import nsga2, summarize_pareto_front
from ribo_switch.turner import TurnerParams

app = FastAPI(title="SMART Riboswitch Design Service")

_params = TurnerParams.turner2004()


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
        if not v or not set(v) <= {'.', '(', ')'}:
            raise ValueError("must be a non-empty dot-bracket string containing only '.', '(', ')'")
        return v


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
