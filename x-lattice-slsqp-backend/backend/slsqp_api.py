"""HTTP API for true SLSQP optimization of the Excel OA-px X-lattice model.

The browser-facing contract accepts physical plate/member thicknesses in mm and
returns both physical responses and the dimensionless objective used by the
workbook reproduction.  The optimizer itself works in normalized thicknesses
x=(tf/Hc, tx/Hc, tp/Hc), with scipy.optimize.minimize(method="SLSQP").

Run locally:
    uvicorn backend.slsqp_api:app --host 0.0.0.0 --port 8000

A static Sites deployment cannot execute this Python module.  Deploy it as a
separate FastAPI service (Render/Railway/Fly/Cloud Run/Azure Functions with an
ASGI adapter) and set the page's API base URL to that service.  CORS is open by
default so a static page can call it; tighten ALLOW_ORIGINS in production.
"""
from __future__ import annotations

import math
import os
import hashlib
import json
import re
from dataclasses import asdict, fields
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, ConfigDict, Field, field_validator
from scipy.optimize import minimize

# Import the single source of truth for OA-px/Excel formulas.
try:  # package import when launched as `uvicorn backend.slsqp_api:app`
    from codex_model.oa_px_theory import (
        DEFAULT_GEOMETRY,
        FCS,
        FFW,
        FFY,
        FIND,
        OAPXGeometry,
        all_responses,
        dimensionless_mass,
        mass,
    )
except ImportError:  # direct invocation from backend/
    from oa_px_theory import (  # type: ignore
        DEFAULT_GEOMETRY,
        FCS,
        FFW,
        FFY,
        FIND,
        OAPXGeometry,
        all_responses,
        dimensionless_mass,
        mass,
    )

F0_DEFAULT_N = 4214.8
VAR_NAMES: Tuple[str, ...] = ("tf", "tx", "tp")
RESPONSE_NAMES: Tuple[str, ...] = ("FY", "FW", "CS", "IND")
RESPONSE_FUNCS = {"FY": FFY, "FW": FFW, "CS": FCS, "IND": FIND}
DEFAULT_LOWER_MM = {"tf": 0.04 * DEFAULT_GEOMETRY.Hc, "tx": 0.05 * DEFAULT_GEOMETRY.Hc, "tp": 0.05 * DEFAULT_GEOMETRY.Hc}

# These are the eight editable formula fields in index.html.  The backend never
# executes formula text as Python/JavaScript.  It validates the expected
# equation symbols and extracts only an explicit numeric multiplier, which makes
# edits such as ``... × 2`` deterministic while rejecting malformed input.
FORMULA_DEFAULTS: Dict[str, str] = {
    "eq1": "rho_mix = rho_x + rho_p(1-rho_x)",
    "eq2": "rho_x = tx*w*(l+b)/(LC*WP*Hc)",
    "eq3": "rho_p = 2*tp*WP*Hc/(LC*WP*Hc)",
    "eq45": "CijH = Cijf + Cijp",
    "eq68": "delta = P*L^3/(48*EIeq) + P*L/(4*GAeq)",
    "eq911": "Pfy = 4*(Hc+tf)*B*tf*sigma_y/L; Pfw = 4*(Hc+tf)*B*tf*sigma_fw/L",
    "eq1219": "PA = 2*B*tf^2*sigma_y/L + 2*B*H*Sigma13star*(1+2*H/L)",
    "eq2024": "Pind = 2*B*tf*sqrt(sigma_y*Sigma33star) + a*B*Sigma33star",
}
FORMULA_ALIASES = {
    **{key: key for key in FORMULA_DEFAULTS},
    **{f"formulaEq{key[2:]}": key for key in FORMULA_DEFAULTS},
    **{f"formulaeq{key[2:]}": key for key in FORMULA_DEFAULTS},
}
FORMULA_REQUIREMENTS = {
    "eq1": ("rho", "=",),
    "eq2": ("rho", "tx", "l", "=",),
    "eq3": ("rho", "tp", "=",),
    "eq45": ("C", "="),
    "eq68": ("EI", "GA", "="),
    "eq911": ("Pfy", "Pfw", "="),
    "eq1219": ("PA", "Sigma", "="),
    "eq2024": ("Pind", "Sigma", "="),
}
FORMULA_ALLOWED_ASCII = {
    "eq1": {"rho", "rho_mix", "rho_x", "rho_p", "x", "p"},
    "eq2": {"rho", "rho_x", "tx", "x", "w", "l", "b", "c", "sin", "cos", "omega", "beta", "LC", "WP", "Hc"},
    "eq3": {"rho", "rho_p", "tp", "p", "WP", "Hc", "LC", "c", "cos", "l"},
    "eq45": {"C", "CijH", "Cijf", "Cijp", "CH"},
    "eq68": {"delta", "P", "L", "EI", "EIeq", "GA", "GAeq", "B", "Hc", "tf", "Es", "Gi", "PL", "eq", "i"},
    "eq911": {"Pfy", "Pfw", "Hc", "tf", "B", "sigma", "sigma_y", "sigma_fw", "L", "kfw", "pi", "Eeq", "L0", "OA", "OB", "Wp", "l", "c", "omega", "beta", "i", "cos", "fw", "s", "y"},
    "eq1219": {"PA", "PB", "B", "tf", "sigma", "sigma_y", "L", "H", "Sigma", "Sigma13star", "Sigmaix", "Sigmaip", "sigma_crx", "kx", "Esx", "tx", "l", "rho", "rho_x", "rho_p", "eta", "sin", "cos", "omega", "beta", "kp", "kpf", "Esp", "tp", "Wp", "Hc", "nup", "sigmay", "c", "i", "x", "p", "y", "yp", "yx", "BH", "Btf", "crp", "crx", "ip", "ix", "b", "OA", "OB"},
    "eq2024": {"Pind", "B", "tf", "sigma", "sigma_y", "Sigma", "Sigma33star", "sigma_x33", "sigma_p33", "sigma_crx", "sigma_crp", "rho", "rho_x", "rho_p", "eta", "sin", "cos", "omega", "a", "tp", "Wp", "Esp", "nup", "sigmay", "sqrt", "x", "p", "y", "yp", "Btf", "aB", "crp", "crx", "p33", "x33"},
}


def _has_symbol(text: str, requirement: str) -> bool:
    """Accept ASCII and the Greek notation used by the editable page."""
    alternatives = {
        "rho": ("rho", "ρ̄"),
        "tx": ("tx",),
        "tp": ("tp",),
        "l": ("l", "ℓ"),
        "C": ("C", "c"),
        "EI": ("EI", "(EI)"),
        "GA": ("GA", "(GA)"),
        "Pfy": ("Pfy",),
        "Pfw": ("Pfw",),
        "PA": ("PA",),
        "Pind": ("Pind",),
        "Sigma": ("Sigma", "Σ"),
        "=": ("=",),
    }
    return any(part in text for part in alternatives.get(requirement, (requirement,)))


def _formula_payload(raw: Mapping[str, str]) -> Tuple[Dict[str, str], Dict[str, float], str]:
    """Validate formula text and return canonical text, factors and SHA-256.

    Only numeric multipliers explicitly written after ``×``, ``*`` or Chinese
    ``乘以`` are applied.  There is no eval/exec path.  A bad submitted formula
    raises 422 from the endpoint; it never silently falls back to old text.
    """
    if not isinstance(raw, Mapping):
        raise ValueError("formula_texts must be an object")
    canonical: Dict[str, str] = dict(FORMULA_DEFAULTS)
    supplied = set()
    for key, value in raw.items():
        key_s = str(key)
        norm = FORMULA_ALIASES.get(key_s, key_s.lower())
        if norm not in FORMULA_DEFAULTS:
            raise ValueError(f"unknown formula field {key!r}; expected eq1, eq2, eq3, eq45, eq68, eq911, eq1219, eq2024")
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{key} must be a non-empty formula string")
        text = value.strip()
        if len(text) > 12000:
            raise ValueError(f"{key} is too long")
        for required in FORMULA_REQUIREMENTS[norm]:
            if not _has_symbol(text, required):
                raise ValueError(f"{key} is invalid: missing required symbol {required!r}")
        # Block control characters and code-like delimiters before hashing.
        if any(ord(ch) < 32 and ch not in "\n\r\t" for ch in text) or any(token in text for token in ("__", "<script", "`", "{{", "}}")):
            raise ValueError(f"{key} contains unsupported characters")
        pairs = {')': '(', ']': '[', '}': '{'}
        stack: List[str] = []
        for char in text:
            if char in '([{':
                stack.append(char)
            elif char in pairs:
                if not stack or stack.pop() != pairs[char]:
                    raise ValueError(f"{key} has unbalanced brackets")
        if stack:
            raise ValueError(f"{key} has unbalanced brackets")
        tokens = set(re.findall(r"[A-Za-z][A-Za-z0-9_]*", text))
        unknown_tokens = tokens - FORMULA_ALLOWED_ASCII[norm]
        if unknown_tokens:
            raise ValueError(f"{key} contains unsupported symbol(s): {sorted(unknown_tokens)}")
        # Every explicit '*' or '×/乘以' multiplier must be numeric.  The
        # ordinary formulas use '*' only as multiplication, so rejecting *x
        # prevents ambiguity and accidental code-like expressions.
        for match in re.finditer(r"(?:×|\*|乘以)", text):
            # In the manuscript `Σ*33` means Sigma-star-33, not "multiply by
            # 33".  Keep this notation as a symbol; explicit editable factors
            # remain `... × 2`, `... *2`, or `... 乘以2` elsewhere.
            if match.group(0) == "*" and match.start() > 0 and text[match.start() - 1] in "Σσ":
                continue
            suffix = text[match.end():].lstrip()
            # `*name` is ordinary formula multiplication.  A numeric suffix is
            # the explicit editable multiplier; an empty/punctuation suffix is
            # invalid.  Unicode × and Chinese 乘以 are reserved for numeric
            # multipliers, avoiding ambiguity.
            numeric = re.match(r"^[0-9]+(?:\.[0-9]+)?", suffix) is not None
            if match.group(0) != "*" and not numeric:
                raise ValueError(f"{key} has a non-numeric multiplier near {text[max(0, match.start()-8):match.end()+8]!r}")
            if match.group(0) == "*" and not suffix:
                raise ValueError(f"{key} has a non-numeric multiplier near {text[max(0, match.start()-8):match.end()+8]!r}")
        numbers = [float(m.group(1)) for m in re.finditer(r"(?:×|\*|乘以)\s*([0-9]+(?:\.[0-9]+)?)", text) if not (m.group(0).startswith("*") and m.start() > 0 and text[m.start() - 1] in "Σσ")]
        factor = 1.0
        for number in numbers:
            if not math.isfinite(number) or number <= 0 or number > 1e6:
                raise ValueError(f"{key} multiplier must be in (0, 1e6]")
            factor *= number
        canonical[norm] = text
        supplied.add(norm)
    # A stable version/hash lets the UI prove exactly which formula text was
    # used.  Defaults are included so omitted fields are explicit in the hash.
    serial = json.dumps(canonical, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(serial.encode("utf-8")).hexdigest()
    factors = {key: 1.0 for key in FORMULA_DEFAULTS}
    for key, text in canonical.items():
        numbers = [float(m.group(1)) for m in re.finditer(r"(?:×|\*|乘以)\s*([0-9]+(?:\.[0-9]+)?)", text) if not (m.group(0).startswith("*") and m.start() > 0 and text[m.start() - 1] in "Σσ")]
        for number in numbers:
            factors[key] *= number
    return canonical, factors, digest


class OptimizeRequest(BaseModel):
    """JSON request accepted at POST /api/slsqp.

    `initial_mm` and bounds use physical mm.  `variables` selects the optimized
    coordinates; coordinates not selected remain fixed at initial_mm.  For
    backward compatibility, `x` may be supplied as normalized thicknesses and
    is converted to mm when initial_mm is omitted.
    """

    model_config = ConfigDict(extra="forbid")
    model: str = Field(default="OA-px", description="Only OA-px is implemented")
    target_load_N: float = Field(default=F0_DEFAULT_N, gt=0)
    variables: List[str] = Field(default_factory=lambda: list(VAR_NAMES), min_length=1)
    initial_mm: Dict[str, float] = Field(default_factory=lambda: {"tf": 0.5, "tx": 0.8, "tp": 0.8})
    x: Optional[List[float]] = Field(default=None, min_length=3, max_length=3, description="Optional normalized [tf/Hc,tx/Hc,tp/Hc]")
    lower_mm: Dict[str, float] = Field(default_factory=lambda: dict(DEFAULT_LOWER_MM))
    upper_mm: Optional[Dict[str, float]] = None
    geometry: Dict[str, float] = Field(default_factory=dict)
    starts_mm: Optional[List[Dict[str, float]]] = None
    formula_texts: Dict[str, str] = Field(default_factory=dict, description="Editable formula fields from the page; aliases formulaEq1... accepted")
    max_iter: int = Field(default=200, ge=1, le=2000)
    ftol: float = Field(default=1e-9, gt=0, le=1)

    @field_validator("model")
    @classmethod
    def model_is_oapx(cls, v: str) -> str:
        if v.lower().replace("_", "-") not in {"oa-px", "oapx"}:
            raise ValueError("only the OA-px model is available")
        return "OA-px"

    @field_validator("variables")
    @classmethod
    def variables_are_known(cls, values: List[str]) -> List[str]:
        if len(set(values)) != len(values):
            raise ValueError("variables must be unique")
        unknown = set(values) - set(VAR_NAMES)
        if unknown:
            raise ValueError(f"unknown variables: {sorted(unknown)}")
        return values


def _geometry(overrides: Mapping[str, float]) -> OAPXGeometry:
    """Construct a validated geometry from Excel defaults plus finite overrides."""
    valid = {f.name for f in fields(OAPXGeometry)}
    unknown = set(overrides) - valid
    if unknown:
        raise ValueError(f"unknown geometry fields: {sorted(unknown)}")
    clean: Dict[str, Any] = {}
    for key, value in overrides.items():
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f"geometry.{key} must be finite")
        clean[key] = value
    return OAPXGeometry(**{**asdict(DEFAULT_GEOMETRY), **clean})


def _finite_map(values: Mapping[str, float], label: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for key, value in values.items():
        if key not in VAR_NAMES:
            raise ValueError(f"{label} contains unknown variable {key!r}")
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f"{label}.{key} must be finite")
        out[key] = value
    return out


def _prepare(req: OptimizeRequest, g: OAPXGeometry):
    initial = _finite_map(req.initial_mm, "initial_mm")
    if req.x is not None and not req.initial_mm:
        initial = dict(zip(VAR_NAMES, np.asarray(req.x, dtype=float) * g.Hc))
    # Pydantic's default initial_mm is populated, so explicitly supplied x wins
    # only when the caller omitted initial_mm in the raw request.  This helper is
    # called after model validation; detect x by checking whether all defaults
    # are still present and x differs from them.
    if req.x is not None:
        x_mm = dict(zip(VAR_NAMES, np.asarray(req.x, dtype=float) * g.Hc))
        if req.initial_mm == {"tf": 0.5, "tx": 0.8, "tp": 0.8}:
            initial = x_mm
    for name in VAR_NAMES:
        initial.setdefault(name, {"tf": 0.5, "tx": 0.8, "tp": 0.8}[name])
    lower = dict(DEFAULT_LOWER_MM)
    lower.update(_finite_map(req.lower_mm, "lower_mm"))
    upper = _finite_map(req.upper_mm or {}, "upper_mm")
    for name in VAR_NAMES:
        if name not in upper:
            upper[name] = math.inf
        if lower[name] <= 0:
            raise ValueError(f"lower_mm.{name} must be > 0")
        if upper[name] <= lower[name]:
            raise ValueError(f"upper_mm.{name} must exceed lower_mm.{name}")
        initial[name] = max(lower[name], min(upper[name], initial[name]))
    return initial, lower, upper


def _evaluate(x_norm: Sequence[float], g: OAPXGeometry, target: float, factors: Optional[Mapping[str, float]] = None) -> Dict[str, Any]:
    """Evaluate the OA-px Excel formulas with validated editable-text factors."""
    x = tuple(float(v) for v in x_norm)
    fm = {key: 1.0 for key in FORMULA_DEFAULTS}
    if factors:
        fm.update({key: float(value) for key, value in factors.items()})
    tf, tx, tp = (value * g.Hc for value in x)
    rho_x = tx * g.w * (4 * g.b + 4 * g.l) / (g.LC * g.WP * g.Hc) * fm["eq2"]
    rho_p = 2 * tp * g.WP * g.Hc / (g.LC * g.WP * g.Hc) * fm["eq3"]
    rho_px = (rho_x + rho_p * (1 - rho_x)) * fm["eq1"]
    Eeq = g.Es / (1 - g.nu**2)
    sy_eq = 2 * g.sigma_y / math.sqrt(3)
    om = math.radians(g.omega)
    be = math.radians(g.beta)
    co = math.cos(om) * math.cos(be)
    ka = 2 * g.Es * g.w * tx / g.l
    kb = 12 * g.Es / (g.l / 2) ** 3 * ((g.w * tp**3 + tp * g.w**3) / 24)
    kd = ka - kb
    c55x = rho_x * ((g.l / 2) / (4 * g.w * tx)) * (kb * (co**2 + math.sin(om) ** 2) + 4 * kd * co**2 * math.sin(be) ** 2)
    c55p = rho_p * Eeq * (tp / g.Hc) ** 2
    c55star = c55x + c55p
    ffy = 4 * (g.Hc + tf) * g.B * tf * sy_eq / g.L
    unsupported = g.LC / 4
    sigma_fw = g.kfw**2 * math.pi**2 * Eeq * tf**2 / (12 * unsupported**2)
    ffw = 4 * (g.Hc + tf) * g.B * tf / g.L * sigma_fw
    eta = g.l / (g.l + g.b)
    sigma_crx = min(g.kx**2 * math.pi**2 * g.Es / 12 * (tx / (0.5 * g.l)) ** 2, g.sigma_y)
    sigma13x = rho_x * eta * sigma_crx * math.sin(om) * ((2 * g.l * math.cos(om) * math.sin(be) + 2 * g.b) / (2 * g.l))
    fcs = 2 * g.B * tf**2 * sy_eq / g.L + 2 * g.B * g.Hc * sigma13x * (1 + 2 * g.H / g.L)
    sigma33x = rho_x * eta * sigma_crx * math.sin(om) ** 2
    sigma33p = min(g.kp * math.pi**2 * Eeq / 12 * (tp / g.WP) ** 2, g.sigma_y) * rho_p
    sigma33star = sigma33x + (1 - rho_x) * sigma33p
    find = 2 * g.B * tf * math.sqrt(max(0.0, sy_eq * sigma33star)) + g.a * g.B * sigma33star
    # Formula-group edits are intentionally explicit and auditable.  Eq.(9–11)
    # controls both face modes; Eq.(12–19) controls OA model-A shear; Eq.(20–24)
    # controls indentation.  Eq.(4–8) are validated and hashed but do not enter
    # these failure loads in the OA-px worksheet objective.
    responses = {
        "FY": float(ffy * fm["eq911"]),
        "FW": float(ffw * fm["eq911"]),
        "CS": float(fcs * fm["eq1219"]),
        "IND": float(find * fm["eq2024"]),
    }
    dimensionless = 2 * tf / g.Hc + rho_px
    mass_kg = g.density_kg_m3 * 1e-9 * (2 * g.L * g.B * tf + rho_px * g.L * g.B * g.Hc)
    margins = {name: responses[name] / target - 1.0 for name in RESPONSE_NAMES}
    return {
        "x_normalized": list(x),
        "x_mm": [float(v * g.Hc) for v in x],
        "objective_mbar": float(dimensionless),
        "mass_kg": float(mass_kg),
        "rhoX": float(rho_x),
        "rhoP": float(rho_p),
        "rhoMix": float(rho_px),
        "responses_N": responses,
        "Pfy": responses["FY"],
        "Pfw": responses["FW"],
        "PA": responses["CS"],
        "Pind": responses["IND"],
        "constraint_margins": margins,
        "feasible": bool(min(margins.values()) >= -1e-7),
    }


def _solve_once(start_mm: Mapping[str, float], req: OptimizeRequest, g: OAPXGeometry, lower: Mapping[str, float], upper: Mapping[str, float], factors: Mapping[str, float]) -> Dict[str, Any]:
    vars_opt = list(req.variables)
    fixed = {name: float(start_mm[name]) / g.Hc for name in VAR_NAMES}
    z0 = np.array([fixed[name] for name in vars_opt], dtype=float)
    bounds = [(float(lower[name] / g.Hc), None if math.isinf(upper[name]) else float(upper[name] / g.Hc)) for name in vars_opt]

    def full(z: Sequence[float]) -> Tuple[float, float, float]:
        out = [fixed[name] for name in VAR_NAMES]
        for idx, name in enumerate(vars_opt):
            out[VAR_NAMES.index(name)] = float(z[idx])
        return tuple(out)  # type: ignore[return-value]

    def objective(z: Sequence[float]) -> float:
        return float(_evaluate(full(z), g, req.target_load_N, factors)["objective_mbar"])

    def constraint(z: Sequence[float], name: str) -> float:
        return float(_evaluate(full(z), g, req.target_load_N, factors)["constraint_margins"][name])

    constraints = [{"type": "ineq", "fun": (lambda z, name=name: constraint(z, name))} for name in RESPONSE_NAMES]
    trace: List[Dict[str, Any]] = []
    trace.append({"iteration": 0, **_evaluate(full(z0), g, req.target_load_N, factors), "accepted": True})

    def callback(z: np.ndarray) -> None:
        trace.append({"iteration": len(trace), **_evaluate(full(z), g, req.target_load_N, factors), "accepted": True})

    result = minimize(
        objective,
        z0,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        callback=callback,
        options={"maxiter": int(req.max_iter), "ftol": float(req.ftol), "disp": False},
    )
    final_eval = _evaluate(full(result.x), g, req.target_load_N, factors)
    # `success` from SciPy can be optimistic near a violated constraint; expose
    # a strict feasibility check so the UI never labels an infeasible point pass.
    success = bool(result.success and final_eval["feasible"])
    active = [name for name, value in final_eval["constraint_margins"].items() if abs(value) <= 1e-5]
    for name in vars_opt:
        margin = (final_eval["x_mm"][VAR_NAMES.index(name)] - lower[name]) / g.Hc
        if abs(margin) <= 1e-5:
            active.append(f"{name}_lb")
    return {
        "success": success,
        "scipy_success": bool(result.success),
        "message": str(result.message),
        "iterations": len(trace) - 1,
        "variables": vars_opt,
        "active_constraints": active,
        "trace": trace,
        **final_eval,
    }


def optimize(req: OptimizeRequest) -> Dict[str, Any]:
    formula_texts, formula_factors, formula_hash = _formula_payload(req.formula_texts)
    g = _geometry(req.geometry)
    initial, lower, upper = _prepare(req, g)
    starts = req.starts_mm or [initial]
    normalized_starts: List[Dict[str, float]] = []
    for start in starts:
        merged = dict(initial)
        merged.update(_finite_map(start, "starts_mm"))
        for name in VAR_NAMES:
            merged[name] = max(lower[name], min(upper[name], merged[name]))
        normalized_starts.append(merged)
    runs = [_solve_once(start, req, g, lower, upper, formula_factors) for start in normalized_starts]
    feasible = [r for r in runs if r["success"]]
    best = min(feasible or runs, key=lambda r: r["objective_mbar"])
    return {
        "solver": "scipy.optimize.minimize(method='SLSQP')",
        "model": "OA-px",
        "cs_model": "A",
        "formula_version": "oa-px-editable-v1",
        "formula_hash": formula_hash,
        "formula_texts": formula_texts,
        "formula_multipliers": formula_factors,
        "target_load_N": req.target_load_N,
        "geometry": asdict(g),
        "bounds_mm": {"lower": lower, "upper": {k: (None if math.isinf(v) else v) for k, v in upper.items()}},
        "success": bool(best["success"]),
        "message": "best feasible multistart result" if feasible else "no feasible multistart result",
        "best": best,
        "runs": runs,
    }


app = FastAPI(title="X-Lattice OA-px SLSQP API", version="1.0.0")
origins = [origin.strip() for origin in os.getenv("ALLOW_ORIGINS", "*").split(",") if origin.strip()]
app.add_middleware(CORSMiddleware, allow_origins=origins, allow_methods=["*"], allow_headers=["*"], allow_credentials=False)


@app.get("/api/health")
def health() -> Dict[str, Any]:
    return {"ok": True, "solver": "scipy.optimize.minimize(method='SLSQP')", "model": "OA-px", "cs_model": "A"}


@app.post("/api/slsqp")
def run_slsqp(req: OptimizeRequest) -> Dict[str, Any]:
    try:
        return optimize(req)
    except (ValueError, TypeError, FloatingPointError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:  # return an actionable API error rather than HTML/traceback
        raise HTTPException(status_code=500, detail=f"optimizer failure: {exc}") from exc
