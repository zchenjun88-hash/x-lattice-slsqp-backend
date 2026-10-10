# X-Lattice true SLSQP backend

`index.html`/Sites is a static app (`.openai/hosting.json` points at `dist`), so it cannot run Python/SciPy. This backend is a separate ASGI service that exposes the same Excel OA-px formulas and calls `scipy.optimize.minimize(..., method="SLSQP")`.

## Run

From the project root:

```bash
python -m pip install -r backend/requirements.txt
uvicorn backend.slsqp_api:app --host 0.0.0.0 --port 8000
```

The browser can call `POST https://<service>/api/slsqp`; CORS is open by default (`ALLOW_ORIGINS=https://<site>` is recommended in production).

## Contract

`POST /api/slsqp` JSON:

```json
{
  "model": "OA-px",
  "target_load_N": 4214.8,
  "variables": ["tf", "tx", "tp"],
  "initial_mm": {"tf": 0.5, "tx": 0.8, "tp": 0.8},
  "lower_mm": {"tf": 0.3864, "tx": 0.483, "tp": 0.483},
  "upper_mm": {"tf": 3.0, "tx": 3.0, "tp": 3.0},
  "geometry": {},
  "formula_texts": {
    "formulaEq1": "ρ̄ = ρ̄x + ρ̄p(1 − ρ̄x)",
    "formulaEq2": "ρ̄x = tx·w·(l+b) / [l·sinω·(l·cosω·sinβ+b)·(l·cosω·cosβ+c)]",
    "formulaEq3": "ρ̄p = tp / (l·cosω·cosβ+c)"
  },
  "max_iter": 200,
  "ftol": 1e-9
}
```

All thicknesses are physical mm. `variables` may be a subset; unselected coordinates stay fixed. `starts_mm` optionally supplies multistart dictionaries. `x` (normalized `[tf/Hc,tx/Hc,tp/Hc]`) is accepted when `initial_mm` is omitted.

The response reports `best.objective_mbar`, `best.mass_kg`, `best.x_mm`, `best.responses_N` (`FY`, `FW`, `CS`, `IND`), normalized `best.constraint_margins`, `best.active_constraints`, `best.trace`, and all multistart `runs`. `formula_texts` accepts all eight editable fields (`formulaEq1`, `formulaEq2`, `formulaEq3`, `formulaEq45`, `formulaEq68`, `formulaEq911`, `formulaEq1219`, `formulaEq2024`) or short keys (`eq1`, ...). The parser is a safe, non-evaluating validator. Explicit numeric factors (`×2`, `*2`, `乘以2`) are applied to the density and failure groups that feed this OA-px mass/constraint objective (eq1–3, eq911, eq1219, eq2024); eq45/eq68 are validated and hashed for traceability but their stiffness outputs are not used by the worksheet failure-load objective. The response returns canonical `formula_texts`, `formula_multipliers`, `formula_version`, and a SHA-256 `formula_hash`. Invalid/unknown symbols return HTTP 422 and never fall back to an earlier formula.

`GET /api/health` is a cheap readiness probe.
