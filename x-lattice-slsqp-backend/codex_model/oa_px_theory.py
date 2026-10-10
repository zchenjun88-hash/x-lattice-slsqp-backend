"""OA-px plate-fin + X-lattice theory model.

The equations below are a direct transcription of the OA-px worksheet in
20260904X-lattice夹芯板三点弯(1).xlsx. Only OA-px is implemented; no SQP
optimization is executed in this stage.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import cos, pi, radians, sin, sqrt
from typing import Dict, Tuple

NormalizedThicknesses = Tuple[float, float, float]


@dataclass(frozen=True)
class OAPXGeometry:
    """Fixed OA-px values copied from the workbook baseline row."""

    L: float = 98.0
    H: float = 14.0
    B: float = 36.0
    LC: float = 14.0
    WP: float = 12.0
    Hc: float = 9.66
    omega: float = 52.3
    beta: float = 60.0
    l: float = 10.8
    w: float = 2.0
    b: float = 2.5
    c: float = 1.5
    kfw: float = 1.42
    kx: float = 2.0
    kp: float = 4.0
    a: float = 5.0
    Es: float = 163000.0
    sigma_y: float = 489.0
    nu: float = 0.3
    density_kg_m3: float = 7900.0


DEFAULT_GEOMETRY = OAPXGeometry()


def _thicknesses(x: NormalizedThicknesses, g: OAPXGeometry):
    if len(x) != 3:
        raise ValueError("x must be (tbar_f, tbar_x, tbar_p)")
    return tuple(float(v) * g.Hc for v in x)


def _core_relative_densities(tx: float, tp: float, g: OAPXGeometry):
    # Q5, S5 and I5 in the OA-px worksheet.
    rho_x = tx * g.w * (4 * g.b + 4 * g.l) / (g.LC * g.WP * g.Hc)
    rho_p = 2 * tp * g.WP * g.Hc / (g.LC * g.WP * g.Hc)
    rho_px = rho_x + rho_p * (1 - rho_x)
    return rho_x, rho_p, rho_px


def _intermediates(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY) -> Dict[str, float]:
    tf, tx, tp = _thicknesses(x, g)
    rho_x, rho_p, rho_px = _core_relative_densities(tx, tp, g)

    # W5 and X5.
    Eeq = g.Es / (1 - g.nu**2)
    sy_eq = 2 * g.sigma_y / sqrt(3)
    om = radians(g.omega)
    be = radians(g.beta)
    co = cos(om) * cos(be)

    # Y5, Z5, AA5, AB5, AC5, AD5.
    ka = 2 * g.Es * g.w * tx / g.l
    kb = 12 * g.Es / (g.l / 2) ** 3 * ((g.w * tp**3 + tp * g.w**3) / 24)
    kd = ka - kb
    C55x = rho_x * ((g.l / 2) / (4 * g.w * tx)) * (
        kb * (co**2 + sin(om) ** 2) + 4 * kd * co**2 * sin(be) ** 2
    )
    C55p = rho_p * Eeq * (tp / g.Hc) ** 2
    C55star = C55x + C55p

    # AE5 and AF5.
    EIeq = g.Es * g.B * tf * (g.Hc + tf) ** 2 / 2
    GAeq = g.B * (g.Hc + tf) ** 2 * C55star / g.Hc

    # AG5: face yielding, kN.
    Ffy = 4 * (g.Hc + tf) * g.B * tf * sy_eq / g.L / 1000

    # AJ5: face wrinkling, kN; AH5 is the unsupported length LC/4.
    unsupported = g.LC / 4
    sigma_fw = g.kfw**2 * pi**2 * Eeq * tf**2 / (12 * unsupported**2)
    Ffw = 4 * (g.Hc + tf) * g.B * tf / g.L * sigma_fw / 1000

    # OA uses AQ5 / model A. AO5 is explicitly zero in the worksheet.
    eta = g.l / (g.l + g.b)
    sigma_crx = min(g.kx**2 * pi**2 * g.Es / 12 * (tx / (0.5 * g.l)) ** 2, g.sigma_y)
    Sigma13x = rho_x * eta * sigma_crx * sin(om) * (
        (2 * g.l * cos(om) * sin(be) + 2 * g.b) / (2 * g.l)
    )
    Sigma13star = Sigma13x
    Fcs = (2 * g.B * tf**2 * sy_eq / g.L + 2 * g.B * g.Hc * Sigma13star * (1 + 2 * g.H / g.L)) / 1000

    # AT5, AV5, AW5, AY5: indentation.
    Sigma33x = rho_x * eta * sigma_crx * sin(om) ** 2
    sigma_33p_member = min(g.kp * pi**2 * Eeq / 12 * (tp / g.WP) ** 2, g.sigma_y)
    Sigma33p = sigma_33p_member * rho_p
    Sigma33star = Sigma33x + (1 - rho_x) * Sigma33p
    Find = (2 * g.B * tf * sqrt(sy_eq * Sigma33star) + g.a * g.B * Sigma33star) / 1000

    return {
        "tf": tf, "tx": tx, "tp": tp,
        "rho_x": rho_x, "rho_p": rho_p, "rho_px": rho_px,
        "Eeq": Eeq, "sy_eq": sy_eq,
        "ka": ka, "kb": kb, "kd": kd,
        "C55x": C55x, "C55p": C55p, "C55star": C55star,
        "EIeq": EIeq, "GAeq": GAeq,
        "Ffy_kN": Ffy, "Ffw_kN": Ffw,
        "eta": eta, "sigma_crx": sigma_crx,
        "Sigma13x": Sigma13x, "Sigma13star": Sigma13star,
        "Fcs_kN": Fcs,
        "Sigma33x": Sigma33x, "sigma_33p_member": sigma_33p_member,
        "Sigma33p": Sigma33p, "Sigma33star": Sigma33star,
        "Find_kN": Find,
    }


def mass(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY) -> float:
    """Actual OA-px material mass in kg: two faces + hybrid core."""
    v = _intermediates(x, g)
    solid_mm3 = 2 * g.L * g.B * v["tf"] + v["rho_px"] * g.L * g.B * g.Hc
    return g.density_kg_m3 * 1e-9 * solid_mm3


def dimensionless_mass(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY) -> float:
    """m/(rho_s*L*B*Hc) = 2*tbar_f + rho_px."""
    v = _intermediates(x, g)
    return 2 * v["tf"] / g.Hc + v["rho_px"]


def FFY(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY) -> float:
    """Face-sheet yielding load, N."""
    return 1000 * _intermediates(x, g)["Ffy_kN"]


def FFW(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY) -> float:
    """Face-sheet wrinkling load, N."""
    return 1000 * _intermediates(x, g)["Ffw_kN"]


def FCS(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY) -> float:
    """OA core-shear load, model A / workbook PA, N."""
    return 1000 * _intermediates(x, g)["Fcs_kN"]


def FIND(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY) -> float:
    """Indentation load, N."""
    return 1000 * _intermediates(x, g)["Find_kN"]


def baseline_design(g: OAPXGeometry = DEFAULT_GEOMETRY) -> NormalizedThicknesses:
    return 0.5 / g.Hc, 0.8 / g.Hc, 0.8 / g.Hc


def all_responses(x: NormalizedThicknesses, g: OAPXGeometry = DEFAULT_GEOMETRY):
    return {
        "m_kg": mass(x, g), "mbar": dimensionless_mass(x, g),
        "FY_N": FFY(x, g), "FW_N": FFW(x, g),
        "CS_N": FCS(x, g), "IND_N": FIND(x, g),
    }


def geometry_as_dict(g: OAPXGeometry = DEFAULT_GEOMETRY):
    return asdict(g)


__all__ = ["OAPXGeometry", "DEFAULT_GEOMETRY", "mass", "dimensionless_mass", "FFY", "FFW", "FCS", "FIND", "baseline_design", "all_responses", "geometry_as_dict"]
