"""
lm1_loading.py
==============
LM1 vertical load and braking force calculations for a buried box culvert.

Metre-strip analysis — 1 m strip along the barrel axis (LL direction).
Vehicles travel transversely across the culvert span (B_ext direction).

References
----------
BS EN 1991-2:2003  Cl. 4.3 (LM1), Cl. 4.4.1 (braking)
UK NA to BS EN 1991-2 (alpha = 1.0 for new structures)
PD6694-1  load dispersion through fill
"""

import math
from dataclasses import dataclass, field
from typing import List

# ── LM1 characteristic values — BS EN 1991-2 Table 4.2 ───────────────────────
# Q_ik  tandem axle load (kN) — 2 axles per lane, each axle has 2 wheels
# q_ik  UDL intensity (kN/m²)
LM1_CHAR: dict = {
    1: {"Q_ik": 300.0, "q_ik": 9.0},
    2: {"Q_ik": 200.0, "q_ik": 2.5},
    3: {"Q_ik": 100.0, "q_ik": 2.5},
    4: {"Q_ik":   0.0, "q_ik": 2.5},   # no tandem in lane 4
}

# α adjustment factors — UK NA, new structures
ALPHA_Q: dict = {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}
ALPHA_q: dict = {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}

# Tandem system (TS) geometry
AXLE_SPACING  = 1.2   # m  between the two tandem axles (travel / B_ext direction)
WHEEL_SPACING = 2.0   # m  wheel centre to wheel centre (transverse / LL direction)
CONTACT_L     = 0.4   # m  wheel contact patch in travel direction
CONTACT_T     = 0.4   # m  wheel contact patch in transverse direction

# Load dispersion through fill — 1 horizontal : 1 vertical (45°)
DISP = 1.0


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class DispersionGeometry:
    """Dispersed footprint of the tandem system at crown level."""
    H_c:         float   # cover depth (m)
    disp_LL:     float   # dispersed width in LL direction (m)
    disp_Bext:   float   # dispersed length in B_ext direction (m)
    Bext_loaded: float   # effective loaded length in B_ext (capped at B_ext)
    axles_merged: bool   # True if two axles merge at crown level


@dataclass
class LaneResult:
    lane:        int
    Q_ik:        float   # design axle load (kN)   [characteristic × alpha]
    q_ik:        float   # design UDL (kN/m²)       [characteristic × alpha]
    udl_per_m:   float   # UDL vertical load on culvert per metre strip (kN/m)
    ts_per_m:    float   # TS vertical load on culvert per metre strip (kN/m)
    total_per_m: float   # combined (kN/m)


@dataclass
class BrakingResult:
    """
    Horizontal braking/acceleration force from Lane 1 (governs).
    Applied at road surface level in the direction of vehicle travel (B_ext).
    """
    Q_lk_raw:   float   # formula result before clamping (kN)
    Q_lk:       float   # clamped to [180, 900] kN per UK NA
    Q_lk_per_m: float   # per metre strip in LL direction (kN/m)


@dataclass
class LM1Result:
    n_lanes:      int
    B_ext:        float
    LL:           float
    H_c:          float
    lane_width:   float
    dispersion:   DispersionGeometry
    lanes:        List[LaneResult]
    braking:      BrakingResult
    max_V_per_m:  float   # maximum total vertical LL per metre strip (kN/m)
    min_V_per_m:  float   # minimum = no live load = 0 kN/m


# ── Core calculation ──────────────────────────────────────────────────────────

def _dispersion(H_c: float, B_ext: float) -> DispersionGeometry:
    """
    Compute dispersed tandem footprint at crown level.

    Transverse (LL) direction — wheels side by side:
        disp_LL = WHEEL_SPACING + CONTACT_T + 2 × DISP × H_c

    Travel (B_ext) direction — two axles fore-aft:
        If dispersion from each axle reaches the other at crown level they
        merge into one block; otherwise use outer-edge to outer-edge extent.
        disp_Bext = AXLE_SPACING + CONTACT_L + 2 × DISP × H_c  (in both cases)
    """
    spread = DISP * H_c

    disp_LL   = WHEEL_SPACING + CONTACT_T + 2 * spread
    disp_Bext = AXLE_SPACING  + CONTACT_L + 2 * spread
    merged    = (AXLE_SPACING - 2 * spread) <= 0

    Bext_loaded = min(disp_Bext, B_ext)

    return DispersionGeometry(
        H_c=H_c,
        disp_LL=disp_LL,
        disp_Bext=disp_Bext,
        Bext_loaded=Bext_loaded,
        axles_merged=merged,
    )


def _ts_per_m(Q_ik: float, dg: DispersionGeometry) -> float:
    """
    Tandem vertical load on culvert per metre strip (kN/m) after dispersion.

    Total tandem force landing on culvert = 2 × Q_ik × (Bext_loaded / disp_Bext).
    This is spread uniformly over (Bext_loaded × disp_LL) at crown level.
    Load per metre strip (1 m in LL) = total_on_culvert / disp_LL.
    """
    if Q_ik <= 0 or dg.disp_LL <= 0:
        return 0.0
    travel_factor = dg.Bext_loaded / dg.disp_Bext if dg.disp_Bext > 0 else 1.0
    return 2 * Q_ik * travel_factor / dg.disp_LL


def compute(
    B_ext:      float,
    LL:         float,
    H_c:        float,
    lane_width: float,
    n_lanes:    int,
) -> LM1Result:
    """
    Compute LM1 vertical loads and braking for a buried box culvert.

    Parameters
    ----------
    B_ext      : External culvert width in vehicle travel direction (m)
    LL         : Barrel length perpendicular to travel (m)
    H_c        : Cover depth — road surface to crown of culvert (m)
    lane_width : Width of each notional lane (m)
    n_lanes    : Number of notional lanes to load (max 4)

    Returns
    -------
    LM1Result with per-lane breakdown, braking, max/min totals
    """
    n  = min(max(1, n_lanes), 4)
    dg = _dispersion(H_c, B_ext)

    lane_results: List[LaneResult] = []
    max_V = 0.0

    for ln in range(1, n + 1):
        Q_ik = LM1_CHAR[ln]["Q_ik"] * ALPHA_Q[ln]
        q_ik = LM1_CHAR[ln]["q_ik"] * ALPHA_q[ln]

        # UDL acts uniformly over full B_ext per 1 m strip
        udl = q_ik * B_ext

        # TS dispersed through fill
        ts = _ts_per_m(Q_ik, dg)

        total = udl + ts
        max_V += total

        lane_results.append(LaneResult(
            lane=ln, Q_ik=Q_ik, q_ik=q_ik,
            udl_per_m=udl, ts_per_m=ts, total_per_m=total,
        ))

    # ── Braking force (BS EN 1991-2 Cl. 4.4.1) ──────────────────────────────
    # Lane 1 governs; loaded length L = B_ext (span in travel direction)
    Q1k = LM1_CHAR[1]["Q_ik"] * ALPHA_Q[1]
    q1k = LM1_CHAR[1]["q_ik"] * ALPHA_q[1]
    w1  = lane_width

    Q_lk_raw = 0.6 * 2 * Q1k + 0.10 * q1k * w1 * B_ext
    Q_lk     = max(180.0, min(900.0, Q_lk_raw))   # UK NA limits

    # Per metre strip: braking distributed uniformly over lane 1 width (LL direction)
    Q_lk_per_m = Q_lk / w1   # kN/m of barrel within lane 1

    braking = BrakingResult(
        Q_lk_raw=Q_lk_raw,
        Q_lk=Q_lk,
        Q_lk_per_m=Q_lk_per_m,
    )

    return LM1Result(
        n_lanes=n,
        B_ext=B_ext,
        LL=LL,
        H_c=H_c,
        lane_width=lane_width,
        dispersion=dg,
        lanes=lane_results,
        braking=braking,
        max_V_per_m=max_V,
        min_V_per_m=0.0,   # minimum = no live load
    )


# ── Text summary ──────────────────────────────────────────────────────────────

def summary(r: LM1Result) -> str:
    dg = r.dispersion
    merged_str = "merged" if dg.axles_merged else "separate"
    lines = [
        "=" * 70,
        f"LM1 Loading — {r.n_lanes} lane(s)  |  B_ext={r.B_ext:.2f}m  LL={r.LL:.2f}m  Hc={r.H_c:.3f}m",
        "-" * 70,
        f"Dispersion (1:1) at crown:",
        f"  LL direction  : {WHEEL_SPACING:.1f} + {CONTACT_T:.1f} + 2×{DISP}×{r.H_c:.3f} = {dg.disp_LL:.3f} m",
        f"  B_ext direction: {AXLE_SPACING:.1f} + {CONTACT_L:.1f} + 2×{DISP}×{r.H_c:.3f} = {dg.disp_Bext:.3f} m  [{merged_str}]",
        f"  Effective loaded length in B_ext: min({dg.disp_Bext:.3f}, {r.B_ext:.2f}) = {dg.Bext_loaded:.3f} m",
        "-" * 70,
        f"{'Lane':<6}{'Q_ik (kN)':<12}{'q_ik (kN/m²)':<15}{'UDL/m (kN/m)':<15}{'TS/m (kN/m)':<14}{'Total/m (kN/m)':<14}",
        "-" * 70,
    ]
    for ln in r.lanes:
        lines.append(
            f"{ln.lane:<6}{ln.Q_ik:<12.1f}{ln.q_ik:<15.1f}{ln.udl_per_m:<15.2f}{ln.ts_per_m:<14.2f}{ln.total_per_m:<14.2f}"
        )
    lines += [
        "-" * 70,
        f"  Max vertical per metre strip (all lanes):  {r.max_V_per_m:>8.2f} kN/m",
        f"  Min vertical per metre strip (no LL):      {r.min_V_per_m:>8.2f} kN/m",
        "=" * 70,
        "Braking / acceleration force (Lane 1 governs):",
        f"  Q_lk = 0.6 × 2 × {LM1_CHAR[1]['Q_ik']:.0f} + 0.10 × {LM1_CHAR[1]['q_ik']:.1f} × {r.lane_width:.2f} × {r.B_ext:.2f}",
        f"       = {r.braking.Q_lk_raw:.1f} kN  ->  clamped [{180},{900}] kN  =  {r.braking.Q_lk:.1f} kN",
        f"  Per metre strip (÷ lane width {r.lane_width:.2f} m): {r.braking.Q_lk_per_m:.2f} kN/m",
        "=" * 70,
    ]
    return "\n".join(lines)


# ── Quick test ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    result = compute(
        B_ext=2.5,    # culvert external width (travel direction)
        LL=6.0,       # barrel length
        H_c=1.0,      # cover depth
        lane_width=3.65,
        n_lanes=2,
    )
    print(summary(result))
