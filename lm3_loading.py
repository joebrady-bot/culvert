"""
lm3_loading.py
==============
LM3 Special Vehicle (SV) loading for a buried box culvert.

!! ALL SV AXLE LOADS AND POSITIONS MUST BE VERIFIED AGAINST UK NA TABLE NA.5 !!

Conventions
-----------
- SV vehicle travels in the B_ext direction (transverse to barrel axis).
- Metre-strip analysis: 1 m strip in LL direction (along barrel axis).
- SV occupies Lane 1; remaining lanes carry LM1 Lane 2 / 3 / 4 TS + UDL.
- Worst-case SV position found by scanning B_ext offsets (0.05 m step).
- 30° load dispersion through cover depth H_c per PD6694-1 Figure 11 (same as LM1).

References
----------
BS EN 1991-2:2003  Cl. 4.3.4 (LM3 Special Vehicles)
UK NA to BS EN 1991-2  Table NA.5 (SV vehicle definitions)
PD6694-1  load dispersion through fill
"""

import math
from dataclasses import dataclass, field
from typing import List, Tuple

import lm1_loading
from lm1_loading import (
    LM1_CHAR, ALPHA_Q, ALPHA_q,
    DISP,
    LaneResult, DispersionGeometry,
)

# ── SV geometry constants (UK NA Table NA.5 / worked example) ────────────────
SV_CONTACT_L     = 0.35   # m — wheel contact patch in B_ext direction
SV_CONTACT_T     = 0.35   # m — wheel contact patch in LL direction
SV_WHEEL_SPACING = 2.65   # m — wheel centre-to-centre in LL direction


# ── SV Vehicle definitions ────────────────────────────────────────────────────
# !! APPROXIMATE DATA — VERIFY ALL ENTRIES AGAINST UK NA TABLE NA.5 !!

@dataclass
class SVVehicle:
    name:       str
    gvw:        float          # Gross vehicle weight (kN)
    axle_loads: List[float]    # Load per axle (kN), front → rear
    axle_pos:   List[float]    # Distance of each axle from the front axle (m)


# !! APPROXIMATE — VERIFY AGAINST UK NA TABLE NA.5 BEFORE USE IN DESIGN !!
SV_VEHICLES: dict = {
    "SV80": SVVehicle(
        name="SV80", gvw=800.0,
        axle_loads=[100.0] * 8,
        # 4 paired groups, 1.5 m within pair, ~5.0 m between pairs
        axle_pos=[0.0, 1.5, 6.5, 8.0, 13.0, 14.5, 19.5, 21.0],
    ),
    "SV100": SVVehicle(
        name="SV100", gvw=1000.0,
        axle_loads=[100.0] * 10,
        # 5 paired groups
        axle_pos=[0.0, 1.5, 6.5, 8.0, 13.0, 14.5, 19.5, 21.0, 26.0, 27.5],
    ),
    "SV150": SVVehicle(
        name="SV150", gvw=1500.0,
        axle_loads=[100.0] * 15,
        # 5 groups of 3 axles, 1.5 m within group, 3.0 m between groups
        axle_pos=[
            0.0,  1.5,  4.5,  6.0,  9.0, 10.5,
            13.5, 15.0, 18.0, 19.5, 22.5, 24.0,
            27.0, 28.5, 31.5,
        ],
    ),
    "SV196": SVVehicle(
        name="SV196",
        # gvw = sum of BASIC (unfactored) axle loads — used for braking Q = 0.25 × GVW
        # 100 + 2×180 + 9×165 = 1945 kN  →  0.25 × 1945 = 486 kN braking
        gvw=1945.0,
        # DAF-factored axle loads (UK NA / worked example D. Childs):
        #   Steer (×1): 100 kN × DAF 1.20 = 120.0 kN
        #   Drive (×2): 180 kN × DAF 1.10 = 198.0 kN
        #   Trailer(×9): 165 kN × DAF 1.12 = 184.8 kN
        axle_loads=[120.0, 198.0, 198.0] + [184.8] * 9,
        # Approximate axle positions (m from steer axle) — verify against UK NA Table NA.5
        # Steer | Drive bogie | king-pin gap | 9 trailer axles at 1.2 m spacing
        axle_pos=[0.0, 3.5, 5.0, 8.0, 9.2, 10.4, 11.6, 12.8, 14.0, 15.2, 16.4, 17.6],
    ),
}


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class SVDispersionGeometry:
    """Per-axle dispersed footprint at crown level."""
    H_c:        float   # cover depth (m)
    disp_LL:    float   # dispersed width in LL direction (m)  [same formula as LM1]
    disp_B:     float   # dispersed length per axle in B_ext direction (m)


@dataclass
class LM3BrakingResult:
    """
    Horizontal braking/acceleration force for an SV vehicle.
    BS EN 1991-2 Cl. 4.4.4: Q_brk = 0.25 × basic GVW (unfactored axle loads).
    For in-situ box structures, distributed over barrel length LL via in-plane rigidity
    (PD6694-1 Cl. 10.2.8.2).
    Applied at crown level; arm = H_ext from base.
    """
    Q_brk_raw:     float   # 0.25 × basic GVW (kN)
    Q_brk_per_m:   float   # per metre strip in LL direction (kN/m) — ÷ LL


@dataclass
class LM3Result:
    vehicle_name:       str
    gvw:                float
    n_lanes:            int
    B_ext:              float
    LL:                 float
    H_c:                float
    lane_width:         float
    dispersion:         SVDispersionGeometry
    sv_load_per_m:      float   # SV contribution per metre strip (kN/m)
    sv_worst_offset:    float   # front-axle offset from culvert leading edge (m)
    sv_n_axles:         int     # number of SV axles contributing at worst position
    secondary_lanes:    List[LaneResult]    # LM1 Lanes 2, 3, 4
    secondary_per_m:    float              # LM1 secondary-lane total per metre strip (kN/m)
    braking:            LM3BrakingResult   # SV horizontal braking force
    max_V_per_m:        float              # SV + secondary lanes (kN/m)
    min_V_per_m:        float              # = 0


# ── Core helpers ──────────────────────────────────────────────────────────────

def _dispersion(H_c: float) -> SVDispersionGeometry:
    """
    Per-axle dispersed footprint at crown level.

    B_ext direction (one axle):
        disp_B = SV_CONTACT_L + 2 * DISP * H_c

    LL direction — SV wheels are 2.65 m apart; each wheel disperses independently.
    When the dispersed patches of adjacent wheels don't overlap (gap >= 0), the worst
    1m strip is under one wheel only:
        load/m = (axle_load/2) / disp_LL_single   [one wheel per strip]
    Encoding this as ax_load / disp_LL by setting disp_LL = 2 × disp_LL_single.
    When patches overlap (gap < 0) they merge:
        disp_LL = SV_WHEEL_SPACING + disp_LL_single
    """
    spread         = DISP * H_c
    disp_LL_single = SV_CONTACT_T + 2 * spread              # one wheel's LL dispersal
    gap_LL         = SV_WHEEL_SPACING - disp_LL_single       # > 0 → patches separate
    if gap_LL >= 0:
        disp_LL = 2 * disp_LL_single                         # non-overlapping: per-wheel
    else:
        disp_LL = SV_WHEEL_SPACING + disp_LL_single          # merged: both wheels
    disp_B = SV_CONTACT_L + 2 * spread
    return SVDispersionGeometry(H_c=H_c, disp_LL=disp_LL, disp_B=disp_B)


def _sv_load_on_culvert(
    sv:      SVVehicle,
    B_ext:   float,
    disp:    SVDispersionGeometry,
    step:    float = 0.05,
) -> Tuple[float, float, int]:
    """
    Scan SV vehicle positions along B_ext and return the worst case.

    For each front-axle offset the load landing on the culvert [0, B_ext]
    from each axle is:

        overlap   = max(0, min(axle_right, B_ext) - max(axle_left, 0))
        fraction  = overlap / disp.disp_B
        load/m    = axle_load × fraction / disp.disp_LL

    Returns
    -------
    (max_load_per_m, worst_offset, n_contributing_axles)
    """
    half_B       = disp.disp_B / 2.0
    total_length = sv.axle_pos[-1] if sv.axle_pos else 0.0

    # Scan range: first axle just before culvert → last axle just past culvert
    offset = -(total_length + half_B)
    end    = B_ext + half_B

    max_load   = 0.0
    worst_off  = offset
    worst_n    = 0

    while offset <= end + 1e-9:
        load_here = 0.0
        n_contrib = 0
        for ax_pos, ax_load in zip(sv.axle_pos, sv.axle_loads):
            centre   = offset + ax_pos
            ax_left  = centre - half_B
            ax_right = centre + half_B
            overlap  = max(0.0, min(ax_right, B_ext) - max(ax_left, 0.0))
            frac     = overlap / disp.disp_B if disp.disp_B > 0 else 0.0
            load_here += ax_load * frac / disp.disp_LL
            if frac > 0:
                n_contrib += 1
        if load_here > max_load:
            max_load  = load_here
            worst_off = offset
            worst_n   = n_contrib
        offset += step

    return max_load, worst_off, worst_n


# ── Main computation ──────────────────────────────────────────────────────────

def compute(
    vehicle_name: str,
    B_ext:        float,
    LL:           float,
    H_c:          float,
    lane_width:   float,
    n_lanes:      int,
) -> LM3Result:
    """
    Compute LM3 vertical loads for a buried box culvert.

    Parameters
    ----------
    vehicle_name : One of 'SV80', 'SV100', 'SV150', 'SV196'
    B_ext        : External culvert width in vehicle travel direction (m)
    LL           : Barrel length perpendicular to travel (m)
    H_c          : Cover depth — road surface to culvert crown (m)
    lane_width   : Width of each notional lane (m)
    n_lanes      : Total number of notional lanes (max 4)

    Returns
    -------
    LM3Result with SV contribution + LM1 secondary-lane breakdown
    """
    sv   = SV_VEHICLES[vehicle_name]
    disp = _dispersion(H_c)
    n    = min(max(1, n_lanes), 4)

    # ── Lane 1: SV vehicle (worst-case position) ──────────────────────────────
    sv_load, sv_offset, sv_n = _sv_load_on_culvert(sv, B_ext, disp)

    # ── Lanes 2+ : LM1 secondary lanes ───────────────────────────────────────
    lm1_dg = lm1_loading._dispersion(H_c, B_ext)   # LM1 tandem dispersion
    secondary: List[LaneResult] = []
    sec_total = 0.0

    for ln in range(2, n + 1):
        Q_ik = LM1_CHAR[ln]["Q_ik"] * ALPHA_Q[ln]
        q_ik = LM1_CHAR[ln]["q_ik"] * ALPHA_q[ln]
        udl  = q_ik * B_ext
        ts   = lm1_loading._ts_per_m(Q_ik, lm1_dg)
        tot  = udl + ts
        sec_total += tot
        secondary.append(LaneResult(
            lane=ln, Q_ik=Q_ik, q_ik=q_ik,
            udl_per_m=udl, ts_per_m=ts, total_per_m=tot,
        ))

    max_V = sv_load + sec_total

    # ── SV braking — BS EN 1991-2 Cl. 4.4.4 / PD6694-1 Cl. 10.2.8.2 ──────────
    # 25% of basic (unfactored) GVW, distributed over barrel length via in-plane rigidity.
    Q_brk_raw   = 0.25 * sv.gvw   # 25% of basic GVW (e.g. 0.25×1945 = 486 kN for SV196)
    Q_brk_per_m = Q_brk_raw / LL  # distribute over barrel length Lj = LL

    braking = LM3BrakingResult(
        Q_brk_raw=Q_brk_raw,
        Q_brk_per_m=Q_brk_per_m,
    )

    return LM3Result(
        vehicle_name=vehicle_name,
        gvw=sv.gvw,
        n_lanes=n,
        B_ext=B_ext,
        LL=LL,
        H_c=H_c,
        lane_width=lane_width,
        dispersion=disp,
        sv_load_per_m=sv_load,
        sv_worst_offset=sv_offset,
        sv_n_axles=sv_n,
        secondary_lanes=secondary,
        secondary_per_m=sec_total,
        braking=braking,
        max_V_per_m=max_V,
        min_V_per_m=0.0,
    )


# ── Text summary ──────────────────────────────────────────────────────────────

def summary(r: LM3Result) -> str:
    dg = r.dispersion
    sv = SV_VEHICLES[r.vehicle_name]
    n_axles = len(sv.axle_loads)
    lines = [
        "=" * 70,
        f"LM3 Loading — {r.vehicle_name} ({r.gvw:.0f} kN GVW, {n_axles} axles)",
        f"  {r.n_lanes} lane(s)  |  B_ext={r.B_ext:.2f} m  LL={r.LL:.2f} m  Hc={r.H_c:.3f} m",
        "-" * 70,
        "!! SV axle data is approximate — verify against UK NA Table NA.5 !!",
        "-" * 70,
        f"Per-axle dispersion (30°, tan30°={DISP:.4f}) at crown (contact patch {SV_CONTACT_L*1000:.0f}×{SV_CONTACT_T*1000:.0f}mm, wheel spacing {SV_WHEEL_SPACING:.2f}m):",
        f"  LL per wheel  : {SV_CONTACT_T:.3f} + 2x{DISP:.4f}x{dg.H_c:.3f} = {SV_CONTACT_T + 2*DISP*dg.H_c:.3f} m  -> disp_LL = 2x = {dg.disp_LL:.3f} m",
        f"  B_ext per axle: {SV_CONTACT_L:.3f} + 2×{DISP:.4f}×{dg.H_c:.3f} = {dg.disp_B:.3f} m",
        "-" * 70,
        f"SV vehicle (Lane 1):",
        f"  Worst offset (front axle from culvert edge): {r.sv_worst_offset:+.3f} m",
        f"  Contributing axles at worst position:        {r.sv_n_axles}",
        f"  SV load per metre strip:                     {r.sv_load_per_m:.2f} kN/m",
        "-" * 70,
    ]
    if r.secondary_lanes:
        lines.append(
            f"{'Lane':<6}{'Q_ik (kN)':<12}{'q_ik (kN/m²)':<15}"
            f"{'UDL/m (kN/m)':<15}{'TS/m (kN/m)':<14}{'Total/m (kN/m)':<14}"
        )
        lines.append("-" * 70)
        for ln in r.secondary_lanes:
            lines.append(
                f"{ln.lane:<6}{ln.Q_ik:<12.1f}{ln.q_ik:<15.1f}"
                f"{ln.udl_per_m:<15.2f}{ln.ts_per_m:<14.2f}{ln.total_per_m:<14.2f}"
            )
        lines.append("-" * 70)
        lines.append(f"  LM1 secondary lanes total:  {r.secondary_per_m:.2f} kN/m")
    else:
        lines.append("  No secondary lanes (1 lane total).")
    brk = r.braking
    lines += [
        "=" * 70,
        f"  Max vertical per metre strip (SV + secondary): {r.max_V_per_m:>8.2f} kN/m",
        f"  Min vertical per metre strip (no LL):          {r.min_V_per_m:>8.2f} kN/m",
        "=" * 70,
        "Braking / acceleration force (BS EN 1991-2 Cl. 4.4.4 / PD6694-1 Cl. 10.2.8.2):",
        f"  Q_brk = 0.25 × {r.gvw:.0f} kN basic GVW = {brk.Q_brk_raw:.1f} kN",
        f"  Distributed over barrel length LL = {r.LL:.2f} m (in-plane rigidity)",
        f"  Per metre strip: {brk.Q_brk_per_m:.2f} kN/m",
        "=" * 70,
    ]
    return "\n".join(lines)


# ── Quick test ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    for vname in SV_VEHICLES:
        result = compute(
            vehicle_name=vname,
            B_ext=2.5,
            LL=6.0,
            H_c=1.0,
            lane_width=3.65,
            n_lanes=3,
        )
        print(summary(result))
        print()
