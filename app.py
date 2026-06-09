import streamlit as st
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import lm1_loading
import lm3_loading

st.set_page_config(page_title="Box Culvert Stability", layout="wide")

st.title("Buried Box Culvert — Stability Checker")
st.caption("Metre-strip · PD6694-1 Annex B Tables B.4/B.5/B.6 · BS EN 1991-2 LM1 · BS EN 1997-1")

# ── PD6694-1 Annex B ──────────────────────────────────────────────────────────
# Ka  : design active coefficient — Tables B.4/B.5/B.6 (includes γM · γSd;K = 1.2)
# Kmax: design max coefficient for restrained side — Tables B.4/B.5 (includes γM · γSd;K)
# Kr  : Rankine passive for B.6 sliding resistance — computed from design φ (not tabulated)
#
# EC7 partial factors
#   gG_u / gG_f: unfavourable / favourable on permanent actions
#   g_phi / g_c: on shear strength parameters (resistance side)
LS = {
    "SLS":        {"Ka": 0.33, "Kmax": 0.60, "gG_u": 1.00, "gG_f": 1.00, "g_phi": 1.00, "g_c": 1.00, "gQ": 1.00},
    "EQU":        {"Ka": 0.44, "Kmax": 0.60, "gG_u": 1.10, "gG_f": 0.90, "g_phi": 1.25, "g_c": 1.25, "gQ": 1.35},
    "STR/GEO C1": {"Ka": 0.40, "Kmax": 0.72, "gG_u": 1.35, "gG_f": 1.00, "g_phi": 1.00, "g_c": 1.00, "gQ": 1.35},
    "STR/GEO C2": {"Ka": 0.49, "Kmax": 0.84, "gG_u": 1.00, "gG_f": 1.00, "g_phi": 1.25, "g_c": 1.25, "gQ": 1.35},
}
LS_NAMES = list(LS.keys())

# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Geometry")
    B   = st.number_input("Internal Width, B (m)",    0.1, 20.0, 2.0,  0.1,  "%.2f")
    H   = st.number_input("Internal Height, H (m)",   0.1, 10.0, 1.5,  0.1,  "%.2f")
    LL  = st.number_input("Overall Length, LL (m)",   0.5, 50.0, 6.0,  0.5,  "%.1f",
                           help="Longitudinal barrel length — used for load dispersion and braking reduction factor η.")
    t_w = st.number_input("Wall Thickness, t_w (m)",  0.05, 2.0, 0.25, 0.05, "%.2f")
    t_s = st.number_input("Slab Thickness, t_s (m)",  0.05, 2.0, 0.30, 0.05, "%.2f")

    st.header("Cover Layers (top to bottom)")
    st.markdown("**Road construction**")
    t_road  = st.number_input("Thickness (m)",       0.0, 2.0,  0.10, 0.01, "%.3f", key="t_road")
    γ_road  = st.number_input("Unit Weight (kN/m³)", 10.0, 30.0, 24.0, 0.5,         key="g_road")
    st.markdown("**Subbase**")
    t_sub   = st.number_input("Thickness (m)",       0.0, 2.0,  0.25, 0.05, "%.3f", key="t_sub")
    γ_sub   = st.number_input("Unit Weight (kN/m³)", 10.0, 25.0, 20.0, 0.5,         key="g_sub")
    st.markdown("**Fill**")
    t_fill  = st.number_input("Thickness (m)",       0.0, 20.0, 0.65, 0.05, "%.3f", key="t_fill")
    γ_fill  = st.number_input("Unit Weight (kN/m³)", 10.0, 25.0, 18.0, 0.5,         key="g_fill")

    st.header("Backfill Properties")
    st.caption("Used for passive (Kr) resistance on culvert walls")
    φ_fill_deg = st.number_input("Friction Angle φ'_fill (°)", 0.0, 45.0, 35.0, 1.0, key="phi_fill")
    c_fill     = st.number_input("Cohesion c'_fill (kPa)",     0.0, 200.0,  0.0, 1.0, key="c_fill")

    st.header("Founding Layer Properties")
    st.caption("Used for base sliding resistance and bearing")
    φ_fnd_deg = st.number_input("Friction Angle φ'_fnd (°)", 0.0, 45.0, 28.0, 1.0, key="phi_fnd")
    c_fnd     = st.number_input("Cohesion c'_fnd (kPa)",     0.0, 200.0,  0.0, 1.0, key="c_fnd")
    q_Rd      = st.number_input(
        "Bearing Resistance q_Rd (kPa)", 50.0, 5000.0, 300.0, 10.0,
        help="Design bearing resistance from ground model (applied uniformly across limit states).")

    st.header("Road Geometry")
    st.caption("Road runs transversely across the culvert (vehicles travel in B_ext direction)")
    cw_width   = st.number_input("Carriageway Width (m)", 1.0, 30.0, 7.30, 0.05, "%.2f", key="cw_width",
                                  help="Total width between kerbs, centred on barrel length LL")
    lane_width = st.number_input("Lane Width (m)",        1.0, 5.0,  3.65, 0.05, "%.2f", key="lane_width",
                                  help="Individual lane width")

    st.header("LM3 Special Vehicle")
    st.caption("!! Verify SV data against UK NA Table NA.5 !!")
    lm3_vehicle = st.selectbox("SV Vehicle", list(lm3_loading.SV_VEHICLES.keys()))

    st.header("Water Table")
    h_wt = st.number_input("Depth below GL (m)", 0.0, 30.0, 5.0, 0.1, key="h_wt",
                            help="Shown on diagram. Affects effective stresses and uplift in B.6 checks.")

    st.header("Material")
    γ_c = st.number_input("Concrete Unit Weight γ_c (kN/m³)", 20.0, 26.0, 24.0, 0.5)

# ── Derived geometry ───────────────────────────────────────────────────────────
B_ext  = B + 2 * t_w
H_ext  = H + 2 * t_s
H_c    = t_road + t_sub + t_fill   # total cover depth to crown
H_inv  = H_c + H_ext

# Road geometry — notional lanes
n_lanes = max(1, int(cw_width / lane_width))

A_conc = B_ext * H_ext - B * H
W_conc = A_conc * γ_c
# Overburden weight: sum of each layer (characteristic, not factored)
W_road = γ_road * t_road * B_ext
W_sub  = γ_sub  * t_sub  * B_ext
W_fill = γ_fill * t_fill * B_ext
W_soil = W_road + W_sub + W_fill
N_G_k  = W_conc + W_soil   # characteristic permanent vertical load (per metre strip)

φ_fill_k = math.radians(φ_fill_deg)   # backfill — governs Kr passive resistance
φ_fnd_k  = math.radians(φ_fnd_deg)   # founding layer — governs base sliding resistance

γ_w = 9.81  # kN/m³

# Total vertical stresses at crown and invert
σ_top = γ_road * t_road + γ_sub * t_sub + γ_fill * t_fill
σ_bot = σ_top + γ_fill * H_ext

# Pore water pressures at crown and invert (zero above water table)
u_top = γ_w * max(0.0, H_c   - h_wt)
u_bot = γ_w * max(0.0, H_inv - h_wt)

# Effective vertical stresses (used for K·σ'v horizontal earth pressure)
σ_eff_top = max(0.0, σ_top - u_top)
σ_eff_bot = max(0.0, σ_bot - u_bot)

# Characteristic uplift on base per unit length
U_k = u_bot * B_ext

def trapz_resultant(K, σ_t, σ_b, h):
    """Resultant of trapezoidal horizontal pressure K·σ over height h."""
    F   = 0.5 * K * (σ_t + σ_b) * h
    arm = h * (2*σ_t + σ_b) / (3*(σ_t + σ_b)) if (σ_t + σ_b) > 0 else h/3
    return F, arm

# ── LM1 vehicle loading (computed first — Q_vk feeds into run()) ───────────────
lm1 = lm1_loading.compute(
    B_ext=B_ext, LL=LL, H_c=H_c,
    lane_width=lane_width, n_lanes=n_lanes,
)
Q_vk = lm1.max_V_per_m   # characteristic LM1 vertical load (all lanes, per metre strip)

# ── LM3 special vehicle loading ───────────────────────────────────────────────
lm3 = lm3_loading.compute(
    vehicle_name=lm3_vehicle,
    B_ext=B_ext, LL=LL, H_c=H_c,
    lane_width=lane_width, n_lanes=n_lanes,
)

# ── PD6694-1 Table 6 — horizontal traffic loads on active wall ────────────────
# Kd = Ka (characteristic, without γM/γSd;K) per Table 6 Note 4, active side only
Ka_char   = math.tan(math.pi/4 - φ_fill_k/2) ** 2

# Note 5 reduction factor for F line loads (buried structure, Hc in metres)
_rF       = (1.0 - H_c / 2.0) ** 2  if H_c < 2.0  else 0.0

# Characteristic horizontal forces per metre strip (LL direction)
Q_h_udl_k = 20.0  * Ka_char * H_ext                        # σh UDL resultant (kN/m)
Q_h_F_k   = 2.0 * 300.0 * Ka_char * _rF / lane_width      # two F line loads (kN/m)
# Moment arms about the base (invert level):
#   σh uniform over H_ext  → arm = H_ext/2
#   F applied at crown     → arm = H_ext
Q_h_M_k   = Q_h_udl_k * (H_ext / 2.0) + Q_h_F_k * H_ext  # driving moment (kNm/m)

# ── Per-limit-state calculations ───────────────────────────────────────────────
def run(p):
    Ka, Kmax   = p["Ka"],   p["Kmax"]
    gG_u, gG_f = p["gG_u"], p["gG_f"]
    gQ         = p["gQ"]
    g_phi, g_c = p["g_phi"],p["g_c"]

    # Design soil parameters
    φ_fill_d = math.atan(math.tan(φ_fill_k) / g_phi)
    c_fill_d = c_fill / g_c
    φ_fnd_d  = math.atan(math.tan(φ_fnd_k)  / g_phi)
    c_fnd_d  = c_fnd  / g_c
    Kp       = math.tan(math.pi/4 + φ_fill_d/2) ** 2

    # ── B.4 / B.5 — total stresses, no water table effect ───────────────────
    # Traffic is unfavourable for B.4 bearing (max vertical); excluded from
    # B.5/overturning resistance where minimum vertical governs.
    V_u = gG_u * N_G_k + gQ * Q_vk   # B.4 bearing: max vertical (Gk + Qk)
    V_f = gG_f * N_G_k                # B.4/B.5 OT/sliding resistance: min vertical (Gk only)

    F_Ka,   arm_Ka   = trapz_resultant(Ka,   σ_top, σ_bot, H_ext)
    F_Kmax, arm_Kmax = trapz_resultant(Kmax, σ_top, σ_bot, H_ext)
    F_net    = F_Kmax - F_Ka
    M_net_B4 = F_Kmax * arm_Kmax - F_Ka * arm_Ka   # earth pressure only (for display)
    R_fric   = math.tan(φ_fnd_d) * V_f + c_fnd_d * B_ext
    M_stb    = V_f * B_ext / 2

    # Table 6 horizontal traffic loads — active side, factored by γQ
    F_h_tr   = gQ * (Q_h_udl_k + Q_h_F_k)   # factored horizontal force (kN/m)
    M_h_tr   = gQ * Q_h_M_k                   # factored driving moment (kNm/m)

    # B.4/B.5 OT and sliding include Table 6 traffic horizontal on active side
    F_drv_B45  = F_net    + F_h_tr
    M_drv_B45  = M_net_B4 + M_h_tr

    q_B4       = V_u / B_ext
    UR_B4_bear = q_B4 / q_Rd if q_Rd > 0 else float("inf")
    UR_B4_ov   = M_drv_B45 / M_stb  if M_stb  > 0 else float("inf")
    UR_B4_sl   = F_drv_B45 / R_fric if R_fric > 0 else float("inf")

    q_B5       = V_f / B_ext
    UR_B5_bear = q_B5 / q_Rd if q_Rd > 0 else float("inf")
    UR_B5_ov   = UR_B4_ov
    UR_B5_sl   = UR_B4_sl

    # ── B.6 — effective stresses + uplift + traffic ──────────────────────────
    V_u_B6 = max(gG_u * N_G_k + gQ * Q_vk - gG_f * U_k, 0.0)
    V_f_B6 = max(gG_f * N_G_k              - gG_u * U_k, 0.0)

    F_Ka_B6,   arm_Ka_B6 = trapz_resultant(Ka, σ_eff_top, σ_eff_bot, H_ext)
    F_Kr,      arm_Kr     = trapz_resultant(Kp, σ_eff_top, σ_eff_bot, H_ext)
    R_fric_B6  = math.tan(φ_fnd_d) * V_f_B6 + c_fnd_d * B_ext
    M_stb_B6   = V_f_B6 * B_ext / 2

    # Table 6 traffic horizontal on B.6 active (Ka) driving side
    F_drv_B6   = F_Ka_B6 + F_h_tr
    M_Ka_B6    = F_Ka_B6 * arm_Ka_B6 + M_h_tr
    M_net_B6   = max(0.0, M_Ka_B6 - F_Kr * arm_Kr)
    R_B6       = F_Kr + R_fric_B6

    q_B6       = V_u_B6 / B_ext
    UR_B6_bear = q_B6  / q_Rd if q_Rd > 0 else float("inf")
    UR_B6_ov   = M_net_B6 / M_stb_B6 if M_stb_B6 > 0 else float("inf")
    UR_B6_sl   = F_drv_B6 / R_B6     if R_B6     > 0 else float("inf")

    return dict(
        Ka=Ka, Kmax=Kmax, gG_u=gG_u, gG_f=gG_f, gQ=gQ, g_phi=g_phi, g_c=g_c,
        V_u=V_u, V_f=V_f, V_u_B6=V_u_B6, V_f_B6=V_f_B6, U_k=U_k,
        u_top=u_top, u_bot=u_bot, σ_eff_top=σ_eff_top, σ_eff_bot=σ_eff_bot,
        φ_fill_d_deg=math.degrees(φ_fill_d), c_fill_d=c_fill_d,
        φ_fnd_d_deg=math.degrees(φ_fnd_d),   c_fnd_d=c_fnd_d,
        F_Ka=F_Ka, F_Ka_B6=F_Ka_B6, F_Kmax=F_Kmax, F_net=F_net,
        arm_Ka=arm_Ka, arm_Ka_B6=arm_Ka_B6, arm_Kmax=arm_Kmax, arm_Kr=arm_Kr,
        Kp=Kp, F_Kr=F_Kr, R_fric=R_fric, R_fric_B6=R_fric_B6, R_B6=R_B6,
        F_h_tr=F_h_tr, M_h_tr=M_h_tr,
        F_drv_B45=F_drv_B45, M_drv_B45=M_drv_B45,
        F_drv_B6=F_drv_B6, M_Ka_B6=M_Ka_B6,
        q_B4=q_B4, q_B5=q_B5, q_B6=q_B6,
        M_stb=M_stb, M_stb_B6=M_stb_B6, M_net_B4=M_net_B4, M_net_B6=M_net_B6,
        UR_B4_bear=UR_B4_bear, UR_B4_ov=UR_B4_ov, UR_B4_sl=UR_B4_sl,
        UR_B5_bear=UR_B5_bear, UR_B5_ov=UR_B5_ov, UR_B5_sl=UR_B5_sl,
        UR_B6_bear=UR_B6_bear, UR_B6_ov=UR_B6_ov, UR_B6_sl=UR_B6_sl,
    )

res = {n: run(LS[n]) for n in LS_NAMES}

# ── Shared drawing helper ──────────────────────────────────────────────────────
def dim_arrow(ax, x, y0, y1, label, side="r", fontsize=6.5):
    ax.annotate("", xy=(x, y1), xytext=(x, y0),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    dx = 0.08 if side == "r" else -0.08
    ax.text(x + dx, (y0 + y1) / 2, label, fontsize=fontsize, va="center",
            ha="left" if side == "r" else "right")

def h_dim_arrow(ax, x0, x1, y, label, fontsize=6.5):
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text((x0 + x1) / 2, y + 0.06, label, fontsize=fontsize,
            ha="center", va="bottom")

legend_patches = []
if t_road > 0:
    legend_patches.append(patches.Patch(fc="#404040", label=f"Road ({t_road*1000:.0f} mm)"))
if t_sub > 0:
    legend_patches.append(patches.Patch(fc="#B0B0B0", label=f"Subbase ({t_sub*1000:.0f} mm)"))
legend_patches.append(patches.Patch(fc="#C8A86E", label=f"Fill ({t_fill:.2f} m)"))
legend_patches.append(patches.Patch(fc="#8B7355", label="Founding layer"))

fnd_depth = 0.5   # founding layer depth shown below invert in both diagrams

# ══ Section 1: Geometry ═══════════════════════════════════════════════════════
with st.expander("Geometry", expanded=True):
    col_xs, col_elev = st.columns(2)

# ── Cross-section diagram ──────────────────────────────────────────────────────
with col_xs:
    st.subheader("Cross-Section")
    pad  = max(B_ext * 0.6, 0.8)
    xlim = (-B_ext/2 - pad, B_ext/2 + pad)
    ylim = (-(H_inv + fnd_depth + 0.2), 0.4)

    fig, ax = plt.subplots(figsize=(4, 5))
    ax.set_facecolor("white")

    # Founding layer — sits at base, below invert
    ax.add_patch(patches.Rectangle(
        (xlim[0], -(H_inv + fnd_depth)), xlim[1]-xlim[0], fnd_depth,
        fc="#8A8278", ec="none", zorder=0))
    ax.text(0, -(H_inv + fnd_depth/2), "Founding Layer",
            ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")

    # Fill layer — from invert level up to GL (alongside culvert walls and above)
    ax.add_patch(patches.Rectangle(
        (xlim[0], -H_inv), xlim[1]-xlim[0], H_inv,
        fc="#C8A86E", ec="none", zorder=1))

    # Subbase layer (overwrites fill above)
    if t_sub > 0:
        ax.add_patch(patches.Rectangle(
            (xlim[0], -(t_road + t_sub)), xlim[1]-xlim[0], t_sub,
            fc="#B0B0B0", ec="none", zorder=1))

    # Road construction layer
    if t_road > 0:
        ax.add_patch(patches.Rectangle(
            (xlim[0], -t_road), xlim[1]-xlim[0], t_road,
            fc="#404040", ec="none", zorder=1))

    # Culvert concrete box
    ax.add_patch(patches.Rectangle(
        (-B_ext/2, -H_inv), B_ext, H_ext, fc="#BBBBBB", ec="#333333", lw=2, zorder=3))
    # Internal void
    ax.add_patch(patches.Rectangle(
        (-B/2, -H_inv+t_s), B, H, fc="#F0F8FF", ec="#777777", lw=1, zorder=4))
    ax.axhline(0, color="#3E1E00", lw=2.5, zorder=5)
    ax.text(xlim[1]-0.05, 0.08, "GL", color="#3E1E00", fontsize=8,
            ha="right", va="bottom", fontweight="bold")

    # Water table
    if ylim[0] < -h_wt < ylim[1]:
        ax.axhline(-h_wt, color="#1A6EBD", lw=1.2, ls="--", zorder=6)
        ax.text(xlim[0] + 0.05, -h_wt + 0.05, "▼ WT",
                color="#1A6EBD", fontsize=6.5, va="bottom", fontweight="bold")

    xr = B_ext/2 + 0.22
    xl = -B_ext/2 - 0.22
    if H_c > 0:
        dim_arrow(ax, xr, -H_c, 0, f"Hc={H_c:.2f} m")
    dim_arrow(ax, xr, -H_inv, -H_c, f"Hext={H_ext:.2f} m")
    dim_arrow(ax, xl, -H_inv, 0, f"Hinv={H_inv:.2f} m", "l")
    ax.annotate("", xy=(B_ext/2, -H_inv-0.32), xytext=(-B_ext/2, -H_inv-0.32),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text(0, -H_inv-0.50, f"Bext = {B_ext:.2f} m", ha="center", fontsize=6.5)

    ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_aspect("equal")
    ax.set_xlabel("(m)", fontsize=8); ax.set_ylabel("Depth below GL (m)", fontsize=8)
    ax.tick_params(labelsize=7); fig.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

# ── Longitudinal elevation diagram ────────────────────────────────────────────
with col_elev:
    st.subheader("Longitudinal Section")

    fig_e, ax_e = plt.subplots(figsize=(5, 4))
    ax_e.set_facecolor("white")

    # Founding layer — greyed warm tone, distinct from engineered fill above
    ax_e.add_patch(patches.Rectangle(
        (0, -(H_inv + fnd_depth)), LL, fnd_depth,
        fc="#8A8278", ec="#4A4240", lw=0.8, zorder=0))
    ax_e.text(LL / 2, -(H_inv + fnd_depth / 2), "Founding Layer",
              ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")

    # Fill (full width, from invert up to GL)
    ax_e.add_patch(patches.Rectangle(
        (0, -H_inv), LL, H_inv, fc="#C8A86E", ec="none", zorder=1))

    # Subbase
    if t_sub > 0:
        ax_e.add_patch(patches.Rectangle(
            (0, -(t_road + t_sub)), LL, t_sub, fc="#B0B0B0", ec="none", zorder=2))

    # Road construction
    if t_road > 0:
        ax_e.add_patch(patches.Rectangle(
            (0, -t_road), LL, t_road, fc="#404040", ec="none", zorder=2))

    # Culvert box (top and bottom slab + side walls visible in section)
    # Top slab
    ax_e.add_patch(patches.Rectangle(
        (0, -H_c), LL, t_s, fc="#BBBBBB", ec="#333333", lw=1.2, zorder=3))
    # Bottom slab
    ax_e.add_patch(patches.Rectangle(
        (0, -H_inv), LL, t_s, fc="#BBBBBB", ec="#333333", lw=1.2, zorder=3))
    # Left headwall hint
    ax_e.add_patch(patches.Rectangle(
        (0, -H_inv), t_w, H_ext, fc="#AAAAAA", ec="#333333", lw=1.2, zorder=3))
    # Right headwall hint
    ax_e.add_patch(patches.Rectangle(
        (LL - t_w, -H_inv), t_w, H_ext, fc="#AAAAAA", ec="#333333", lw=1.2, zorder=3))
    # Interior void
    ax_e.add_patch(patches.Rectangle(
        (t_w, -H_inv + t_s), LL - 2*t_w, H, fc="#F0F8FF", ec="#777777", lw=0.8, zorder=4))

    # Carriageway — centred on LL, width = cw_width
    cw_x0    = LL / 2 - cw_width / 2
    cw_x1    = LL / 2 + cw_width / 2
    cw_y_top = 0.22
    kerb_w   = max(0.05, cw_width * 0.015)

    # Road surface
    ax_e.add_patch(patches.Rectangle(
        (cw_x0, 0), cw_width, cw_y_top, fc="#3A3A3A", ec="none", zorder=5))

    # Dashed lane dividers
    for i in range(1, n_lanes):
        x_div = cw_x0 + i * lane_width
        if cw_x0 < x_div < cw_x1:
            ax_e.plot([x_div, x_div], [0, cw_y_top], color="white", lw=0.7, ls="--", zorder=6)

    # Lane labels
    for i in range(n_lanes):
        x_mid = cw_x0 + (i + 0.5) * lane_width
        if cw_x0 < x_mid < cw_x1:
            ax_e.text(x_mid, cw_y_top / 2, f"L{i+1}",
                      ha="center", va="center", fontsize=5, color="white")

    # Kerbs — full height of carriageway strip, at each edge
    for kx in [cw_x0 - kerb_w, cw_x1]:
        ax_e.add_patch(patches.Rectangle(
            (kx, 0), kerb_w, cw_y_top, fc="#888888", ec="#444444", lw=0.5, zorder=7))

    # Width annotation
    ax_e.annotate("", xy=(cw_x1, cw_y_top + 0.07), xytext=(cw_x0, cw_y_top + 0.07),
                  arrowprops=dict(arrowstyle="<->", color="#CCCCCC", lw=0.9))
    ax_e.text(LL / 2, cw_y_top + 0.13,
              f"CW={cw_width:.2f}m  {n_lanes}ln  lw={lane_width:.2f}m",
              ha="center", fontsize=5.5, color="#CCCCCC")

    # GL
    ax_e.axhline(0, color="#3E1E00", lw=2, zorder=5)
    ax_e.text(LL + 0.05, 0.04, "GL", color="#3E1E00", fontsize=7,
              ha="left", va="bottom", fontweight="bold")

    # Water table
    elev_ylim_bot = -(H_inv + fnd_depth + 0.35)
    if elev_ylim_bot < -h_wt < 0.45:
        ax_e.axhline(-h_wt, color="#1A6EBD", lw=1.2, ls="--", zorder=6)
        ax_e.text(-0.05, -h_wt + 0.03, "▼ WT",
                  color="#1A6EBD", fontsize=6, va="bottom", ha="right", fontweight="bold")

    # Dimensions
    xr_e = LL + 0.15
    if H_c > 0:
        dim_arrow(ax_e, xr_e, -H_c, 0,      f"Hc={H_c:.2f}m")
    dim_arrow(ax_e, xr_e, -H_inv, -H_c,     f"Hext={H_ext:.2f}m")
    dim_arrow(ax_e, -0.15, -H_inv, 0,        f"Hinv={H_inv:.2f}m", "l")
    h_dim_arrow(ax_e, 0, LL, -(H_inv + fnd_depth + 0.15), f"LL = {LL:.1f} m")

    ax_e.set_xlim(-0.5, LL + 0.55)
    ax_e.set_ylim(-(H_inv + fnd_depth + 0.35), 0.65)
    ax_e.set_xlabel("Length (m)", fontsize=8)
    ax_e.set_ylabel("Depth below GL (m)", fontsize=8)
    ax_e.tick_params(labelsize=7)
    fig_e.tight_layout()
    st.pyplot(fig_e, use_container_width=True)
    plt.close()

# ══ Section 2: LM1 Vehicle Loading ════════════════════════════════════════════
with st.expander("LM1 Vehicle Loading", expanded=True):
    _lm1_left, _lm1_right = st.columns([1, 1])

    with _lm1_left:
        st.markdown("**LM1 characteristic loads per metre strip** (BS EN 1991-2, UK NA α=1.0)")
        dg = lm1.dispersion
        lm1_rows = []
        for ln in lm1.lanes:
            lm1_rows.append({
                "Lane": f"Lane {ln.lane}",
                "Q_ik (kN/axle)": f"{ln.Q_ik:.0f}",
                "q_ik (kN/m²)":   f"{ln.q_ik:.1f}",
                "UDL/m (kN/m)":   f"{ln.udl_per_m:.2f}",
                "TS/m (kN/m)":    f"{ln.ts_per_m:.2f}",
                "Total/m (kN/m)": f"{ln.total_per_m:.2f}",
            })
        st.dataframe(pd.DataFrame(lm1_rows), hide_index=True, use_container_width=True)

        st.markdown(f"""
**Dispersion through {H_c:.3f} m cover (1:1):**
- LL direction: {dg.disp_LL:.3f} m &nbsp;·&nbsp; B_ext direction: {dg.disp_Bext:.3f} m ({'merged' if dg.axles_merged else 'separate'})
- Effective loaded length: min({dg.disp_Bext:.3f}, {B_ext:.3f}) = {dg.Bext_loaded:.3f} m

**Totals:**
- Max vertical (all lanes): **{lm1.max_V_per_m:.2f} kN/m**
- Min vertical (no LL): **0.00 kN/m**

**Braking (Lane 1 governs — BS EN 1991-2 Cl. 4.4.1):**
- Q_lk = 0.6×2×{lm1_loading.LM1_CHAR[1]['Q_ik']:.0f} + 0.10×{lm1_loading.LM1_CHAR[1]['q_ik']:.1f}×{lane_width:.2f}×{B_ext:.2f} = {lm1.braking.Q_lk_raw:.1f} kN
- Clamped [180, 900] kN → **{lm1.braking.Q_lk:.1f} kN**
- Per metre strip (÷ {lane_width:.2f} m lane width): **{lm1.braking.Q_lk_per_m:.2f} kN/m**
""")

    with _lm1_right:
        st.markdown("**Worst-case tandem position — Lane 1 (travel / B_ext direction)**")

        # --- Vehicle position diagram (in travel direction, B_ext on x-axis) ---
        fig_v, ax_v = plt.subplots(figsize=(5, 4))
        ax_v.set_facecolor("white")

        pad_v = max(B_ext * 0.35, 0.6)
        x0_v  = -pad_v
        x1_v  =  B_ext + pad_v
        y_bot = -(H_c + 0.3)    # show a little below crown

        # Fill layers behind (in travel direction cross-section)
        ax_v.add_patch(patches.Rectangle((x0_v, y_bot), x1_v - x0_v, -y_bot,
                                          fc="#C8A86E", ec="none", zorder=0))
        if t_sub > 0:
            ax_v.add_patch(patches.Rectangle((x0_v, -(t_road + t_sub)), x1_v - x0_v, t_sub,
                                              fc="#B0B0B0", ec="none", zorder=1))
        if t_road > 0:
            ax_v.add_patch(patches.Rectangle((x0_v, -t_road), x1_v - x0_v, t_road,
                                              fc="#404040", ec="none", zorder=1))

        # Culvert top slab (just the top edge in this section)
        ax_v.add_patch(patches.Rectangle((0, -H_c), B_ext, t_s,
                                          fc="#BBBBBB", ec="#333333", lw=1.5, zorder=2))

        # GL line
        ax_v.axhline(0, color="#3E1E00", lw=2, zorder=5)
        ax_v.text(x1_v - 0.05, 0.04, "GL", color="#3E1E00", fontsize=7,
                  ha="right", va="bottom", fontweight="bold")

        # Crown line
        ax_v.axhline(-H_c, color="#333333", lw=0.8, ls=":", zorder=3)
        ax_v.text(x1_v - 0.05, -H_c + 0.02, "Crown", color="#333333",
                  fontsize=6, ha="right", va="bottom")

        # Tandem worst position: centred over culvert in B_ext direction
        ctr  = B_ext / 2
        ax1  = ctr - lm1_loading.AXLE_SPACING / 2   # front axle x-centre
        ax2  = ctr + lm1_loading.AXLE_SPACING / 2   # rear  axle x-centre
        cl   = lm1_loading.CONTACT_L                 # contact patch length
        wh   = 0.08                                   # wheel drawing height (visual)

        for ax_x in [ax1, ax2]:
            # Wheel contact patch
            ax_v.add_patch(patches.Rectangle((ax_x - cl/2, 0), cl, wh,
                                              fc="#222222", ec="#111111", lw=0.8, zorder=6))
            # Dispersion lines at 1:1 through H_c
            left_edge  = ax_x - cl/2
            right_edge = ax_x + cl/2
            ax_v.plot([left_edge,  left_edge  - H_c], [0, -H_c], "r--", lw=0.8, zorder=4)
            ax_v.plot([right_edge, right_edge + H_c], [0, -H_c], "r--", lw=0.8, zorder=4)

        # Dispersed footprint at crown level
        left_disp  = ctr - dg.disp_Bext / 2
        right_disp = ctr + dg.disp_Bext / 2
        ax_v.annotate("", xy=(min(right_disp, B_ext), -H_c - 0.06),
                      xytext=(max(left_disp, 0), -H_c - 0.06),
                      arrowprops=dict(arrowstyle="<->", color="red", lw=1.0))
        ax_v.text(ctr, -H_c - 0.14, f"disp={dg.disp_Bext:.2f}m",
                  ha="center", fontsize=6, color="red")

        # Axle load labels — Q1k is per axle; tandem = 2 axles = 2 × Q1k
        Q1k = lm1_loading.LM1_CHAR[1]["Q_ik"]
        for label, ax_x in zip(["Axle 1", "Axle 2"], [ax1, ax2]):
            ax_v.text(ax_x, wh + 0.04, f"{label}\n{Q1k:.0f} kN",
                      ha="center", va="bottom", fontsize=5.5, color="#222222", fontweight="bold",
                      linespacing=1.3)
        # Total tandem label centred between axles
        ax_v.text(ctr, wh + 0.32, f"Tandem: 2 × {Q1k:.0f} = {2*Q1k:.0f} kN total",
                  ha="center", fontsize=6, color="#222222",
                  bbox=dict(fc="lightyellow", ec="#AAAAAA", lw=0.5, pad=2))

        # Axle spacing dimension (below wheel patches)
        ax_v.annotate("", xy=(ax2, -0.06), xytext=(ax1, -0.06),
                      arrowprops=dict(arrowstyle="<->", color="#555555", lw=0.8))
        ax_v.text(ctr, -0.11, f"{lm1_loading.AXLE_SPACING:.1f} m",
                  ha="center", fontsize=6, color="#555555")

        # B_ext span arrow
        ax_v.annotate("", xy=(B_ext, y_bot + 0.05), xytext=(0, y_bot + 0.05),
                      arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
        ax_v.text(ctr, y_bot + 0.10, f"B_ext={B_ext:.2f}m", ha="center", fontsize=6.5)

        # H_c depth label
        ax_v.annotate("", xy=(x0_v + 0.15, -H_c), xytext=(x0_v + 0.15, 0),
                      arrowprops=dict(arrowstyle="<->", color="black", lw=0.8))
        ax_v.text(x0_v + 0.22, -H_c / 2, f"Hc={H_c:.3f}m", fontsize=6, va="center")

        ax_v.set_xlim(x0_v, x1_v)
        ax_v.set_ylim(y_bot, wh + 0.50)
        ax_v.set_xlabel("B_ext direction (m)", fontsize=8)
        ax_v.set_ylabel("Depth (m)", fontsize=8)
        ax_v.tick_params(labelsize=7)
        fig_v.tight_layout()
        st.pyplot(fig_v, use_container_width=True)
        plt.close()

    # ── LL direction diagram ──────────────────────────────────────────────────
    st.markdown("---")
    st.markdown("**Lane distribution — LL direction** (barrel axis; vehicles travel into page)")

    # Dynamic bottom margin to fit all staggered dispersion arrows
    _arr_space = 0.10 * len(lm1.lanes) + 0.18
    _y_bot_ll  = -(H_c + _arr_space + 0.22)
    _ll_cw_x0  = LL / 2 - cw_width / 2   # carriageway left edge in LL

    fig_ll, ax_ll = plt.subplots(figsize=(10, 3.8))
    ax_ll.set_facecolor("white")

    # Background fill layers
    ax_ll.add_patch(patches.Rectangle((0, _y_bot_ll), LL, abs(_y_bot_ll),
                                      fc="#C8A86E", ec="none", zorder=0))
    if t_sub > 0:
        ax_ll.add_patch(patches.Rectangle((0, -(t_road + t_sub)), LL, t_sub,
                                          fc="#B0B0B0", ec="none", zorder=1))
    if t_road > 0:
        ax_ll.add_patch(patches.Rectangle((0, -t_road), LL, t_road,
                                          fc="#404040", ec="none", zorder=1))

    # Culvert top slab (full LL width)
    ax_ll.add_patch(patches.Rectangle((0, -H_c), LL, t_s,
                                      fc="#BBBBBB", ec="#333333", lw=1.2, zorder=2))

    # GL and crown reference lines
    ax_ll.axhline(0,    color="#3E1E00", lw=2,   zorder=5)
    ax_ll.axhline(-H_c, color="#333333", lw=0.8, ls=":", zorder=3)
    ax_ll.text(LL + 0.06, 0.03,   "GL",    color="#3E1E00", fontsize=7, ha="left", va="bottom")
    ax_ll.text(LL + 0.06, -H_c + 0.02, "Crown", color="#333333", fontsize=6, ha="left", va="bottom")

    _wh_ll  = 0.07   # wheel patch visual height (not structural)
    _lclrs  = ["#5A8FD0", "#5CA06A", "#C86060", "#9A70C0"]

    for _i, _ln in enumerate(lm1.lanes):
        _lctr  = _ll_cw_x0 + (_i + 0.5) * lane_width
        _lleft = _ll_cw_x0 + _i * lane_width
        _lc    = _lclrs[_i % len(_lclrs)]

        # Lane strip above GL (colour tint so lanes are distinguishable)
        ax_ll.add_patch(patches.Rectangle((_lleft, 0), lane_width, 0.11,
                                          fc=_lc, ec="white", lw=0.5, alpha=0.75, zorder=4))
        ax_ll.text(_lctr, 0.055, f"Lane {_ln.lane}  —  {_ln.Q_ik:.0f} kN/axle",
                   ha="center", va="center", fontsize=5.5, color="white",
                   fontweight="bold", zorder=7)

        # Lane boundary
        ax_ll.axvline(_lleft, color="#999999", lw=0.5, ls=":", zorder=3)

        # Wheel patches and dispersion for this lane
        _w_l = _lctr - lm1_loading.WHEEL_SPACING / 2
        _w_r = _lctr + lm1_loading.WHEEL_SPACING / 2
        for _wx in [_w_l, _w_r]:
            _wx0 = _wx - lm1_loading.CONTACT_T / 2
            _wx1 = _wx + lm1_loading.CONTACT_T / 2
            ax_ll.add_patch(patches.Rectangle((_wx0, 0), lm1_loading.CONTACT_T, _wh_ll,
                                              fc="#222222", ec="#111111", lw=0.6, zorder=6))
            ax_ll.plot([_wx0, _wx0 - H_c], [0, -H_c], color=_lc, lw=0.9, ls="--",
                       alpha=0.80, zorder=5)
            ax_ll.plot([_wx1, _wx1 + H_c], [0, -H_c], color=_lc, lw=0.9, ls="--",
                       alpha=0.80, zorder=5)

        # Staggered dispersion arrow at crown level
        _dsp_l = _lctr - dg.disp_LL / 2
        _dsp_r = _lctr + dg.disp_LL / 2
        _y_arr = -H_c - 0.09 - _i * 0.10
        ax_ll.annotate("", xy=(_dsp_r, _y_arr), xytext=(_dsp_l, _y_arr),
                        arrowprops=dict(arrowstyle="<->", color=_lc, lw=1.0))
        ax_ll.text(_lctr, _y_arr - 0.05,
                   f"disp_LL = {dg.disp_LL:.2f} m",
                   ha="center", fontsize=5, color=_lc)

    # Right carriageway edge
    ax_ll.axvline(_ll_cw_x0 + cw_width, color="#999999", lw=0.5, ls=":", zorder=3)

    # Wheel-to-wheel spacing annotation (Lane 1 only, above wheels)
    _ln1_ctr = _ll_cw_x0 + 0.5 * lane_width
    ax_ll.annotate("",
                   xy  =(_ln1_ctr + lm1_loading.WHEEL_SPACING / 2, _wh_ll + 0.04),
                   xytext=(_ln1_ctr - lm1_loading.WHEEL_SPACING / 2, _wh_ll + 0.04),
                   arrowprops=dict(arrowstyle="<->", color="#333333", lw=0.9))
    ax_ll.text(_ln1_ctr, _wh_ll + 0.08,
               f"Wheel spacing = {lm1_loading.WHEEL_SPACING:.1f} m",
               ha="center", fontsize=5.5, color="#333333")

    # 1m strip indicator (centred on Lane 1 tandem centroid)
    ax_ll.add_patch(patches.Rectangle(
        (_ln1_ctr - 0.5, _y_bot_ll), 1.0, abs(_y_bot_ll) + 0.13,
        fc="none", ec="darkorange", lw=1.3, ls="-", alpha=0.8, zorder=8))
    ax_ll.text(_ln1_ctr, _y_bot_ll + 0.05,
               "1 m strip (worst)", ha="center", fontsize=5.5, color="darkorange")

    # H_c depth annotation
    ax_ll.annotate("", xy=(-0.3, -H_c), xytext=(-0.3, 0),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=0.8))
    ax_ll.text(-0.35, -H_c / 2, f"Hc\n{H_c:.3f} m",
               fontsize=5.5, va="center", ha="right", linespacing=1.2)

    # Carriageway width annotation
    ax_ll.annotate("",
                   xy  =(_ll_cw_x0 + cw_width, 0.18),
                   xytext=(_ll_cw_x0, 0.18),
                   arrowprops=dict(arrowstyle="<->", color="#CCCCCC", lw=0.9))
    ax_ll.text(LL / 2, 0.21,
               f"Carriageway = {cw_width:.2f} m  ({n_lanes} notional lanes × {lane_width:.2f} m)",
               ha="center", fontsize=6, color="#CCCCCC")

    # LL total span annotation
    ax_ll.annotate("", xy=(LL, _y_bot_ll + 0.07), xytext=(0, _y_bot_ll + 0.07),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax_ll.text(LL / 2, _y_bot_ll + 0.14, f"LL = {LL:.1f} m", ha="center", fontsize=6.5)

    ax_ll.set_xlim(-0.7, LL + 0.7)
    ax_ll.set_ylim(_y_bot_ll, 0.28)
    ax_ll.set_xlabel("LL direction — barrel axis (m)", fontsize=8)
    ax_ll.set_ylabel("Depth (m)", fontsize=8)
    ax_ll.tick_params(labelsize=7)
    fig_ll.tight_layout()
    st.pyplot(fig_ll, use_container_width=True)
    plt.close()

# ══ Section 2b: LM3 Special Vehicle Loading ═══════════════════════════════════
with st.expander("LM3 Special Vehicle Loading", expanded=False):
    sv = lm3_loading.SV_VEHICLES[lm3.vehicle_name]
    st.warning("SV axle data is **approximate** — verify all values against UK NA Table NA.5 before use in design.")

    col_l3, col_r3 = st.columns([1, 1])

    with col_l3:
        st.markdown(f"""
**{lm3.vehicle_name}** — {lm3.gvw:.0f} kN GVW · {len(sv.axle_loads)} axles · {len(sv.axle_loads)} × {sv.axle_loads[0]:.0f} kN/axle

**Per-axle dispersion through {H_c:.3f} m cover (1:1):**
- LL direction:   {lm3.dispersion.disp_LL:.3f} m &nbsp;=&nbsp; {lm1_loading.WHEEL_SPACING:.1f} + {lm1_loading.CONTACT_T:.1f} + 2×{H_c:.3f}
- B_ext per axle: {lm3.dispersion.disp_B:.3f} m &nbsp;=&nbsp; {lm1_loading.CONTACT_L:.1f} + 2×{H_c:.3f}

**Worst-case position scan (0.05 m step):**
- Front-axle offset from culvert edge: **{lm3.sv_worst_offset:+.3f} m**
- Axles contributing at worst position: **{lm3.sv_n_axles}**
- SV load (Lane 1) per metre strip: &nbsp;**{lm3.sv_load_per_m:.2f} kN/m**
""")

        if lm3.secondary_lanes:
            sec_rows = []
            for ln in lm3.secondary_lanes:
                sec_rows.append({
                    "Lane": f"Lane {ln.lane}",
                    "Q_ik (kN)":      f"{ln.Q_ik:.0f}",
                    "q_ik (kN/m²)":   f"{ln.q_ik:.1f}",
                    "UDL/m (kN/m)":   f"{ln.udl_per_m:.2f}",
                    "TS/m (kN/m)":    f"{ln.ts_per_m:.2f}",
                    "Total/m (kN/m)": f"{ln.total_per_m:.2f}",
                })
            st.markdown("**Secondary lanes (LM1 Lane 2 / 3 / 4 TS + UDL):**")
            st.table(pd.DataFrame(sec_rows).set_index("Lane"))
        else:
            st.markdown("*No secondary lanes (only 1 notional lane).*")

        st.markdown(f"""
---
- LM1 secondary lanes total: &nbsp;**{lm3.secondary_per_m:.2f} kN/m**
- **Max vertical (SV + secondary): {lm3.max_V_per_m:.2f} kN/m**
- Min vertical (no LL): **0.00 kN/m**
""")

    with col_r3:
        # Diagram: SV vehicle worst-case position in B_ext direction
        st.markdown(f"**Worst-case axle positions — {lm3.vehicle_name} (B_ext direction)**")

        fig_sv, ax_sv = plt.subplots(figsize=(5, 4))
        ax_sv.set_facecolor("white")

        pad_sv  = max(B_ext * 0.35, 0.6)
        x0_sv   = -pad_sv
        x1_sv   = B_ext + pad_sv
        y_bot_sv = -(H_c + 0.3)

        # Fill layers
        ax_sv.add_patch(patches.Rectangle((x0_sv, y_bot_sv), x1_sv - x0_sv, abs(y_bot_sv),
                                          fc="#C8A86E", ec="none", zorder=0))
        if t_sub > 0:
            ax_sv.add_patch(patches.Rectangle((x0_sv, -(t_road + t_sub)), x1_sv - x0_sv, t_sub,
                                              fc="#B0B0B0", ec="none", zorder=1))
        if t_road > 0:
            ax_sv.add_patch(patches.Rectangle((x0_sv, -t_road), x1_sv - x0_sv, t_road,
                                              fc="#404040", ec="none", zorder=1))

        # Culvert top slab
        ax_sv.add_patch(patches.Rectangle((0, -H_c), B_ext, t_s,
                                          fc="#BBBBBB", ec="#333333", lw=1.5, zorder=2))

        # GL and crown
        ax_sv.axhline(0,    color="#3E1E00", lw=2,   zorder=5)
        ax_sv.axhline(-H_c, color="#333333", lw=0.8, ls=":", zorder=3)
        ax_sv.text(x1_sv - 0.05, 0.04,     "GL",    color="#3E1E00", fontsize=7, ha="right", va="bottom")
        ax_sv.text(x1_sv - 0.05, -H_c+0.02,"Crown", color="#333333", fontsize=6,  ha="right", va="bottom")

        wh_sv = 0.07
        half_B_sv = lm3.dispersion.disp_B / 2.0
        offset = lm3.sv_worst_offset

        for ax_pos, ax_load in zip(sv.axle_pos, sv.axle_loads):
            centre = offset + ax_pos
            fp_l   = centre - half_B_sv
            fp_r   = centre + half_B_sv
            overlap = max(0.0, min(fp_r, B_ext) - max(fp_l, 0.0))
            on_culvert = overlap > 0.0

            colour = "#CC3300" if on_culvert else "#888888"
            # Wheel patch at GL level
            ax_sv.add_patch(patches.Rectangle((centre - lm1_loading.CONTACT_L/2, 0),
                                              lm1_loading.CONTACT_L, wh_sv,
                                              fc=colour, ec="#333333", lw=0.6, zorder=6))
            # Dispersion lines: start at contact-patch edges (GL), spread to crown footprint edges
            cl = lm1_loading.CONTACT_L
            ax_sv.plot([centre - cl/2, fp_l], [0, -H_c], color=colour, lw=0.7, ls="--", alpha=0.7, zorder=4)
            ax_sv.plot([centre + cl/2, fp_r], [0, -H_c], color=colour, lw=0.7, ls="--", alpha=0.7, zorder=4)

        # B_ext span arrow
        ax_sv.annotate("", xy=(B_ext, y_bot_sv + 0.05), xytext=(0, y_bot_sv + 0.05),
                       arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
        ax_sv.text(B_ext/2, y_bot_sv + 0.10, f"B_ext={B_ext:.2f}m", ha="center", fontsize=6.5)

        # H_c depth
        ax_sv.annotate("", xy=(x0_sv + 0.15, -H_c), xytext=(x0_sv + 0.15, 0),
                       arrowprops=dict(arrowstyle="<->", color="black", lw=0.8))
        ax_sv.text(x0_sv + 0.22, -H_c/2, f"Hc={H_c:.3f}m", fontsize=6, va="center")

        ax_sv.set_xlim(x0_sv, x1_sv)
        ax_sv.set_ylim(y_bot_sv, wh_sv + 0.45)
        ax_sv.set_xlabel("B_ext direction (m)", fontsize=8)
        ax_sv.set_ylabel("Depth (m)", fontsize=8)
        ax_sv.tick_params(labelsize=7)
        ax_sv.set_title(
            f"{lm3.vehicle_name} worst position (offset={lm3.sv_worst_offset:+.3f} m) — "
            f"red=contributing, grey=outside",
            fontsize=6.5, pad=4
        )
        fig_sv.tight_layout()
        st.pyplot(fig_sv, use_container_width=True)
        plt.close()

    # ── LL direction diagram ──────────────────────────────────────────────────
    st.markdown("---")
    st.markdown(
        f"**Lane distribution — LL direction** (barrel axis; {lm3.vehicle_name} travels into page)"
    )

    _lm3_n_disp   = lm3.n_lanes           # total lanes shown (SV + secondary)
    _lm3_arr_sp   = 0.10 * _lm3_n_disp + 0.18
    _lm3_y_bot    = -(H_c + _lm3_arr_sp + 0.22)
    _lm3_cw_x0    = LL / 2 - cw_width / 2
    _lm3_disp_LL  = lm3.dispersion.disp_LL   # same formula as LM1

    fig_ll3, ax_ll3 = plt.subplots(figsize=(10, 3.8))
    ax_ll3.set_facecolor("white")

    # Background fill layers
    ax_ll3.add_patch(patches.Rectangle(
        (0, _lm3_y_bot), LL, abs(_lm3_y_bot), fc="#C8A86E", ec="none", zorder=0))
    if t_sub > 0:
        ax_ll3.add_patch(patches.Rectangle(
            (0, -(t_road + t_sub)), LL, t_sub, fc="#B0B0B0", ec="none", zorder=1))
    if t_road > 0:
        ax_ll3.add_patch(patches.Rectangle(
            (0, -t_road), LL, t_road, fc="#404040", ec="none", zorder=1))

    # Culvert top slab
    ax_ll3.add_patch(patches.Rectangle(
        (0, -H_c), LL, t_s, fc="#BBBBBB", ec="#333333", lw=1.2, zorder=2))

    ax_ll3.axhline(0,    color="#3E1E00", lw=2,   zorder=5)
    ax_ll3.axhline(-H_c, color="#333333", lw=0.8, ls=":", zorder=3)
    ax_ll3.text(LL + 0.06, 0.03,      "GL",    color="#3E1E00", fontsize=7, ha="left", va="bottom")
    ax_ll3.text(LL + 0.06, -H_c+0.02, "Crown", color="#333333", fontsize=6,  ha="left", va="bottom")

    _wh_ll3 = 0.07
    _lclrs3 = ["#5A8FD0", "#5CA06A", "#C86060", "#9A70C0"]

    # Lane 1 — SV vehicle (i = 0)
    _sv_axle_kn  = sv.axle_loads[0]
    _ln1_ctr3    = _lm3_cw_x0 + 0.5 * lane_width
    _lc1_3       = _lclrs3[0]

    ax_ll3.add_patch(patches.Rectangle(
        (_lm3_cw_x0, 0), lane_width, 0.11,
        fc=_lc1_3, ec="white", lw=0.5, alpha=0.75, zorder=4))
    ax_ll3.text(_ln1_ctr3, 0.055,
                f"Lane 1 ({lm3.vehicle_name} — {_sv_axle_kn:.0f} kN/axle)",
                ha="center", va="center", fontsize=5.5, color="white",
                fontweight="bold", zorder=7)
    ax_ll3.axvline(_lm3_cw_x0, color="#999999", lw=0.5, ls=":", zorder=3)

    for _wx3 in [_ln1_ctr3 - lm1_loading.WHEEL_SPACING/2,
                 _ln1_ctr3 + lm1_loading.WHEEL_SPACING/2]:
        _wx0_3 = _wx3 - lm1_loading.CONTACT_T / 2
        _wx1_3 = _wx3 + lm1_loading.CONTACT_T / 2
        ax_ll3.add_patch(patches.Rectangle(
            (_wx0_3, 0), lm1_loading.CONTACT_T, _wh_ll3,
            fc="#222222", ec="#111111", lw=0.6, zorder=6))
        ax_ll3.plot([_wx0_3, _wx0_3 - H_c], [0, -H_c],
                    color=_lc1_3, lw=0.9, ls="--", alpha=0.80, zorder=5)
        ax_ll3.plot([_wx1_3, _wx1_3 + H_c], [0, -H_c],
                    color=_lc1_3, lw=0.9, ls="--", alpha=0.80, zorder=5)

    # Lane 1 dispersion arrow (i = 0, no stagger offset)
    _y_arr3_1 = -H_c - 0.09
    ax_ll3.annotate("",
                    xy=(_ln1_ctr3 + _lm3_disp_LL/2, _y_arr3_1),
                    xytext=(_ln1_ctr3 - _lm3_disp_LL/2, _y_arr3_1),
                    arrowprops=dict(arrowstyle="<->", color=_lc1_3, lw=1.0))
    ax_ll3.text(_ln1_ctr3, _y_arr3_1 - 0.05,
                f"disp_LL = {_lm3_disp_LL:.2f} m",
                ha="center", fontsize=5, color=_lc1_3)

    # Secondary lanes (LM1 Lane 2 / 3 / 4) — loop index i starts at 1
    for _ln3 in lm3.secondary_lanes:
        _i3    = _ln3.lane - 1   # 0-based (lane 2 → 1, lane 3 → 2, …)
        _lctr3 = _lm3_cw_x0 + (_i3 + 0.5) * lane_width
        _lft3  = _lm3_cw_x0 + _i3 * lane_width
        _lc3   = _lclrs3[_i3 % len(_lclrs3)]

        ax_ll3.add_patch(patches.Rectangle(
            (_lft3, 0), lane_width, 0.11,
            fc=_lc3, ec="white", lw=0.5, alpha=0.75, zorder=4))
        ax_ll3.text(_lctr3, 0.055,
                    f"Lane {_ln3.lane} (LM1 — {_ln3.Q_ik:.0f} kN/axle)",
                    ha="center", va="center", fontsize=5.5, color="white",
                    fontweight="bold", zorder=7)
        ax_ll3.axvline(_lft3, color="#999999", lw=0.5, ls=":", zorder=3)

        for _wxs in [_lctr3 - lm1_loading.WHEEL_SPACING/2,
                     _lctr3 + lm1_loading.WHEEL_SPACING/2]:
            _wx0_s = _wxs - lm1_loading.CONTACT_T / 2
            _wx1_s = _wxs + lm1_loading.CONTACT_T / 2
            ax_ll3.add_patch(patches.Rectangle(
                (_wx0_s, 0), lm1_loading.CONTACT_T, _wh_ll3,
                fc="#222222", ec="#111111", lw=0.6, zorder=6))
            ax_ll3.plot([_wx0_s, _wx0_s - H_c], [0, -H_c],
                        color=_lc3, lw=0.9, ls="--", alpha=0.80, zorder=5)
            ax_ll3.plot([_wx1_s, _wx1_s + H_c], [0, -H_c],
                        color=_lc3, lw=0.9, ls="--", alpha=0.80, zorder=5)

        _y_arr3_s = -H_c - 0.09 - _i3 * 0.10
        ax_ll3.annotate("",
                        xy=(_lctr3 + _lm3_disp_LL/2, _y_arr3_s),
                        xytext=(_lctr3 - _lm3_disp_LL/2, _y_arr3_s),
                        arrowprops=dict(arrowstyle="<->", color=_lc3, lw=1.0))
        ax_ll3.text(_lctr3, _y_arr3_s - 0.05,
                    f"disp_LL = {_lm3_disp_LL:.2f} m",
                    ha="center", fontsize=5, color=_lc3)

    # Right carriageway edge
    ax_ll3.axvline(_lm3_cw_x0 + cw_width, color="#999999", lw=0.5, ls=":", zorder=3)

    # Wheel spacing annotation — Lane 1
    ax_ll3.annotate("",
                    xy  =(_ln1_ctr3 + lm1_loading.WHEEL_SPACING/2, _wh_ll3 + 0.04),
                    xytext=(_ln1_ctr3 - lm1_loading.WHEEL_SPACING/2, _wh_ll3 + 0.04),
                    arrowprops=dict(arrowstyle="<->", color="#333333", lw=0.9))
    ax_ll3.text(_ln1_ctr3, _wh_ll3 + 0.08,
                f"Wheel spacing = {lm1_loading.WHEEL_SPACING:.1f} m",
                ha="center", fontsize=5.5, color="#333333")

    # 1m strip indicator — centred on Lane 1
    ax_ll3.add_patch(patches.Rectangle(
        (_ln1_ctr3 - 0.5, _lm3_y_bot), 1.0, abs(_lm3_y_bot) + 0.13,
        fc="none", ec="darkorange", lw=1.3, ls="-", alpha=0.8, zorder=8))
    ax_ll3.text(_ln1_ctr3, _lm3_y_bot + 0.05,
                "1 m strip (worst)", ha="center", fontsize=5.5, color="darkorange")

    # H_c depth annotation
    ax_ll3.annotate("", xy=(-0.3, -H_c), xytext=(-0.3, 0),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=0.8))
    ax_ll3.text(-0.35, -H_c/2, f"Hc\n{H_c:.3f} m",
                fontsize=5.5, va="center", ha="right", linespacing=1.2)

    # Carriageway width annotation
    ax_ll3.annotate("",
                    xy  =(_lm3_cw_x0 + cw_width, 0.18),
                    xytext=(_lm3_cw_x0, 0.18),
                    arrowprops=dict(arrowstyle="<->", color="#CCCCCC", lw=0.9))
    ax_ll3.text(LL/2, 0.21,
                f"Carriageway = {cw_width:.2f} m  ({n_lanes} notional lanes × {lane_width:.2f} m)",
                ha="center", fontsize=6, color="#CCCCCC")

    # LL span annotation
    ax_ll3.annotate("", xy=(LL, _lm3_y_bot + 0.07), xytext=(0, _lm3_y_bot + 0.07),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax_ll3.text(LL/2, _lm3_y_bot + 0.14, f"LL = {LL:.1f} m", ha="center", fontsize=6.5)

    ax_ll3.set_xlim(-0.7, LL + 0.7)
    ax_ll3.set_ylim(_lm3_y_bot, 0.28)
    ax_ll3.set_xlabel("LL direction — barrel axis (m)", fontsize=8)
    ax_ll3.set_ylabel("Depth (m)", fontsize=8)
    ax_ll3.tick_params(labelsize=7)
    fig_ll3.tight_layout()
    st.pyplot(fig_ll3, use_container_width=True)
    plt.close()

# ══ Section 3: Stability check results ════════════════════════════════════════
with st.expander("Stability Check Results", expanded=True):
    col_tbl_l, col_tbl_r = st.columns(2)

with col_tbl_l:
    # ── K-value parameter table ────────────────────────────────────────────────
    st.subheader("Limit State Parameters — PD6694-1 Annex B")
    param_df = pd.DataFrame({
        "Limit State":  LS_NAMES,
        "Ka (active)":  [LS[n]["Ka"]   for n in LS_NAMES],
        "Kmax (restrained)": [LS[n]["Kmax"] for n in LS_NAMES],
        "γG unfav.":    [LS[n]["gG_u"] for n in LS_NAMES],
        "γG fav.":      [LS[n]["gG_f"] for n in LS_NAMES],
        "γQ (traffic)": [LS[n]["gQ"]   for n in LS_NAMES],
        "γφ (resist.)": [LS[n]["g_phi"] for n in LS_NAMES],
    })
    st.dataframe(param_df, hide_index=True, use_container_width=True)

with col_tbl_r:
    # ── Utilisation ratio summary ──────────────────────────────────────────────
    st.subheader("Utilisation Ratios  Ed / Rd  (≤ 1.00 = PASS)")

    def fmt(ur):
        s = f"{ur:.3f}" if ur < 1e6 else "∞"
        return f"{s} {'✅' if ur <= 1.0 else '❌'}"

    rows = [
        ("Bearing",             "B.4 (Gk+Qk)",   "UR_B4_bear"),
        ("Overturning",         "B.4/B.5 (Gk)",  "UR_B4_ov"),
        ("Sliding (fric only)", "B.4/B.5 (Gk)",  "UR_B4_sl"),
        ("Bearing",             "B.5 (Gk only)",  "UR_B5_bear"),
        ("Bearing",             "B.6 (Gk+Qk)",   "UR_B6_bear"),
        ("Overturning",         "B.6 (Gk)",      "UR_B6_ov"),
        ("Sliding (fric+Kr)",   "B.6 (Gk)",      "UR_B6_sl"),
    ]
    summary_data = {
        "Check":     [r[0] for r in rows],
        "Case":      [r[1] for r in rows],
        **{n: [fmt(res[n][r[2]]) for r in rows] for n in LS_NAMES},
    }
    st.dataframe(pd.DataFrame(summary_data), hide_index=True, use_container_width=True)

    # Overall verdict
    all_keys = ("UR_B4_bear", "UR_B4_ov", "UR_B4_sl",
                "UR_B5_bear", "UR_B6_bear", "UR_B6_ov", "UR_B6_sl")
    all_urs = [res[n][k] for n in LS_NAMES for k in all_keys]
    if all(u <= 1.0 for u in all_urs if u < 1e6):
        st.success("**PASS** — All checks satisfied across all four limit states.")
    else:
        fails = []
        check_map = [("Bearing B.4",    "UR_B4_bear"), ("Overturning B.4/B.5", "UR_B4_ov"),
                     ("Sliding B.4/B.5","UR_B4_sl"),   ("Bearing B.5",         "UR_B5_bear"),
                     ("Bearing B.6",    "UR_B6_bear"), ("Overturning B.6",     "UR_B6_ov"),
                     ("Sliding B.6",   "UR_B6_sl")]
        for n in LS_NAMES:
            for label, key in check_map:
                if res[n][key] > 1.0:
                    fails.append(f"{label} [{n}]")
        st.error("**FAIL** — " + " · ".join(fails))

# ══ Section 4: Detailed calculations ══════════════════════════════════════════
with st.expander("Detailed Calculations", expanded=False):
    tabs = st.tabs(LS_NAMES)
for tab, name in zip(tabs, LS_NAMES):
    r = res[name]
    with tab:
        st.markdown(f"""
**Limit state parameters** · Ka = {r['Ka']} · Kmax = {r['Kmax']} · γG_u = {r['gG_u']:.2f} · γG_f = {r['gG_f']:.2f} · γQ = {r['gQ']:.2f} · γφ = {r['g_phi']:.2f}

**Characteristic vertical loads**
- N_G,k (permanent) = {N_G_k:.2f} kN/m &nbsp;·&nbsp; Q_v,k (LM1 traffic, all lanes) = {Q_vk:.2f} kN/m

**Factored vertical loads**
- V_u (B.4 bearing, max) = {r['gG_u']:.2f}×{N_G_k:.2f} + {r['gQ']:.2f}×{Q_vk:.2f} = **{r['V_u']:.2f} kN/m**
- V_f (B.5 resist., min) = {r['gG_f']:.2f}×{N_G_k:.2f} + 0×{Q_vk:.2f} &nbsp;= **{r['V_f']:.2f} kN/m** *(traffic excluded — variable action not favourable)*

**Horizontal earth pressure** (σ_top = {σ_top:.2f} kPa · σ_bot = {σ_bot:.2f} kPa)
- F_Ka   = ½ × {r['Ka']} × ({σ_top:.2f}+{σ_bot:.2f}) × {H_ext:.3f} = **{r['F_Ka']:.2f} kN/m** (active)
- F_Kmax = ½ × {r['Kmax']} × ({σ_top:.2f}+{σ_bot:.2f}) × {H_ext:.3f} = **{r['F_Kmax']:.2f} kN/m** (restrained)
- F_net (earth only) = {r['F_Kmax']:.2f} − {r['F_Ka']:.2f} = **{r['F_net']:.2f} kN/m**

**Table 6 horizontal traffic — active side (PD6694-1)**
- Ka_char = tan²(45−{φ_fill_deg:.1f}°/2) = **{Ka_char:.4f}**
- σh = 20 × {Ka_char:.4f} × {H_ext:.3f} m = **{Q_h_udl_k:.2f} kN/m** (arm {H_ext/2:.3f} m)
- F reduction (Note 5, Hc={H_c:.3f} m): r_F = **{_rF:.4f}** {"(Hc≥2m — F=0)" if H_c >= 2.0 else ""}
- F_line = 2×300×{Ka_char:.4f}×{_rF:.4f}/{lane_width:.2f} = **{Q_h_F_k:.2f} kN/m** (arm {H_ext:.3f} m)
- γQ × (σh + F_line) = {r['gQ']:.2f} × ({Q_h_udl_k:.2f}+{Q_h_F_k:.2f}) = **{r['F_h_tr']:.2f} kN/m** (factored)
- γQ × M_h = {r['gQ']:.2f} × {Q_h_M_k:.2f} = **{r['M_h_tr']:.2f} kNm/m** (factored driving moment)
""")

        c1, c2, c3 = st.columns(3)

        with c1:
            st.markdown(f"""
**Table B.4 — Max vertical (Gk + Qk)**

*Bearing*
- V_u = {r['gG_u']:.2f}×{N_G_k:.2f} + {r['gQ']:.2f}×{Q_vk:.2f} = **{r['V_u']:.2f} kN/m**
- q = V_u/B_ext = {r['V_u']:.2f}/{B_ext:.3f} = **{r['q_B4']:.2f} kPa**
- q_Rd = {q_Rd:.1f} kPa → **UR = {r['UR_B4_bear']:.3f}** {"✅" if r["UR_B4_bear"] <= 1.0 else "❌"}

*Overturning (min V = Gk only)*
- M_earth = {r['M_net_B4']:.2f} kNm/m &nbsp;·&nbsp; M_traffic = {r['M_h_tr']:.2f} kNm/m
- M_dst = {r['M_net_B4']:.2f} + {r['M_h_tr']:.2f} = **{r['M_drv_B45']:.2f} kNm/m**
- M_stb = {r['V_f']:.2f}×{B_ext/2:.3f} = **{r['M_stb']:.2f} kNm/m**
- **UR = {r['UR_B4_ov']:.3f}** {"✅" if r["UR_B4_ov"] <= 1.0 else "❌"}

*Sliding (min V = Gk only)*
- F_earth = {r['F_net']:.2f} kN/m &nbsp;·&nbsp; F_traffic = {r['F_h_tr']:.2f} kN/m
- F_drv = {r['F_net']:.2f} + {r['F_h_tr']:.2f} = **{r['F_drv_B45']:.2f} kN/m**
- R_fric = tan({r['φ_fnd_d_deg']:.2f}°)×{r['V_f']:.2f} + {r['c_fnd_d']:.2f}×{B_ext:.3f} = **{r['R_fric']:.2f} kN/m**
- **UR = {r['UR_B4_sl']:.3f}** {"✅" if r["UR_B4_sl"] <= 1.0 else "❌"}
""")

        with c2:
            st.markdown(f"""
**Table B.5 — Min vertical (Gk only)**

*Bearing*
- V_f = {r['gG_f']:.2f}×{N_G_k:.2f} = **{r['V_f']:.2f} kN/m** *(traffic excluded)*
- q = V_f/B_ext = {r['V_f']:.2f}/{B_ext:.3f} = **{r['q_B5']:.2f} kPa**
- q_Rd = {q_Rd:.1f} kPa → **UR = {r['UR_B5_bear']:.3f}** {"✅" if r["UR_B5_bear"] <= 1.0 else "❌"}

*Overturning*
- M_dst = {r['M_net_B4']:.2f} kN·m/m *(same K as B.4)*
- M_stb = **{r['M_stb']:.2f} kN·m/m**
- **UR = {r['UR_B5_ov']:.3f}** {"✅" if r["UR_B5_ov"] <= 1.0 else "❌"}

*Sliding (friction only)*
- Driving F_net = **{r['F_net']:.2f} kN/m**
- R_fric = **{r['R_fric']:.2f} kN/m**
- **UR = {r['UR_B5_sl']:.3f}** {"✅" if r["UR_B5_sl"] <= 1.0 else "❌"}
""")

        with c3:
            st.markdown(f"""
**Table B.6 — Ka driving / Kr passive (WT: {h_wt:.2f} m BGL)**

*Eff. stresses:* σ'top={r['σ_eff_top']:.2f} kPa, σ'bot={r['σ_eff_bot']:.2f} kPa · Uplift Uk={r['U_k']:.2f} kN/m

*Bearing (max V − uplift)*
- Vu B6 = {r['gG_u']:.2f}×{N_G_k:.2f} + {r['gQ']:.2f}×{Q_vk:.2f} − {r['gG_f']:.2f}×{r['U_k']:.2f} = **{r['V_u_B6']:.2f} kN/m**
- q = Vu B6 / Bext = **{r['q_B6']:.2f} kPa** → **UR = {r['UR_B6_bear']:.3f}** {"✅" if r["UR_B6_bear"] <= 1.0 else "❌"}

*Resistance vertical (min V − uplift)*
- Vf B6 = {r['gG_f']:.2f}×{N_G_k:.2f} − {r['gG_u']:.2f}×{r['U_k']:.2f} = **{r['V_f_B6']:.2f} kN/m**

*Overturning*
- M_Ka (earth) = {r['F_Ka_B6']:.2f}×{r['arm_Ka_B6']:.3f} = **{r['F_Ka_B6']*r['arm_Ka_B6']:.2f} kNm/m**
- M_traffic = **{r['M_h_tr']:.2f} kNm/m**
- M_Ka total = {r['F_Ka_B6']*r['arm_Ka_B6']:.2f}+{r['M_h_tr']:.2f} = **{r['M_Ka_B6']:.2f} kNm/m**
- M_Kr = {r['F_Kr']:.2f}×{r['arm_Kr']:.3f} = **{r['F_Kr']*r['arm_Kr']:.2f} kNm/m**
- M_net = max(0, {r['M_Ka_B6']:.2f}−{r['F_Kr']*r['arm_Kr']:.2f}) = **{r['M_net_B6']:.2f} kNm/m**
- M_stb = {r['V_f_B6']:.2f}×{B_ext/2:.3f} = **{r['M_stb_B6']:.2f} kNm/m**
- **UR = {r['UR_B6_ov']:.3f}** {"✅" if r["UR_B6_ov"] <= 1.0 else "❌"}

*Sliding (friction + passive Kr)*
- F_Ka (earth) = **{r['F_Ka_B6']:.2f} kN/m** · F_traffic = **{r['F_h_tr']:.2f} kN/m**
- F_drv = {r['F_Ka_B6']:.2f}+{r['F_h_tr']:.2f} = **{r['F_drv_B6']:.2f} kN/m**
- φ_fill_d = **{r['φ_fill_d_deg']:.2f}°** → Kp_d = **{r['Kp']:.3f}**
- F_Kr = **{r['F_Kr']:.2f} kN/m**
- R_fric_B6 = tan({r['φ_fnd_d_deg']:.2f}°)×{r['V_f_B6']:.2f} + {r['c_fnd_d']:.2f}×{B_ext:.3f} = **{r['R_fric_B6']:.2f} kN/m**
- R_B6 = {r['F_Kr']:.2f}+{r['R_fric_B6']:.2f} = **{r['R_B6']:.2f} kN/m**
- **UR = {r['UR_B6_sl']:.3f}** {"✅" if r["UR_B6_sl"] <= 1.0 else "❌"}
""")

# ── Geometry / assumptions note ────────────────────────────────────────────────
with st.expander("Geometry & assumptions"):
    st.markdown(f"""
| Parameter | Value |
|---|---|
| B_ext | {B:.2f} + 2×{t_w:.2f} = **{B_ext:.3f} m** |
| H_ext | {H:.2f} + 2×{t_s:.2f} = **{H_ext:.3f} m** |
| H_c | {t_road:.3f} + {t_sub:.3f} + {t_fill:.3f} = **{H_c:.3f} m** |
| H_inv | {H_c:.3f} + {H_ext:.3f} = **{H_inv:.3f} m** |
| Concrete area | **{A_conc:.4f} m²/m** |
| W_concrete | **{W_conc:.2f} kN/m** |
| W_road (cover) | {γ_road:.1f} × {t_road:.3f} × {B_ext:.3f} = **{W_road:.2f} kN/m** |
| W_subbase (cover) | {γ_sub:.1f} × {t_sub:.3f} × {B_ext:.3f} = **{W_sub:.2f} kN/m** |
| W_fill (cover) | {γ_fill:.1f} × {t_fill:.3f} × {B_ext:.3f} = **{W_fill:.2f} kN/m** |
| W_soil total | **{W_soil:.2f} kN/m** |
| N_G,k (permanent) | **{N_G_k:.2f} kN/m** |
| Q_v,k (LM1 all lanes) | **{Q_vk:.2f} kN/m** |
| σ_v at crown | {γ_road:.1f}×{t_road:.3f} + {γ_sub:.1f}×{t_sub:.3f} + {γ_fill:.1f}×{t_fill:.3f} = **{σ_top:.2f} kPa** |
| σ_v at invert | {σ_top:.2f} + {γ_fill:.1f}×{H_ext:.3f} = **{σ_bot:.2f} kPa** |

**Load case framework (PD6694-1 Annex B):**
- **B.4 bearing** — Max vertical (γG_u × N_G,k + γQ × Q_v,k); Ka on active wall → governs bearing pressure.
- **B.4/B.5 overturning/sliding** — Min vertical (γG_f × N_G,k, traffic excluded); Kmax driving → governs overturning and sliding.
- **B.5 bearing** — Min vertical (γG_f × N_G,k, traffic excluded); same horizontal → governs when uplift reduces base pressure.
- **B.6** — Effective stresses + water table uplift; Ka driving, Kr (Rankine passive) + base friction resisting.

**Traffic load integration:**
- Q_v,k = LM1 max vertical per metre strip (all lanes loaded simultaneously) — BS EN 1991-2.
- γQ = 1.35 at ULS (STR/GEO, EQU); 1.00 at SLS — BS EN 1990 Table A2.4.
- Traffic vertical is unfavourable for bearing (adds to V_u); excluded from OT/sliding resistance.
- **Horizontal traffic (PD6694-1 Table 6, Figure 2) — active side only:**
  - σh = 20 × Ka_char = {20*Ka_char:.2f} kN/m² over H_ext → resultant {Q_h_udl_k:.2f} kN/m (arm H_ext/2)
  - F line loads: 2 × 300 × Ka_char × r_F / lane_width = {Q_h_F_k:.2f} kN/m at crown (arm H_ext); r_F = {_rF:.3f} {"(Hc≥2m, F=0)" if H_c >= 2.0 else f"(Note 5, Hc={H_c:.2f}m)"}
  - Factored by γQ; added to driving force and moment in B.4/B.5 and B.6 OT/sliding checks.

**Other notes:**
- Ka and Kmax from PD6694-1 Tables B.4/B.5 include γM and γSd;K = 1.2 (Classes 6N/6P backfill assumed).
- Vertical stresses σ_v use characteristic γ (soil unit weight is not a partial-factored load in EC7).
- Kr (passive) computed from design φ_d = arctan(tan φ'_k / γφ); no γSd;K applied to resistance side.
- q_Rd applied uniformly across all limit states — user to verify suitability for each LS.
""")
