import streamlit as st
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

st.set_page_config(page_title="Box Culvert Stability", layout="wide")

st.title("Buried Box Culvert — Stability Checker")
st.caption("Per-unit-length · PD6694-1 Annex B Tables B.4, B.5, B.6 · BS EN 1997-1 · No traffic loading")

# ── PD6694-1 Annex B ──────────────────────────────────────────────────────────
# Ka  : design active coefficient — Tables B.4/B.5/B.6 (includes γM · γSd;K = 1.2)
# Kmax: design max coefficient for restrained side — Tables B.4/B.5 (includes γM · γSd;K)
# Kr  : Rankine passive for B.6 sliding resistance — computed from design φ (not tabulated)
#
# EC7 partial factors
#   gG_u / gG_f: unfavourable / favourable on permanent actions
#   g_phi / g_c: on shear strength parameters (resistance side)
LS = {
    "SLS":        {"Ka": 0.33, "Kmax": 0.60, "gG_u": 1.00, "gG_f": 1.00, "g_phi": 1.00, "g_c": 1.00},
    "EQU":        {"Ka": 0.44, "Kmax": 0.60, "gG_u": 1.10, "gG_f": 0.90, "g_phi": 1.25, "g_c": 1.25},
    "STR/GEO C1": {"Ka": 0.40, "Kmax": 0.72, "gG_u": 1.35, "gG_f": 1.00, "g_phi": 1.00, "g_c": 1.00},
    "STR/GEO C2": {"Ka": 0.49, "Kmax": 0.84, "gG_u": 1.00, "gG_f": 1.00, "g_phi": 1.25, "g_c": 1.25},
}
LS_NAMES = list(LS.keys())

# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Geometry")
    B   = st.number_input("Internal Width, B (m)",   0.1, 20.0, 2.0,  0.1,  "%.2f")
    H   = st.number_input("Internal Height, H (m)",  0.1, 10.0, 1.5,  0.1,  "%.2f")
    t_w = st.number_input("Wall Thickness, t_w (m)", 0.05, 2.0, 0.25, 0.05, "%.2f")
    t_s = st.number_input("Slab Thickness, t_s (m)", 0.05, 2.0, 0.30, 0.05, "%.2f")

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

    st.header("Material")
    γ_c = st.number_input("Concrete Unit Weight γ_c (kN/m³)", 20.0, 26.0, 24.0, 0.5)

# ── Derived geometry ───────────────────────────────────────────────────────────
B_ext  = B + 2 * t_w
H_ext  = H + 2 * t_s
H_c    = t_road + t_sub + t_fill   # total cover depth to crown
H_inv  = H_c + H_ext

A_conc = B_ext * H_ext - B * H
W_conc = A_conc * γ_c
# Overburden weight: sum of each layer (characteristic, not factored)
W_road = γ_road * t_road * B_ext
W_sub  = γ_sub  * t_sub  * B_ext
W_fill = γ_fill * t_fill * B_ext
W_soil = W_road + W_sub + W_fill
N_v_k  = W_conc + W_soil

φ_fill_k = math.radians(φ_fill_deg)   # backfill — governs Kr passive resistance
φ_fnd_k  = math.radians(φ_fnd_deg)   # founding layer — governs base sliding resistance

# Vertical stress at crown = sum of overburden layers (used for horizontal earth pressure)
σ_top = γ_road * t_road + γ_sub * t_sub + γ_fill * t_fill
# Vertical stress at invert = crown stress + fill unit weight × structural height
# (material alongside culvert walls is fill)
σ_bot = σ_top + γ_fill * H_ext

def trapz_resultant(K, σ_t, σ_b, h):
    """Resultant of trapezoidal horizontal pressure K·σ over height h."""
    F   = 0.5 * K * (σ_t + σ_b) * h
    arm = h * (2*σ_t + σ_b) / (3*(σ_t + σ_b)) if (σ_t + σ_b) > 0 else h/3
    return F, arm

# ── Per-limit-state calculations ───────────────────────────────────────────────
def run(p):
    Ka, Kmax   = p["Ka"],   p["Kmax"]
    gG_u, gG_f = p["gG_u"], p["gG_f"]
    g_phi, g_c = p["g_phi"],p["g_c"]

    # Design vertical loads
    V_u = gG_u * N_v_k    # unfavourable — for bearing
    V_f = gG_f * N_v_k    # favourable  — for sliding/overturning resistance

    # Design soil parameters — backfill (passive wall Kr) and founding layer (base friction)
    φ_fill_d = math.atan(math.tan(φ_fill_k) / g_phi)
    c_fill_d = c_fill / g_c
    φ_fnd_d  = math.atan(math.tan(φ_fnd_k)  / g_phi)
    c_fnd_d  = c_fnd  / g_c

    # ── Horizontal force components ────────────────────────────────────────────
    F_Ka,   arm_Ka   = trapz_resultant(Ka,   σ_top, σ_bot, H_ext)
    F_Kmax, arm_Kmax = trapz_resultant(Kmax, σ_top, σ_bot, H_ext)

    F_net    = F_Kmax - F_Ka
    M_net_B4 = F_Kmax * arm_Kmax - F_Ka * arm_Ka

    # B.6 passive — Rankine Kp from design backfill φ
    Kp      = math.tan(math.pi/4 + φ_fill_d/2) ** 2
    F_Kr, _ = trapz_resultant(Kp, σ_top, σ_bot, H_ext)

    # Base friction — founding layer design parameters
    R_fric = math.tan(φ_fnd_d) * V_f + c_fnd_d * B_ext

    # ── Table B.4 — max vertical: Bearing (primary), Overturning, Sliding ────────
    q_d        = V_u / B_ext
    UR_B4_bear = q_d / q_Rd if q_Rd > 0 else float("inf")

    M_stb_B4   = V_f * B_ext / 2       # self-weight always favourable for OT stability
    UR_B4_ov   = M_net_B4 / M_stb_B4 if M_stb_B4 > 0 else float("inf")

    UR_B4_sl   = F_net / R_fric if R_fric > 0 else float("inf")

    # ── Table B.5 — min vertical: Overturning, Sliding (same K, different V for driving check)
    # K values identical to B.4; without traffic, stability results are the same
    # (min vertical only changes the resistance, which already uses V_f in B.4 above)
    UR_B5_ov   = UR_B4_ov
    UR_B5_sl   = UR_B4_sl

    # ── Table B.6 — sliding with passive soil resistance ──────────────────────
    R_B6     = F_Kr + R_fric           # passive soil + base friction
    UR_B6_sl = F_Ka / R_B6 if R_B6 > 0 else float("inf")

    return dict(
        Ka=Ka, Kmax=Kmax, gG_u=gG_u, gG_f=gG_f, g_phi=g_phi, g_c=g_c,
        V_u=V_u, V_f=V_f,
        φ_fill_d_deg=math.degrees(φ_fill_d), c_fill_d=c_fill_d,
        φ_fnd_d_deg=math.degrees(φ_fnd_d),   c_fnd_d=c_fnd_d,
        F_Ka=F_Ka, F_Kmax=F_Kmax, F_net=F_net,
        Kp=Kp, F_Kr=F_Kr, R_fric=R_fric, R_B6=R_B6,
        q_d=q_d,
        M_net_B4=M_net_B4, M_stb_B4=M_stb_B4,
        UR_B4_bear=UR_B4_bear, UR_B4_ov=UR_B4_ov, UR_B4_sl=UR_B4_sl,
        UR_B5_ov=UR_B5_ov,  UR_B5_sl=UR_B5_sl,
        UR_B6_sl=UR_B6_sl,
    )

res = {n: run(LS[n]) for n in LS_NAMES}

# ── Layout ─────────────────────────────────────────────────────────────────────
col_l, col_r = st.columns([1, 1.8])

# Cross-section diagram
with col_l:
    st.subheader("Cross-Section")
    pad  = max(B_ext * 0.6, 0.8)
    xlim = (-B_ext/2 - pad, B_ext/2 + pad)
    ylim = (-(H_inv + 0.7), 0.6)

    fig, ax = plt.subplots(figsize=(4, 5))
    ax.set_facecolor("white")

    # Fill layer (bottom of cover layers down to invert, and alongside culvert walls)
    ax.add_patch(patches.Rectangle(
        (xlim[0], ylim[0]), xlim[1]-xlim[0], -ylim[0],
        fc="#C8A86E", ec="none", zorder=0))

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

    def dim_arrow(ax, x, y0, y1, label, side="r"):
        ax.annotate("", xy=(x, y1), xytext=(x, y0),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
        dx = 0.10 if side == "r" else -0.10
        ax.text(x+dx, (y0+y1)/2, label, fontsize=6.5, va="center",
                ha="left" if side == "r" else "right")

    xr = B_ext/2 + 0.22
    xl = -B_ext/2 - 0.22
    if H_c > 0:
        dim_arrow(ax, xr, -H_c, 0, f"Hc={H_c:.2f} m")
    dim_arrow(ax, xr, -H_inv, -H_c, f"Hext={H_ext:.2f} m")
    dim_arrow(ax, xl, -H_inv, 0, f"Hinv={H_inv:.2f} m", "l")

    # Layer legend
    legend_patches = []
    if t_road > 0:
        legend_patches.append(patches.Patch(fc="#404040", label=f"Road ({t_road*1000:.0f} mm)"))
    if t_sub > 0:
        legend_patches.append(patches.Patch(fc="#B0B0B0", label=f"Subbase ({t_sub*1000:.0f} mm)"))
    legend_patches.append(patches.Patch(fc="#C8A86E", label=f"Fill ({t_fill:.2f} m)"))
    if legend_patches:
        ax.legend(handles=legend_patches, fontsize=6, loc="lower right")
    ax.annotate("", xy=(B_ext/2, -H_inv-0.32), xytext=(-B_ext/2, -H_inv-0.32),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text(0, -H_inv-0.50, f"Bext = {B_ext:.2f} m", ha="center", fontsize=6.5)

    ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_aspect("equal")
    ax.set_xlabel("(m)", fontsize=8); ax.set_ylabel("Depth below GL (m)", fontsize=8)
    ax.tick_params(labelsize=7); fig.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

# Results panel
with col_r:
    # ── K-value parameter table ────────────────────────────────────────────────
    st.subheader("Limit State Parameters — PD6694-1 Annex B")
    param_df = pd.DataFrame({
        "Limit State":  LS_NAMES,
        "Ka (active)":  [LS[n]["Ka"]   for n in LS_NAMES],
        "Kmax (restrained)": [LS[n]["Kmax"] for n in LS_NAMES],
        "γG unfav.":    [LS[n]["gG_u"] for n in LS_NAMES],
        "γG fav.":      [LS[n]["gG_f"] for n in LS_NAMES],
        "γφ (resist.)": [LS[n]["g_phi"] for n in LS_NAMES],
    })
    st.dataframe(param_df, hide_index=True, use_container_width=True)

    st.divider()

    # ── Utilisation ratio summary ──────────────────────────────────────────────
    st.subheader("Utilisation Ratios  Ed / Rd  (≤ 1.00 = PASS)")

    def fmt(ur):
        s = f"{ur:.3f}" if ur < 1e6 else "∞"
        return f"{s} {'✅' if ur <= 1.0 else '❌'}"

    rows = [
        ("Bearing",                    "B.4 (max V)",  "UR_B4_bear"),
        ("Overturning",                "B.4/B.5",      "UR_B4_ov"),
        ("Sliding — friction only",    "B.5",          "UR_B5_sl"),
        ("Sliding — friction + Kr",    "B.6",          "UR_B6_sl"),
    ]
    summary_data = {
        "Check":        [r[0] for r in rows],
        "Load case":    [r[1] for r in rows],
        **{n: [fmt(res[n][r[2]]) for r in rows] for n in LS_NAMES},
    }
    st.dataframe(pd.DataFrame(summary_data), hide_index=True, use_container_width=True)

    # Overall verdict
    all_urs = [res[n][key] for n in LS_NAMES
               for key in ("UR_B4_bear", "UR_B4_ov", "UR_B5_sl", "UR_B6_sl")]
    if all(u <= 1.0 for u in all_urs if u < 1e6):
        st.success("**PASS** — All checks satisfied across all four limit states.")
    else:
        fails = []
        check_map = [("Bearing (B.4)", "UR_B4_bear"), ("Overturning (B.4/B.5)", "UR_B4_ov"),
                     ("Sliding B.5",   "UR_B5_sl"),   ("Sliding B.6",           "UR_B6_sl")]
        for n in LS_NAMES:
            for label, key in check_map:
                if res[n][key] > 1.0:
                    fails.append(f"{label} [{n}]")
        st.error("**FAIL** — " + " · ".join(fails))

# ── Detailed calculations per limit state ──────────────────────────────────────
st.subheader("Detailed Calculations")
tabs = st.tabs(LS_NAMES)
for tab, name in zip(tabs, LS_NAMES):
    r = res[name]
    with tab:
        st.markdown(f"""
**Limit state parameters** · Ka = {r['Ka']} · Kmax = {r['Kmax']} · γG_u = {r['gG_u']:.2f} · γG_f = {r['gG_f']:.2f} · γφ = {r['g_phi']:.2f}

**Vertical loads** (characteristic N_v_k = {N_v_k:.2f} kN/m)
- V_d unfav. = {r['gG_u']:.2f} × {N_v_k:.2f} = **{r['V_u']:.2f} kN/m**
- V_d fav.   = {r['gG_f']:.2f} × {N_v_k:.2f} = **{r['V_f']:.2f} kN/m**

**Horizontal forces** (σ_top = {σ_top:.2f} kPa · σ_bot = {σ_bot:.2f} kPa)
- F_Ka   = ½ × {r['Ka']} × ({σ_top:.2f}+{σ_bot:.2f}) × {H_ext:.3f} = **{r['F_Ka']:.2f} kN/m** (active side)
- F_Kmax = ½ × {r['Kmax']} × ({σ_top:.2f}+{σ_bot:.2f}) × {H_ext:.3f} = **{r['F_Kmax']:.2f} kN/m** (restrained side)
- **F_net = F_Kmax − F_Ka = {r['F_net']:.2f} kN/m** (B.4/B.5 net driving force)
""")

        c1, c2, c3, c4 = st.columns(4)

        with c1:
            st.markdown(f"""
**Bearing (B.4)**
- q_d = V_u / B_ext
  = {r['V_u']:.2f} / {B_ext:.3f}
  = **{r['q_d']:.2f} kPa**
- q_Rd = {q_Rd:.1f} kPa
- **UR = {r['UR_B4_bear']:.3f}** {"✅" if r["UR_B4_bear"] <= 1.0 else "❌"}
""")

        with c2:
            st.markdown(f"""
**Overturning (B.4/B.5)**
- M_dst = {r['M_net_B4']:.2f} kN·m/m
  *(net Kmax−Ka moment)*
- M_stb = V_f × B/2
  = {r['V_f']:.2f} × {B_ext/2:.3f}
  = **{r['M_stb_B4']:.2f} kN·m/m**
- **UR = {r['UR_B4_ov']:.3f}** {"✅" if r["UR_B4_ov"] <= 1.0 else "❌"}
""")

        with c3:
            st.markdown(f"""
**Sliding — B.5**
*(base friction only)*
- Driving F_net = **{r['F_net']:.2f} kN/m**
- φ_fnd_d = arctan(tan{φ_fnd_deg:.0f}°/{r['g_phi']:.2f}) = **{r['φ_fnd_d_deg']:.2f}°**
- c_fnd_d = {c_fnd:.1f}/{r['g_c']:.2f} = **{r['c_fnd_d']:.2f} kPa**
- R_fric = tan(φ_fnd_d)×V_f + c_fnd_d×B_ext
  = **{r['R_fric']:.2f} kN/m**
- **UR = {r['UR_B5_sl']:.3f}** {"✅" if r["UR_B5_sl"] <= 1.0 else "❌"}
""")

        with c4:
            st.markdown(f"""
**Sliding — B.6**
*(base friction + passive Kr)*
- Driving F_Ka = **{r['F_Ka']:.2f} kN/m**
- φ_fill_d = arctan(tan{φ_fill_deg:.0f}°/{r['g_phi']:.2f}) = **{r['φ_fill_d_deg']:.2f}°**
- Kp_d = tan²(45+φ_fill_d/2) = **{r['Kp']:.3f}**
- F_Kr = **{r['F_Kr']:.2f} kN/m**
- R_B6 = F_Kr + R_fric
  = {r['F_Kr']:.2f} + {r['R_fric']:.2f}
  = **{r['R_B6']:.2f} kN/m**
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
| N_v,k | **{N_v_k:.2f} kN/m** |
| σ_v at crown | {γ_road:.1f}×{t_road:.3f} + {γ_sub:.1f}×{t_sub:.3f} + {γ_fill:.1f}×{t_fill:.3f} = **{σ_top:.2f} kPa** |
| σ_v at invert | {σ_top:.2f} + {γ_fill:.1f}×{H_ext:.3f} = **{σ_bot:.2f} kPa** |

**Load case framework (PD6694-1 Annex B):**
- **B.4** — Maximum vertical load; Ka (active) on one wall, Kmax (restrained) on other → Bearing (primary) + overturning/sliding.
- **B.5** — Minimum vertical load; same Ka/Kmax as B.4 → Overturning and sliding (no traffic: results same as B.4 stability).
- **B.6** — Ka (active) on one wall, Kr (Rankine passive from design φ) on other + base friction → Sliding stability.

**Notes:**
- Ka and Kmax from PD6694-1 Tables B.4/B.5 include γM and γSd;K = 1.2 (Classes 6N/6P backfill assumed).
- Horizontal traffic surcharge omitted (no traffic loading in this analysis).
- Vertical stresses σ_v use characteristic γ (soil unit weight is not a partial-factored load in EC7).
- Kr (passive) computed from design φ_d = arctan(tan φ'_k / γφ); no γSd;K applied to resistance side.
- q_Rd applied uniformly across all limit states — user to verify suitability for each LS.
""")
