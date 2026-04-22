# app_starter.py
# This is the 'Day 1' starter file, to be ran until `app.py` is complete
import streamlit as st
import numpy as np
import plotly.graph_objects as go
import pandas as pd

from src.elements import compute_D
from src.mesh import generate_rect_mesh, generate_plate_with_hole_mesh

st.set_page_config(page_title="CST FEA Solver", layout="wide")
st.title("2D Plane Stress / Plane Strain FEA Solver")

with st.sidebar:
    st.header("Problem Setup")
    mode = st.selectbox("Analysis Mode", ["Plane Stress", "Plane Strain"])
    E = st.number_input("Young's Modulus E (Pa)", value=200e9, format="%.2e")
    nu = st.number_input("Poisson's Ratio ν", value=0.25, min_value=0.0,
                         max_value=0.4999, format="%.4f")
    if nu >= 0.499:
        st.warning("Near-incompressible — results may be unreliable with CST elements. "
                   "See the Locking tab.")
    t = st.number_input("Thickness (m)", value=0.01)
    # L = st.number_input("Plate Length L (m)", value=1.0)
    # h = st.number_input("Plate Height h (m)", value=0.25)
    # P = st.number_input("Tip Load P (N)", value=6000.0)
    # nx = st.slider("Elements in x", 2, 32, 8)
    # ny = st.slider("Elements in y", 2, 16, 4)
    # solve = st.button("Solve", type="primary")
mode_key = "plane_stress" if mode == "Plane Stress" else "plane_strain"

# ──────────────────────────────────────────────
# Tabs
# ──────────────────────────────────────────────
tab1, tab2, tab3 = st.tabs(["🔧 Cantilever Beam", "🕳️ Plate with Hole",
                             "📊 Convergence & Locking"])

# ══════════════════════════════════════════════
# Helper: plot mesh
# ══════════════════════════════════════════════
def plot_mesh(nodes, elements, title="Mesh Preview"):
    fig = go.Figure()
    for tri in elements:
        pts = nodes[list(tri) + [tri[0]]]
        fig.add_trace(go.Scatter(x=pts[:, 0], y=pts[:, 1],
                                 mode='lines', line=dict(color='#2138ab', width=0.5),
                                 showlegend=False))
    fig.update_layout(yaxis_scaleanchor="x", title=title,
                      height=350, margin=dict(l=0, r=0, t=30, b=0))
    return fig


def plot_deformed_mesh(nodes, elements, u, stresses, scale, title="Deformed Mesh"):
    vm = compute_von_mises(stresses)
    u_disp = u.reshape(-1, 2)
    nodes_def = nodes + scale * u_disp

    fig = go.Figure(go.Mesh3d(
        x=nodes_def[:, 0], y=nodes_def[:, 1], z=np.zeros(len(nodes_def)),
        i=elements[:, 0], j=elements[:, 1], k=elements[:, 2],
        intensity=vm, colorscale="Jet", showscale=True,
        colorbar=dict(title="σ_vm (Pa)"),
        flatshading=True,
    ))
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title="x (m)", yaxis_title="y (m)",
            zaxis=dict(visible=False),
            aspectmode='data',
            camera=dict(eye=dict(x=0, y=0, z=2)),
        ),
        height=500, margin=dict(l=0, r=0, t=30, b=0),
    )
    return fig

# ══════════════════════════════════════════════
# TAB 1: Cantilever Beam
# ══════════════════════════════════════════════
with tab1:
    st.subheader("Cantilever Plate — Parabolic Tip Shear")

    col_input, col_mesh = st.columns([1, 2])
    with col_input:
        L = st.number_input("Plate Length L (m)", value=1.0, key="cant_L")
        h = st.number_input("Plate Height h (m)", value=0.25, key="cant_h")
        P = st.number_input("Tip Load P (N)", value=6000.0, key="cant_P")
        nx = st.slider("Elements in x", 2, 32, 8, key="cant_nx")
        ny = st.slider("Elements in y", 2, 16, 4, key="cant_ny")
        solve_cant = st.button("Solve Cantilever", type="primary")

    # Mesh preview
    nodes_c, elems_c, tags_c = generate_rect_mesh(L, h, nx, ny)
    with col_mesh:
        st.plotly_chart(plot_mesh(nodes_c, elems_c,
                        f"Mesh: {len(nodes_c)} nodes · {len(elems_c)} elements"),
                        use_container_width=True)

# ══════════════════════════════════════════════
# TAB 2: Plate with Hole
# ══════════════════════════════════════════════
with tab2:
    st.subheader("Plate with Circular Hole — Kirsch Validation")

    col_input2, col_mesh2 = st.columns([1, 2])
    with col_input2:
        W = st.number_input("Half-width W (m)", value=5.0, key="hole_W")
        H = st.number_input("Half-height H (m)", value=5.0, key="hole_H")
        R = st.number_input("Hole radius R (m)", value=1.0, key="hole_R")
        sigma_inf = st.number_input("Far-field stress σ∞ (Pa)", value=1e6, format="%.2e",
                                     key="hole_sigma")
        n_rad = st.slider("Radial divisions", 4, 24, 10, key="hole_nrad")
        n_ang = st.slider("Angular divisions", 4, 24, 12, key="hole_nang")
        rho_mh = st.number_input(
            "Mesh density hear hole ρ_mh", value=1.0, key="rho_mh")
        solve_hole = st.button("Solve Plate with Hole", type="primary")

    # Mesh preview
    nodes_h, elems_h, tags_h = generate_plate_with_hole_mesh(
        W, H, R, n_rad, n_ang, rho_mh)
    with col_mesh2:
        st.plotly_chart(plot_mesh(nodes_h, elems_h,
                        f"Mesh: {len(nodes_h)} nodes · {len(elems_h)} elements"),
                        use_container_width=True)

# ══════════════════════════════════════════════
# TAB 3: Convergence & Locking
# ══════════════════════════════════════════════
with tab3:
    st.subheader("h-Refinement Convergence")

st.info("Configure inputs in the sidebar and click Solve.")
