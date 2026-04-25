# src/mesh.py
import numpy as np
import os


def generate_rect_mesh(L, h, nx, ny):
    """
    Generate a structured triangular mesh over a rectangular domain [0, L] x [-h/2, h/2].

    Each rectangular cell is split into 2 triangles. The split direction must be
    consistent across the mesh (e.g., always along the same diagonal).

    Parameters
    ----------
    L : float
        Length of the plate (x-direction).
    h : float
        Height of the plate (y-direction).
    nx : int
        Number of element divisions in x.
    ny : int
        Number of element divisions in y.

    Returns
    -------
    nodes : ndarray, shape (n_nodes, 2)
        Nodal coordinates.
    elements : ndarray, shape (n_elems, 3)
        Element connectivity (node indices per triangle).
    boundary_tags : dict
        Must contain at minimum:
        - 'fixed': list of node indices on the x=0 boundary (cantilever root)
        - 'loaded': list of node indices on the x=L boundary (cantilever tip)
    """
    nx_nodes = nx + 1
    ny_nodes = ny + 1
    n_nodes = nx_nodes * ny_nodes

    dx = L / nx
    dy = h / ny

    nodes = np.array([
        (i_x*dx, -h/2 + j_y*dy)
        for j_y in range(ny_nodes) for i_x in range(nx_nodes)])

    elements = []
    for row in range(ny_nodes-1):
        for col in range(nx_nodes-1):
            i = nx_nodes*row + col
            j = i + 1
            k = i + nx_nodes
            l = k + 1

            elements += [[i, j, k], [l, k, j]]
    elements = np.array(elements)

    boundary_tags = {
        "fixed" : [n for n in range(n_nodes) if nodes[n][0] == 0],
        "loaded" : [n for n in range(n_nodes) if nodes[n][0] == L]}

    return (nodes, elements, boundary_tags)


def generate_plate_with_hole_mesh(W, H, R, n_radial, n_angular, rho_mh=1.0):
    """
    Generate a triangular mesh for a rectangular plate [0, W] x [0, H]
    with a circular hole of radius R centered at the origin.

    Due to symmetry, only the quarter-plate (first quadrant) needs to be meshed.
    The hole boundary and plate edges must be tagged for boundary conditions.

    Parameters
    ----------
    W : float
        Half-width of the plate (x-direction extent from center).
    H : float
        Half-height of the plate (y-direction extent from center).
    R : float
        Radius of the circular hole.
    n_radial : int
        Number of element divisions in the radial direction (from hole to edge).
    n_angular : int
        Number of element divisions in the angular direction (quarter circle).
    rho_mh : float
        Mesh density factor near the hole (1 = not denser near the hole, >1
        denser near the hole).

    Returns
    -------
    nodes : ndarray, shape (n_nodes, 2)
        Nodal coordinates.
    elements : ndarray, shape (n_elems, 3)
        Element connectivity.
    boundary_tags : dict
        Must contain:
        - 'hole': list of node indices on the hole boundary
        - 'right': list of node indices on the x=W boundary (applied tension)
        - 'sym_x': list of node indices on the y=0 boundary (symmetry: v=0)
        - 'sym_y': list of node indices on the x=0 boundary (symmetry: u=0)

    Notes
    -----
    Option A (manual): Create a structured mesh in (r, theta) space for
    r in [R, outer] and theta in [0, pi/2], then map to (x, y) using
    x = r cos(theta), y = r sin(theta). The outer boundary must conform to the
    rectangular plate edges.

    Option B (Gmsh): Use the gmsh Python API to define the geometry
    (rectangle minus circle), set mesh sizes, generate the mesh, and
    extract nodes, elements, and physical groups for boundary tagging.
    See https://gmsh.info/doc/textures/gmsh_api.html for the Python API docs.
    If you use Gmsh, add 'gmsh' to your requirements.txt.

    Option C (fallback): Use load_fallback_hole_mesh() below to load the
    pre-generated mesh from data/plate_with_hole_mesh.npz. This lets you
    proceed with the rest of the project while you work on your own mesher.
    """

    nr_nodes = n_radial + 1
    na_nodes = n_angular + 1
    n_nodes = nr_nodes * na_nodes

    theta_0 = 0.0
    theta_diag = np.arctan(H/W)
    theta_f = np.pi/2

    # Split the angular divisions relative to the size of H vs. W
    na_right = max(1, np.round(n_angular*H/(W+H)))
    na_top = n_angular - na_right
    dtheta_height = (theta_diag - theta_0)/na_right
    dtheta_width = (theta_f - theta_diag)/na_top

    theta_array = np.array([
        (dtheta_height*j_theta if j_theta <= na_right
         else theta_diag + dtheta_width*(j_theta - na_right))
        for j_theta in range(na_nodes)])

    nodes_polar = np.array([(
        (R + ((i_r/n_radial)**rho_mh)*(
            -R + (
                W/np.cos(theta_j) if theta_j < theta_diag
                else H/np.sin(theta_j)))),
        theta_j) for theta_j in theta_array for i_r in range(nr_nodes)])

    # Map from polar to Cartesian
    nodes = np.array([
        ((0.0 if theta_j == theta_f else r_i*np.cos(theta_j)), r_i*np.sin(theta_j))
        for r_i, theta_j in zip(nodes_polar[:,0], nodes_polar[:,1])])

    elements = []
    for row in range(na_nodes-1):
        for col in range(nr_nodes-1):
            i = nr_nodes*row + col
            j = i + 1
            k = i + nr_nodes
            l = k + 1
            elements += [[i, j, k], [l, k, j]]
    elements = np.array(elements)

    tol = 1e-9*max(H, W)
    boundary_tags = {
        "hole" : [n for n in range(n_nodes) if abs(nodes_polar[n,0]-R) < tol],
        "right" : [n for n in range(n_nodes) if abs(nodes[n,0]-W) < tol],
        "sym_x" : [n for n in range(n_nodes) if abs(nodes[n,1]) < tol],
        "sym_y" : [n for n in range(n_nodes) if abs(nodes[n,0]) < tol],
    }

    return (nodes, elements, boundary_tags)

def load_fallback_hole_mesh(filepath=None):
    """
    Load the pre-generated plate-with-hole mesh from the .npz file.

    Parameters
    ----------
    filepath : str, optional
        Path to the .npz file. Defaults to data/plate_with_hole_mesh.npz
        relative to the project root.

    Returns
    -------
    nodes : ndarray, shape (n_nodes, 2)
    elements : ndarray, shape (n_elems, 3)
    boundary_tags : dict with keys 'hole', 'right', 'sym_x', 'sym_y'

    Notes
    -----
    This is a fallback so you can work on the solver and validation plots
    without being blocked by mesh generation. Your final submission should
    include your own mesh generator (Option A or B above).
    """
    if filepath is None:
        filepath = os.path.join(os.path.dirname(__file__), '..', 'data', 'plate_with_hole_mesh.npz')
    data = np.load(filepath)
    nodes = data['nodes']
    elements = data['elements']
    boundary_tags = {
        'hole': data['hole'].tolist(),
        'right': data['right'].tolist(),
        'sym_x': data['sym_x'].tolist(),
        'sym_y': data['sym_y'].tolist(),
    }
    return nodes, elements, boundary_tags
