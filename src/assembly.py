# src/assembly.py
import numpy as np
from scipy.sparse import lil_matrix
from src.elements import compute_B, compute_D, compute_k


def assemble_K(nodes, elements, D, thickness):
    """
    Assemble the global stiffness matrix from element contributions.

    Parameters
    ----------
    nodes : ndarray, shape (n_nodes, 2)
    elements : ndarray, shape (n_elems, 3)
    D : ndarray, shape (3, 3)
    thickness : float

    Returns
    -------
    K : scipy.sparse.csr_matrix, shape (n_dof, n_dof)
        where n_dof = 2 * n_nodes.

    Notes
    -----
    Use scipy.sparse.lil_matrix for assembly (efficient for incremental insertion),
    then convert to CSR before returning (efficient for solving).
    The DOF ordering convention is: [u0, v0, u1, v1, ...] -- interleaved.
    For each element, extract the 6 global DOF indices from the 3 node indices,
    then scatter the 6x6 element stiffness into the global matrix.
    """
    n_nodes = nodes.shape[0]
    n_dof = 2*n_nodes
    K = lil_matrix((n_dof, n_dof))

    for element in elements:
        coords = nodes[element]
        i, j, k = element

        k_e = compute_k(coords, D, thickness)
        gl = {
            0 : 2*i,
            1 : 2*i+1,
            2 : 2*j,
            3 : 2*j+1,
            4 : 2*k,
            5 : 2*k+1}
        for el_r in range(k_e.shape[0]):
            for el_c in range(k_e.shape[1]):
                K[gl[el_r], gl[el_c]] += k_e[el_r, el_c]

    return K.tocsr()

def assemble_R_parabolic_shear(nodes, loaded_nodes, P, h):
    """
    Assemble the global load vector for a parabolic shear traction at the cantilever tip.

    The traction distribution along the tip edge (x = L) is:
        t_y(y) = (3P / 2h) * (1 - 4y^2/h^2)

    This must be integrated consistently using shape functions along each edge
    segment between adjacent loaded nodes (not applied as point loads).

    Parameters
    ----------
    nodes : ndarray, shape (n_nodes, 2)
    loaded_nodes : list of int
        Node indices along the loaded edge (x = L), to be sorted by y-coordinate.
    P : float
        Total applied tip shear force.
    h : float
        Plate height.

    Returns
    -------
    R : ndarray, shape (n_dof,)

    Notes
    -----
    For each edge segment between two adjacent loaded nodes, use at least 2-point
    Gauss quadrature to integrate t_y(y) * N_a(y) dy and t_y(y) * N_b(y) dy,
    where N_a and N_b are the linear (1D) shape functions along the edge.
    Verify: R.sum() should equal P (global force equilibrium).
    """

    loaded_nodes = np.array(loaded_nodes)
    loaded_nodes = [
        int(ind) for ind in list(
            loaded_nodes[nodes[loaded_nodes][:,1].argsort()[::-1]])]

    n_dof = 2*nodes.shape[0]

    R = np.zeros(shape=(n_dof,))

    t_y = lambda y : (3*P/(2*h))*(1-4*(y**2)/(h**2))
    xi_1, xi_2 = np.array([1.0, -1.0])/np.sqrt(3.0)
    R = np.zeros(shape=(n_dof,))

    for i, k in zip(loaded_nodes[:-1], loaded_nodes[1:]):
        y_i, y_k = nodes[[i, k],1]
        y = lambda xi : 0.5*((y_i + y_k) + (y_i - y_k)*xi)

        f = lambda xi : np.array([
            (1 + xi)*t_y(y(xi)),
            (1 - xi)*t_y(y(xi))])

        f_syi, f_syk = (0.25)*(y_i-y_k)*(f(xi_1) + f(xi_2))

        R[2*i + 1] += f_syi
        R[2*k + 1] += f_syk

    return R


def assemble_R_uniform_tension(nodes, loaded_nodes, sigma_inf, thickness):
    """
    Assemble the global load vector for uniform tension applied
    to a set of boundary nodes (used for the plate-with-hole problem).

    The traction is: t_x = sigma_inf applied along the loaded edge.

    Parameters
    ----------
    nodes : ndarray, shape (n_nodes, 2)
    loaded_nodes : list of int
    sigma_inf : float
    thickness : float

    Returns
    -------
    R : ndarray, shape (n_dof,)
    """
    loaded_nodes = np.array(loaded_nodes)
    loaded_nodes = [
        int(ind) for ind in list(
            loaded_nodes[nodes[loaded_nodes][:,1].argsort()[::-1]])]

    n_dof = 2*nodes.shape[0]

    R = np.zeros(shape=(n_dof,))

    for i, k in zip(loaded_nodes[:-1], loaded_nodes[1:]):
        y_i, y_k = nodes[[i, k],1]
        f_sxi, f_sxk = (thickness*(y_i-y_k)*sigma_inf/2)*np.array([1, 1])
        R[2*i + 1] += f_sxi
        R[2*k + 1] += f_sxk

    return R
