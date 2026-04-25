# src/postprocess.py
import numpy as np
from src.elements import compute_B


def compute_stresses(nodes, elements, u, D):
    """
    Recover element stresses from the displacement solution.

    Parameters
    ----------
    nodes : ndarray, shape (n_nodes, 2)
    elements : ndarray, shape (n_elems, 3)
    u : ndarray, shape (n_dof,)
    D : ndarray, shape (3, 3)

    Returns
    -------
    stresses : ndarray, shape (n_elems, 3)
        Stress components [sigma_xx, sigma_yy, tau_xy] at each element centroid.
    """

    stresses = []
    for element in elements:
        i, j, k = element
        B = compute_B(nodes[element])

        u_e = u[[2*i, 2*i+1, 2*j, 2*j+1, 2*k, 2*k+1]]
        stresses += [np.matmul(np.matmul(D, B), u_e)]
    stresses = np.array(stresses)

    return stresses


def compute_von_mises(stresses):
    """Von Mises stress from [sigma_xx, sigma_yy, tau_xy] per element."""

    sigma_VM = lambda sigma_xx, sigma_yy, tau_xy : np.sqrt(
        sigma_xx**2 + sigma_yy**2 - sigma_xx*sigma_yy + 3*(tau_xy**2))

    return sigma_VM(stresses[:,0], stresses[:,1], stresses[:,2])


def strain_energy(K, u):
    """Return 0.5 * u^T K u."""
    return 0.5*(u.T @ K @ u)
