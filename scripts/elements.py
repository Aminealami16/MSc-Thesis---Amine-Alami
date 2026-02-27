import numpy as np
import sys
import os
sys.path.append(os.path.abspath('..'))
from models import timoshenko_model as tm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D





def plot_elements(elements):
    
    fig = plt.figure(figsize=(25, 25))
    ax = fig.add_subplot(111, projection='3d')

    for i, element in enumerate(elements):
        n1, n2, _, _ = element
        x = [n1[0], n2[0]]
        y = [n1[1], n2[1]]
        z = [n1[2], n2[2]]
        ax.plot(x, y, z, 'b-')
        
        # Plot nodes
        ax.scatter(*n1, color='red', s=50)
        ax.scatter(*n2, color='red', s=50)
        
        # Add element number at midpoint
        mid_x = (n1[0] + n2[0]) / 2
        mid_y = (n1[1] + n2[1]) / 2
        mid_z = (n1[2] + n2[2]) / 2
        ax.text(mid_x, mid_y, mid_z, str(i), fontsize=15, color='black')
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_zlim(0, 100)
    plt.show()
    
def plot_elements2d(elements):
    """
    Plot elements in 2D using x and y coordinates.
    """
    fig, ax = plt.subplots(figsize=(20, 20))
    
    for i, element in enumerate(elements):
        n1, n2, _, _ = element
        x = [n1[0], n2[0]]
        y = [n1[1], n2[1]]
        ax.plot(x, y, 'b-')
        
        # Plot nodes
        ax.scatter(*n1[:2], color='red', s=50)
        ax.scatter(*n2[:2], color='red', s=50)
        
        # Add element number at midpoint
        mid_x = (n1[0] + n2[0]) / 2
        mid_y = (n1[1] + n2[1]) / 2
        ax.text(mid_x, mid_y, str(i), fontsize=15, color='black')
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_aspect('equal')
    plt.show()


def elements(n1, n2, ep_K, ep_M):
    """
    Create an element connecting two nodes.

    Parameters:
    n1 (np.ndarray): Coordinates of the first node.
    n2 (np.ndarray): Coordinates of the second node.

    Returns:
    tuple: A tuple containing the coordinates of the two nodes.
    """
    ex = [n1[0], n2[0]]
    ey = [n1[1], n2[1]]
    ez = [n1[2], n2[2]]
    eo = [0, -1, 0]

    M_e, K_e = tm.T_element(ex, ey, ez, eo, ep_K, ep_M)

    return n1, n2, M_e, K_e


def elements_added_mass(n1, n2, ep_K, ep_M):
    """
    Create an element connecting two nodes.

    Parameters:
    n1 (np.ndarray): Coordinates of the first node.
    n2 (np.ndarray): Coordinates of the second node.

    Returns:
    tuple: A tuple containing the coordinates of the two nodes.
    """
    ex = [n1[0], n2[0]]
    ey = [n1[1], n2[1]]
    ez = [n1[2], n2[2]]
    eo = [0, -1, 0]

    M_e, K_e = tm.T_element_added_mass_retaining_wall(ex, ey, ez, eo, ep_K, ep_M)

    return n1, n2, M_e, K_e


def hermite_shape_functions(xi, L):
    """
    Returns the 4 Hermite cubic shape functions evaluated at xi.
    
    Parameters:
    xi : float or array
        Normalized coordinate along element (0 <= xi <= 1)
    L : float
        Element length (for scaling slope terms)
    
    Returns:
    H : array of shape (4,)
        Hermite cubic functions [H1, H2, H3, H4]
    """
    H1 = 1 - 3*xi**2 + 2*xi**3
    H2 = xi - 2*xi**2 + xi**3
    H3 = 3*xi**2 - 2*xi**3
    H4 = -xi**2 + xi**3
    
    return np.array([H1, H2*L, H3, H4*L])


def expand_eigenvectors(eigvecs_reduced, keep_dofs, total_dofs):
    """
    Expand reduced eigenvectors to full DOF size by inserting zeros
    at constrained DOFs.
    
    Parameters:
    -----------
    eigvecs_reduced : array, shape (num_free_dofs, num_modes)
        Eigenvectors from reduced system
    keep_dofs : list or array
        Indices of free DOFs kept in reduced system
    total_dofs : int
        Total number of DOFs in the full system
    
    Returns:
    --------
    eigvecs_full : array, shape (total_dofs, num_modes)
        Eigenvectors in full DOF space, zeros at constrained DOFs
    """
    num_modes = eigvecs_reduced.shape[1]
    eigvecs_full = np.zeros((total_dofs, num_modes))
    
    for j, dof in enumerate(keep_dofs):
        eigvecs_full[dof, :] = eigvecs_reduced[j, :]
    
    return eigvecs_full



def extract_displacement(arr, keep=3, skip=3):
    """
    arr: 2D numpy array where each column is an eigenvector
    keep: number of rows to keep in each cycle
    skip: number of rows to skip in each cycle
    """
    n_rows = arr.shape[0]
    step = keep + skip

    # Generate the indices once
    idx = np.hstack([
        np.arange(i, min(i + keep, n_rows))
        for i in range(0, n_rows, step)
    ])

    # Apply same row indices to ALL columns → returns a compact matrix
    return arr[idx, :]


def remove_close_frequencies(frequencies, threshold=1e-3):
    unique_freqs = []
    for freq in frequencies:
        if all(abs(freq - uf) > threshold for uf in unique_freqs):
            unique_freqs.append(freq)
    return np.array(unique_freqs)