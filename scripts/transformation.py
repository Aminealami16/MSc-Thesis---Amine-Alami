import numpy as np

def transformation_matrix(n1, n2):
    """Create an element connecting two nodes. """
    nodes = [n1, n2]

    # Extract coordinates of the nodes
    x1, y1, z1 = n1
    x2, y2, z2 = n2

    # Reference vector for local axis orientation
    v = np.array([1, 0, 0], dtype=float)  # Reference vector along the x-axis

    # Element vector
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    L = np.sqrt(dx**2 + dy**2 + dz**2)  # Length of the element

    # Local axis orientation
    lx = dx / L
    ly = dy / L
    lz = dz / L
    x_local = np.array([lx, ly, lz], dtype=float)

    # Local z-axis = x × v
    z_local = np.cross(x_local, v)
    z_norm = np.linalg.norm(z_local)
    if z_norm < 1e-12:
        raise ValueError("Reference vector is parallel to element axis. Choose another reference vector.")
    z_local /= z_norm

    # Local y-axis = z × x
    y_local = np.cross(z_local, x_local)

    # Rotation matrix (3×3)
    Lambda = np.vstack((x_local, y_local, z_local))

    # Build full 12×12 transformation matrix
    T = np.zeros((12, 12))
    T[0:3, 0:3] = Lambda
    T[3:6, 3:6] = Lambda
    T[6:9, 6:9] = Lambda
    T[9:12, 9:12] = Lambda


    return T


