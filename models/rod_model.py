import numpy as np

def rod_model(EA, L, rho, A, omega):

    """Generate element dynamic stiffness matrix for a 1D rod element."""

    beta = omega * np.sqrt(rho * A / EA)
    
    # Damping coefficient (proportional damping example)
    alpha = 0.0  # mass proportional damping
    beta_damp = 0 # stiffness proportional damping (adjust as needed)

    # Element dynamic stiffness matrix (undamped)
    k_e = EA * beta * np.array([[1 / np.tan(beta*L), -1 / np.sin(beta*L)],
                                [-1 / np.sin(beta*L), 1 / np.tan(beta*L)]])

    # Element mass matrix (consistent mass matrix for a rod)
    m_e = (rho* A * L / 6) * np.array([[2, 1],
                                     [1, 2]])

    # Add damping (Rayleigh damping: C = alpha*M + beta*K)
    c_e = alpha * m_e + beta_damp * k_e

    # Dynamic stiffness with damping: K_dyn = K + i*omega*C
    k_e = k_e + 1j * omega * c_e
    
    return k_e

    