def timoshenko_model(EI, GA, L, rhoA, kappa=5/6):
    """Generate element stiffness and mass matrices for a Timoshenko beam element."""

    # Stiffness matrix
    k_e = [[(12*EI)/(L**3), (6*EI)/(L**2), -(12*EI)/(L**3), (6*EI)/(L**2)],
           [(6*EI)/(L**2), (GA*kappa)/L + (4*EI)/L, -(6*EI)/(L**2), (2*EI)/L - (GA*kappa)/L],
           [-(12*EI)/(L**3), -(6*EI)/(L**2), (12*EI)/(L**3), -(6*EI)/(L**2)],
           [(6*EI)/(L**2), (2*EI)/L - (GA*kappa)/L, -(6*EI)/(L**2), (GA*kappa)/L + (4*EI)/L]]

    # Mass matrix
    m_e = (rhoA * L / 420) * [[156, 22*L, 54, -13*L],
                               [22*L, 4*L**2, 13*L, -3*L**2],
                               [54, 13*L, 156, -22*L],
                               [-13*L, -3*L**2, -22*L, 4*L**2]]

    return k_e, m_e