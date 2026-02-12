from matplotlib.pylab import beta
import numpy as np

def euler_bernoulli_model2D(E, I, rho, A, L, omega):

    beta = (rho * A * omega**2 / (E * I))**(1/4)

    # First column of dynamic stiffness matrix

    k11 = (
    ((np.cos(beta*L) + np.sin(beta*L)) * np.exp(2*beta*L)
     - np.cos(beta*L) + np.sin(beta*L))
    * beta**3 * E * I
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    ) 

    k21 = (
    -np.sin(beta*L) * (np.exp(2*beta*L) - 1) * E * I * beta**2
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k31 = (
    E * I * beta**3
    * (2*np.exp(beta*L)*np.sin(beta*L) + np.exp(2*beta*L) - 1)
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k41 = (
        E * I * beta**2
        * (np.exp(2*beta*L) - 2*np.exp(beta*L)*np.cos(beta*L) + 1)
        / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    # Second column of dynamic stiffness matrix

    k12 = (
    -np.sin(beta*L) * (np.exp(2*beta*L) - 1) * E * I * beta**2
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k_22 = (
    -E * I * beta
    * ((np.cos(beta*L) - np.sin(beta*L)) * np.exp(2*beta*L) - np.sin(beta*L) - np.cos(beta*L))
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k32 = (
    E * I * beta**2
    * (2*np.exp(beta*L)*np.cos(beta*L) - np.exp(2*beta*L) - 1)
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k42 = (
    E * I * beta
    * (-np.exp(2*beta*L) + 2*np.exp(beta*L)*np.sin(beta*L) + 1)
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    # Third column of dynamic stiffness matrix


    k13 = (
        -(np.exp(2*beta*L) + 2*np.exp(beta*L)*np.sin(beta*L) - 1) * E * I * beta**3
        / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k23 = (
    E * I * beta**2
    * (np.exp(2*beta*L) - 2*np.exp(beta*L)*np.cos(beta*L) + 1)
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k33 = (
    -E * I * ((np.cos(beta*L) + np.sin(beta*L)) * np.exp(2*beta*L)
              - np.cos(beta*L) + np.sin(beta*L)) * beta**3
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k43 = (
    -np.sin(beta*L) * (np.exp(2*beta*L) - 1) * E * I * beta**2
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    # Fourth column of dynamic stiffness matrix

    k14 = (
    E * I * ((np.cos(beta*L) + np.sin(beta*L)) * np.exp(2*beta*L)
             - np.cos(beta*L) + np.sin(beta*L)) * beta**3
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k24 = (
    -np.sin(beta*L) * (np.exp(2*beta*L) - 1) * E * I * beta**2
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )

    k34 = (
    (np.exp(2*beta*L) + 2*np.exp(beta*L)*np.sin(beta*L) - 1) * E * I * beta**3
    / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )   

    k44 = (
        E * I * beta**2
        * (np.exp(2*beta*L) - 2*np.exp(beta*L)*np.cos(beta*L) + 1)
        / (np.exp(2*beta*L)*np.cos(beta*L) - 2*np.exp(beta*L) + np.cos(beta*L))
    )


    k_dyn = np.array([[k11, k12, k13, k14],
                                         [k21, k_22, k23, k24], 
                                            [k31, k32, k33, k34],
                                            [k41, k42, k43, k44]])
    
    return k_dyn


def euler_bernoulli_model3D(E, I, rho, A, L, omega):
















    

    
