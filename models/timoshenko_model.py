import numpy as np

def timoshenko_model(G, E, I, rho, A, omega, L, kappa=5/6):
    """Generate element stiffness and mass matrices for a Timoshenko beam element."""

    # Dynamic stiffness matrix


    # Wave numbers for bending and shear deformation:
    a = (kappa* G * A * rho * I * omega**2  - E * I * rho * A * omega**2)  / (E* I * kappa * G * A)
    b = ((kappa * G * A - rho * I * omega**2) * rho * A * omega**2) / (E * I * kappa * G * A)
    lambda_1 = np.emath.sqrt(a + np.emath.sqrt(a**2 + 4*b) / 2)
    lambda_2 = np.emath.sqrt(a - np.emath.sqrt(a**2 + 4*b) / 2)

    # Factor needed for the rotation function

    alpha_1 = -(kappa * G * A * lambda_1) / (kappa * G * A - rho * I * omega**2 - E * I * lambda_1**2)
    alpha_2 = -(kappa * G * A * lambda_2) / (kappa * G * A - rho * I * omega**2 - E * I * lambda_2**2)



    numerator11 = -kappa * np.exp(2 * lambda_1 * L) * (-alpha_1 * lambda_2 + alpha_2 * lambda_1) * (
    np.exp(L * (lambda_1 - lambda_2)) + np.exp(L * (lambda_1 + lambda_2)) -
    np.exp(L * (3 * lambda_1 - lambda_2)) - np.exp(L * (3 * lambda_1 + lambda_2)))

    denominator11 = ((alpha_1 + alpha_2) * np.exp(L * (3 * lambda_1 - lambda_2)) +
                     (-alpha_1 + alpha_2) * np.exp(L * (5 * lambda_1 - lambda_2)) -
                     2 * alpha_2 * np.exp(2 * L * (2 * lambda_1 + lambda_2)) +
                     (-alpha_1 + 3 * alpha_2) * np.exp(L * (3 * lambda_1 + lambda_2)) +
                     (alpha_1 + 3 * alpha_2) * np.exp(L * (5 * lambda_1 + lambda_2)) -
                     2 * alpha_2 * (np.exp(2 * lambda_1 * L) + np.exp(4 * lambda_1 * L) + np.exp(6 * lambda_1 * L)))
    

    numerator12 = kappa * (
              alpha_1 * (lambda_1 + lambda_2 + alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2)) -
              alpha_1 * (lambda_1 - lambda_2 + alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2)) -
              2 * alpha_2 * (lambda_1 + alpha_1) * np.exp(2*L*(2*lambda_1 + lambda_2)) +
              (-alpha_1**2 + (-lambda_1 + lambda_2 + 3*alpha_2)*alpha_1 + 2*alpha_2*lambda_1) * np.exp(L*(3*lambda_1 + lambda_2)) +
              (alpha_1**2 + (lambda_1 + lambda_2 + 3*alpha_2)*alpha_1 + 2*alpha_2*lambda_1) * np.exp(L*(5*lambda_1 + lambda_2)) -
              2 * alpha_1 * (lambda_2 + alpha_2) * np.exp(2*lambda_1*L) -
              2 * alpha_2 * (lambda_1 + alpha_1) * np.exp(4*lambda_1*L) -
              2 * alpha_1 * (lambda_2 + alpha_2) * np.exp(6*lambda_1*L))

    denominator12 = (
              ( (alpha_1 + alpha_2) * np.exp(L*(3*lambda_1 - lambda_2)) +
              (-alpha_1 + alpha_2) * np.exp(L*(5*lambda_1 - lambda_2)) -
              2 * alpha_2 * np.exp(2*L*(2*lambda_1 + lambda_2)) +
              (-alpha_1 + 3*alpha_2) * np.exp(L*(3*lambda_1 + lambda_2)) +
              (alpha_1 + 3*alpha_2) * np.exp(L*(5*lambda_1 + lambda_2)) -
              2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
              ) * alpha_1)


    numerator13 = 2 * kappa * np.exp(3 * lambda_1 * L) * (np.exp(2 * lambda_1 * L) - 1) * (-alpha_1 * lambda_2 + alpha_2 * lambda_1)

    denominator13 = (
              (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2)) +
              (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2)) +
              2 * alpha_2 * np.exp(2*L*(2*lambda_1 + lambda_2)) +
              (alpha_1 - 3*alpha_2) * np.exp(L*(3*lambda_1 + lambda_2)) +
              (-alpha_1 - 3*alpha_2) * np.exp(L*(5*lambda_1 + lambda_2)) +
              2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L)))

    numerator14 = -kappa * np.exp(2 * lambda_1 * L) * (-alpha_1 * lambda_2 + alpha_2 * lambda_1) * (
              np.exp(L * (lambda_1 - lambda_2)) + np.exp(L * (lambda_1 + lambda_2)) -
              np.exp(L * (3 * lambda_1 - lambda_2)) - np.exp(L * (3 * lambda_1 + lambda_2)))

    denominator14 = (
              (alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2)) +
              (-alpha_1 + alpha_2) * np.exp(L * (5*lambda_1 - lambda_2)) -
              2 * alpha_2 * np.exp(2*L*(2*lambda_1 + lambda_2)) +
              (-alpha_1 + 3*alpha_2) * np.exp(L*(3*lambda_1 + lambda_2)) +
              (alpha_1 + 3*alpha_2) * np.exp(L*(5*lambda_1 + lambda_2)) -
              2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L)))

 

    # --- Numerator (numerator21) ---
    numerator21 = (
       -E * I * np.exp(2 * lambda_1 * L) * alpha_1 * alpha_2 *
       (
              (lambda_1 + lambda_2) * np.exp(L * (3*lambda_1 - lambda_2))
              + (lambda_1 + lambda_2) * np.exp(L * (3*lambda_1 + lambda_2))
              + (lambda_1 - lambda_2) * np.exp(L * (lambda_1 - lambda_2))
              - 2 * lambda_1 * np.exp(2 * L * (lambda_1 + lambda_2))
              + (lambda_1 - lambda_2) * np.exp(L * (lambda_1 + lambda_2))
              - 2 * lambda_1 * np.exp(2 * lambda_1 * L)))

# --- Denominator (denominator21) ---
    denominator21 = (
       (alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (-alpha_1 + alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       - 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (-alpha_1 + 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (alpha_1 + 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       - 2 * alpha_2 * (
              np.exp(2 * lambda_1 * L)
              + np.exp(4 * lambda_1 * L)
              + np.exp(6 * lambda_1 * L)))



    k11 = numerator11 / denominator11
    k12 = numerator12 / denominator12
    k13 = numerator13 / denominator13
    k14 = numerator14 / denominator14

    k21 = numerator21 / denominator21
# ========= k22 =========
    k22 = (
       (E * I) * (
              (alpha_1 * lambda_1 + alpha_2 * lambda_2) * np.exp(L * (3*lambda_1 - lambda_2))
              + (alpha_1 * lambda_1 + alpha_2 * lambda_2) * np.exp(L * (5*lambda_1 - lambda_2))
              + (((-2*lambda_1 + lambda_2) * alpha_2) - alpha_1 * lambda_1) * np.exp(L * (3*lambda_1 + lambda_2))
              + (((2*lambda_1 + lambda_2) * alpha_2) - alpha_1 * lambda_1) * np.exp(L * (5*lambda_1 + lambda_2))
              - 2 * lambda_2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(6*lambda_1*L))
       )
       ) / (
       (alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (-alpha_1 + alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       - 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (-alpha_1 + 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (alpha_1 + 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       - 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k23 =========
    k23 = (
       2 * E * I * alpha_1 * alpha_2 * np.exp(3*lambda_1*L) *
       (
              -lambda_1 * np.exp(2*lambda_1*L)
              - lambda_2 * np.exp(2*lambda_1*L)
              + 2 * lambda_1 * np.exp(L*(lambda_1 + lambda_2))
              - lambda_1
              + lambda_2
       )
       ) / (
       (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       + 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (alpha_1 - 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (-alpha_1 - 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       + 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k24 =========
    k24 = -(
       (E * I) * alpha_2 * alpha_1 * np.exp(2*lambda_1*L) * (
              (lambda_1 + lambda_2) * np.exp(L * (3*lambda_1 - lambda_2))
              + (lambda_1 + lambda_2) * np.exp(L * (3*lambda_1 + lambda_2))
              + (lambda_1 - lambda_2) * np.exp(L * (lambda_1 - lambda_2))
              - 2 * lambda_1 * np.exp(2 * L * (lambda_1 + lambda_2))
              + (lambda_1 - lambda_2) * np.exp(L * (lambda_1 + lambda_2))
              - 2 * lambda_1 * np.exp(2*lambda_1*L)
       )
       ) / (
       (alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (-alpha_1 + alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       - 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (-alpha_1 + 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (alpha_1 + 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       - 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k31 =========
    k31 = (
       kappa * G * A * (
              alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 - lambda_2))
              - alpha_2 * (lambda_1 + lambda_2) * np.exp(L * (3*lambda_1 + 2*lambda_2))
              + alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (5*lambda_1 + 2*lambda_2))
              - alpha_2 * (lambda_1 + lambda_2) * np.exp(L * (6*lambda_1 - lambda_2))
              + alpha_2 * (lambda_1 + lambda_2) * np.exp(L * (2*lambda_1 + lambda_2))
              - alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (6*lambda_1 + lambda_2))
              + (((-lambda_1 + lambda_2) * alpha_2) - 2 * alpha_1 * lambda_2) * np.exp(3*lambda_1*L)
              + (alpha_2 * (lambda_1 + lambda_2) + 2*alpha_1*lambda_2) * np.exp(5*lambda_1*L)
       )
       ) / (
       (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       + 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (alpha_1 - 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (-alpha_1 - 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       + 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k32 =========
    k32 = (
       kappa * G * A * (
                     -alpha_1 * (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 - lambda_2))
              + alpha_2 * (lambda_1 + lambda_2) * np.exp(L * (3*lambda_1 + 2*lambda_2))
              + alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (5*lambda_1 + 2*lambda_2))
              + alpha_1 * (lambda_1 + lambda_2) * np.exp(L * (6*lambda_1 - lambda_2))
              + alpha_1 * (lambda_1 + lambda_2) * np.exp(L * (2*lambda_1 + lambda_2))
              - alpha_1 * (lambda_1 - lambda_2) * np.exp(L * (6*lambda_1 + lambda_2))
              - 4 * lambda_1 * alpha_2 * np.exp(L * (4*lambda_1 + lambda_2))
              + (alpha_2 * (lambda_1 - lambda_2) - 2 * alpha_1 * lambda_2) * np.exp(3*lambda_1*L)
              + (alpha_2 * (lambda_1 + lambda_2) - 2 * alpha_1 * lambda_2) * np.exp(5*lambda_1*L)
       )
       ) / (
       (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       + 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (alpha_1 - 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (-alpha_1 - 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       + 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k33 =========
    k33 = -(
       kappa * G * A * (
              lambda_2 * (alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
              - lambda_2 * (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
              + 2 * lambda_2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
              + ((2*lambda_1 - lambda_2) * alpha_2 + alpha_1 * lambda_2)
              * np.exp(L * (3*lambda_1 + lambda_2))
              + ((-2*lambda_1 - lambda_2) * alpha_2 - alpha_1 * lambda_2)
              * np.exp(L * (5*lambda_1 + lambda_2))
              - 2 * alpha_2 * (
              lambda_1 * np.exp(2*lambda_1*L)
              - lambda_1 * np.exp(6*lambda_1*L)
              + lambda_2 * np.exp(4*lambda_1*L)
                     )
              )
              ) / (
              (alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
              + (-alpha_1 + alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
              - 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
              + (-alpha_1 + 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
              + (alpha_1 + 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
              - 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
              )


# ========= k34 (same structure as k31) =========
    k34 = (
              kappa * G * A * (
                     alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 - lambda_2))
                     - alpha_2 * (lambda_1 + lambda_2) * np.exp(L * (3*lambda_1 + 2*lambda_2))
                     + alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (5*lambda_1 + 2*lambda_2))
                     - alpha_2 * (lambda_1 + lambda_2) * np.exp(L * (6*lambda_1 - lambda_2))
                     + alpha_2 * (lambda_1 + lambda_2) * np.exp(L * (2*lambda_1 + lambda_2))
                     - alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (6*lambda_1 + lambda_2))
                     + (((-lambda_1 + lambda_2) * alpha_2) - 2 * alpha_1 * lambda_2) * np.exp(3*lambda_1*L)
                     + (alpha_2 * (lambda_1 + lambda_2) + 2*alpha_1*lambda_2) * np.exp(5*lambda_1*L)
              )
              ) / (
              (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
              + (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
              + 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
              + (alpha_1 - 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
              + (-alpha_1 - 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
              + 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
              )


# ========= k41 =========
    k41 = -(
       np.exp(L * (3*lambda_1 - lambda_2)) * E * I * alpha_2 * (
              (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 + 3*lambda_2))
              + (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 + lambda_2))
              - 2 * lambda_1 * np.exp(L * (lambda_1 + 2*lambda_2))
              + (lambda_1 + lambda_2) * np.exp(3*lambda_2*L)
              + (lambda_1 + lambda_2) * np.exp(lambda_2*L)
              - 2 * lambda_1 * np.exp(lambda_1*L)
       ) * alpha_1
       ) / (
       (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       + 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (alpha_1 - 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (-alpha_1 - 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       + 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k42 =========
    k42 = -(
       (
              -alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 + 3*lambda_2))
              + 2 * alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (3*lambda_1 + 2*lambda_2))
              - 2 * alpha_2 * (lambda_1 + lambda_2) * np.exp(-L * (lambda_1 - 2*lambda_2))
              - alpha_2 * (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 + lambda_2))
              - 2 * lambda_1 * alpha_1 * np.exp(L * (lambda_1 + 2*lambda_2))
              + alpha_2 * (lambda_1 + lambda_2) * np.exp(3*lambda_2 * L)
              + alpha_2 * (lambda_1 + lambda_2) * np.exp(lambda_2 * L)
              + 2 * lambda_1 * alpha_1 * np.exp(lambda_1 * L)
       ) * np.exp(L * (3*lambda_1 - lambda_2)) * E * I
       ) / (
       (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       + 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (alpha_1 - 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (-alpha_1 - 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       + 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k43 =========
    k43 = (
       2 * E * I * alpha_1 * alpha_2 * np.exp(L * (3*lambda_1 + lambda_2)) *
       (
              -lambda_1 * np.exp(2*lambda_1*L)
              + lambda_2 * np.exp(2*lambda_1*L)
              + 2 * lambda_1 * np.exp(L * (lambda_1 - lambda_2))
              - lambda_1
              - lambda_2
       )
       ) / (
       (alpha_1 + alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (-alpha_1 + alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       - 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (-alpha_1 + 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (alpha_1 + 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       - 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )


# ========= k44 =========
    k44 = -(
       np.exp(L * (3*lambda_1 - lambda_2)) * E * I * alpha_2 * (
              (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 + 3*lambda_2))
              + (lambda_1 - lambda_2) * np.exp(L * (2*lambda_1 + lambda_2))
              - 2 * lambda_1 * np.exp(L * (lambda_1 + 2*lambda_2))
              + (lambda_1 + lambda_2) * np.exp(3*lambda_2 * L)
              + (lambda_1 + lambda_2) * np.exp(lambda_2 * L)
              - 2 * lambda_1 * np.exp(lambda_1 * L)
       ) * alpha_1
       ) / (
       (-alpha_1 - alpha_2) * np.exp(L * (3*lambda_1 - lambda_2))
       + (alpha_1 - alpha_2) * np.exp(L * (5*lambda_1 - lambda_2))
       + 2 * alpha_2 * np.exp(2 * L * (2*lambda_1 + lambda_2))
       + (alpha_1 - 3*alpha_2) * np.exp(L * (3*lambda_1 + lambda_2))
       + (-alpha_1 - 3*alpha_2) * np.exp(L * (5*lambda_1 + lambda_2))
       + 2 * alpha_2 * (np.exp(2*lambda_1*L) + np.exp(4*lambda_1*L) + np.exp(6*lambda_1*L))
       )
       

    k_e = np.array([[k11, k12, k13, k14], 
                             [k21, k22, k23, k24], 
                             [k31, k32, k33, k34], 
                             [k41, k42, k43, k44]])
    return k_e