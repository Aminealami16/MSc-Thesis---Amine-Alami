import numpy as np

def effective_truss_stiffness(d1, d2,t1, t2, L, h): 
    '''Calculates equivalent moment of inertia and area to be used in a Timoshenko beam model

    d1: Diameter of the diagonal truss elements
    d2: Diameter of the horizontal truss elements
    t1: Thickness of the diagonal truss elements
    t2: Thickness of the horizontal truss elements
    L: Length of each truss section
    h: Height of the truss
    z_NC: Distance from the neutral axis to the centroid of the truss elements

    returns:
    A_eq: Effective area of the truss
    I_eqy: Effective moment of inertia around the y-axis (strong axis)
    I_eqz: Effective moment of inertia around the z-axis (weak axis)
    b_eq: Effective width of the truss for shear deformation calculations
    '''
    # Calculate the neutral axis based on the top and bottom rules of the truss

    z_NCy = (0.25 * np.pi * (d2 - 2 * t2)**2  * h) / (0.25 * np.pi * (d2 - 2 * t2)**2  * 3)
    y_NC = (0.25 * np.pi * (d2 - 2 * t2)**2  * L + 0.25 * np.pi * (d2 - 2 * t2)**2  * 0.5*L) / (0.25 * np.pi * (d2 - 2 * t2)**2  * 3)

    # Calculate the area of the diagonal and horizontal truss elements
    A1 = np.pi * (d1/2)**2 - np.pi * ((d1/2) - t1)**2
    A2 = np.pi * (d2/2)**2 - np.pi * ((d2/2) - t2)**2

    # Calculate the moment of inertia for the diagonal and horizontal truss elements around the y-axis (strong axis)
    I1y = (np.pi / 4) * ((d1/2)**4 - ((d1/2) - t1)**4)
    I2y = (np.pi / 4) * ((d2/2)**4 - ((d2/2) - t2)**4)
    I2z = I2y

    # Calculate the effective area and moment of inertia for the Timoshenko beam model
    A_eq = A1 + A2 + 2 * A2 # Adding the contribution of the horizontal truss elements and he diagonal truss elements
    I_eqy  = 2 * (I2y + A2 * z_NCy**2) + I2y + A2 * (h - z_NCy)**2 # Moment of inertia is only based on the top and bottom rules
    I_eqz = I2z + A2 * y_NC**2 + I2z + A2 * (L - y_NC)**2 + I2z + A2 * (0.5*L - y_NC)**2 # Moment of inertia is based on all horizontal truss elements
    b_eq = A_eq / h 
    return A_eq, I_eqy, I_eqz, b_eq


    
def effective_retaining_wall_stiffness():

    '''Calculates the effective stiffness of a retaining wall to be used in a Timoshenko beam model

    This function is a placeholder and should be implemented based on the specific geometry and material properties of the retaining wall.
    '''
    pass


