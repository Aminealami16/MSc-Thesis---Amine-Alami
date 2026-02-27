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
    A1 = np.pi * (d1/2)**2 - np.pi * ((d1/2) - 2 * t1)**2
    A2 = np.pi * (d2/2)**2 - np.pi * ((d2/2) - 2 * t2)**2

    # Calculate the moment of inertia for the diagonal and horizontal truss elements around the y-axis (strong axis)
    I1y = (np.pi / 4) * ((d1/2)**4 - ((d1/2) - t1)**4)
    I2y = (np.pi / 4) * ((d2/2)**4 - ((d2/2) - t2)**4)
    I2z = I2y

    # Calculate the effective area and moment of inertia for the Timoshenko beam model
    A_eq = A1 * 4 + A2 + 2 * A2 # Adding the contribution of the horizontal truss elements and he diagonal truss elements
    I_eqy  = 2 * (I2y + A2 * z_NCy**2) + I2y + A2 * (h - z_NCy)**2 # Moment of inertia is only based on the top and bottom rules
    I_eqz = I2z + A2 * y_NC**2 + I2z + A2 * (L - y_NC)**2 + I2z + A2 * (0.5*L - y_NC)**2 # Moment of inertia is based on all horizontal truss elements
    b_eq = A_eq / h 
    return A_eq, I_eqy, I_eqz, b_eq


    
def effective_retaining_wall_stiffness(t):

    '''Calculates the effective stiffness of a retaining wall to be used in a Timoshenko beam model
    '''
    h = 22 #m
    h1 = 14.5 #m
    h2 = 7.5 #m
    b1 = 8 #m
    b2 = 12 #m



    z_NC = (t * b1 * 0 + h * t * 0.5 * h + h1 * t * 0.5 * h1 + (b2 - b1) * t * h1 + h2 * t * (h2 * 0.5 + h1) + b2 * t * h) / (2 * t * h + b1 * t + b2 * t + (b2 - b1) * t)
    y_NC = (h * t * 0 + t * b1 * 0.5 * b1 + t * b2 * 0.5 * b2 + h1 * t * b1 + h2 * t * b2 + (b2 - b1) * t * (b1 + 0.5 * (b2 - b1))) / (2 * t * h + b1 * t + b2 * t + (b2- b1) * t)
    I_eqy = (1/12) * b1 * t**3 + (1/12) * t * h**3 + (1/12) * t * h1**3 + (1/12) * t * h2**3 + (1/12) * (b2 - b1) * t**3 + (1/12) * b2 * t**3 + b1 * t * z_NC**2 + t * h * (z_NC - 0.5 * h)**2 + h1 * t * (0.5*h1 - z_NC)**2 + h2 * t * (h1 + 0.5 * h2 - z_NC)**2 + b2 * t * (h - z_NC)**2 + (b2 - b1) * t * (z_NC - h1)**2
    I_eqz = (1/12) * h * t ** 3 + (1/12) * t * b2**3 + (1/12) * t * b1**3 + (1/12) * h1 * t**3 + (1/12) * t * (b2 - b1) **3  + (1/12) * h2 * t**3 + h * t * y_NC**2 + b1 * t * (y_NC - 0.5 * b1)**2 + b2 * t * (y_NC - 0.5 * b2)**2 + h1 * t * (b1 - y_NC)**2 + h2 * t * (b2 - y_NC)**2 + (b2 - b1) * t * (b1+ 0.5 * (b2 - b1) - y_NC)**2

    A_eq = 2 * t * h + b1 * t + b2 * t + (b2 - b1) * t
    b_eq = A_eq / h

    return A_eq, I_eqy, I_eqz, b_eq




def stiffness_fenders():

    '''Calculates the stiffness of the fenders to be used in a Timoshenko beam model

    '''

    # K_spring_linear = 4 * (250 * 1000) / (90 / 1000) # retrieved from report Handboek 2: Kerende wand en vakwerkarmen
    K_spring_linear = ((1200+1750) / 2) * 10**6  # Retrieved from HVR engineering report
    return K_spring_linear
   


def stiffness_connecting_beams():

    '''Calculates the stiffness of the connecting beams to be used in a Timoshenko beam model
    '''
    h = 18
    d_out = 1.8 
    d_in = 1.8  - (2 * 80) / 1000
    A = np.pi * (d_out/2)**2 - np.pi * (d_in/2)**2
    Iy = 2 * ( (np.pi / 4) * ((d_out/2)**4 - (d_in/2)**4) + A * (h/2)**2) 
    Iz = 2 * ( (np.pi / 4) * ((d_out/2)**4 - (d_in/2)**4) )
    Ip = Iy + Iz    
    It = Ip
    

    return Iy, Iz, Ip, It, A

