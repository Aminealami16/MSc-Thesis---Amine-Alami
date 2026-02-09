import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def nodes(x, y, z):
    """
    Create a 3D grid of nodes.

    Parameters:
    x (int): Number of nodes in the x-direction.
    y (int): Number of nodes in the y-direction.
    z (int): Number of nodes in the z-direction.

    Returns:
    np.ndarray: A 3D array of shape (x, y, z) containing the coordinates of the nodes.
    """
    return np.array([x, y, z])

def plot_nodes(nodes):
    """
    Plot the nodes in a 3D space.

    Parameters:
    nodes (np.ndarray): A 3D array containing the coordinates of the nodes.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(nodes)):
        x, y, z = nodes[i]
        ax.scatter(x, y, z)
        ax.text(x, y, z, f'  {i+1}', fontsize=8)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
    plt.show()



def plot_elements(elements):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for element in elements:
        n1, n2 = element
        x = [n1[0], n2[0]]
        y = [n1[1], n2[1]]
        z = [n1[2], n2[2]]
        ax.plot(x, y, z, 'b-')
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show()
    
