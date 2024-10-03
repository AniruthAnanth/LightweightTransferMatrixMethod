import numpy as np
import matplotlib.pyplot as plt

def transfer_matrix(k_0, n, d, theta):
    k_z = k_0 * n * np.cos(theta) # Calculate longitudanal K

    # Reduce redundant calculations in construction T_i
    q_1 = np.cos(k_z * d) 
    q_2 = 1j * np.sin(k_z * d)
    n_cos_th = n * np.cos(theta)

    # Transfer matrix calculation derived from Maxwell's equations
    return np.array([
        [q_1, q_2 / n_cos_th],
        [n_cos_th * q_2, q_1]
    ])

def solve_tmm(layers, wavelength, theta):
    # Scales wave propagation
    k_0 = (2 * np.pi) / wavelength

    # Global transfer matrix initialization
    M = np.eye(2)

    # 
    n_0 = layers[0][0]
    n_l = layers[-1][0]
    
    for layer in layers[1:-1]:
        # Get layer parameters
        n = layer[0](wavelength) # Refractive index as a function
        d = layer[1]

        # Calculate transfer matrix
        theta_i = np.arcsin(n_0 * np.sin(theta) / n)
        T_i = transfer_matrix(k_0, n, d, theta_i)

        # Update global transfer matrix with dot product
        M = np.dot(M, T_i)

    # Calculate M_in and M_out
    q_1 = n_0 * np.cos(theta)
    theta_l = np.arcsin(n_0 * np.sin(theta) / n_l)
    q_2 = n_l * np.cos(theta_l)
                   
    M_in = np.array([
            [1, 1],
            [q_1, -q_1],
        ])
    
    M_out = np.array([
            [1, 1],
            [q_2, -q_2],
        ])

    # Calculate M_global final
    M_total = np.dot(np.linalg.pinv(M_in), np.dot(M, M_out))

    # Extract reflection and transmission
    r = M_total[1, 0] / M_total[0, 0]
    t = 1 / M_total[0, 0]
    R = np.abs(r) ** 2
    T = np.abs(t) ** 2 * (n_l * np.cos(theta_l)) / (n_0 * np.cos(theta))

    return R, T
