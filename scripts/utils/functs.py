import numpy as np


# FUNCTIONS
def pendol_funct(t: float, pos : np.ndarray):
    return np.array([pos[1], -np.sin(pos[0])])


def s1_ex4_functA(t: float, pos : np.ndarray):
    return np.array([2*pos[0]/t])

def s1_ex4_functB(t: float, pos : np.ndarray):
    return np.array([pos[1], - pos[0]])

def lorentz_funct(t : float, pos : np.ndarray,
                  sigma_ : float = 3, p_ : float = 26.5, b_ : float = 1.):
    x, y, z = pos
    return np.array([sigma_*(y - x), - x*z + p_*x - y, x*y - b_*z])

def rtbp_funct(t : float, pos : np.ndarray, mu_ : float):

    # pos z, y, z, vel u, v, w
    x, y, z, u, v, w = pos

    # Computed elements for a
    a1 = 2*v + x - (mu_ - x)*(mu_ - 1)/((mu_ - x)**2 + y**2 + z**2)**(3/2) + \
        (mu_ - x - 1)*mu_/((mu_ - x - 1)**2 + y**2 + z**2)**(3/2)
    a2 = -2*u + y + (mu_ - 1)*y/((mu_ - x)**2 + y**2 + z**2)**(3/2) - \
        mu_*y/((mu_ - x - 1)**2 + y**2 + z**2)**(3/2)
    a3 = (mu_ - 1)*z/((mu_ - x)**2 + y**2 + z**2)**(3/2) - mu_*z/((mu_ - x - 1)**2 + \
        y**2 + z**2)**(3/2)

    return np.array([u, v, w] + [a1, a2, a3])

def s2_funct(t : float, pos : np.ndarray, alpha_ : float):
    r2 = sum(pos**2)
    return np.array([alpha_*(1 - r2)*pos[0] - pos[1], pos[0] + alpha_*(1 - r2)*pos[1]])

def s2_dfunct(t : float, pos : np.ndarray, alpha_ : float):
    x, y = pos
    return np.array([[alpha_*(1 - y**2 - 3*x**2), - (1 + 2*alpha_*x*y)],
                     [1 - 2*alpha_*x*y, alpha_*(1 - x**2 - 3*y**2)]])