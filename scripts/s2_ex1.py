import argparse
import numpy as np
from functools import partial
from utils.functs import s2_funct, s2_dfunct
from utils.include_diffs import include_diffs
from utils.rk78 import rk78

def args():
    parser = argparse.ArgumentParser()

    return parser.parse_args()

def main(args):

    # Create a rk78 instance
    numeric_algorithm = rk78()

    # Hyperparameters ------------------------
    h_min, h_max = 0.0001, 0.005
    tol = 0.0001
    max_steps = 10000
    h = 0.001
    # ----------------------------------------
    # Initial conditions ---------------------
    t0 = 0
    init_pos = np.array([1, 0] + [1, 0, 0, 1])
    # ----------------------------------------
    # Time predictions -----------------------
    t_pred = 0.5
    # ----------------------------------------
    # Parameters -----------------------------
    alpha = 0.5
    # ----------------------------------------

    # Build system function
    inc_diffs = include_diffs(partial(s2_funct, alpha_=alpha),
                              partial(s2_dfunct, alpha_=alpha))
    funct = inc_diffs.camp

    # Predict
    output = numeric_algorithm.flux(t_pred, t0, init_pos,
                                    h, h_min, h_max, tol, max_steps, funct)
    
    # Predict delta
    delta = 0.00001
    delta_x0 = np.array([delta, 0])
    phi_pdelta_x = numeric_algorithm.flux(t_pred, t0, init_pos[:2] + delta_x0,
                                         h, h_min, h_max, tol, max_steps,
                                         partial(s2_funct, alpha_=alpha))['pos']
    delta_y0 = np.array([0, delta])
    phi_pdelta_y = numeric_algorithm.flux(t_pred, t0, init_pos[:2] + delta_y0,
                                         h, h_min, h_max, tol, max_steps,
                                         partial(s2_funct, alpha_=alpha))['pos']
    
    # Predict -delta
    delta = -0.00001
    delta_x0 = np.array([delta, 0])
    phi_mdelta_x = numeric_algorithm.flux(t_pred, t0, init_pos[:2] + delta_x0,
                                         h, h_min, h_max, tol, max_steps,
                                         partial(s2_funct, alpha_=alpha))['pos']
    delta_y0 = np.array([0, delta])
    phi_mdelta_y = numeric_algorithm.flux(t_pred, t0, init_pos[:2] + delta_y0,
                                         h, h_min, h_max, tol, max_steps,
                                         partial(s2_funct, alpha_=alpha))['pos']

    # Build D_xPhi(t=0.5)
    dphi = np.concatenate([((phi_pdelta_x - phi_mdelta_x)/(2*delta))[None],
                           ((phi_pdelta_y - phi_mdelta_y)/(2*delta))[None]], axis=0).T
    
    # Compare with the one computed first
    print("Error between the computed matrices:")
    print(dphi - output['pos'][2:].reshape(2, 2))

if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)