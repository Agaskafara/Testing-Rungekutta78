import argparse
import numpy as np
from utils.functs import s1_ex4_functA, s1_ex4_functB
from utils.rk78 import rk78

def args():
    parser = argparse.ArgumentParser()

    return parser.parse_args()

def main(args):

    # Create a rk78 instance
    numeric_algorithm = rk78()

    # Hyperparameters ------------------------
    h_min, h_max = 0.01, 0.05
    tol = 0.0001
    max_steps = 1000
    h = 0.01
    # ----------------------------------------
    # Initial conditions ---------------------
    t0 = {"A": 1, "B": 0}
    init_pos = {"A": np.array([1]),
                "B": np.array([1, 0])}
    # ----------------------------------------
    # time to predict ------------------------
    t_pred = {"A": 3.2506, "B": 3.2506}
    # ----------------------------------------

    # Predict A
    outputA = numeric_algorithm.flux(t_pred["A"], t0["A"], init_pos["A"],
                                     h, h_min, h_max, tol, max_steps, s1_ex4_functA)

    # Show prediction and check if tol is satisfied
    print(outputA)
    if outputA['success']:
        print('Does the estimation satisify the tolerance?:',
              abs(outputA['pos'] - np.array([t_pred["A"]])**2) < tol)


    # Predict B
    outputB = numeric_algorithm.flux(t_pred["B"], t0["B"], init_pos["B"],
                                    h, h_min, h_max, tol, max_steps, s1_ex4_functB)

    # Show prediction and check if tol is satisfied
    print(outputB)
    if outputB['success']:
        print('Does the estimation satisify the tolerance?:',
              abs(outputB['pos'] - np.array([np.cos(t_pred["B"]), - np.sin(t_pred["B"])])) < tol)

if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)