import argparse
import numpy as np
from utils.system_solver import solve_system


def args():
    parser = argparse.ArgumentParser()
    return parser.parse_args()

def main(args):

    # Init system matrices
    A = np.array([[0, -4], [0, 0], [5, -2]])
    x = np.array([0.3, -0.25])
    b = np.array([1, 3, 2])

    # Solve system
    sol = solve_system(A, b)

    # Show results
    print("Check solution:")
    print('Real solution :', x)
    print('Found solution:', sol)

if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)