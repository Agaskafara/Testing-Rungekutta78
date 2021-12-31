import argparse
import numpy as np
from utils.system_solver import solve_system


def args():
    parser = argparse.ArgumentParser()
    return parser.parse_args()

def main(args):

    # Create a list with random integers
    n_list = np.random.randint(3, 50, 10)
    print('\nMax element error in each of 10 random systems of random dimension:\n')

    for n in n_list:

        # Init a random matrix
        A = np.random.random((n, n))
        x = np.ones(n)

        # Find b such that x is solution to Ax = b
        b = A.dot(x)
        
        # solve system
        sol = solve_system(A, b)

        # Show results
        print("n = ", n)
        print("Max element error | x_i - 1|:")
        print(max(abs(sol - x)), "\n")

if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)