import os
import argparse
import numpy as np
from utils.rk78 import rk78
from utils.functs import pendol_funct

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--init_conditions', required=True,
                        help='txt file with initial conditions.')
    parser.add_argument('--output_folder', required=True,
                        help='Folder where file is stored.')
    return parser.parse_args()

def main(args):

    # Create a rk78 algorithm instance
    numeric_algorithm = rk78()

    # Load initial conditions
    initial_conditions = np.loadtxt(args.init_conditions)

    # Hyperparameters ------------------------    
    h_min, h_max = 0.01, 0.05
    tol = 0.0001
    # ----------------------------------------
    # Initial conditions ---------------------
    t0 = 0
    # ----------------------------------------

    # Prepare output folder and file path
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    file_path = os.path.join(output_folder, 'orbites.txt')

    # For each initial condition
    for initial_condition in initial_conditions:

        # State initial position, step size, number of steps
        init_pos = initial_condition[:2]
        h, steps = initial_condition[2:]
        
        # Compute rk78 multistep
        vars_track = numeric_algorithm.rk78_multistep(t0, init_pos, h, h_min, h_max,
                                                      tol, pendol_funct, int(steps))

        # Store data
        with open(file_path, "a") as file_id:
            np.savetxt(file_id, vars_track[:, 1:])  # , fmt="%d %lf %lf %lf %lf")
            file_id.write("\n")
        file_id.close()

if __name__ == "__main__":
    ARGS = args()
    main(ARGS)

