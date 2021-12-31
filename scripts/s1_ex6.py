import argparse
import os
from functools import partial
import numpy as np
from utils.functs import rtbp_funct
from utils.rk78 import rk78

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--halos_input', required=True,
                        help='Halos file txt with parameters.')
    parser.add_argument('--output_folder', required=True,
                        help='Folder where file is stored.')
    return parser.parse_args()

def main(args):

    # Create rk78 instance
    numeric_algorithm = rk78()

    # Hyperparameters ------------------------
    h_min, h_max = 0.0001, 0.005
    tol = 0.000001
    # ----------------------------------------
    # Parameters -----------------------------
    mu = 0.01215058560962404
    # ----------------------------------------
    # Initial conditions ---------------------
    t0 = 0
    input_vars = np.loadtxt(args.halos_input)
    # ----------------------------------------

    # Prepare output file
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Files path
    paths = {'orbites': os.path.join(output_folder, 'orbites.txt'),
             'points': os.path.join(output_folder, 'points.txt')}
    
    # If the file already exists, erase it before storing data
    if os.path.exists(paths['orbites']):
        os.remove(paths['orbites'])

    # For each initial condition
    for indx_par, halo_params in enumerate(input_vars):

        init_pos = halo_params[:6]
        h, num_evals = halo_params[6:]

        # Predict
        output = numeric_algorithm.rk78_multistep(t0, init_pos, h, h_min, h_max, tol,
                                              partial(rtbp_funct, mu_=mu), steps=int(num_evals))

        # Store values
        # File Orbites
        with open(paths['orbites'], "a") as file_id:
            np.savetxt(file_id, output[:, 2:-1])
            file_id.write("\n\n")
        file_id.close()

        # File Points
        if indx_par == 0:
            with open(paths['points'], "w") as file_id:
                np.savetxt(file_id, [[mu - 1 , 0, 0], [mu, 0, 0]])
            file_id.close()


if __name__ == "__main__":
    ARGS = args()
    main(ARGS)