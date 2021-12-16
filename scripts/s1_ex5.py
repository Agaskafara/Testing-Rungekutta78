import argparse
import os
from functools import partial
import numpy as np
from utils.functs import lorentz_funct
from utils.rk78 import rk78

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', required=True,
                        help='Folder where file is stored.')
    return parser.parse_args()

def main(args):

    # Create rk78 instance
    numeric_algorithm = rk78()

    # Hyperparameters ------------------------
    h, h_min, h_max = 0.01, 0.01, 0.05
    max_steps = 1000
    tol = 0.0001
    # ----------------------------------------
    # Initial conditions and time prediction -
    t0, t_pred = 0, 100
    init_pos = np.array([- 0.416,
                         - 0.909,
                         0.014])
    num_evals = 10000
    # ----------------------------------------
    # Function parameters --------------------
    sigma, p, b = 3, 26.5, 1.0
    # ----------------------------------------

    # Prepare output folder and filepath
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    file_path = os.path.join(output_folder, 'orbites.txt')

    # Predict
    output = numeric_algorithm.flux_multistep(t_pred, num_evals, t0, init_pos,
                                              h, h_min, h_max, tol, max_steps,
                                              partial(lorentz_funct, sigma_=sigma, p_=p, b_=b))

    # Store values
    with open(file_path, "w") as file_id:
        np.savetxt(file_id, np.concatenate((output['time_track'].reshape(-1, 1),
                                           np.array(output['pos_track'])), axis=1))
    file_id.close()


if __name__ == "__main__":
    ARGS = args()
    main(ARGS)