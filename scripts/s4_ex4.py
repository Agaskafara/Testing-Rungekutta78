import argparse
import os
import numpy as np
from functools import partial
from utils.functs import rtbp_funct, rtbp_dfunct
from utils.maneuvers import Maneuvers


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--init_conditions', required=True,
                        help='txt file with initial conditions.')
    parser.add_argument('--output_folder', required=True,
                        help='Folder where file is stored.')
    return parser.parse_args()

def main(args):

    # Prepare output file
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    file_path = os.path.join(output_folder, "boosts.txt")

    # If the file already exists, erase it before storing data
    if os.path.exists(file_path):
        os.remove(file_path)

    # rk78 paramaters ----------------
    h, h_min, h_max = 0.0005, 0.00005, 0.001
    max_steps = 100000
    tol = 1e-15
    # --------------------------------
    # Function parameters ------------
    mu = 0.01215058560962404
    # --------------------------------
    # Initial values -----------------
    inputs = np.loadtxt(args.init_conditions)
    step_time = inputs[0]
    init_pos = inputs[1:7]
    end_pos = inputs[7:]
    boosts = np.array([[0, 0, 0], [0, 0, 0]])
    # --------------------------------

    # Create a Maneuvers instance
    man = Maneuvers(h, h_min, h_max, tol, max_steps)
    max_it = 20
    
    # Find right boosts
    output = man.cmani(step_time, init_pos, end_pos, boosts,
                       partial(rtbp_funct, mu_=mu), partial(rtbp_dfunct, mu_=mu), max_it)
    
    # Save right boosts
    with open(file_path, "a") as file_id:
        np.savetxt(file_id, output)
    file_id.close()


if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)