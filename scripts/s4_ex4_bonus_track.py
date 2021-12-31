import argparse
import os
from functools import partial
import numpy as np
from utils.functs import rtbp_funct
from utils.rk78 import rk78

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--init_conditions', required=True,
                        help='txt file with initial conditions.')
    parser.add_argument('--boosts', required=True,
                        help='txt file with boosts.')
    parser.add_argument('--output_folder', required=True,
                        help='Folder where file is stored.')
    return parser.parse_args()

def main(args):

    # Load inputs
    inputs = np.loadtxt(args.init_conditions)
    boosts = np.loadtxt(args.boosts)

    # Create rk78 instance
    numeric_algorithm = rk78()

    # Hyperparameters ------------------------
    h, h_min, h_max = 0.0005, 0.0001, 0.001
    max_steps = 10000
    tol = 1e-15
    # ----------------------------------------
    # Parameters -----------------------------
    mu = 0.01215058560962404
    # ----------------------------------------
    # Initial conditions and time prediction -
    t0, t_pred = 0, inputs[0]
    # ----------------------------------------
    # Manual choice of initial position/velocities
    init_pos = inputs[1:7]
    init_pos[3:] += boosts[0]  # Add boosts
    # ----------------------------------------
    # Number of evaluations ------------------
    num_evals = 1000
    # ----------------------------------------

    # Prepare output folder
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Files path
    paths = {'orbites': os.path.join(output_folder, 'orbites.txt'),
             'points': os.path.join(output_folder, 'points.txt')}

    # Predict first boost
    output_boost1 = numeric_algorithm.flux_multistep(t_pred/2, num_evals, t0, init_pos,
                                                     h, h_min, h_max, tol, max_steps,
                                                     partial(rtbp_funct, mu_=mu))
    middle_pos = output_boost1['pos_track'][-1].copy()
    middle_pos[3:] += boosts[1]

    # Predict second boost
    output_boost2 = numeric_algorithm.flux_multistep(t_pred, num_evals, t_pred/2., middle_pos,
                                                     h, h_min, h_max, tol, max_steps,
                                                     partial(rtbp_funct, mu_=mu))
    # Add the following satellite trajectory
    output_extra = numeric_algorithm.flux_multistep(t_pred*3, num_evals, t_pred, 
                                                    output_boost2['pos_track'][-1].copy(),
                                                    h, h_min, h_max, tol, max_steps,
                                                    partial(rtbp_funct, mu_=mu))

    # Extract the position coordinates
    pos_track_boost1 = np.concatenate((output_boost1['time_track'].reshape(-1, 1),
                                       np.array(output_boost1['pos_track'])[:, :3]), axis=1)
    pos_track_boost2 = np.concatenate((output_boost2['time_track'].reshape(-1, 1),
                                       np.array(output_boost2['pos_track'])[:, :3]), axis=1)
    pos_track_extra = np.concatenate((output_extra['time_track'].reshape(-1, 1),
                                      np.array(output_extra['pos_track'])[:, :3]), axis=1)

    # Store values
    # File Orbites
    with open(paths['orbites'], "w") as file_id:
        np.savetxt(file_id, pos_track_boost1)
        file_id.write("\n\n")
        np.savetxt(file_id, pos_track_boost2)
        file_id.write("\n\n")
        np.savetxt(file_id, pos_track_extra)
    file_id.close()

    # File Points
    with open(paths['points'], "w") as file_id:
        np.savetxt(file_id, [[mu - 1 , 0, 0], [mu, 0, 0]])
    file_id.close()


if __name__ == "__main__":
    ARGS = args()
    main(ARGS)