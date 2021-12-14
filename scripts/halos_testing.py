import argparse
import os
from functools import partial
import numpy as np
from utils.functs import rtbp_funct
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
    h, h_min, h_max = 0.0005, 0.0001, 0.001
    max_steps = 10000
    tol = 0.00001
    # ----------------------------------------

    # Parameters -----------------------------
    mu = 0.01215058560962404
    # ----------------------------------------

    # Initial conditions and time prediction -
    t0, t_pred = 0, 1.5*np.pi
    # ----------------------------------------

    # Manual choice of initial position/velocities
    rot = np.array([[np.cos(t0), -np.sin(t0), 0], [np.sin(t0), np.cos(t0), 0], [0, 0, 1]])
    inv_rot = np.array([[np.cos(t0), np.sin(t0), 0], [-np.sin(t0), np.cos(t0), 0], [0, 0, 1]])
    dinv_rot = np.array([[- np.sin(t0), np.cos(t0), 0], [-np.cos(t0), -np.sin(t0), 0], [0, 0, 1]])
    orig_vel = np.array([0, 0, 0])
    sinodic_pos = [20, 0, 0]
    sinodic_vel = dinv_rot.dot(rot).dot(sinodic_pos) + inv_rot.dot(orig_vel)
    init_pos = np.array(sinodic_pos + list(sinodic_vel))

    num_evals = 1000
    # ----------------------------------------

    # Prepare output folder
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Files path
    paths = {'orbites': os.path.join(output_folder, 'orbites_test.txt'),
             'points': os.path.join(output_folder, 'points_test.txt')}

    # Predict
    output = numeric_algorithm.flux_multistep(t_pred, num_evals, t0, init_pos,
                                              h, h_min, h_max, tol, max_steps,
                                              partial(rtbp_funct, mu_=mu))

    # Store values
    # File Orbites
    with open(paths['orbites'], "w") as file_id:
        np.savetxt(file_id, np.concatenate((output['time_track'].reshape(-1, 1),
                                            np.array(output['pos_track'])[:, :3]), axis=1))
    file_id.close()

    # File Points
    with open(paths['points'], "w") as file_id:
        np.savetxt(file_id, [[mu - 1 , 0, 0], [mu, 0, 0]])
    file_id.close()


if __name__ == "__main__":
    ARGS = args()
    main(ARGS)