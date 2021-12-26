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
    t0, t_pred = 0, 2.746177778348061
    # ----------------------------------------

    # Manual choice of initial position/velocities
    init_pos = np.array([-0.8457719086638686, 0.05934028672713117, 0,
                         0.03769526738828564, 0.01515062570613701, 0.04577799834681743])
    # Add boosts
    init_pos[3:] += np.array([-9.35186676e+04, 5.40396394e+04, -9.18213888e-02])

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
    output_boost1 = numeric_algorithm.flux_multistep(t_pred/2, num_evals, t0, init_pos,
                                                     h, h_min, h_max, tol, max_steps,
                                                     partial(rtbp_funct, mu_=mu))
    middle_pos = output_boost1['pos_track'][-1].copy()
    middle_pos[3:] += np.array([-3.14638200e+04, -1.08432204e+05, 9.24594824e-02])

    output_boost2 = numeric_algorithm.flux_multistep(t_pred, num_evals, t_pred/2., middle_pos,
                                                     h, h_min, h_max, tol, max_steps,
                                                     partial(rtbp_funct, mu_=mu))
    
    pos_track_boost1 = np.concatenate((output_boost1['time_track'].reshape(-1, 1),
                                       np.array(output_boost1['pos_track'])[:, :3]), axis=1)
    pos_track_boost2 = np.concatenate((output_boost2['time_track'].reshape(-1, 1),
                                       np.array(output_boost2['pos_track'])[:, :3]), axis=1)
    full_pos_track = np.concatenate([pos_track_boost1, pos_track_boost2])

    # Store values
    # File Orbites
    with open(paths['orbites'], "w") as file_id:
        np.savetxt(file_id, full_pos_track)
    file_id.close()

    # File Points
    with open(paths['points'], "w") as file_id:
        np.savetxt(file_id, [[mu - 1 , 0, 0], [mu, 0, 0]])
    file_id.close()


if __name__ == "__main__":
    ARGS = args()
    main(ARGS)