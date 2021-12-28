import argparse
import numpy as np
from functools import partial
from utils.functs import rtbp_funct, rtbp_dfunct
from utils.maneuvers import Maneuvers


def args():
    parser = argparse.ArgumentParser()

    return parser.parse_args()

def main(args):

    # rk78 paramaters ----------------
    h, h_min, h_max = 0.005, 0.001, 0.01
    max_steps = 100000
    tol = 1e-12
    # --------------------------------
    # Function parameters ------------
    mu = 0.01215058560962404
    # --------------------------------
    # Initial values -----------------
    step_time = 2.746177778348061
    init_pos = np.array([-0.8457719086638686, 0.05934028672713117, 0,
                         0.03769526738828564, 0.01515062570613701, 0.04577799834681743])
    end_pos = np.array([-0.8518851062423166, 0.06132143191581388, 0.0005120725158794939,
                        0.0206128163364461, 0.01844596878578819, 0.04641615142301939])
    boosts = np.array([[0, 0, 0], [0, 0, 0]])
    # --------------------------------

    # Create a Maneuvers instance
    man = Maneuvers(h, h_min, h_max, tol, max_steps)
    max_it = 20
    
    # Find right boosts
    output = man.cmani(step_time, init_pos, end_pos, boosts,
                       partial(rtbp_funct, mu_=mu), partial(rtbp_dfunct, mu_=mu), max_it)
    print(output)


if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)