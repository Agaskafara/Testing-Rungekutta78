import argparse
import numpy as np
from utils.functs import pendol_funct, pendol_dfunct
from utils.maneuvers import Maneuvers


def args():
    parser = argparse.ArgumentParser()

    return parser.parse_args()

def main(args):

    # rk78 paramaters ----------------
    h, h_min, h_max = 0.0005, 0.0001, 0.001
    max_steps = 10000
    tol = 0.0001
    # --------------------------------
    # --------------------------------
    # Initial values -----------------
    step_time = np.pi/2.
    init_pos = np.array([1, 0])
    end_pos = np.array([0, -np.sqrt(2*(1 - np.cos(1.)))])
    boosts = np.array([[0], [0]])
    # --------------------------------

    # Create a Maneuvers instance
    man = Maneuvers(h, h_min, h_max, tol, max_steps)
    max_it = 4
    
    # Find right boosts
    output = man.cmani(step_time, init_pos, end_pos, boosts,
                       pendol_funct, pendol_dfunct, max_it)
    print(output)


if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)