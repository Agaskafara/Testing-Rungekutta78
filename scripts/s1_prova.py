import argparse
import os
from functools import partial
import numpy as np
from utils.rk78 import rk78



def sample_funct(t : float, pos : np.ndarray):
    x, y, z = pos
    return np.array([-y+np.cos(10*z), x+np.sin(10*z), 0.1])

numeric_algorithm = rk78()
outse = numeric_algorithm.flux_multistep(t_pred_lim=10, num_evals=1000, t0=0,
                            init_pos=np.array([0,0,0]), h=0.01, 
                            h_min=0.01, h_max=0.05,tol=0.0001, max_steps=1000, funct=sample_funct)
file_path = "./orbites.txt"
# Store values
# File 1
with open(file_path, "a") as file_id:
    np.savetxt(file_id, np.concatenate((outse['time_track'].reshape(-1, 1),
                                        np.array(outse['pos_track'])[:, :3]), axis=1))
file_id.close()