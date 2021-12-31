import argparse
import os
import numpy as np
from utils.system_solver import qrres


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder', required=True,
                        help='Folder where file is stored.')
    return parser.parse_args()

def main(args):

    # Prepare output file
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    file_path = os.path.join(output_folder, "ratio.txt")

    # If the file already exists, erase it before storing data
    if os.path.exists(file_path):
        os.remove(file_path)

    # Create an increasing list of integers
    n_list = np.arange(40, 200)

    for n in n_list:

        # Init a random matrix
        A = np.random.random((n, n))
        x = np.ones(n)

        # Find b such that x is solution to Ax = b
        b = A.dot(x)
        
        # solve system
        success, _, _, qr_time = qrres(A, b, count_time=True)

        # Save the ratio in a file to plot
        with open(file_path, "a") as file_id:
            np.savetxt(file_id, np.array([[n + 1, 20000*qr_time/((4.*n**3)/3)]]))
        file_id.close()

        # Show results
        print("n = ", n)
        print("QR successful: ", success)
        print("Compute time:")
        print(qr_time, "\n")

if __name__ == "__main__":
              
    ARGS = args()
    main(ARGS)