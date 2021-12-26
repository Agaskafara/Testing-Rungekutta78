import numpy as np
import time

   


# Numpy function returns 0 for sign(0), and we need sign(0) = 1
# and we don't use an "if" since we are working with floats
def _sign(val: float):
    numpy_sign = np.sign(val)
    return 1. + numpy_sign*(1 - numpy_sign)

# QR factorization
def qrres(A : np.ndarray, b : np.ndarray, tol : float = 0.0005, count_time : bool = False):
    
    # Matrix dimensions
    m, n = A.shape

    # Check dimensions of input matrices
    if m != b.shape[0] or m < n:
        raise ValueError("Wrong dimensions")

    # Workbench variables
    A_k = A.copy().astype(float)
    b_k = b.copy().astype(float)

    # For now let's say there is success
    success = True

    # Start counting time
    start_time = time.time()

    # For each iteration, build A_k and b_k
    for k in range(n):

        # Compute gamma
        u_k = A_k[k:, k].copy()
        gamma_k = np.sqrt(np.sum(u_k**2))

        # Compute u
        u_k[0] += _sign(u_k[0])*gamma_k
        
        # if u_k is zero, more iterations won't change A_k and b_k
        if u_k.dot(u_k) < tol:
            print("break at step " + str(k + 1))
            success = False
            break

        # Compute P
        P_k = np.identity(m - k) - (2/u_k.dot(u_k))*u_k[None].T.dot(u_k[None])
        
        # Update A_k
        A_k[k:, k:] = P_k.dot(A_k[k:, k:])

        # Update b_k
        if not count_time:
            b_k[k:] = P_k.dot(b_k[k:])

    # Stop counting
    time_end = time.time()

    return [success, A_k, b_k] + ([time_end - start_time] if count_time else [0])

# System solver from QR factorization
def solve_qr(Rst : np.ndarray, Qtb : np.ndarray, checkQR : bool = False, tol : float = 0.0005):
    
    # Matrix Rst dimensions
    m, n = Rst.shape

    # check the system is QR factorized
    if checkQR:
        # Check dimensions
        if m != Qtb.shape[0] or m < n:
            raise ValueError("Wrong dimensions")
        
        # Check R is triangular
        for col in range(n):
            if Rst[(col+1):, col].dot(Rst[(col+1):, col]) > tol:
                raise ValueError("R is not triangular")
        if True in (abs(np.diag(Rst)) < tol):
            raise ValueError("R is not triangular")

    # Backward algorithm

    # Initialize the solution
    partial_sol = np.zeros(n, dtype=float)

    # For every row, starting from n row
    for row in range(n -1, -1, -1):
        # Update the solution rows
        partial_sol[row] = (Qtb[row] - Rst[row].dot(partial_sol))/Rst[row, row]

    # Return the solution or partial solution
    return partial_sol

# System solver (Global)
def solve_system(A : np.ndarray, b : np.ndarray, tol : float = 0.0005):
    
    # transform system into a QR system
    success, Rst, Qtb, _ = qrres(A, b, tol=tol,count_time=False)
    
    # If not possible, return no solution
    if not success:
        print("No solution available")
        pass

    # Else, solve the QR system without checking inputs
    # and return the solution
    return solve_qr(Rst, Qtb, checkQR=False, tol=tol)

