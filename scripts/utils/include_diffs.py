import numpy as np

class include_diffs:
    def __init__(self, funct, dfunct):
        self.funct = funct
        self.dfunct = dfunct

    def camp(self, t : float, pos : np.ndarray):
        # Find function dimension.
        # len(pos) should be n*(n+1) where n is the function dimension.
        # Therefore,
        n = (-1 + np.sqrt(1 + 4*len(pos)))/2

        # If the computed n is not an integer, raise an exception
        if n - int(n) > 0.000001:
            raise ValueError("Variable \"pos\": Wrong dimension")
        n = int(n)
        
        # Given pos, compute the A matrix
        A = pos[n:].reshape(n, n)

        # Compute and return the camp funct value
        return np.concatenate([self.funct(t, pos[:n]),
                               self.dfunct(t, pos[:n]).dot(A).reshape(-1)])
