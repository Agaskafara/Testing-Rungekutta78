import numpy as np
from utils.rk78 import rk78
from utils.include_diffs import include_diffs
from utils.system_solver import solve_system

class Maneuvers:

    def __init__(self, h  : float, h_min : float, h_max : float, tol : float, max_steps : int):
        # rk78 parameters
        self.h = h
        self.h_min = h_min
        self.h_max = h_max
        self.tol = tol
        self.max_steps = max_steps
        self._check_inputs()

    def _check_inputs(self):
        # All
        if not isinstance(self.h, float) or not isinstance(self.h_min, float) \
            or not isinstance(self.h_max, float) or not isinstance(self.tol, float) \
            or not isinstance(self.max_steps, int):
            raise ValueError("Input error")
        # h_min and h_max
        if self.h_min > self.h_max or self.h_min <= 0:
            raise ValueError("Wrong h_min, h_max values")
        # tol
        if self.tol <= 0:
            raise ValueError("Wrong tol value")
        # max_steps
        if self.max_steps <= 0:
            raise ValueError("Wrong max_steps value")

    def cmani_gdg(self, step_time : float, start_pos : np.ndarray,
                  end_pos : np.ndarray, boosts: np.ndarray, funct, dfunct, check_inputs : bool = False):

        # If needed, check inputs
        if check_inputs:
            # Check right inputs
            # - Time
            if not isinstance(step_time, float) and not isinstance(step_time, int):
                raise ValueError("Wrong step_time input")
            if step_time <= 0:
                raise ValueError("step_time must be > 0")
            # - Pos
            if not isinstance(start_pos, np.ndarray) or not isinstance(end_pos, np.ndarray):
                raise ValueError("Wrong start_pos or end_pos input")
            try:
                funct(0, start_pos)
                funct(0, end_pos)
                dfunct(0, start_pos)
            except ValueError:
                raise ValueError("Wrong start_pos input")
            # boosts
            if not isinstance(boosts, np.ndarray):
                raise ValueError("Wrong boosts input")
            if boosts.shape[0] != 2:
                raise ValueError("boosts length must be 2")
            if boosts.shape[1] > len(start_pos):
                raise ValueError("Each boost length must be within pos dimensions")
        
        # Init rk78
        numeric_integrator = rk78()

        # Init include_diffs
        inc_diffs = include_diffs(funct, dfunct)
        camp = inc_diffs.camp

        # Compute function dimension and boost dimension
        funct_dim = len(start_pos)
        boost_dim = boosts.shape[1]

        # Prepare init_pos
        init_pos = np.zeros(funct_dim*(funct_dim + 1))
        init_pos[:funct_dim] = start_pos.copy()
        init_pos[boost_dim:funct_dim] += boosts[0]
        init_pos[funct_dim:] = np.identity(funct_dim).reshape(funct_dim*funct_dim)

        # Compute solution at time step_time/2
        phi_1 = numeric_integrator.flux(step_time/2., 0, init_pos,
                                        self.h, self.h_min, self.h_max,
                                        self.tol, self.max_steps, camp)
        if not phi_1['success']:
            raise Exception("r78 flux not successful")

        # Prepare second boost and compute solution
        pos_boost = np.zeros_like(phi_1['pos'])
        pos_boost[boost_dim:funct_dim] += boosts[1]
        phi_2 = numeric_integrator.flux(step_time, step_time/2., phi_1['pos'] + pos_boost,
                                        self.h, self.h_min, self.h_max,
                                        self.tol, self.max_steps, camp)
        
        if not phi_2['success']:
            raise Exception("r78 flux not successful")
        
        # Init diffG matrix
        diffG = np.zeros((funct_dim, 2*boost_dim))
        
        # Store first boost diff
        diffG[:, :boost_dim] = phi_2['pos'][funct_dim:].reshape(funct_dim, funct_dim).dot( \
            phi_1['pos'][funct_dim:].reshape(funct_dim, funct_dim)[:, boost_dim:])
        # Store second boost diff
        diffG[:, boost_dim:] = phi_2['pos'][funct_dim:].reshape(funct_dim, funct_dim)[:, boost_dim:]
        
        return {'G' : phi_2['pos'][:funct_dim] - end_pos, 'DG' : diffG}


    def cmani(self, step_time : float, start_pos : np.ndarray,
              end_pos : np.ndarray, init_boosts: np.ndarray, funct, dfunct, max_it : int):
        # Init working variable
        work_boosts = init_boosts.copy().astype(float)
        boosts_shape = work_boosts.shape

        # Apply Newton's method max_it times
        for _ in range(max_it) :
            # Compute G and DG
            output = self.cmani_gdg(step_time, start_pos, end_pos, work_boosts, funct, dfunct)

            # Solve system and find increase on work_boosts
            work_boosts += solve_system(output['DG'], - output['G'], tol=self.tol).reshape(boosts_shape)
            
            if output['G'].dot(output['G']) < self.tol:
                break
        
        # Return the expected root for G
        return work_boosts
