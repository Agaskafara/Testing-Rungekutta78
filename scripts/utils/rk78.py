import numpy as np

class rk78:
    def __init__(self):
        # step partitions
        self._alpha = np.array([0., 2./27., 1./9., 1./6., 5./12., 0.5,
                               5./6., 1./6., 2./3., 1./3., 1., 0., 1.])

        # k_j coefficients within f
        self._beta_vect = np.array([2./27., 1./36., 1./12., 1./24., 0., 1./8.,
                              5./12.,0., -25./16., 25./16., 0.05,0.,0., 0.25,
                              0.2, -25./108.,0.,0., 125./108., -65./27.,
                              2.*(125./108.), 31./300., 0.,0., 0., 61./225.,
                              -2./9., 13./900.,2., 0.,0.,-53./6., 704./45.,
                              -107./9., 67./90.,3.,  -91./108.,0.,0., 23./108.,
                              -976./135., 311./54., -19./60.,17./6., -1./12.,
                              2383./4100., 0.,0., -341./164., 4496./1025.,
                              -301./82., 2133./4100., 45./82., 45./164., 18./41.,
                              3./205.,0.,0.,0., 0.,-6./41., -3./205.,-3./41.,
                              3./41., 6./41.,0., -1777./4100., 0.,0., -341./164.,
                              4496./1025., -289./82., 2193./4100., 51./82.,
                              33./164.,12./41.,0.,1.])

        # Compute beta as matrix (add extra column full of zeros to allow matrix product)
        self._beta = np.zeros((13 * 13))
        beta_nonzero_indices = np.array([ [i * 13 + j] for i in range(13) for j in range(i)])
        self._beta[beta_nonzero_indices] = self._beta_vect.reshape(78, 1)
        self._beta = self._beta.reshape(13, 13)

        # k_j coefficients in the last sum for rk7
        self._c = np.array([41./840., 0., 0., 0., 0., 34./105., 9./35., 9./35.,
                            9./280., 9./280., 41./840.])

        # k_j coefficients in the last sum for rk8
        self._cp = np.array([0., 0., 0., 0., 0., 34./105., 9./35., 9./35.,
                            9./280., 9./280., 0., 41./840., 41./840.])

    def _compute_algorithm(self, funct):
        """Compute x(t+h) given the parameters specified in params."""

        # Contain rk7 and rk8 k_i values for 1 <= i <= 11 (7), 13 (8)
        k_values = np.zeros((13, self._params['dimension']))

        # Compute within step partitions
        h = self._params['step']
        partial_time = self._params['t'] + self._alpha * h

        # For each step partition
        for s in range(13):

            # compute the partial position
            partial_pos = (self._params['pos'] + h * self._beta[s, None].dot(k_values)).reshape(-1)

            # Evaluate f in t + a_s*h, x + h * (b_(s,i)*k_i) and update k value
            k_values[s, :] = funct(t=partial_time[s], pos=partial_pos)

        # Compute rk7 and rk8 estimations
        rk7_output = self._params['pos'] + h * self._c.dot(k_values[:11, :])
        rk8_output = self._params['pos'] + h * self._cp.dot(k_values)

        # return both values
        return [rk7_output, rk8_output]

    def _compute_condition_values(self, rk7_output_, rk8_output_):
        d = sum(abs(rk7_output_ - rk8_output_))/self._params['dimension']
        dd = sum(abs(rk8_output_))
        e3 = self._params['tol'] * (1 + dd/100)
        return {'d' : d, 'dd' : dd, 'e3' : e3}

    def rk78(self, t0 : float, init_pos : np.ndarray, h : float,
             h_min : float, h_max : float, tol : float, funct):
        
        # Define useful parameters
        self._params = {'t' : t0, 'pos' : init_pos, 'step' : h,
                       'dimension' : len(init_pos), 'tol' : tol}
        # Find 
        h_sign = (np.sign(h) if h != 0 else 1)

        # Until conditions are not fulfilled, reintegrate
        while True:
            
            # Compute rk7 and rk8 algorithms for current parameters
            rk7_output, rk8_output = self._compute_algorithm(funct)

            # Compute condition values
            cond_vals = self._compute_condition_values(rk7_output, rk8_output)

            # Are conditions fulfilled?
            if cond_vals['d'] < cond_vals['e3'] or abs(self._params['step']) <= h_min:
                
                # End loop 
                break
            
            # Else, modify h value and compute algorithm again
            new_h = abs(self._params['step']) * 0.9 * (cond_vals['e3']/cond_vals['d'])**(0.125)
            self._params['step'] = h_sign * max(new_h, h_min)
        
        # Once the procedure is done and we get a good enough prediction, we correct h
        cond_vals['d'] = max(cond_vals['d'], cond_vals['e3']/256.)

        # Fehlberg correction version (7.2.5.16)
        new_h = abs(self._params['step']) * 0.9 * (cond_vals['e3']/cond_vals['d'])**(0.125)
        new_h = h_sign * min(max(new_h, h_min), h_max)

        # Return values
        return {'t_new' : t0 + self._params['step'],
                'pos_new' : rk8_output, 'new_step' : new_h}

    def rk78_multistep(self, t0 : float, init_pos : np.ndarray, h : float,
                       h_min : float, h_max : float, tol : float, funct, steps: int):
        """Apply multiple times rk78."""

        # Pack variables in a dictionary
        self._params = {'t' : t0, 'pos' : init_pos, 'step' : h,
                       'dimension' : len(init_pos), 'tol' : tol}

        # Create a list with the tracked solutions
        vars_track = [[0, t0] + init_pos.tolist() + [h]]

        # For each step
        for iter in range(steps):

            # Retrieve current variables
            tn = vars_track[iter][1]
            h_n = vars_track[iter][-1]
            pos_n = vars_track[iter][2:-1]

            # Compute a single step prediction
            output = self.rk78(tn, np.array(pos_n), h_n, h_min, h_max, tol, funct)

            # Update track list
            vars_track.append([iter + 1, output['t_new']]+
                               output['pos_new'].tolist() + [output['new_step']])
        
        # Return list with tracked variables
        return np.array(vars_track)


    def flux(self, t_pred : float, t0 : float, init_pos : np.ndarray, h : float,
            h_min : int, h_max : int, tol : float, max_steps : int, funct):
        """Apply rk78 until reaching t_pred time."""
        # If t = t0, immediat return
        if t_pred == t0:
            return {'success': False}

        # Else
        # Right step sign
        h0 = (t_pred - t0)/abs(t_pred - t0) * abs(h)
        
        # Create dictionary with non-static variables
        params = {'tn' : t0, 'pos' : init_pos, 'step_size' : h0, 'last_values' : {}}

        # Compute rk78 until overcome t_pred or max_steps reached
        step_iter = 0
        while abs(params['tn'] - t0) < abs(t_pred - t0) and step_iter < max_steps:
            step_iter += 1

            # Predict next position
            output = self.rk78(params['tn'], params['pos'], params['step_size'],
                                            h_min, h_max, tol, funct)

            # Store last step values
            params['last_values'] = {'tn' : params['tn'], 'pos' : params['pos'],
                                     'step_size' : params['step_size']}

            # Update current values
            params['tn'] = output['t_new']
            params['pos'] = output['pos_new']
            params['step_size'] = output['new_step']
        
        # If max steps reached and tn doesnt overcome t_pred, return last value
        if step_iter >= max_steps:
            return {'success' : False}

        else:
            # Retrieve last iteration values
            last_params = params['last_values']

            # Compute new h value to reach t_pred
            h_dif = t_pred - last_params['tn']

            # Recompute the last step with h_dif step
            output = self.rk78(last_params['tn'], last_params['pos'], h_dif,
                                            h_dif, h_max, tol, funct)
            # Update this step values
            params['tn'] = output['t_new']
            params['pos'] = output['pos_new']
            params['step_size'] = output['new_step']

            # Return wether the procedure is succesful and pos obtained.
            return {'success' : True, 'pos' : output['pos_new']}
    
    def flux_multistep(self, t_pred_lim : float, num_evals : int,
                       t0 : float, init_pos : np.ndarray,
                       h : float, h_min : int, h_max : int,
                       tol : float, max_steps : int, funct):
        """Apply rk78 to evaluate multiple times between the given time interval."""

        # Correct step sign
        h0 = (t_pred_lim - t0)/abs(t_pred_lim - t0) * abs(h)

        # Create an equidistant time list
        eval_steps = np.linspace(t0, t_pred_lim, 1 + num_evals)

        # Create a list where solutions for evaluated times are stored
        pos_track = [init_pos]

        # For each evalued time
        for indx, t_target in enumerate(eval_steps):

            # Ignore the first step
            if indx == 0:
                continue

            # Else, compute the flux
            output = self.flux(t_target, eval_steps[indx - 1], pos_track[indx - 1],
                               h0, h_min, h_max, tol, max_steps, funct)
            # If the result was not succesful, explicit it and end procedure
            if not output['success']:
                return {'success': False}

            # Else, update tracking list
            pos_track.append(output['pos'])

        # Return whether it was succesful and the time and positions tracked
        return {'success': True, 'time_track' : eval_steps, 'pos_track' : pos_track}
        