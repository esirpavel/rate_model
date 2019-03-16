import matplotlib.pyplot as pl
import numpy as np

class RateNetwork:
   
    @staticmethod
    def gFun(x, alpha=1.5, cutoff=10.):
        if type(x) == np.ndarray:
            res = x.copy()
            res[np.nonzero(x < cutoff)] = alpha*np.log1p(np.exp(x[np.nonzero(x < cutoff)]/alpha))
            return res
        else:
            return alpha*np.log1p(np.exp(x/alpha)) if x < cutoff else x

    @staticmethod
    def get_angle(m, axis=-1):
        pos = np.linspace(-np.pi, np.pi, m.shape[axis], endpoint=False)
        return np.angle(np.sum(np.exp(1j*pos)*m, axis=axis))

    def __init__(self, sim_time, dt, sampl_dt, N):
        self.sim_time = sim_time
        self.dt = dt
        self.sampl_dt = sampl_dt
        self.n_istep = int(sim_time/dt) # number of integration time steps
        
        self.N = N
        
        self.pos = np.linspace(-np.pi, np.pi, N, endpoint=False)
        X, Y = np.meshgrid(self.pos, self.pos)
        self.dx = X - Y
        self.dx = np.mod(self.dx + np.pi, 2*np.pi) - np.pi
        del X, Y
        
        self.x = np.ones(N)
        self.u = np.ones(N)
        self.hE = np.zeros(N)
        
        n_rstep = int((sim_time)/sampl_dt)
        self.ActU = np.zeros((n_rstep, N), dtype='float32')
        self.ActX = np.zeros((n_rstep, N), dtype='float32')
        self.ActHE = np.zeros((n_rstep, N), dtype='float32')
        self.ActHI = np.zeros(n_rstep, dtype='float32')
        self.ActRE = np.zeros((n_rstep, N), dtype='float32')
        self.ActRI = np.zeros(n_rstep, dtype='float32')
        self.tm = np.linspace(0, sim_time, n_rstep)
    
    @staticmethod
    def trunc_cos(x, J0, J1):
        indices = np.nonzero(np.abs(x) < np.arccos(J0/J1))
        res = x.copy()
        res[:] = J0
        res[indices] = J1*np.cos(x[indices])
        return res    
    
    def set_weights(self, J0, J1, J_EI, J_IE, eps=0., conn_width=1.3, 
                    conn_type='gauss', seed=0):
        self.J0 = J0
        self.J1 = J1
        self.J_EI = J_EI
        self.J_IE = J_IE
        
        np.random.seed(seed)
        W_heter = eps*np.random.randn(self.N, self.N)/np.sqrt(self.N)
        
        if conn_type == 'cos':
            self.W = (J0 + J1*np.cos(self.dx))/self.N + W_heter
        elif conn_type == 'gauss':
            W =  np.exp(-self.dx**2/(2*conn_width**2))
            self.W = (self.J0 + self.J1*W/np.amax(W))/self.N + W_heter
        elif conn_type == 'trunc_cos':
            self.W = self.trunc_cos(self.dx/conn_width, J0, J1)/self.N + W_heter
        else:
            raise RuntimeError('Non-existent connectivity kernel type')
    
    def set_initial_values(self, hE=0., hI=0., x=1., u=None):
        self.hE_iv = hE
        self.hI_iv = hI
        self.x_iv = x
        if u is None:
            self.u_iv = self.U
        else:
            self.u_iv = u
    
    def set_params(self, U, I0, tau_d, tau_f, tau, alpha):
        self.U = U
        self.I0 = I0
        self.tau_d = tau_d
        self.tau_f = tau_f
        self.tau = tau
        self.alpha = alpha
        
        self.u[:] = self.U
    
    def plot_simul(self, figsize = (10, 8)):
        pl.figure(figsize=figsize)
        pl.pcolormesh(self.tm, np.degrees(self.pos), self.ActRE.T)
        pl.plot(self.tm, np.degrees(self.get_angle(self.ActU)), lw=3., c='C3')
        pl.plot(self.tm, np.degrees(self.get_angle(self.ActRE)), lw=3., c='C1')
        pl.xlim((0, self.sim_time))
        pl.ylim((-180, 180))
        pl.xlabel('Time (s)')
        pl.ylabel(r'$\theta$')
        pl.colorbar()
    
    def process_stumulus(self, t):
        if self.stim_idx < len(self.stim_duration):
            if ((not self.is_stimulating) and 
                    t == int(self.stim_start[self.stim_idx]/self.dt)):
                dx1 = np.mod(self.pos - np.radians(self.stim_pos[self.stim_idx]) + np.pi, 2*np.pi) - np.pi
                self.is_stimulating = True
                if self.stim_type[self.stim_idx] == 'cos':
                    self.Iext[:] = self.stim_ampl[self.stim_idx]*np.cos(dx1/self.stim_width[self.stim_idx])
                elif self.stim_type[self.stim_idx] == 'gauss':
                    self.Iext[:] = self.stim_ampl[self.stim_idx]*np.exp(-(dx1/self.stim_width[self.stim_idx])**2/2)
                elif self.stim_type[self.stim_idx] == 'trunc_cos':
                    self.Iext[:] = self.stim_ampl[self.stim_idx]*self.trunc_cos(dx1/self.stim_width, 0, 1)
            elif (self.is_stimulating and 
                    t == int((self.stim_start[self.stim_idx] + 
                              self.stim_duration[self.stim_idx])/self.dt)):
                self.is_stimulating = False
                self.Iext[:] = 0.
                self.stim_idx += 1
    
    def set_stimuli(self, stim_start, stim_duration, stim_ampl, stim_pos, 
                    stim_width, stim_type):
        self.Iext = np.zeros(self.N)
        self.stim_idx = 0
        self.is_stimulating = False
        
        self.stim_start = stim_start
        self.stim_duration = stim_duration
        self.stim_ampl = stim_ampl
        self.stim_pos = stim_pos
        self.stim_width = stim_width
        self.stim_type = stim_type
        
    def simulate_facil(self, backend='python'):
        self.x[:], self.u[:], self.hE[:], self.hI = self.x_iv, self.u_iv, self.hE_iv, self.hI_iv
        
        if backend == 'c':
            import cycover_ring as rn
            rn.set_calc_params(self.N, self.n_istep, self.dt, int(self.sampl_dt/self.dt))
            rn.set_params(self.U, self.J_IE, self.J_EI, self.tau, self.tau_d, 
                          self.tau_f, self.I0, self.alpha)
            rn.init_arrays(self.x, self.u, self.hE, self.hI, self.W, self.ActX, 
                           self.ActU, self.ActHE, self.ActHI)
            
            # @TODO: maybe it would be better to move this dictionary out of here?!
            map_dict = {'cos': 1, 'gauss': 2, 'trunc_cos': 3, 'mask': 4}
            st_type = np.array([*map(map_dict.get, self.stim_type)]).astype('uint32')
            rn.set_stimuli(
                (np.array(self.stim_start)/self.dt).astype('uint32'), 
                (np.array(self.stim_duration)/self.dt).astype('uint32'), 
                np.array(self.stim_ampl), 
                np.radians(self.stim_pos), 
                np.array(self.stim_width), 
                st_type,
                len(self.stim_start))
            
            rn.integrate()
            
            self.ActRE = self.gFun(self.ActHE)
            self.ActRI = self.gFun(self.ActHI)
        elif backend == 'python':
            for i in range(0, self.n_istep):
                self.process_stumulus(i)
                
                rE = self.gFun(self.hE, self.alpha)
                rI = self.gFun(self.hI, self.alpha)
                
                dx = (1.0 - self.x)/self.tau_d - self.u*self.x*rE
                du = (self.U - self.u)/self.tau_f + self.U*(1 - self.u)*rE
                dhE = (-self.hE + np.dot(self.W, self.u*self.x*rE) - self.J_EI*rI + self.I0 + self.Iext)/self.tau
                dhI = (-self.hI + self.J_IE*np.sum(rE)/self.N)/self.tau
                
                self.x[:] = self.x + dx*self.dt
                self.u[:] = self.u + du*self.dt
                self.hE[:] = self.hE + dhE*self.dt
                self.hI = self.hI + dhI*self.dt
                
                if i % int(self.sampl_dt/self.dt) == 0:
                    rec_i = int(i/(self.sampl_dt/self.dt))
                    self.ActU[rec_i] = self.u
                    self.ActX[rec_i] = self.x
                    self.ActHE[rec_i] = self.hE
                    self.ActHI[rec_i] = self.hI
                    self.ActRE[rec_i] = rE
                    self.ActRI[rec_i] = rI
        else:
            raise RuntimeError('Non-existent backend')

    @staticmethod
    def init_all_params(sim_time, dt, sampl_dt, N, J0, J1, J_EI, J_IE, 
                        eps, conn_width, conn_type, seed, U, I0, tau_d, tau_f, tau, alpha):

        self = RateNetwork(sim_time, dt, sampl_dt, N)
        self.set_weights(J0, J1, J_EI, J_IE, eps, conn_width, conn_type, seed)
        self.set_params(U, I0, tau_d, tau_f, tau, alpha)
        self.set_initial_values(u=U)
        
        return self
