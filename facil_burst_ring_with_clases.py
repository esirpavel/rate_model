import matplotlib.pyplot as pl
import numpy as np

class RateNetwork:

    params_arr = [
        'sim_time', 'dt', 'sampl_dt', 'N', 'J0', 'J1', 'J_EI', 'J_IE', 
        'eps', 'conn_width', 'seed', 'U', 'I0', 'tau_d', 'tau_f', 'tau', 'alpha', 
        'stim_width', 'stim_ampl', 'stim_duration', 'stim_pos'
    ]
          
    def gFun(self, x, alpha=1.5, cutoff=10.):
        res = x.copy()
        res[np.nonzero(x < cutoff)] = alpha*np.log1p(np.exp(x[np.nonzero(x < cutoff)]/alpha))
        return res

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
        
        self.Iext = np.zeros(N)
        self.x = np.ones(N)
        self.u = np.ones(N)
        self.hE = np.zeros(N)
        
        n_rstep = int((sim_time + sampl_dt)/sampl_dt)
        self.ActU = np.zeros((n_rstep, N), dtype='float32')
        self.ActX = np.zeros((n_rstep, N), dtype='float32')
        self.ActHE = np.zeros((n_rstep, N), dtype='float32')
        self.ActHI = np.zeros(n_rstep, dtype='float32')
        self.ActRE = np.zeros((n_rstep, N), dtype='float32')
        self.ActRI = np.zeros(n_rstep, dtype='float32')
        self.tm = np.linspace(0, sim_time, n_rstep)
    
    def set_weights(self, J0, J1, J_EI, J_IE, eps=0., conn_width=1.3, conn_type='gauss', seed=0):
        self.J0 = J0
        self.J1 = J1
        self.J_EI = J_EI
        self.J_IE = J_IE
        
        np.random.seed(seed)
        W_heter = eps*np.random.randn(self.N, self.N)/np.sqrt(self.N)
        
        if conn_type == 'cos':
            self.W = (J0 + J1*np.cos(self.dx))/self.N + W_heter
        elif conn_type == 'trunc_cos':
            conn_kernel = lambda x, J0=J0, J1=J1:  J1*np.cos(x) if np.abs(x) < np.arccos(J0/J1) else J0
            conn_kernel_vect = np.vectorize(conn_kernel)
            self.W = conn_kernel_vect(conn_width*self.dx)/self.N + W_heter
        elif conn_type == 'gauss':
            W =  np.exp(-self.dx**2/(2*conn_width**2))
            self.W = (self.J0 + self.J1*W/np.amax(W))/self.N + W_heter
        else:
            raise RuntimeError('Non-existent connectivity kernel type')
    
    def set_initial_values(self, hE=0., hI=0., x=1., u=1.):
        self.hE_iv = hE
        self.hI_iv = hI
        self.x_iv = x
        self.u_iv = u
    
    def set_params(self, U, I0, tau_d, tau_f, tau, alpha):
        self.U = U
        self.I0 = I0
        self.tau_d = tau_d
        self.tau_f = tau_f
        self.tau = tau
        self.alpha = alpha
        
        self.u[:] = self.U
        
    def set_stim_params(self, stim_width, stim_ampl, stim_duration, stim_pos):
        self.stim_width = stim_width
        self.stim_ampl = stim_ampl
        self.stim_duration = stim_duration
        self.stim_pos = stim_pos
    
    
    def plot_simul(self):
        pl.figure(figsize = (10, 8))
        pl.pcolormesh(self.tm, np.degrees(self.pos), self.ActRE.T)
        pl.plot(self.tm, np.degrees(self.get_angle(self.ActRE)), lw=3., c='C3')
        pl.xlim((0, self.sim_time))
        pl.xlabel('Time (s)')
        pl.ylabel(r'$\theta$')
        pl.colorbar()
        
    def save_results(self, fname):
        np.savez(fname, u=self.ActU, x=self.ActX, hE=self.ActHE, hI=self.ActHI, 
                 rE=self.ActRE, rI=self.ActRI)
        
    def simulate_facil(self, backend='python'):
        self.x[:], self.u[:], self.hE[:], self.hI = self.x_iv, self.u_iv, self.hE_iv, self.hI_iv
        
        if backend == 'c':
            import cycover_ring as rn
            rn.set_calc_params(self.N, self.n_istep, self.dt, int(self.sampl_dt/self.dt))
            rn.set_params(self.U, self.J_IE, self.J_EI, self.tau, self.tau_d, self.tau_f, self.I0, self.alpha)
            rn.init_arrays(self.x, self.u, self.hE, self.hI, self.W, self.ActX, self.ActU, self.ActHE, self.ActHI)
            rn.integrate()
            
            self.ActRE = self.gFun(self.ActHE)
            self.ActRI = self.gFun(self.ActHI)
        elif backend == 'python':
            for i in range(0, self.n_istep):
                self.Iext[:] = 0.
                if i < self.stim_duration/self.dt:
                    dx1 = np.mod(self.pos - np.radians(self.stim_pos) + np.pi, 2*np.pi) - np.pi
    #                self.Iext[:] = self.stim_ampl*np.exp(-(dx1/self.stim_width)**2/2)
                    self.Iext[:] = self.stim_ampl*np.cos(dx1)
                
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
        print(np.mean(self.ActRE))

    @staticmethod
    def init_all_params(sim_time, dt, sampl_dt, N, J0, J1, J_EI, J_IE, 
                        eps, conn_width, conn_type, seed, U, I0, tau_d, tau_f, tau, alpha, 
                        stim_width, stim_ampl, stim_duration, stim_pos):

        self = RateNetwork(sim_time, dt, sampl_dt, N)
        self.set_weights(J0, J1, J_EI, J_IE, eps, conn_width, conn_type, seed)
        self.set_params(U, I0, tau_d, tau_f, tau, alpha)
        self.set_stim_params(stim_width, stim_ampl, stim_duration, stim_pos)
        self.set_initial_values(u=U)
        
        return self
    
    @staticmethod
    def get_hash_from_params(params_dict):
        import hashlib
        name = '_'.join([key + '_' + str(params_dict[key]) for key in RateNetwork.params_arr])
        return hashlib.md5(name.encode('utf-8')).hexdigest()
    
    @staticmethod
    def load_results(fname):
        with np.load(fname) as data:
            u = data['u']
            x = data['x']
            hE = data['hE']
            hI = data['hI']
            rE = data['rE']
            rI = data['rI']
        return u, x, hE, hI, rE, rI

if __name__ == '__main__':
    params_dict_Itskov = {
        # main params
        'sim_time': 10.,
        'dt': 0.001,
        'sampl_dt': 0.001,
        'N': 90,
    
        # connectivity params
        'J0': -10.0,
        'J1': 16.0,
        'J_EI': 0.0,
        'J_IE': 0.0,  # @TODO
        'eps': 0.0,
        'conn_width': 1.3, # @TODO
        'conn_type' : 'cos',
        'seed': 0,
    
        # actvity params
        'U': 0.05,
        'I0': 10.0,
        'tau_d': 0.1,
        'tau_f': 1.,
        'tau': 0.01,
        'alpha': 1.5,
        
        # stumulation params
        'stim_width': 0.2, # @TODO
        'stim_ampl': 0.0,
        'stim_duration': 0.05,
        'stim_pos': 0.0
    }
    
    params_dict_Tsodyks = {
        # main params
        'sim_time': 2.,
        'dt': 0.001,
        'sampl_dt': 0.001,
        'N': 90,
    
        # connectivity params
        'J0': -10.0,
        'J1': 40.0,
        'J_EI': 6.0,
        'J_IE': 1.5,
        'eps': 0.0,
        'conn_width': 1.3,
        'conn_type' : 'gauss',
        'seed': 0,
    
        # actvity params
        'U': 0.3,
        'I0': 5.0,
        'tau_d': 0.3,
        'tau_f': 1.5,
        'tau': 0.01,
        'alpha': 1.5,
        
        # stumulation params
        'stim_width': 0.2,
        'stim_ampl': 65.0,
        'stim_duration': 0.03,
        'stim_pos': 0.0 
    }

    rate_network = RateNetwork.init_all_params(**params_dict_Itskov)
    rate_network.set_initial_values(u=params_dict_Itskov['U'], hE=10*np.cos(rate_network.pos))
    
    rate_network.simulate_facil(backend = 'c')
    rate_network.plot_simul()
