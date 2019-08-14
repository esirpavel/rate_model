"""
Created on Thu Feb 28 19:39:36 2019

@author: Pavel Esir
"""
import numpy as np
import matplotlib.pylab as pl
import rate_model as rm
from rate_model.config import *

sim_params = {
    'sim_time': 150.,
    'dt': 0.001,
    'sampl_dt': 0.001,
    'N': 90,
}

stim_stationary = {
    'stim_start': [.0],
    'stim_duration': [.05],
    'stim_ampl': [5.0],
    'stim_pos': [0.0],
    'stim_width': [1.],
    'stim_type': ['cos']
}

stim_burst = {
    'stim_start': [.0],
    'stim_duration': [.05],
    'stim_ampl': [390.0],
    'stim_pos': [0.0],
    'stim_width': [1/2.2],
    'stim_type': ['trunc_cos']
}

if __name__ == '__main__':
    # stim_dict = stim_stationary
    stim_dict = stim_burst
    param_dict = params_dict_burst.copy()
    param_dict.update(sim_params)
    param_dict.update({'eps': 0.5})
    
    rate_network = rm.RateNetwork.init_all_params(**param_dict)
    for st_pos in range(-180, 180, 30):
        stim_dict['stim_pos'] = [st_pos]
        rate_network.set_initial_values(hE=0*np.cos(rate_network.pos))
        rate_network.set_stimuli(**stim_dict)
        
        rate_network.simulate_facil(backend = 'c')
        pl.figure(2, figsize=(10, 7))
        pl.title('eps = {}'.format(param_dict['eps']))
        pl.plot(rate_network.tm[2:], np.degrees(rate_network.get_angle(rate_network.ActU))[2:], lw=3.)
        pl.xlabel("Time (s)")
        pl.ylabel(r'$\theta\, ({deg}^\circ)$')
    
    pl.show()

