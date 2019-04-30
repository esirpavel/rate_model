"""
Created on Mon Mar 11 09:10:10 2019

@author: Pavel Esir
"""
import numpy as np
import numpy.testing as npt
import rate_model as rm
import sqlite_routines
import pytest

params_dict_stationary = {
    # main params
    'sim_time': 10.,
    'dt': 0.001,
    'sampl_dt': 0.01,
    'N': 90,

    # connectivity params
    'J0': -10.0,
    'J1': 16.0,
    'J_EI': 0.0,
    'J_IE': 0.0,
    'eps': 0.5,
    'conn_width': 1.0,
    'conn_type' : 'cos',
    'seed': 0,

    # actvity params
    'U': 0.05,
    'I0': 10.0,
    'tau_d': 0.1,
    'tau_f': 1.,
    'tau': 0.01,
    'alpha': 1.5,
}
   
params_dict_burst = {
    # main params
    'sim_time': 10.,
    'dt': 0.0005,
    'sampl_dt': 0.0005,
    'N': 90,
    
    'J0': -1.*2*np.pi,
    'J1': 12.*2*np.pi,
    'J_EI': 1.9,
    'J_IE': 1.8*2*np.pi,
    'eps': 0.5,
    'conn_width': 1/2.2,
    'conn_type' : 'trunc_cos',
    'seed': 0,
    
    'U': 0.3,
    'I0': -0.1,
    'tau_d': 0.3,
    'tau_f': 4.0,
    'tau': 0.01,
    'alpha': 1.5,
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

def simulate(params, stim_params, backend='python'):
    rate_network = rm.RateNetwork.init_all_params(**params)
    
    rate_network.set_stimuli(**stim_params)
    rate_network.set_initial_values()
    
    rate_network.simulate_facil(backend=backend)
    return rate_network.ActX, rate_network.ActU, rate_network.ActHE, rate_network.ActHI

@pytest.mark.parametrize('backend', ['c', 'python'])
@pytest.mark.parametrize('params, stim_params', 
                          [(params_dict_stationary, stim_stationary), 
                           (params_dict_burst, stim_burst)])
def test_sim(params, stim_params, backend):
    x, u, hE, hI = simulate(params, stim_params, backend)
    x_, u_, hE_, hI_ = sqlite_routines.get_results(params, stim_params)
    npt.assert_allclose(x, x_, atol=0.00001)
    npt.assert_allclose(u, u_, atol=0.00001)
    npt.assert_allclose(hE, hE_, atol=0.00001)
    npt.assert_allclose(hI, hI_, atol=0.00001)

def fill_tables():
    x, u, hE, hI = simulate(params_dict_stationary, stim_stationary)
    sqlite_routines.create_table()
    sqlite_routines.save_results(x, u, hE, hI, params_dict_stationary, stim_stationary)
    
    x, u, hE, hI = simulate(params_dict_burst, stim_burst)
    sqlite_routines.save_results(x, u, hE, hI, params_dict_burst, stim_burst)
    
if __name__ == '__main__':
    fill_tables()
