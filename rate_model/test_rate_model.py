"""
Created on Mon Mar 11 09:10:10 2019

@author: Pavel Esir
"""
import numpy as np
import numpy.testing as npt
import rate_model as rm
import sqlite_routines
import pytest
from config import *

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
    params_dict = {}
    params_dict.update(simulation_params)
    params_dict.update(params)
    
    x, u, hE, hI = simulate(params_dict, stim_params, backend)
    x_, u_, hE_, hI_ = sqlite_routines.get_results(params_dict, stim_params)
    npt.assert_allclose(x, x_, atol=0.00001)
    npt.assert_allclose(u, u_, atol=0.00001)
    npt.assert_allclose(hE, hE_, atol=0.00001)
    npt.assert_allclose(hI, hI_, atol=0.00001)

def fill_tables():
    sqlite_routines.create_table()
    
    params_dict = {}
    params_dict.update(simulation_params)
    
    params_dict.update(params_dict_stationary)
    x, u, hE, hI = simulate(params_dict, stim_stationary)
    sqlite_routines.save_results(x, u, hE, hI, params_dict, stim_stationary)
    
    params_dict.update(params_dict_burst)
    x, u, hE, hI = simulate(params_dict, stim_burst)
    sqlite_routines.save_results(x, u, hE, hI, params_dict, stim_burst)

if __name__ == '__main__':
    fill_tables()
