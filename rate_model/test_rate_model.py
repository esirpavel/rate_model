"""
Created on Mon Mar 11 09:10:10 2019

@author: Pavel Esir
"""
import numpy as np
import numpy.testing as npt
import rate_model as rm
import sqlite_routines
import pytest
from params import *

def simulate(params, stim_params, backend='python'):
    rate_network = rm.RateNetwork.init_all_params(**params)
    
    rate_network.set_stimuli(**stim_params)
    rate_network.set_initial_values()
    
    rate_network.simulate_facil(backend=backend)
    return rate_network.ActX, rate_network.ActU, rate_network.ActHE, rate_network.ActHI

@pytest.mark.parametrize('backend', ['c', 'python'])
@pytest.mark.parametrize('params, stim_params', 
                          [(params_nops, stim_nops), 
                           (params_ps, stim_ps)])
def test_sim(params, stim_params, backend):
    params_dict = {}
    params_dict.update(simulation_params)
    params_dict.update(params)
    params_dict.update({'eps': 0.5})
    
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
    
    params_dict.update(params_nops)
    params_dict.update({'eps': 0.5})
    x, u, hE, hI = simulate(params_dict, stim_nops)
    sqlite_routines.save_results(x, u, hE, hI, params_dict, stim_nops)
    
    params_dict.update(params_ps)
    params_dict.update({'eps': 0.5})
    x, u, hE, hI = simulate(params_dict, stim_ps)
    sqlite_routines.save_results(x, u, hE, hI, params_dict, stim_ps)

if __name__ == '__main__':
    fill_tables()
