# config.py

import numpy as np

test_results_db_path = "../data/"

db_fname = 'test_results.db'
db_tname = 'test_sim_results'

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
    'dt': 0.001,
    'sampl_dt': 0.001,
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

# stimulation params for Tsodyks
stim_Tsodyks = {
    'stim_start': [.0],
    'stim_duration': [.05],
    'stim_ampl': [65.0],
    'stim_pos': [0.0],
    'stim_width': [.2],
    'stim_type': ['gauss']
}

