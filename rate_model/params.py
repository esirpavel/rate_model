# params.py

import numpy as np

simulation_params = {
    'sim_time': 10.,
    'dt': 0.001,
    'sampl_dt': 0.001,
    'N': 90
}

# dynamical noise params
noise_params = {
    'D': 0.,
    'tau_n': 0.1,
    'seed': 0
}

params_nops = {
    # connectivity params
    'J0': -10.0,
    'J1': 16.0,
    'J_EI': 0.0,
    'J_IE': 0.0,
    'eps': 0.0,
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

params_ps = {
    'J0': -1.*2*np.pi,
    'J1': 12.*2*np.pi,
    'J_EI': 1.9,
    'J_IE': 1.8*2*np.pi,
    'eps': 0.0,
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

params_tsodyks = {
    'J0': -10.0,
    'J1': 40.0,
    'J_EI': 6.0,
    'J_IE': 1.5,
    'eps': 0.0,
    'conn_width': 1.3,
    'conn_type' : 'gauss',
    'seed': 0,
    
    'U': 0.3,
    'I0': 5.0,
    'tau_d': 0.3,
    'tau_f': 1.5,
    'tau': 0.01,
    'alpha': 1.5,
}


stim_nops = {
    'stim_start': [.0],
    'stim_duration': [.1],
    'stim_ampl': [20.0],
    'stim_pos': [0.0],
    'stim_width': [1.],
    'stim_type': ['cos']
}

stim_ps = {
    'stim_start': [.0],
    'stim_duration': [.05],
    'stim_ampl': [390.0],
    'stim_pos': [0.0],
    'stim_width': [1/2.2],
    'stim_type': ['trunc_cos']
}

stim_tsodyks = {
    'stim_start': [.0],
    'stim_duration': [.05],
    'stim_ampl': [70.0],
    'stim_pos': [0.0],
    'stim_width': [.5],
    'stim_type': ['gauss']
}

