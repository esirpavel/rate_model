"""
Created on Thu Feb 28 19:39:36 2019

@author: Pavel Esir
"""
import numpy as np
import matplotlib.pylab as pl
import rate_model as rm

if __name__ == '__main__':
    params_dict_Itskov = {
        # main params
        'sim_time': 200.,
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
    
    # Parameters from Yuanyuan
    params_dict_Yuanyuan = {
        # main params
        'sim_time': 20.,
        'dt': 0.001,
        'sampl_dt': 0.01,
        'N': 90,
        
        'J0': -1.*2*np.pi,
        'J1': 12.*2*np.pi,
        'J_EI': 1.9,
        'J_IE': 1.8*2*np.pi,
        'eps': 2.,
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

    
    # stimulating params for Itskov
    stim_Itskov = {
        'stim_start': [.0],
        'stim_duration': [.05],
        'stim_ampl': [5.0],
        'stim_pos': [0.0],
        'stim_width': [1.],
        'stim_type': ['cos']
    }

    # stimulating params for Tsodyks
    stim_Tsodyks = {
        'stim_start': [.0],
        'stim_duration': [.05],
        'stim_ampl': [65.0],
        'stim_pos': [0.0],
        'stim_width': [.2],
        'stim_type': ['gauss']
    }

    # stimulating params for Yuanyuan
    stim_Yuanyuan = {
        'stim_start': [.0],
        'stim_duration': [.05],
        'stim_ampl': [390.0],
        'stim_pos': [30.0],
        'stim_width': [1/2.2],
        'stim_type': ['trunc_cos']
    }
    
    #%%
    rate_network = rm.RateNetwork.init_all_params(**params_dict_Yuanyuan)
    rate_network.set_initial_values(hE=0*np.cos(rate_network.pos))
    rate_network.set_stimuli(**stim_Yuanyuan)
    
    rate_network.simulate_facil(backend = 'c')
    rate_network.plot_simul()
    
    #%%
    stim_dict = stim_Yuanyuan
    param_dict = params_dict_Yuanyuan
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
    # rate_network.plot_simul()
