'''
Created on May 6, 2016

@author: Pavel Esir
'''
cimport numpy as np
import numpy as np
import cython

cdef extern from "ring.h": 
    void setCalcParams(unsigned N, unsigned Tsim, double h, unsigned recH)

    void setParams(double U, double J_IE, double J_EI, double tau, double tau_d, double tau_f, double I0, double alpha)

    void setNoiseParams(double D, double tau_n, unsigned seed);

    void initArrays(double* x, double* u, double* hE, double hI, double* W, float* xRes, float* uRes, float* hERes, float* hIRes)

    void setStimuli(unsigned* stim_start, unsigned* stim_duration, double* stim_ampl, 
                    double* stim_pos, double* stim_width, unsigned* stim_type, unsigned num_stim)

    void c_integrate()

def set_calc_params(unsigned N, Tsim, h, recH):
    setCalcParams(N, Tsim, h, recH)

def set_params(double U, double J_IE, double J_EI, double tau, double tau_d, double tau_f, double I0, double alpha):
    setParams(U, J_IE, J_EI, tau, tau_d, tau_f, I0, alpha)

def set_noise_params(double D, double tau_n, unsigned seed):
    setNoiseParams(D, tau_n, seed)

def init_arrays(np.ndarray[np.float64_t, ndim=1] x, 
                np.ndarray[np.float64_t, ndim=1] u, 
                np.ndarray[np.float64_t, ndim=1] hE, 
                double hI, 
                np.ndarray[np.float64_t, ndim=2] W, 
                np.ndarray[np.float32_t, ndim=2] xRes, 
                np.ndarray[np.float32_t, ndim=2] uRes, 
                np.ndarray[np.float32_t, ndim=2] hERes, 
                np.ndarray[np.float32_t, ndim=1] hIRes):
    
    initArrays(<double*> x.data, 
               <double*> u.data, 
               <double*> hE.data, 
               hI,
               <double*> W.data, 
               <float*> xRes.data, 
               <float*> uRes.data,
               <float*> hERes.data,
               <float*> hIRes.data)

def set_stimuli(np.ndarray[np.uint32_t, ndim=1] stim_start, 
                np.ndarray[np.uint32_t, ndim=1] stim_duration, 
                np.ndarray[np.float64_t, ndim=1] stim_ampl, 
                np.ndarray[np.float64_t, ndim=1] stim_pos, 
                np.ndarray[np.float64_t, ndim=1] stim_width, 
                np.ndarray[np.uint32_t, ndim=1] stim_type, 
                num_stim):
    setStimuli(<unsigned*> stim_start.data, 
                <unsigned*> stim_duration.data, 
                <double*> stim_ampl.data, 
                <double*> stim_pos.data, 
                <double*> stim_width.data, 
                <unsigned*> stim_type.data, 
                num_stim)

def integrate():
    c_integrate()
