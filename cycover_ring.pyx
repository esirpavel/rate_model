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

    void initArrays(double* x, double* u, double* hE, double hI, double* W, float* xRes, float* uRes, float* hERes, float* hIRes)
    
    void c_integrate()
    
def set_calc_params(unsigned N, Tsim, h, recH):
    setCalcParams(N, Tsim, h, recH)
    
def set_params(double U, double J_IE, double J_EI, double tau, double tau_d, double tau_f, double I0, double alpha):
    setParams(U, J_IE, J_EI, tau, tau_d, tau_f, I0, alpha)

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
    
def integrate():
    c_integrate()
