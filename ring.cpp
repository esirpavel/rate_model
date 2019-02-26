/*
 * ring.cpp
 *
 *  Created on: May 5, 2016
 *      Author: pavel
 */
#include <cmath>
#include "ring.h"
#include <cblas.h>
#include <iostream>

namespace ring{
    double* wArr; // 2D matrix of weights stored in row major order
    
    double* xArr;
    double* uArr;
    double* hEArr;
    double hI;
    
    float* xRes;
    float* uRes;
    float* hERes; 
    float* hIRes;
    
    double* rEArr;
    double rI;
    double* mRight; // multiplication of mArr to xArr
    double* mxOld;  // multiplication of mArr to xArr
    double* Iext;

    unsigned stim_idx = 0;
    bool is_stimulating = false;
    
    unsigned* stim_start;
    unsigned* stim_duration;
    double* stim_ampl;
    double* stim_pos;
    double* stim_width;
    unsigned* stim_type;
    unsigned num_stim;
    
    unsigned N;
    double tau_d;
    double tau_f;
    double tau;
    
    double U;
    double I0;
    double J_EI;
    double J_IE;

    double alpha = 1.5;

    unsigned Tsim;
    double h;
    unsigned recH;
}

void setCalcParams(unsigned N, unsigned Tsim, double h, unsigned recH){
    ring::N = N;
    ring::Tsim = Tsim;
    ring::h = h;
    ring::recH = recH;
}

void setParams(double U, double J_IE, double J_EI, 
               double tau, double tau_d, double tau_f, 
               double I0, double alpha){
    ring::U = U;
    ring::tau = tau;
    ring::tau_d = tau_d;
    ring::tau_f = tau_f;
    ring::I0 = I0;
    ring::alpha = alpha;
    ring::J_IE = J_IE;
    ring::J_EI = J_EI;
}

void initArrays(double* x, double* u, double* hE, double hI,
                double* W, float* xRes, float* uRes, float* hERes, float* hIRes){
    ring::xArr = x;
    ring::uArr = u;
    ring::hEArr = hE;
    ring::hI = hI;
    
    ring::wArr = W;
    ring::xRes = xRes;
    ring::uRes = uRes;
    ring::hERes = hERes;
    ring::hIRes = hIRes;
    
    ring::rEArr = new double[ring::N];
    ring::mxOld = new double[ring::N];
    ring::mRight = new double[ring::N];
}

void setStimuli(unsigned* stim_start, unsigned* stim_duration, double* stim_ampl, 
                double* stim_pos, double* stim_width, unsigned* stim_type, unsigned num_stim){
    ring::Iext = new double[ring::N]();
    
    ring::stim_idx = 0;
    ring::is_stimulating = false;
    
    ring::stim_start = stim_start;
    ring::stim_duration = stim_duration;
    ring::stim_ampl = stim_ampl;
    ring::stim_pos = stim_pos;
    ring::stim_width = stim_width;
    ring::stim_type = stim_type;
    ring::num_stim = num_stim;
}

using namespace ring;

double stimFunc(double x, double width, unsigned type){
    switch (type){
        case 1:
            return cos(x/width);
        case 2:
            return exp(-x*x/(width*width*2.0));
        case 3:
            // @TODO here should be truncated cosyne
            break;
        case 4: 
            // @TODO uniform mask
            break;
    } 
}

void process_stumulus(unsigned t){
    if (stim_idx < num_stim) {
        if (!is_stimulating && (t == stim_start[stim_idx])) {
            for (unsigned n = 0; n < N; n++){
                double stim_phase_pos = fmod(2.0*M_PI*n/N - stim_pos[stim_idx], 2.0*M_PI)  - M_PI;
                Iext[n] = stim_ampl[stim_idx]*stimFunc(stim_phase_pos, stim_width[stim_idx], stim_type[stim_idx]);
            }
            is_stimulating = true;
        } else if ( is_stimulating && (t == stim_start[stim_idx] + stim_duration[stim_idx]) ) {
            for (unsigned n = 0; n < N; n++){
                Iext[n] = 0.0;
            }
            is_stimulating = false;
            stim_idx++;
        }
    }
}

double gFun(double x){
    if (x < 10.0){
        return alpha*log1p(exp(x/alpha));
    } else {
        return x;
    }
}

void c_integrate(){
    double I_inh_exc; // global excitatory current to inhibitory pool
    
    for (unsigned t = 0; t < Tsim; t++){
        I_inh_exc = 0.0;
        for (unsigned i = 0; i < N; i++){
            rEArr[i] = gFun(hEArr[i]);
            I_inh_exc += rEArr[i];
            mxOld[i] = rEArr[i]*xArr[i]*uArr[i];
        }
        rI = gFun(hI);
        hI += (-hI + J_IE*I_inh_exc/N)*h/tau;

        cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, wArr, N, mxOld, 1, 0.0, mRight, 1);

        process_stumulus(t);
        
        for (unsigned i = 0; i < N; i++){
            hEArr[i] += (-hEArr[i] + mRight[i] + I0 + Iext[i] - J_EI*rI)*h/tau;
            
            if (tau_d == 0.0){
                xArr[i] = 1.0;
            } else{
                xArr[i] += ((1.0 - xArr[i])/tau_d - mxOld[i])*h;
            }
            if (tau_f == 0.0){
                uArr[i] = U;
            } else{
                uArr[i] += ((U - uArr[i])/tau_f + U*rEArr[i]*(1.0 - uArr[i]))*h;
            }
            
            if ((t % recH) == 0){
                xRes[N*t/recH + i] = xArr[i];
                uRes[N*t/recH + i] = uArr[i];
                hERes[N*t/recH + i] = hEArr[i];

                hIRes[t/recH] = hI;
            }
       }
    }
}
