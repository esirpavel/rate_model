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

using namespace ring;

double gFun(double x){
    // if (x < 0.0)
    //     return 0.0;
    // else
    //     return x;
    return alpha*log(1.0 + exp(x/alpha));
}

void c_integrate(){
    double I_inh_exc; // global excitatory current to inhibitory pool
    
    for (unsigned t = 0; t < Tsim; t++){
        for (unsigned i = 0; i < N; i++){
            rEArr[i] = gFun(hEArr[i]);
            I_inh_exc += J_IE*rEArr[i];
            mxOld[i] = rEArr[i]*xArr[i]*uArr[i];
        }
        rI = gFun(hI);
        hI = (-hI + I_inh_exc/N)*h/tau;

        cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, wArr, N, mxOld, 1, 0.0, mRight, 1);

        for (unsigned i = 0; i < N; i++){
            hEArr[i] += (-hEArr[i] + mRight[i] + I0 - J_EI*rI)*h/tau;
            
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
