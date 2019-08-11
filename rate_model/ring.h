/*
 * ring.h
 *
 *  Created on: May 5, 2016
 *      Author: pavel
 */

#ifndef RING_H_
#define RING_H_

void setCalcParams(unsigned N, unsigned Tsim, double h, unsigned recH);

void setParams(double U, double J_IE, double J_EI, double tau, double tau_d, double tau_f, double I0, double alpha);

void setNoiseParams(double D, double tau_n, unsigned seed);

void initArrays(double* x, double* u, double* hE, double hI,
                double* W, float* xRes, float* uRes, float* hERes, float* hIRes);

void setStimuli(unsigned* stim_start, unsigned* stim_duration, double* stim_ampl, 
                double* stim_pos, double* stim_width, unsigned* stim_type, unsigned num_stim);

void c_integrate();

#endif /* RING_H_ */
