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

void initArrays(double* x, double* u, double* hE, double hI,
                double* W, float* xRes, float* uRes, float* hERes, float* hIRes);

void c_integrate();

#endif /* RING_H_ */
