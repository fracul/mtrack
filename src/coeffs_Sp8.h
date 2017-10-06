//
//  feedback.h
//  mbtrack_mpi
//
//  Created by Jack Borthwick on 23/07/2014.
//  Copyright (c) 2014 synchrotron soleil. All rights reserved.
//

#ifndef mbtrack_mpi_feedback_h
#define mbtrack_mpi_feedback_h
#include <complex.h>

void FIR_coeff(int n, double Qbeta, double ** po, double ** qo, double ** p1, double ** q1);

//Algorithms for matrix inversion: source Numerical Recipes
void gaussj(double **a, int n);

#endif
