//
//  feedback.c
//  mbtrack_mpi
//
//  Created by Jack Borthwick on 23/07/2014.
//  Copyright (c) 2014 synchrotron soleil. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physic.h"
#include "coeffs_Sp8.h"
#include "nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}



void FIR_coeff(int n, double Qbeta, double ** po, double ** qo, double **p1, double **q1){
    unsigned int i,j,k;
    double **alpha=malloc(sizeof(double*)*n);
    double **temp = malloc(sizeof(double*)*5); // Matrice à inverser
    double **temp_r=malloc(sizeof(double*)*5);
    double psi=2*M_PI*Qbeta;
    
    //allocate memory for alpha and fill lines
    for(i=0;i<n;i++){
        alpha[i]=calloc(5,sizeof(double));
        alpha[i][0]=cos(psi*i);
        alpha[i][1]=-psi*i*sin(psi*i);
        alpha[i][2]=sin(psi*i);
        alpha[i][3]=psi*i*cos(psi*i);
        alpha[i][4]= 1;
    }
    
    for(i=0;i<5;i++){
        temp_r[i]=calloc(n,sizeof(double)); // matrice de taille 5xn
        temp[i]=calloc(5,sizeof(double)); //t_alpha*alpha matrice carrée d'ordre 5
    }
    
    //calcul de t_alpha*alpha
    for(i=0;i<5;i++)
    {
        for(j=0;j<5;j++){
            temp[i][j]=0.0; // Verify that the temp table is correctly initialized
            for(k=0;k<n;k++){
                temp[i][j]+=alpha[k][i]*alpha[k][j];
            }
        }
    }
    
    gaussj(temp,5); // compute (t_alpha*alpha)^{-1}
    
    //calcul de (t_alpha*alpha)^{-1} * t_alpha
    for(i=0;i<5;i++)
    {
        for(j=0;j<n;j++){
            temp_r[i][j]=0.0;
            for(k=0;k<5;k++){
                temp_r[i][j]+=temp[i][k]*alpha[j][k];
            }
        }
    }
    
    (*po)=malloc(sizeof(double)*n);
    (*qo)=malloc(sizeof(double)*n);
    (*q1)=malloc(sizeof(double)*n);
    (*p1)=malloc(sizeof(double)*n);
    for(i=0;i<n;i++){
        (*po)[i]=temp_r[0][i]; // FIR coefficients needed for computation of P_0
        (*qo)[i]=temp_r[2][i]; // FIR coefficients needed for computation of Q_0
        (*p1)[i]=temp_r[1][i];
        (*q1)[i]=temp_r[3][i];
    }
    
    
}

//Gauss-Jacobi inversion: Source Numerical recipes
void gaussj(double **a, int n)
{
    int *indxc,*indxr,*ipiv;
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv,temp;
    
    indxc=ivector(1,n);
    indxr=ivector(1,n);
    ipiv=ivector(1,n);
    for (j=0;j<n;j++) ipiv[j+1]=0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv[j+1] != 1)
                for (k=0;k<n;k++) {
                    if (ipiv[k+1] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big=fabs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    } else if (ipiv[k+1] > 1)
                        nrerror("gaussj: Singular Matrix-1");
                }
        ++(ipiv[icol+1]);
        if (irow != icol) {
            for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
                
                }
        indxr[i+1]=irow;
        indxc[i+1]=icol;
        if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
        pivinv=1.0/a[icol][icol];
        a[icol][icol]=1.0;
        for (l=0;l<n;l++) a[icol][l] *= pivinv;
        for (ll=0;ll<n;ll++)
            if (ll != icol) {
                dum=a[ll][icol];
                a[ll][icol]=0.0;
                for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
            }
    }
    for (l=n;l>=1;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<n;k++)
                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free_ivector(ipiv,1,n);
    free_ivector(indxr,1,n);
    free_ivector(indxc,1,n);
}



