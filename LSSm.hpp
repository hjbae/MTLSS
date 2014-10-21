#ifndef LSSm_hpp
#define LSSm_hpp
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpcasa.h"
#include "petscvec.h"
#include "stdlib.h"
#include <iostream>

class lssSolver
{

    private:
        double* t;
        double* s;
        double* zz;
        double* u;
        double* dt;
        double* ds;
        double* uMidf;
        double* uMidc;
        double* dudtc;
        double* dudtf;
        double alpha;
        int* fine;
        int* coarse;
        int* mapfine;
        int* mapcoarse;
        int f;
        int c;
        int p;
        int ntdim;
        int nsdim;
        int nn;
        int rank, size;
        Mat B;
        Mat E;
        Mat wBinv;
        Mat wEinv;
        Mat Pinv;


    public:
        lssSolver(double* t, double* u0, double* zz, int p, int c, int f, int* mapcoarse, int* mapfine, int ntdim, int nsdim, double alpha);
        Vec APinvx(Vec x);
        void CGSolve(Vec b, Vec* w);
        double* lss(int maxIter=8,double atol=1e-8,double rtol=1e-7);


};
#endif
