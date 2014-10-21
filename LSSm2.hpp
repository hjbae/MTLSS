#ifndef LSSm2_hpp
#define LSSm2_hpp
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
        int f, fg, fo;
        int c, cg, co;
        int p;
        int ntdim;
        int ntdimg;
        int nsdim;
        int nsdimg;
        int nn;
        int nng;
        int rank, size;
        int trank, tsize, srank, ssize;
        Mat B;
        Mat E;
        Mat wBinv;
        Mat wEinv;
        Mat Pinv;
        MPI_Comm comm;

    public:
        lssSolver(double* t, double* u0, int p, int c, int cg, int co, int f, int fg, int fo, int* mapcoarse, int* mapfine, int ntdim, int ntdimg, int trank, int tsize, int nsdim, int nsdimg, int srank, int ssize, double alpha, MPI_Comm comm);
        Vec APinvx(Vec x);
        void CGSolve(Vec b, Vec* w);
        double* lss(int maxIter=8,double atol=1e-7,double rtol=1e-4);


};
#endif
