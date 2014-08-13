#ifndef sburger_hpp
#define sburger_hpp

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>

double* F(double* uMid, double p, int ntdim, int nsdim, double* w)
{

	double dx[nsdim+1];
	for (int i = 0; i < 42; i++)
	{
		dx[i] = 1./1024;
	}
	for (int i = 42; i < 512; i ++)
	{
		dx[i] = 1./512;
	}
	for (int i = 512; i < 554; i++)
	{
		dx[i] = 1./1024;
	}

	double Re = 1500;
	double* dudt = new double[(ntdim-1)*nsdim];

	for (int t = 0; t < ntdim-1; t++)
	{
		dudt[nsdim*t] = -(uMid[nsdim*t+1]* uMid[nsdim*t+1])/(2*(dx[1]+dx[0])) + 1/Re* (uMid[nsdim*t+1] - uMid[nsdim*t]*2 + 0) / (dx[1] * dx[0]) + w[nsdim*t];
		dudt[nsdim*(t+1)-1] = (uMid[nsdim*(t+1)-2]*uMid[nsdim*(t+1)-2])/(2*(dx[nsdim]+dx[nsdim-1])) + 1/Re*(0-uMid[nsdim*(t+1)-1]*2 + uMid[nsdim*(t+1)-2]) / (dx[nsdim] * dx[nsdim-1]) + w[nsdim*(t+1)-1];
		for (int i = 1; i < nsdim-1; i++)
		{
			dudt[nsdim*t+i] = -(uMid[nsdim*t+i+1]* uMid[nsdim*t+i+1] - uMid[nsdim*t+i-1]*uMid[nsdim*t+i-1])/(2*(dx[i+1]+dx[i])) + 1/Re * (uMid[nsdim*t+i+1]-uMid[nsdim*t+i]*2+uMid[nsdim*t+i-1])/(dx[i+1]*dx[i])+w[nsdim*t+i];
		}
	}



	return dudt;
}




void dfdu(double* uMid, double* w,int nsdim, int t, double*& val, int*& ia, int*& ja, int& nnz)
{
    double* L = new double[nsdim * nsdim];
    double eps = 1e-5;
    for (unsigned int i = 0; i<nsdim; i++)
    {
        uMid[t*nsdim+i] -= eps;
        double* Fm = new double[nsdim];
        double* Fmm = F(&uMid[t*nsdim],1,2,nsdim, &w[t*nsdim]);
        memcpy (Fm, Fmm, nsdim*sizeof(double));
        uMid[t*nsdim+i] += eps*2;

        double* Fp = new double[nsdim];
        double* Fpp = F(&uMid[t*nsdim],1,2,nsdim, &w[t*nsdim]);
        memcpy (Fp, Fpp, nsdim*sizeof(double));
        uMid[t*nsdim+i] -= eps;

        for (unsigned int j=0; j<nsdim; j++)
        {
            L[j*nsdim+i] = -0.5*(Fp[j] - Fm[j])/(2*eps);
        }

        delete Fp;
        delete Fm;
        delete Fpp;
        delete Fmm;

        }

    nnz=0;
    for (int i = 0; i < nsdim*nsdim; i++)
    {
        if (L[i] != 0)
        {
            nnz++;
        }
    }

    if (val != 0) delete val;
    val = new double[nnz];
    if (ia != 0) delete ia;
    ia = new int[nnz];
    if (ja != 0) delete ja;
    ja = new int[nnz];
    nnz = 0;
    for (int i = 0; i < nsdim; i++)
    {
        for (int j = 0; j < nsdim; j++)
        {
            if (L[i*nsdim+j] != 0)
            {
                val[nnz] = L[i*nsdim+j];
                ia[nnz] = i;
                ja[nnz] = j;
                nnz ++ ;
            }
        }
    }
    delete L;

}



#endif
