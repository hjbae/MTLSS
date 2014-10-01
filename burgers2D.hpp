#ifndef burgers2d_hpp
#define burgers2d_hpp

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>


double* shiftupx(double *y, int nsdim)
{

	double* yupx = new double[nsdim]; 
	for (int i = 0; i < 7; i++)
	{
		for (int j = 1; j < 21; j++)
		{
			yupx[i*21+j-1] = y[i*21+j];
		}
		yupx[i*21+20] = y[i*21];
	}

	for (int i = 0; i < 25; i++)
	{
		for (int j = 1; j < 81; j++)
		{
			yupx[7*21+i*81+j-1] = y[7*21+i*81+j];
		}
		yupx[7*21+i*81+80] = y[7*21+i*81];
	}

	for (int i = 0; i < 7; i ++)
	{
		for (int j = 1; j < 21; j++)
		{
			yupx[7*21+25*81+i*21+j-1] = y[7*21+25*81+i*21+j];
		}
		 yupx[7*21+25*81+i*21+20] = y[7*21+25*81+i*21];
	}

	return yupx;
}


double* shiftdownx(double *y, int nsdim)
{
	double* yupx = new double[nsdim];
        for (int i = 0; i < 7; i++)
        {
                for (int j = 0; j < 20; j++)
                {
                        yupx[i*21+j+1] = y[i*21+j];
                }
                yupx[i*21] = y[i*21+20];
        }

        for (int i = 0; i < 25; i++)
        {
                for (int j = 0; j < 80; j++)
                {
                        yupx[7*21+i*81+j+1] = y[7*21+i*81+j];
                }
                yupx[7*21+i*81] = y[7*21+i*81+80];
        }       

        for (int i = 0; i < 7; i ++)
        {
                for (int j = 0; j < 20; j++)
                {
                        yupx[7*21+25*81+i*21+j+1] = y[7*21+25*81+i*21+j];
                }
                 yupx[7*21+25*81+i*21] = y[7*21+25*81+i*21+20];
        }

        return yupx;

}



double* shiftupy(double *y, int nsdim)
{

	double* yy = new double[nsdim];
	for (int i = 0; i < 6*21; i++)
	{
		yy[i] = y[21+i];
	}
	for (int i = 0; i < 21; i++)
	{
		yy[6*21+i] = y[7*21+i*4];
	}
	for (int i = 0; i < 81*24; i++)
	{
		yy[7*21+i] = y[7*21+81+i];
	}
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			yy[7*21+24*81+i*4+j] = (1-0.25*j)*y[7*21+25*81+i]+(0.25*j)*y[7*21+25*81+i+1];
		}
	}
	yy[7*21+24*81+80] = y[7*21+25*81+20];
	for (int i = 0; i < 6*21; i++)
	{
		yy[7*21+25*81+i] = y[7*21+25*81+21+i];
	}
	for (int i = 0; i < 21; i++)
	{
		yy[13*21+25*81+i] = y[i];
	}

	return yy;
}




double* shiftdowny(double *y, int nsdim)
{

        double* yy = new double[nsdim];
	for (int i = 0; i <21 ; i++)
	{
		yy[i] = y[13*21+25*81+i];
	}
        for (int i = 0; i < 6*21; i++)
        {
                yy[21+i] = y[i];
        }
        for (int i = 0; i < 20; i++)
        {
                for (int j = 0; j < 4; j++)
                {
                        yy[7*21+i*4+j] = (1-0.25*j)*y[6*21+i]+(0.25*j)*y[6*21+i+1];
                }
        }
        yy[7*21+80] = y[6*21+20];
        for (int i = 0; i < 81*24; i++)
        {
                yy[7*21+81+i] = y[7*21+i];
        }
	for (int i = 0; i < 21; i++)
	{
		yy[7*21+81*25+i] = y[7*21+24*81+i*4];
	}
        for (int i = 0; i < 6*21; i++)
        {
                yy[7*21+25*81+21+i] = y[7*21+81*25+i];
        }

        return yy;
}


double* F(double* uMid, double p, int ntdim, int nsdim)
{
	double dx[nsdim];
        for (int i = 0; i < 7*21; i++)
	{
		dx[i] = 0.4;
	}
	for (int i = 7*21; i < 7*21+81*25; i++)
	{
		dx[i] = 0.1;
	}
	for (int i = 7*21+81*25; i < 14*21+81*25; i++)
	{
		dx[i] = 0.4;
	}


	double c = 0.025;
	double* dudt = new double[(ntdim-1)*nsdim];
	for (int t = 0; t < ntdim-1; t++)
	{
		double* bxup = shiftupx(&uMid[t*nsdim],nsdim);
		double* bxdown = shiftdownx(&uMid[t*nsdim],nsdim);
		double* byup = shiftupy(&uMid[t*nsdim],nsdim);
		double* bydown = shiftdowny(&uMid[t*nsdim],nsdim);
		double* dxup = shiftupx(&dx[0],nsdim);
		double* dxdown = shiftdownx(&dx[0], nsdim);
		double* dyup = shiftupy(&dx[0],nsdim);
		double* dydown = shiftdowny(&dx[0],nsdim);

		for (int i = 0; i < nsdim; i++)
		{
			dudt[t*nsdim+i] = -(bxup[i]*bxup[i]-bxdown[i]*bxdown[i])/(dxup[i]+dxdown[i]) + c*((bxup[i] - 2*uMid[t*nsdim+i] + bxdown[i])/(dxdown[i]*dxup[i]) + (byup[i] - 2*uMid[t*nsdim+i] + bydown[i])/(dyup[i]*dydown[i]));
		}

		delete bxup;
		delete bxdown;
		delete byup;
		delete bydown;
		delete dxup;
		delete dxdown;
		delete dyup;
		delete dydown;

	}
	return dudt;
}




void dfdu(double* uMid, double* s,int nsdim, int t, double*& val, int*& ia, int*& ja, int& nnz)
{
    double* L = new double[nsdim * nsdim];
    double eps = 1e-5;
    for (unsigned int i = 0; i<nsdim; i++)
    {
        uMid[t*nsdim+i] -= eps;
        double* Fm = new double[nsdim];
        double* Fmm = F(&uMid[t*nsdim],1,2,nsdim);
        memcpy (Fm, Fmm, nsdim*sizeof(double));
        uMid[t*nsdim+i] += eps*2;

        double* Fp = new double[nsdim];
        double* Fpp = F(&uMid[t*nsdim],1,2,nsdim);
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
