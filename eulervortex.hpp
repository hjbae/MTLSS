#ifndef eulervortex_hpp
#define eulervortex_hpp

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
double* shiftup(double* y, int nsdim)
{

    double* z = new double[nsdim/4];

    for (int i = 0; i<21*11; i++)
    {
        z[i] = y[i+21];
    }
    for (int i = 21*11; i<21*11+13; i++)
    {
        z[i] = y[i+21];
    }
    for (int i = 0; i < 6; i++)
    {
        z[21*11+13+i] = y[12*21+13+2*i+1];
    }
    z[21*11+19] = y[21*12+25];
    z[21*11+20] = y[21*12+26];


    for (int i = 0; i<6; i++)
    {
        for (int j = 0; j < 12; j++)
        {
            z[12*21+i*40+j] = y[12*21+(i+1)*40+j];
        }
        for (int j = 12; j<25; j++)
        {
            z[12*21+i*40+j] = y[12*21+i*40+j+15];
        }
        for (int j = 25; j<27; j++)
        {
                        z[12*21+i*40+j] = y[12*21+(i+1)*40+j];
                }
        for (int j = 27; j<40; j++)
        {
            z[12*21+i*40+j] = y[12*21+(i+1)*40+j-15];
        }
    }

    for (int i = 0; i<13; i++)
    {
        z[12*21+6*40+i] = y[12*21+6*40+27+i];
    }
    for (int i = 0; i<6; i++)
    {
        z[12*21+6*40+12+2*i+1] = (y[12*21+6*40+27+12+i] + y[12*21+6*40+27+12+i+1])/2;
        z[12*21+6*40+13+2*i+1] = y[12*21+6*40+27+12+i+1];
    }
    for (int i = 25; i<27; i++)
    {
                z[12*21+6*40+i] = y[12*21+6*40+21+i];
        }

    for (int i = 12*21+6*40+27; i < 13*21+6*40+27; i++)
    {
        z[i] = y[i+21];
    }
    for (int i = 13*21+6*40+27; i < 14*21+6*40+27; i++)
    {
        z[i] = y[i-(13*21+6*40+27)];
    }

    return z;
}


double* shiftdown(double* y, int nsdim)
{

    double* z = new double[nsdim/4];
        for (int i = 0; i<21; i++)
    {
        z[i] = y[nsdim/4-21+i];
    }
    for (int i = 21; i<21*12; i++)
        {
                z[i] = y[i-21];
        }

    for (int i = 0; i<13; i++)
        {
                z[12*21+i] = y[11*21+i];
        }
        for (int i = 0; i<6; i++)
        {
                z[12*21+12+2*i+1] = (y[11*21+12+i] + y[11*21+12+i+1])/2;
                z[12*21+13+2*i+1] = y[11*21+12+i+1];
        }
        for (int i = 25; i<27; i++)
        {
                z[12*21+i] = y[11*21+i-6];
        }


        for (int i = 0; i<6; i++)
        {
                for (int j = 0; j < 12; j++)
                {
                        z[12*21+(i+1)*40+j] = y[12*21+i*40+j];
                }
                for (int j = 12; j<25; j++)
                {
                        z[12*21+(i+1)*40+j] = y[12*21+i*40+j+15];
                }
                for (int j = 25; j<27; j++)
                {
                        z[12*21+(i+1)*40+j] = y[12*21+i*40+j];
                }
                for (int j = 27; j<40; j++)
                {
                        z[12*21+i*40+j] = y[12*21+i*40+j-15];
                }
        }


    for (int i = 21*12+40*6+27; i<21*12+40*6+27+12; i++)
        {
                z[i] = y[i-27];
        }
        for (int i = 0; i < 7; i++)
        {
                z[21*12+40*6+27+12+i] = y[12*21+40*6+12+2*i];
        }
        z[21*12+6*40+27+19] = y[21*12+6*40+25];
        z[21*12+6*40+27+20] = y[21*12+6*40+26];

    for (int i = 21*12+6*40+27+21; i < 21*12+6*40+27+2*21; i++)
    {
        z[i] = y[i-21];
    }


    return z;

}


double* changeaxis(double* y, int nsdim)
{

    double* z = new double[nsdim/4];

    for (int i = 0; i < 12; i++)
    {
        for (int j = 0; j < 12; j++)
        {
            z[i*21+j] = y[j*21+i];
        }
        for (int j = 12; j<19; j++)
        {
            z[i*21+j] = y[12*21+(j-12)*40+i];
        }
        for (int j = 19; j <21; j++)
        {
            z[i*21+j] = y[12*21+6*40+27+(j-19)*21+i];
        }
    }

    for (int i=12; i< 19; i++)
    {
        for (int j = 0; j<12; j++)
        {
            z[12*21+(i-12)*40+j] = y[j*21 + i];
        }
        for (int j = 12; j<25; j++)
        {
            if (j % 2 == 0)    
            {
                z[12*21+(i-12)*40+j] = y[12*21 + (j-12)*20 + (i-12)*2 + 12];
            }
            else
            {
                z[12*21+(i-12)*40+j] = y[12*21 + (j-13)*20 + (i-12)*2 + 27];
            }
        }
        for (int j = 25; j<27; j++)
        {
            z[12*21+(i-12)*40+j] = y[12*21+6*40+27+(j-25)*21+i];
        }
        for (int j = 27; j < 40; j++)
        {
            if (j % 2 == 1)
            {
                z[12*21+(i-12)*40+j] = y[12*21 + (j-27)*20+(i-12)*2 + 13]; 
            }
            else
            {
                z[12*21+(i-12)*40+j] = y[12*21 + (j-28)*20 + 28 + (i-12)*2];
            }
        }
    }

    for (int i = 19; i<21; i++)
    {
        for (int j = 0; j < 12; j++)
        {
            z[12*21+6*40+27+(i-19)*21+j] = y[j*21+i];
        }
        for (int j = 12; j < 19; j++)
        {
            z[12*21+6*40+27+(i-19)*21+j] = y[12*21+(j-12)*40+i+6];
        }
        for (int j = 19; j < 21; j++)
        {
            z[12*21+6*40+27+(i-19)*21+j] = y[12*21+6*40+27+(j-19)*21+i];
        }
    }
    return z;
}



double* shiftupx(double *u, int nsdim)
{
    double *temp1, *temp2, *output;
    temp1 = changeaxis(u,nsdim);
    temp2 = shiftup(temp1, nsdim);
    delete temp1;
    output = changeaxis(temp2,nsdim);
    delete temp2;

    return output;
}



double* shiftdownx(double *u, int nsdim)
{
    double *temp1, *temp2, *output;
    temp1 = changeaxis(u,nsdim);
    temp2 = shiftdown(temp1, nsdim);
    delete temp1;
    output = changeaxis(temp2,nsdim);
    delete temp2;

    return output;
}








double* F(double* uMid, double p, int ntdim, int nsdim)
{

    double dx[nsdim/4];
    for (int i = 0; i < nsdim/4; i++)
    {
        dx[i] = 0.5;
    }
    for (int i = 0; i < 6; i++)
    {
        for (int j = 12; j < 25; j++)
        {
            dx[12*21+i*40+j] = 0.25;
        }
        for (int j = 27; j < 40; j++)
        {
            dx[12*21+i*40+j] = 0.25;
        }
    }
    for (int j = 12; j<25; j++)
    {
        dx[12*21+6*40+j] = 0.25;
    }


    double gamma = 1.67;
    double* dudt = new double[(int)(nsdim*floor((ntdim-1)/p))];

    for (int t=0; t<(ntdim-1); t+=p)
    {
        double* drhoudx1 = shiftupx(&uMid[t*nsdim+nsdim/4],nsdim);
        double* drhoudx2 = shiftdownx(&uMid[t*nsdim+nsdim/4],nsdim);
        double* drhovdy1 = shiftup(&uMid[t*nsdim+nsdim/2],nsdim);
        double* drhovdy2 = shiftdown(&uMid[t*nsdim+nsdim/2],nsdim);
        for (int i = 0; i<nsdim/4; i++)
        {
            dudt[t*nsdim+i] = -(drhoudx1[i]-drhoudx2[i])/(dx[i]*2) - (drhovdy1[i]-drhovdy2[i])/(dx[i]*2);
        }

        delete drhoudx1;
        delete drhoudx2;
        delete drhovdy1;
        delete drhovdy2;
        
	double* rhoSqupp = new double[nsdim/4];
        double* rhouv = new double[nsdim/4];
        for (int i = 0; i<nsdim/4; i++)
        {
            rhoSqupp[i] = uMid[t*nsdim+nsdim/4+i]*uMid[t*nsdim+nsdim/4+i]/uMid[t*nsdim+i]+pow(uMid[t*nsdim+i],gamma)/gamma;
            rhouv[i] = uMid[t*nsdim+nsdim/4+i]*uMid[t*nsdim+nsdim/2+i]/uMid[t*nsdim+i];
        }
        double* drhoSquppdx1 = shiftupx(rhoSqupp,nsdim);
        double* drhoSquppdx2 = shiftdownx(rhoSqupp,nsdim);
        double* drhouvdy1 = shiftup(rhouv,nsdim);
        double* drhouvdy2 = shiftdown(rhouv,nsdim);

        for (int i = 0; i<nsdim/4; i++)
        {
            dudt[t*nsdim+nsdim/4+i] = -(drhoSquppdx1[i]-drhoSquppdx2[i])/(dx[i]*2) - (drhouvdy1[i]-drhouvdy2[i])/(dx[i]*2);
        }

        delete rhoSqupp;
        delete drhoSquppdx1;
        delete drhoSquppdx2;
        delete drhouvdy1;
        delete drhouvdy2;

        double* rhoSqvpp = new double[nsdim/4];
        for (int i = 0; i<nsdim/4; i++)
        {
            rhoSqvpp[i] = uMid[t*nsdim+nsdim/2+i]*uMid[t*nsdim+nsdim/2+i]/uMid[t*nsdim+i]+pow(uMid[t*nsdim+i],gamma)/gamma;
        }
        double* drhouvdx1 = shiftupx(rhouv,nsdim);
        double* drhouvdx2 = shiftdownx(rhouv,nsdim);
        double* drhoSqvppdy1 = shiftup(rhoSqvpp,nsdim);
        double* drhoSqvppdy2 = shiftdown(rhoSqvpp,nsdim);

        for (int i = 0; i<nsdim/4; i++)
        {
            dudt[t*nsdim+nsdim/2+i] = -(drhouvdx1[i]-drhouvdx2[i])/(dx[i]*2) - (drhoSqvppdy1[i] - drhoSqvppdy2[i])/(dx[i]*2);
        }


        delete rhoSqvpp;
        delete rhouv;
        delete drhouvdx1;
        delete drhouvdx2;
        delete drhoSqvppdy1;
        delete drhoSqvppdy2;

        double* rhoeupp = new double[nsdim/4];
        double* rhoevpp = new double[nsdim/4];

        for (int i =0; i<nsdim/4; i++)
        {
            rhoeupp[i] = (uMid[t*nsdim+3*nsdim/4+i]+pow(uMid[t*nsdim+i],gamma)/gamma)/uMid[t*nsdim+i]*uMid[t*nsdim+nsdim/4+i];
            rhoevpp[i] = (uMid[t*nsdim+3*nsdim/4+i]+pow(uMid[t*nsdim+i],gamma)/gamma)/uMid[t*nsdim+i]*uMid[t*nsdim+nsdim/2+i];
        }
        
        double* drhoeuppdx1 = shiftupx(rhoeupp,nsdim);
        double* drhoeuppdx2 = shiftdownx(rhoeupp,nsdim);
        double* drhoevppdy1 = shiftup(rhoevpp,nsdim);
        double* drhoevppdy2 = shiftdown(rhoevpp,nsdim);

        for (int i =0; i<nsdim/4; i++)
        {
            dudt[t*nsdim+3*nsdim/4+i] = -(drhoeuppdx1[i]-drhoeuppdx2[i])/(dx[i]*2) - (drhoevppdy1[i]-drhoevppdy2[i])/(dx[i]*2);
        }

        delete rhoeupp;
        delete rhoevpp;
        delete drhoeuppdx1;
        delete drhoeuppdx2;
        delete drhoevppdy1;
        delete drhoevppdy2;



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
