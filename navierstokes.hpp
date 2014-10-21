#ifndef ns_hpp
#define ns_hpp

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>


class nssolver
{
    private:
        int Nx, Ny, Nz;
        int bs;
        int ntdim, nsdim, nsdim_x, nsdim_y, nsdim_z;
        int srank, ssize, srank_x, srank_y, srank_z, ssize_x, ssize_y, ssize_z;
        MPI_COMM comm;
        double* u;

    public:
        nssolver(int ntdim, int nsdim, int nsdim_x, int nsdim_y, int nsdim_z, int bs, MPI_COMM comm);
        int ID(int i, int j, int k);
        double* F(double* uMid, double p)
        void dfdu(double* uMid, double* s, int t, double*& val, int*& ia, int*& ja, int& nnz)
}


nssolver::nssolver(int ntdim, int nsdim, MPI_COMM comm);
{
    comm = this->comm;
    ntdim = this->ntdim;
    nsdim = this->nsdim;

    MPI_Comm_Rank(comm, &srank);
    MPI_Comm_Size(comm, &ssize);

    int i = int(pow((nsdim*nsdim_x*nsdim_x)/(nsdim_y*nsdimz),1./3.))+1;
    for (int k = i; k > 0; k--)
    {
        if (nsdim % k == 0)
        {
            ssize_x = k;
            break;
        }
    }
    i = int(pow(nsdim*nsdim_y/(ssize_x*nsdim_z),1./2.))+1;
    for (int k = i; k > 0; k--)
    {
        if (nsdim % (ssize_x*k) == 0)
        {
            ssize_y = k;
            break;
        }
    }
    ssize_z = nsdim / (ssize_x * ssize_y);

    srank_x = srank % (ssize_y * ssize_z);
    srank_y = ((srank-srank_x)/ssize_y) % ssize_z;
    srank_z = (((srank-srank_x)/ssize_y)-srank_y) / ssize_z; 


    nsdim_x = 





    u = new(double*)[Nx*Ny*Nz];

}



int ID(int i, int j, int k)
{

    return ( i * Ny * Nz + j * Nz + k );

}

double* F(double* uMid, double p, int ntdim, int nsdim)
{

    double dx[nsdim/5];
    for (int i = 0; i < nsdim/5; i++)
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
        double* rhou_xu = shiftupx(&uMid[t*nsdim+nsdim/5],nsdim);
        double* rhou_xd = shiftdownx(&uMid[t*nsdim+nsdim/5],nsdim);
        double* rhov_yu = shiftup(&uMid[t*nsdim+2*nsdim/5],nsdim);
        double* rhov_yd = shiftdown(&uMid[t*nsdim+2*nsdim/5],nsdim);
        double* rhow_zu = shiftup(&uMid[t*nsdim+3*nsdim/5],nsdim);
        double* rhow_zd = shiftdown(&uMid[t*nsdim+3*nsdim/5],nsdim);
        for (int i = 0; i<nsdim/5; i++)
        {
            dudt[t*nsdim+i] = -(rhou_xu[i]-rhou_xd[i])/(dx[i]*2) - (rhov_yu[i]-rhov_yd[i])/(dx[i]*2);
        }

        delete rhou_xu;
        delete rhou_xd;
        delete rhov_yu;
        delete rhov_yd;
        delete rhow_zu;
        delete rhow_zd;
        
	double* drho = new double[nsdim/5];
        double* rhouv = new double[nsdim/5];
        for (int i = 0; i<nsdim/4; i++)
        {
            rhoSqupp[i] = uMid[t*nsdim+nsdim/4+i]*uMid[t*nsdim+nsdim/4+i]/uMid[t*nsdim+i]+pow(uMid[t*nsdim+i],gamma)/gamma;
            rhouv[i] = uMid[t*nsdim+nsdim/4+i]*uMid[t*nsdim+nsdim/2+i]/uMid[t*nsdim+i];
        }
        double* drhox1 = shiftupx(rhoSqupp,nsdim);
        double* drhox2 = shiftdownx(rhoSqupp,nsdim);
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
