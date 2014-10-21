#include <math.h>
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpcasa.h"
#include "petscvec.h"
#include "burgers2D.hpp"
//#include "eulervortex.hpp"
#include "stdlib.h"
#include <iostream>
#include "LSSm2.hpp"
#include <stdio.h>
#include <fstream>
        lssSolver::lssSolver(double* t, double* u0, int p, int c, int cg, int co, int f, int fg, int fo, int* mapcoarse, int* mapfine, int ntdim, int ntdimg, int trank, int tsize, int nsdim, int nsdimg, int srank, int ssize, double alpha, MPI_Comm comm)
        {
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            this->t = t;
            u = u0;
            this->t = t;
            this->p = p;
            this->ntdim = ntdim; 
            this->ntdimg = ntdimg;
            nn = (ntdim-1)/p+1; 
            nng = (ntdimg-1)/p+1;
            this->nsdim = nsdim;
            this->nsdimg = nsdimg;
            this->trank = trank;
            this->tsize = tsize;
            this->srank = srank;
            this->ssize = ssize;
            this->f = f;
            this->c = c;
            this->fg = fg;
            this->cg = cg;
            this->fo = fo;
            this->co = co;
            this->comm = comm; 
// Define dt
            dt = new double[ntdim-1];
            for (int i =0; i<ntdim-1; i++)
            {
                dt[i] = t[i+1]-t[i];
            }
// Define s, the larger timestep, and ds       
            s = new double[nn];
            for (int i = 0; i<nn; i++)
            {
                s[i] = t[i*p];
            }
            ds = new double[nn-1];
            for (int i = 0; i<nn-1; i++)
            {
                ds[i] = s[i+1] - s[i];
            }
            
// Find the linear interpolation of the coarse grid values
            for (int i = 0; i < nn-1; i++)
            {
                for ( int k = 0 ; k < c; k ++)
                {
                    double sub = 0;
                    for (int j =1; j < p; j ++)
                    {
                        sub += dt[i*p+j-1];
                        u[(i*p+j)*nsdim+mapcoarse[k]] = (1-sub/ds[i])*u[i*p*nsdim+mapcoarse[k]] + (sub/ds[i])*u[(i+1)*p*nsdim+mapcoarse[k]];
                    }
                }
            }

// Calculate mid values of fine grid values
            uMidf = new double[nsdim*(ntdim-1)];
            for (int i = 0; i<ntdim-1; i++)
            {
                for (int j = 0; j<nsdim; j++)
                {
                    //uMidf[i*nsdim+j] = 0.5*(u[(i+1)*nsdim+j] + u[i*nsdim+j]);
                    uMidf[i*nsdim+j] = u[i*nsdim+j];
                }
            }
           
 // Calculate mid values of coarse grid values
            uMidc = new double[nsdim*(nn-1)];
            for (int i = 0; i < nn-1; i++)
            {
                for (int j = 0; j < nsdim; j++)
                {
                    //uMidc[i*nsdim+j] = 0.5*(u[(i+1)*p*nsdim+j] + u[i*p*nsdim+j]);
                    uMidc[i*nsdim+j] = u[i*p*nsdim+j];
                }
            }


// Calculate dudt in fine and coarse grid points
            dudtf = new double[(ntdim-1)*f];
            for (int i = 0; i<ntdim-1; i++)
            {
                for (int j = 0; j< f; j++)
                {
                    dudtf[i*f+j]  = (u[(i+1)*nsdim+mapfine[j]] - u[i*nsdim+mapfine[j]])/dt[i];
		}
            }

            dudtc = new double[(nn-1)*c];
            for (int i = 0; i<nn-1; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    dudtc[i*c+j] = (u[(i+1)*p*nsdim+mapcoarse[j]] - u[i*p*nsdim+mapcoarse[j]])/ds[i];
                }
            }

            this->mapcoarse = mapcoarse;
            this->mapfine = mapfine;
//Inverse mat of coarse and fine
            coarse = new int[nsdim];
            for (int i = 0; i<nsdim; i++)
            {
                coarse[i] = -1;
            }
            for (int i = 0; i < c; i++)
            {
                coarse[mapcoarse[i]] = i;
            }

            fine = new int[nsdim];
            for (int i = 0; i< nsdim; i++)
            {
                fine[i] = -1;
            }
            for (int i = 0; i < f; i++)
            {
                fine[mapfine[i]] = i;
            }

            this->alpha = alpha;
        };

// Apply diagonal preconditioner
        Vec lssSolver::APinvx(Vec x)
        {
            Vec y;
            VecDuplicate(x,&y);
            VecAssemblyBegin(y);
            VecAssemblyEnd(y);


            Vec temp1;
            VecCreate(PETSC_COMM_WORLD, &temp1);
            VecSetSizes(temp1,(nn-1)*(p+1),(nng-1)*(p+1));
            VecSetFromOptions(temp1);
            VecAssemblyBegin(temp1);
            VecAssemblyEnd(temp1);

            Vec temp2;
            VecCreate(PETSC_COMM_WORLD, &temp2);
            VecDuplicate(temp1, &temp2);
            VecAssemblyBegin(temp2);
            VecAssemblyEnd(temp2);
/*
            Vec temp3;
            VecCreate(PETSC_COMM_WORLD, &temp3);
            VecSetSizes(temp3, (nn-1)*(p*f+c), (nng-1)*(p*fg+cg));
            VecSetFromOptions(temp3);
            VecAssemblyBegin(temp3);
            VecAssemblyEnd(temp3);
*/        
            MatMultTranspose(E, x, temp1);
            MatMult(wEinv, temp1, temp2);
            MatMult(E, temp2, y);

            //MatMult(E, temp2, temp3);
            //MatMult(Pinv,temp3, y);

            VecDestroy(&temp1);
            VecDestroy(&temp2);
            //VecDestroy(&temp3);

            VecCreate(PETSC_COMM_WORLD, &temp1);
            if (trank == tsize-1)
            {
                VecSetSizes(temp1,(nn-1)*(p*f+c)+f+c,(nng-1)*(p*fg+cg)+fg+cg);
            }
            else
            {
                VecSetSizes(temp1,(nn-1)*(p*f+c), (nng-1)*(p*fg+cg)+fg+cg);
            }
            VecSetFromOptions(temp1);
            VecAssemblyBegin(temp1);
            VecAssemblyEnd(temp1);

            VecCreate(PETSC_COMM_WORLD, &temp2);
            VecDuplicate(temp1, &temp2);
            VecAssemblyBegin(temp2);
            VecAssemblyEnd(temp2);
/*
            VecCreate(PETSC_COMM_WORLD, &temp3);
            VecSetSizes(temp3, (nn-1)*(p*f+c),(nng-1)*(p*fg+cg));
            VecSetFromOptions(temp3);
            VecAssemblyBegin(temp3);
            VecAssemblyEnd(temp3);
*/

            Vec Ax;
            VecCreate(PETSC_COMM_WORLD, &Ax);
            VecSetSizes(Ax,(nn-1)*(p*f+c),(nng-1)*(p*fg+cg));
            VecSetFromOptions(Ax);
            VecAssemblyBegin(Ax);
            VecAssemblyEnd(Ax);

            MatMultTranspose(B, x, temp1);
            MatMult(wBinv, temp1, temp2);
            MatMult(B, temp2, Ax);
            
            //MatMult(B, temp2, temp3);
            //MatMult(Pinv, temp3, Ax);

            VecDestroy(&temp1);
            VecDestroy(&temp2);
            //VecDestroy(&temp3);

            VecAXPY(Ax, 1, y);

            VecDestroy(&y);
            return Ax;
        }
// CG solve
        void lssSolver::CGSolve(Vec b, Vec* w)
        {

            double Tol = 1.e-7;
            int maxiter = 1000;
            Vec r, p;
            
            VecDuplicate(b,&r);
            VecCopy(b, r);
            Vec Ax0 = APinvx(*w); 

            double norm0;

            VecAXPY(r, -1, Ax0);
            VecDestroy(&Ax0);

            VecDuplicate(r, &p);
            VecCopy(r, p);
            VecNorm(r,NORM_2,&norm0);
            if (rank == 0)
            std::cout << "Norm of r0 is " << norm0 << std::endl;

            int k = 0;

            while (1)
            {

                PetscScalar rr,pAp;
                VecTDot(r,r,&rr);
                Vec Ap = APinvx(p);
                VecTDot(p, Ap, &pAp);

                PetscScalar alpha = rr/pAp;
                
                VecAXPY(*w, alpha, p);

                VecAXPY(r, -alpha, Ap);

                PetscScalar norm;
                VecNorm(r,NORM_2,&norm);
                if ((norm < Tol*norm0) || (k > maxiter))                
                {
                    if (k>maxiter)
                    {
                        if (rank == 0)
                        std::cout << "exceeded maxiter" << std::endl;
                    }

                    VecDestroy(&r);
                    VecDestroy(&p);
                    VecDestroy(&Ap);
                    if (rank == 0)
                    std::cout << "converged in " << k << " iterations, with norm of " << norm << std::endl;
                    break;
                }

                PetscScalar rr2;
                VecTDot(r,r,&rr2);
                PetscScalar beta = rr2 / rr;
                VecAYPX(p,beta,r);
                k ++;
            } 
        }
// main solver of MTLSS
        double* lssSolver::lss(int maxIter,double atol,double rtol)
        {

            Vec b,w,v,eta,temp4;

            VecCreate(PETSC_COMM_WORLD, &b);
            VecSetSizes(b, (nn-1)*(p*f+c), (nng-1)*(p*fg+cg));
            VecSetFromOptions(b);
            double* funcvalf = F(uMidf, 1, ntdim, nsdim);
            double* funcvalc = F(uMidc, 1, nn, nsdim);

            for (int i=0; i< nn-1; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    VecSetValue(b, (trank*(nn-1)+i)*(p*fg+cg)+co+j, dudtc[i*c+j] - funcvalc[i*nsdim+mapcoarse[j]], INSERT_VALUES);
		}
                for (int j = 0; j < p; j++)
                {
                    for (int k = 0; k < f; k++)
                    {
                        VecSetValue(b, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, dudtf[(i*p+j)*f+k] - funcvalf[(i*p+j)*nsdim+mapfine[k]], INSERT_VALUES);
		    }
                }
            }
            

            VecAssemblyBegin(b);
            VecAssemblyEnd(b);
            PetscReal norm_b0;
            VecNorm(b,NORM_2, &norm_b0);
            if (rank == 0)
            {
                std::cout << "Norm of b is " << norm_b0 << std::endl;        
            }

            for (int iNewton=0; iNewton<maxIter; iNewton++)
            {

                MatCreate(PETSC_COMM_WORLD, &B);
                if (trank == tsize-1)
                {
                    MatSetSizes(B, (nn-1)*(p*f+c),(nn-1)*(p*f+c)+f+c, (nng-1)*(p*fg+cg),(nng-1)*(p*fg+cg)+fg+cg);
                }
                else
                {
                    MatSetSizes(B, (nn-1)*(p*f+c),(nn-1)*(p*f+c), (nng-1)*(p*f+c),(nng-1)*(p*fg+cg)+fg+cg);
                }
                MatSetType(B, MATMPIAIJ);
                MatSetFromOptions(B);
                MatSetUp(B);
                MatMPIAIJSetPreallocation(B,2*nsdim,NULL,2*nsdim,NULL);
		double* val = 0;
		int* ia = 0;
		int* ja = 0;
		int nnz;
                for (int i = 0; i < nn-1; i++)
                {
                    dfdu(uMidc,s,nsdim,i,val,ia,ja,nnz);
                    for( int k = 0; k < nnz; k++)
                    {
                        if ((coarse[ia[k]]!=-1)&&(coarse[ja[k]]!=-1))
                        {
                            MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+co+coarse[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+co+coarse[ja[k]], val[k], INSERT_VALUES);
                            MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+co+coarse[ia[k]], (trank*(nn-1)+i+1)*(p*f+c)+co+coarse[ja[k]], val[k], INSERT_VALUES);
                        }
                        else { 
                            if ((coarse[ia[k]]!=-1)&&(fine[ja[k]]!=-1))
                            {
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+co+coarse[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+cg+fo+fine[ja[k]], val[k], INSERT_VALUES);
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+co+coarse[ia[k]], (trank*(nn-1)+i+1)*(p*fg+cg)+cg+fo+fine[ja[k]], val[k], INSERT_VALUES);
                            }
                        }

                    }

                    dfdu(uMidf,s,nsdim,i*p,val,ia,ja,nnz);
                    for (int k =0; k<nnz; k++)
                    {
                        if ((fine[ia[k]]!=-1)&&(fine[ja[k]]!=-1))
                        {
                            MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+fo+fine[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+cg+fo+fine[ja[k]], val[k], INSERT_VALUES);
                            MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+fo+fine[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+cg+fg+fo+fine[ja[k]], val[k], INSERT_VALUES);
                        }
                        
                        else { 
                            if ((fine[ia[k]]!=-1)&&(coarse[ja[k]]!=-1))
                            {
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+fo+fine[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+co+coarse[ja[k]], (2-dt[i*p]/ds[i])*val[k], INSERT_VALUES);
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+fo+fine[ia[k]], (trank*(nn-1)+i+1)*(p*fg+cg)+co+coarse[ja[k]], dt[i*p]/ds[i]*val[k], INSERT_VALUES);
                            }
                        }
                    }
                    double sub=0;
                    for (int j = 1; j < p; j++)
                    {
                        sub += dt[i*p+j-1];
                        dfdu(uMidf,s,nsdim,i*p+j,val,ia,ja,nnz);    
                        for (int k = 0; k < nnz; k++)
                        {
                            if ((fine[ia[k]]!=-1)&&(fine[ja[k]]!=-1))
                            {
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+fine[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+fine[ja[k]], val[k], INSERT_VALUES);
                                if (j != p-1)
                                {
                                    MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+fine[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+cg+(j+1)*fg+fo+fine[ja[k]], val[k], INSERT_VALUES);
                                }
                                else
                                {
                                    MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+fine[ia[k]], (trank*(nn-1)+i+1)*(p*fg+cg)+cg+fo+fine[ja[k]], val[k], INSERT_VALUES);
                                }
                            }

                            else if ((fine[ia[k]]!=-1)&&(coarse[ja[k]]!=-1))
                            {
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+fine[ia[k]], (trank*(nn-1)+i)*(p*fg+cg)+co+coarse[ja[k]], (2-(2*sub+dt[i*p+j])/ds[i])*val[k], INSERT_VALUES);
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+fine[ia[k]], (trank*(nn-1)+i+1)*(p*fg+cg)+co+coarse[ja[k]], (2*sub+dt[i*p+j])/ds[i]*val[k], INSERT_VALUES);
                            }
                        }
                    }
                }
		delete val;
		delete ia;
		delete ja;


                MatAssemblyBegin(B, MAT_FLUSH_ASSEMBLY);
                MatAssemblyEnd(B, MAT_FLUSH_ASSEMBLY);

                for (int i = 0; i < nn-1; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+co+j, (trank*(nn-1)+i)*(p*fg+cg)+co+j, -1/ds[i], ADD_VALUES);
                        MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+co+j, (trank*(nn-1)+i+1)*(p*fg+cg)+co+j, 1/ds[i], ADD_VALUES);
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, -1/dt[i*p+j], ADD_VALUES);
                            if (j != p-1)
                            {
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (trank*(nn-1)+i)*(p*fg+cg)+cg+(j+1)*fg+fo+k, 1/dt[i*p+j], ADD_VALUES);
                            }
                            else {
                                MatSetValue(B, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (trank*(nn-1)+i+1)*(p*fg+cg)+cg+fo+k, 1/dt[i*p+j], ADD_VALUES);
                            }
                        }
                    }
                }
                MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

                MatCreate(PETSC_COMM_WORLD, &E);
                MatSetSizes(E, (nn-1)*(p*f+c), (nn-1)*(p+1), (nng-1)*(p*fg+cg), (nng-1)*(p+1));
                MatSetType(E, MATMPIAIJ);
                MatSetFromOptions(E);
                MatSetUp(E);
                MatMPIAIJSetPreallocation(E,1,NULL,1,NULL);
                for (int i = 0; i < (nn-1); i++)
                {
                    for (int j = 0; j< c; j++)
                    {
                        MatSetValue(E, (trank*(nn-1)+i)*(p*fg+cg)+co+j, (trank*(nn-1)+i)*(p+1), dudtc[i*c+j], INSERT_VALUES);
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            MatSetValue(E, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (trank*(nn-1)+i)*(p+1)+1+j, dudtf[i*p*f+j*f+k], INSERT_VALUES);
                        }
                    }
                }

                MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY);


                MatCreate(PETSC_COMM_WORLD, &wEinv);
                if (srank == 0)
                    MatSetSizes(wEinv, (nn-1)*(p+1), (nn-1)*(p+1), (nng-1)*(p+1), (nng-1)*(p+1));
                else
                    MatSetSizes(wEinv, 0, 0, (nng-1)*(p+1), (nng-1)*(p+1));
                MatSetType(wEinv, MATMPIAIJ);
                MatSetFromOptions(wEinv);
                MatSetUp(wEinv);
                MatMPIAIJSetPreallocation(wEinv,1,NULL,1,NULL);

                if (srank == 0)
                {
                    for (int i = 0; i< (nn-1); i++)
                    {
                        MatSetValue(wEinv, (trank*(nn-1)+i)*(p+1), (trank*(nn-1)+i)*(p+1), (((double) tsize/(trank+1))*t[ntdim-1])/(ds[i]*pow(alpha,2)),INSERT_VALUES);
                        for (int j = 0; j<p; j++)
                        {
                            MatSetValue(wEinv, (trank*(nn-1)+i)*(p+1)+1+j, (trank*(nn-1)+i)*(p+1)+1+j, (((double) tsize/(trank+1))*t[ntdim-1])/(dt[i*p+j]*pow(alpha,2)), INSERT_VALUES);
                        }
                    }
                }
                    
                MatAssemblyBegin(wEinv, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(wEinv, MAT_FINAL_ASSEMBLY);

                MatCreate(PETSC_COMM_WORLD, &wBinv);
		if (trank != tsize-1)
		{
                    MatSetSizes(wBinv, (nn-1)*(p*f+c), (nn-1)*(p*f+c), (nng-1)*(p*fg+cg)+fg+cg, (nng-1)*(p*fg+cg)+fg+cg);
                }
                else 
                {
                    MatSetSizes(wBinv, (nn-1)*(p*f+c)+f+c, (nn-1)*(p*f+c)+f+c, (nng-1)*(p*fg+cg)+fg+cg, (nng-1)*(p*fg+cg)+fg+cg);
                }
                MatSetType(wBinv, MATMPIAIJ);
                MatSetFromOptions(wBinv);
                MatSetUp(wBinv);
                MatMPIAIJSetPreallocation(wBinv,1,NULL,1,NULL);
                for (int i = 0; i < nn-1; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        if ((trank != 0) || (i != 0))
                        {    
                            MatSetValue(wBinv, (trank*(nn-1)+i)*(p*fg+cg)+co+j,(trank*(nn-1)+i)*(p*fg+cg)+co+j, (((double) tsize/(trank+1))*t[ntdim-1])/((ds[i])), INSERT_VALUES);
                        }
                        else {
                            MatSetValue(wBinv, (trank*(nn-1)+i)*(p*fg+cg)+co+j,(trank*(nn-1)+i)*(p*fg+cg)+co+j, (((double) tsize/(trank+1))*t[ntdim-1])/(ds[i]/2), INSERT_VALUES);
                        }
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k ++)
                        {
                            if ((trank!= 0) || ((i+j)!=0))
                            {
                                MatSetValue(wBinv, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (((double) tsize/(trank+1))*t[ntdim-1])/((dt[i*p+j])), INSERT_VALUES);
                            }
                            else
                            {
                                MatSetValue(wBinv, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, (((double) tsize/(trank+1))*t[ntdim-1])/(dt[i*p+j]/2), INSERT_VALUES);
                            }
                        }
                    }
                }
		if (trank == tsize-1)
                {
                    for (int j = 0; j < c; j++)
                    {
                        MatSetValue(wBinv, (nng-1)*(p*fg+cg)+co+j, (nng-1)*(p*fg+cg)+co+j, (((double) tsize/(trank+1))*t[ntdim-1])/(ds[nn-2]/2), INSERT_VALUES);
                    }
                    for (int j = 0; j < f; j++)
                    {
                        MatSetValue(wBinv, (nng-1)*(p*fg+cg)+cg+fo+j, (nng-1)*(p*fg+cg)+cg+fo+j, (((double) tsize/(trank+1))*t[ntdim-1])/(dt[ntdim-2]/2), INSERT_VALUES);
                    }
                }
                MatAssemblyBegin(wBinv, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(wBinv, MAT_FINAL_ASSEMBLY);

                VecDuplicate(b,&w);
                VecSet(w, 0.);
                VecAssemblyBegin(w);
                VecAssemblyEnd(w);
/*            
                Vec b1;

                VecDuplicate(w, &b1);
                VecAssemblyBegin(b1);
                VecAssemblyEnd(b1);

                MatCreate(PETSC_COMM_WORLD, &Pinv);
                MatSetType(Pinv, MATMPIAIJ);
                MatSetSizes(Pinv, (nn-1)*(p*f+c), (nn-1)*(p*f+c), (nng-1)*(p*fg+cg), (nng-1)*(p*fg+cg));
                MatSetFromOptions(Pinv);
                MatSetUp(Pinv);
                MatMPIAIJSetPreallocation(Pinv,1,NULL,1,NULL);

                for (int i = 0; i < (nn-1)*(p*f+c); i++)
                {
*/
/*
                    PetscInt ncols;
                    const PetscInt *cols;
                    const PetscScalar *vals;
                    MatGetRow(B,(trank*(nn-1))*(p*f+c)+i,&ncols,&cols,&vals);

                    double sum = 0;
                    for (int j = 0; j <ncols; j++)
                    {
                        int blah = cols[j];
                        PetscScalar wBinvval;
                        MatGetValues(wBinv,1,&blah,1,&blah,&wBinvval);
                        sum += vals[j]*vals[j]*wBinvval;
                    }


                    MatRestoreRow(B,(trank*(nn-1))*(p*f+c)+i,&ncols,&cols,&vals);

                    MatGetRow(E,(trank*(nn-1))*(p*f+c)+i,&ncols,&cols,&vals);


                    for (int j = 0; j <ncols; j++)
                    {
                        int blah = cols[j];
                        PetscScalar wEinvval;
                        MatGetValues(wEinv,1,&blah,1,&blah,&wEinvval);
                        sum += vals[j]*vals[j]*wEinvval;
                    }

                    MatRestoreRow(E, (trank*(nn-1))*(p*f+c)+i, &ncols, &cols, &vals);
*/
/*
// this part needs work
                    double sum = 1;
                    MatSetValue(Pinv, trank*(nn-1)*(p*fg+cg)+i, trank*(nn-1)*(p*fg+cg)+i, 1./sum, INSERT_VALUES);
                }

                MatAssemblyBegin(Pinv,MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(Pinv,MAT_FINAL_ASSEMBLY);

                MatMult(Pinv,b,b1);
*/
                //CGSolve(b1,&w);
                CGSolve(b,&w);


                //VecDestroy(&b1);
                VecDestroy(&b);
                //MatDestroy(&Pinv);

                VecCreate(PETSC_COMM_WORLD, &temp4);
                if (trank == tsize-1)
                {
                    VecSetSizes(temp4,(nn-1)*(p*f+c)+f+c,(nng-1)*(p*fg+cg)+fg+cg);
                }
                else
                {
                    VecSetSizes(temp4,(nn-1)*(p*f+c), (nng-1)*(p*fg+cg)+fg+cg);
                }
                VecSetFromOptions(temp4);
                VecSet(temp4, 0.);
                VecAssemblyBegin(temp4);
                VecAssemblyEnd(temp4);


                VecDuplicate(temp4,&v);
                VecAssemblyBegin(v);
                VecAssemblyEnd(v);

                MatMultTranspose(B, w, temp4);
                MatMult(wBinv, temp4, v);
                VecDestroy(&temp4);

                VecCreate(PETSC_COMM_WORLD, &temp4);
                VecSetSizes(temp4, (nn-1)*(p+1), (nng-1)*(p+1));
                VecSetFromOptions(temp4);
                VecSet(temp4, 0.);
                VecAssemblyBegin(temp4);
                VecAssemblyEnd(temp4);

                VecDuplicate(temp4,&eta);
                VecAssemblyBegin(eta);
                VecAssemblyEnd(eta);

                MatMultTranspose(E, w, temp4);
                MatMult(wEinv, temp4, eta);
                VecDestroy(&temp4);

                MatDestroy(&B);
                MatDestroy(&E);
                MatDestroy(&wEinv);
                MatDestroy(&wBinv);

                PetscScalar y;


                for (int i = 0; i < nn-1 ; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        int loc = (trank*(nn-1)+i)*(p*fg+cg)+co+j;
                        VecGetValues(v, 1 , &loc, &y);
                        u[i*p*nsdim+mapcoarse[j]] = u[i*p*nsdim+mapcoarse[j]]-y;
                    }
                    for (int j = 0; j <p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            int loc = (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k;
                            VecGetValues(v, 1, &loc, &y);
                            u[(i*p+j)*nsdim + mapfine[k]] = u[(i*p+j)*nsdim + mapfine[k]] - y;
                        }
                    }
                }
                if (trank == size-1)
                {
                    for (int j = 0; j < c; j++)
                    {
                        int loc = (nng-1)*(p*fg+cg)+co+j;
                        VecGetValues(v,1,&loc,&y);
                        u[(nn-1)*p*nsdim+mapcoarse[j]] = u[(nn-1)*p*nsdim+mapcoarse[j]] - y;
                    }
                    for (int j = 0; j < f; j++)
                    {
                        int loc = (nng-1)*(p*fg+cg)+cg+fo+j;
                        VecGetValues(v,1,&loc,&y);
                        u[(nn-1)*p*nsdim+mapfine[j]] = u[(nn-1)*p*nsdim+mapfine[j]] - y;
                    }
                }
                MPI_Request request_send, request_recv;
                MPI_Status status;
                
                double* rbuf = new double[nsdim];
                double* sbuf = new double[nsdim];

                for (int bb = 0; bb < nsdim; bb++)
                {
                    sbuf[bb] = u[bb];
                }
                if (trank != 0)
                    MPI_Send(&sbuf[0], nsdim, MPI_DOUBLE, (trank*ssize + srank)-1, 123, MPI_COMM_WORLD);

                if (trank != tsize-1)
                    MPI_Recv(&rbuf[0], nsdim, MPI_DOUBLE, (trank*ssize + srank)+1, 123, MPI_COMM_WORLD,&status);

                if (trank != tsize-1)
                {
                    for (int bb = 0; bb < nsdim; bb++)
                    {
                        u[(nn-1)*p*nsdim+bb] = rbuf[bb];
                    }
                }


                for (int i = 0; i < nn-1; i++)    
                {
                    for (int j =0; j < p; j++)
                    {
                        int loc1 = (trank*(nn-1)+i)*(p+1)+1+j;
                        int loc2 = (trank*(nn-1)+i)*(p+1);
                        VecGetValues(eta,1,&loc1, &y);
                        dt[i*p+j] *= exp(y);
                        VecGetValues(eta,1,&loc2, &y);
                        dt[i*p+j] *= exp(y*dt[i*p+j]/ds[i]);
                    }
                }

                t[0] = 0;
                for (int i =1; i<ntdim; i++)
                {
                    t[i] = t[i-1]+dt[i-1];
                }
    
                for (int i = 0; i<nn; i++)
                {
                    s[i] = t[i*p];
                }
                for (int i = 0; i<nn-1; i++)
                {
                    ds[i] = s[i+1] - s[i];
                }

                VecDestroy(&eta);
                VecDestroy(&v);
                VecDestroy(&w);


                for (int i = 0; i < nn-1; i++)
                {
                    for ( int k = 0 ; k < c; k ++)
                    {
                        double sub = 0;
                        for (int j =1; j < p; j ++)
                        {
                            sub = t[i*p+j]-s[i];
                            u[(i*p+j)*nsdim+mapcoarse[k]] = (1-sub/ds[i])*u[i*p*nsdim+mapcoarse[k]] + (sub/ds[i])*u[(i+1)*p*nsdim+mapcoarse[k]];
                        }
                    }
                }

                for (int i = 0; i<ntdim - 1; i++)
                {
                    for (int j = 0; j< nsdim; j++)
                    { 
                        //uMidf[i*nsdim+j] = 0.5*(u[(i+1)*nsdim+j] + u[i*nsdim+j]);
                        uMidf[i*nsdim+j] = u[i*nsdim+j];
                    }
                }
    
                for (int i = 0; i < nn-1; i++)
                {
                    for (int j = 0; j < nsdim; j++)
                    {
                        //uMidc[i*nsdim+j] = 0.5*(u[(i+1)*p*nsdim+j] + u[i*p*nsdim+j]);
                        uMidc[i*nsdim+j] = u[i*p*nsdim+j];
                    }
                }


                for (int i = 0; i<ntdim-1; i++)
                {
                    for (int j = 0; j< f; j++)
                    {
                        dudtf[i*f+j]  = (u[(i+1)*nsdim+mapfine[j]] - u[i*nsdim+mapfine[j]])/dt[i];
                    }
                }


                for (int i = 0; i<nn-1; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        dudtc[i*c+j] = (u[(i+1)*p*nsdim+mapcoarse[j]] - u[i*p*nsdim+mapcoarse[j]])/ds[i];
                    }
                }


                VecCreate(PETSC_COMM_WORLD, &b);
                VecSetSizes(b, (nn-1)*(p*f+c),(nng-1)*(p*fg+cg));
                VecSetFromOptions(b);

                funcvalf = F(uMidf, 1, ntdim, nsdim);
                funcvalc = F(uMidc, 1, nn, nsdim);


                for (int i=0; i< nn-1; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        VecSetValue(b, (trank*(nn-1)+i)*(p*fg+cg)+co+j, dudtc[i*c+j] - funcvalc[i*nsdim+mapcoarse[j]], INSERT_VALUES);
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            VecSetValue(b, (trank*(nn-1)+i)*(p*fg+cg)+cg+j*fg+fo+k, dudtf[(i*p+j)*f+k] - funcvalf[(i*p+j)*nsdim+mapfine[k]], INSERT_VALUES);
                        }
                    }    
                }


                VecAssemblyBegin(b);
                VecAssemblyEnd(b);


                PetscReal norm_b;
                VecNorm(b,NORM_2, &norm_b);

		if (rank == 0)
                std::cout << "Norm of b is " << norm_b << std::endl;
                if ((norm_b < atol) || (norm_b < rtol*norm_b0))
                {

			std::ofstream outfile ("dt.txt");
                        outfile.precision(15);
    			if(outfile.is_open())
			{
        			for(int i = 0 ; i < ntdim-1; i++)
                		{
					outfile << dt[i] ;
            				outfile <<  std::endl;
				}
        		}
    	
                    return u;
                }


            }

            std::cout << "Failed in iteration " << maxIter << std::endl;
            return u;

        };
