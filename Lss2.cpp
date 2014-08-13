#ifndef lsssolver_hpp
#define lsssolver_hpp
#include <cstdio>
#include <ctime>
#include <math.h>
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpcasa.h"
#include "petscvec.h"
#include "eulervortex.hpp"
#include "stdlib.h"
#include <iostream>

class lssSolver
{

    private:
        double* t;
        double* s;
        double* u;
        double* dt;
        double* uMid;
        double* dudt;
        double alpha;
        int ntdim;
        int nsdim;
        Mat B;
        Mat E;
        Mat wBinv;
        Mat wEinv;


    public:
        lssSolver(double* t, double* u0, double* s, int ntdim, int nsdim, double alpha)
        {

            this->t = t;
            u = u0;
            this->s = s;
            this->t = t;
            this->ntdim = ntdim;
            this->nsdim = nsdim;

            dt = (double*)malloc((ntdim-1)*sizeof(double));
            for (int i =0; i<ntdim-1; i++)
            {
                    dt[i] = t[i+1]-t[i];
            }
                
            uMid = (double*)malloc(nsdim*(ntdim-1)*sizeof(double));
            for (int i = 0; i<ntdim-1; i++)
            {

                for (int j = 0; j<nsdim; j++)
                {
                    uMid[i*nsdim+j] = 0.5*(u[(i+1)*nsdim+j] + u[i*nsdim+j]);
                }
            }


            dudt = (double*)malloc((ntdim-1)*nsdim*sizeof(double));
            for (int i = 0; i<ntdim-1; i++)
            {
                for (int j = 0; j< nsdim; j++)
                {
                    dudt[i*nsdim+j]  = (u[(i+1)*nsdim+j] - u[i*nsdim+j])/dt[i];
                }
            }

            this->alpha = alpha;
        };


        double* lss(int maxIter=8,double atol=1e-7,double rtol=1e-4)
        {

            Mat Smat;
            Vec b,w,v,eta,temp4;

            VecCreate(MPI_COMM_WORLD, &b);
            VecSetSizes(b, PETSC_DECIDE, nsdim*(ntdim-1));
            VecSetFromOptions(b);

            double* funcval = F(uMid, s, ntdim, nsdim);

            for (int i=0; i<(ntdim-1)*nsdim; i++)
            {
                VecSetValue(b, i, dudt[i] - funcval[i], INSERT_VALUES);
            }
            
            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            PetscReal norm_b0;
            VecNorm(b,NORM_2, &norm_b0);

            std::cout << "Norm of b is " << norm_b0 << std::endl;        


            for (int iNewton=0; iNewton<maxIter; iNewton++)
            {

                MatCreate(MPI_COMM_WORLD, &B);
                            MatSetType(B, MATSEQAIJ);
                            MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, (ntdim-1)*nsdim, ntdim*nsdim);
                            MatSetFromOptions(B);
                            MatSetUp(B);
                            MatSeqAIJSetPreallocation(B,2*nsdim,NULL);

                            for (int i = 0; i < ntdim-1; i++)
                            {

                                    double* val;
                                    int* ia;
                                    int* ja;
                                    int nnz;
                                    dfdu(uMid,s,nsdim,i,val,ia,ja,nnz);

                            for (int k = 0; k < nsdim; k++)
                                        {
                                                MatSetValue(B, i*nsdim+k, i*nsdim+k, 0, INSERT_VALUES);
                                                MatSetValue(B, i*nsdim+k, (i+1)*nsdim+k, 0, INSERT_VALUES);
                                        }
                            for( int k = 0; k < nnz; k++)
                                    {
                                            MatSetValue(B, i*nsdim+ia[k], i*nsdim+ja[k], val[k], INSERT_VALUES);
                                            MatSetValue(B, i*nsdim+ia[k], (i+1)*nsdim+ja[k], val[k], INSERT_VALUES);
                                    }
                            }
                            MatAssemblyBegin(B, MAT_FLUSH_ASSEMBLY);
                MatAssemblyEnd(B, MAT_FLUSH_ASSEMBLY);


                            for (int i = 0; i < ntdim-1; i++)
                            {

                                    for (int k = 0; k < nsdim; k++)
                                    {
                                            MatSetValue(B, i*nsdim+k, i*nsdim+k, -1/dt[i], ADD_VALUES);
                                            MatSetValue(B, i*nsdim+k, (i+1)*nsdim+k, 1/dt[i], ADD_VALUES);
                                    }

                            }
    
                            MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
                            MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);




                            MatCreate(MPI_COMM_WORLD, &E);
                            MatSetType(E, MATSEQAIJ);
                            MatSetSizes(E, PETSC_DECIDE, PETSC_DECIDE, (ntdim-1)*nsdim, (ntdim-1));
                            MatSetFromOptions(E);
                            MatSetUp(E);
                            MatSeqAIJSetPreallocation(E,1,NULL);
                            for (int i = 0; i < (ntdim-1); i++)
                            {
                                    for (int j = 0; j < nsdim; j++)
                                    {
                                            MatSetValue(E, i*nsdim+j, i, dudt[i*nsdim+j], INSERT_VALUES);
                                    }
                            }

                            MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY);
                            MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY);
    

                MatCreate(MPI_COMM_WORLD, &wEinv);
                               MatSetType(wEinv, MATSEQAIJ);
                            MatSetSizes(wEinv, PETSC_DECIDE, PETSC_DECIDE, (ntdim-1), (ntdim-1));
                            MatSetFromOptions(wEinv);
                            MatSetUp(wEinv);
                            MatSeqAIJSetPreallocation(wEinv,1,NULL);



                            for (int i = 0; i< (ntdim-1); i++)
                            {
                                    MatSetValue(wEinv, i, i, (t[ntdim-1]-t[0])/(dt[i]*pow(alpha,2)), INSERT_VALUES);
                }

                            MatAssemblyBegin(wEinv, MAT_FINAL_ASSEMBLY);
                            MatAssemblyEnd(wEinv, MAT_FINAL_ASSEMBLY);

                            MatCreate(MPI_COMM_WORLD, &wBinv);
                            MatSetType(wBinv, MATSEQAIJ);
                            MatSetSizes(wBinv, PETSC_DECIDE, PETSC_DECIDE, nsdim*ntdim, nsdim*ntdim);
                            MatSetFromOptions(wBinv);
                            MatSetUp(wBinv);
                            MatSeqAIJSetPreallocation(wBinv,1,NULL);

                            for (int j = 0; j < nsdim; j++)
                            {
                                    MatSetValue(wBinv, j,j, (t[ntdim-1]-t[0])/(dt[0]/2), INSERT_VALUES);
                                    for (int i = 1; i< (ntdim-1); i++)
                                    {
                                            MatSetValue(wBinv, i*nsdim+j, i*nsdim+j, (t[ntdim-1]-t[0])/((dt[i]+dt[i-1])/2), INSERT_VALUES);
                                    }
                                    MatSetValue(wBinv, (ntdim-1)*nsdim+j, (ntdim-1)*nsdim+j, (t[ntdim-1]-t[0])/(dt[ntdim-2]/2), INSERT_VALUES);
                            }

                            MatAssemblyBegin(wBinv, MAT_FINAL_ASSEMBLY);
                            MatAssemblyEnd(wBinv, MAT_FINAL_ASSEMBLY);

                            Mat temp1, temp2, temp3;
    
                            std::cout << "Perform Matrix Multiplication" << std::endl;
                            MatMatMult(B, wBinv, MAT_INITIAL_MATRIX, 1., &temp1);
                            MatMatTransposeMult(temp1,B, MAT_INITIAL_MATRIX, 2., &temp2);
                            MatMatMult(E, wEinv, MAT_INITIAL_MATRIX, 1., &temp3);
                            MatMatTransposeMult(temp3, E, MAT_INITIAL_MATRIX, 2., &Smat);

                            std::cout << "Performing Matrix Addition" << std::endl;
                            MatAXPY(Smat, 1., temp2, DIFFERENT_NONZERO_PATTERN);

                MatSetOption(Smat, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
                MatAssemblyBegin(Smat, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(Smat, MAT_FINAL_ASSEMBLY);

                            MatDestroy(&temp1);
                            MatDestroy(&temp2);
                            MatDestroy(&temp3);




                std::cout << "Computed Smat" << std::endl;


		std::clock_t start;
		double duration;
		start = std::clock();


                KSP solver;
                PC pc;
                KSPConvergedReason reason;
                PetscInt niter = -1;

                KSPCreate(MPI_COMM_WORLD, &solver);
        
                KSPSetType(solver, KSPCG);    
                KSPGetPC(solver,&pc);
                PCSetType(pc,PCJACOBI);
                
                KSPSetOperators(solver, Smat, Smat, DIFFERENT_NONZERO_PATTERN);
                KSPMonitorSet(solver, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
                KSPSetTolerances(solver,1.e-10, 1.e-10, 1.e2, 100);
                KSPSetUp(solver);
                std::cout << "Initiate Solver" << std::endl;


                VecDuplicate(b,&w);
                VecAssemblyBegin(w);
                VecAssemblyEnd(w);


                KSPSolve(solver, b, w);

                KSPGetConvergedReason(solver, &reason);
                KSPGetIterationNumber(solver, &niter);

                if (reason >=0)
                {
                    PetscPrintf(MPI_COMM_WORLD, "Solver converged with %d iterations \n", niter);
                }
                else
                {
                    PetscPrintf(MPI_COMM_WORLD, "Solver did not converge! \n");
                }

                KSPDestroy(&solver);
                MatDestroy(&Smat);
		duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		std::cout<<"timer: "<< duration << std::endl;


                VecDestroy(&b);
                



                VecCreate(MPI_COMM_WORLD, &temp4);
                VecSetSizes(temp4,PETSC_DECIDE,ntdim*nsdim);
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


                VecCreate(MPI_COMM_WORLD, &temp4);
                VecSetSizes(temp4, PETSC_DECIDE, ntdim-1);
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

                for (int i = 0; i < nsdim*ntdim; i++)
                {
                    VecGetValues(v, 1, &i, &y);
                    u[i] = u[i] - y;
                }
                for (int i = 0; i < ntdim-1; i++)    
                {
                    VecGetValues(eta,1,&i,&y);
                    dt[i] = dt[i] * exp(y);
                    t[i+1] = t[i] + dt[i];
                }



                VecDestroy(&eta);
                VecDestroy(&v);

                for (int i = 0; i<ntdim - 1; i++)
                {
                    for (int j = 0; j< nsdim; j++)
                    { 
                        uMid[i*nsdim+j] = 0.5*(u[(i+1)*nsdim+j] + u[i*nsdim+j]);
                    }
                }



                for (int i = 0; i<ntdim-1; i++)
                {
                    for (int j = 0; j<nsdim; j++)
                    {
                        dudt[i*nsdim+j]  = (u[(i+1)*nsdim+j] - u[i*nsdim+j])/dt[i];
                    }
                }


                VecCreate(MPI_COMM_WORLD, &b);
                VecSetSizes(b, PETSC_DECIDE, nsdim*(ntdim-1));
                VecSetFromOptions(b);

                funcval = F(uMid, s, ntdim, nsdim);

                for (int i=0; i<((ntdim-1)*nsdim); i++)
                {
                    VecSetValue(b, i, dudt[i] - funcval[i], INSERT_VALUES);
                }
                VecAssemblyBegin(b);
                VecAssemblyEnd(b);

                PetscReal norm_b;
                VecNorm(b,NORM_2, &norm_b);


                std::cout << "Norm of b is " << norm_b << std::endl;

                if ((norm_b < atol) || (norm_b < rtol*norm_b0))
                {
                    return u;
                }




            }

            std::cout << "Failed in iteration " << maxIter << std::endl;
            return u;

        };


};


#endif
