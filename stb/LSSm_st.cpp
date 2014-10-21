#include <math.h>
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpcasa.h"
#include "petscvec.h"
//#include "burgers2D.hpp"
//#include "eulervortex.hpp"
#include "stochasticburgers.hpp"
#include "stdlib.h"
#include <iostream>
#include "LSSm_st.hpp"
#include <stdio.h>

        lssSolver::lssSolver(double* t, double* u0, double* w, int p, int c, int f, int* mapcoarse, int* mapfine, int ntdim, int nsdim, double alpha)
        {

            this->t = t;
            u = u0;
            this->w = w;
            this->t = t;
            this->p = p;
            this->ntdim = ntdim;
            nn = (ntdim-1)/p+1;
            this->nsdim = nsdim;
            this->f = f;
            this->c = c;

            dt = new double[ntdim-1];
            for (int i =0; i<ntdim-1; i++)
            {
                dt[i] = t[i+1]-t[i];
            }
            
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


            uMidf = new double[nsdim*(ntdim-1)];
            for (int i = 0; i<ntdim-1; i++)
            {
                for (int j = 0; j<nsdim; j++)
                {
                    uMidf[i*nsdim+j] = 0.5*(u[(i+1)*nsdim+j] + u[i*nsdim+j]);
                }
            }
           
 
            uMidc = new double[nsdim*(nn-1)];
            for (int i = 0; i < nn-1; i++)
            {
                for (int j = 0; j < nsdim; j++)
                {
                    uMidc[i*nsdim+j] = 0.5*(u[(i+1)*p*nsdim+j] + u[i*p*nsdim+j]);
                }
            }

            wf = new double[nsdim*(ntdim-1)];
            for (int i = 0; i<ntdim-1; i++)
            {
                for (int j = 0; j<nsdim; j++)
                {
                    wf[i*nsdim+j] = 0.5*(w[(i+1)*nsdim+j] + w[i*nsdim+j]);
                }
            }


            wc = new double[nsdim*(nn-1)];
            for (int i = 0; i < nn-1; i++)
            {
                for (int j = 0; j < nsdim; j++)
                {
                    wc[i*nsdim+j] = 0.5*(w[(i+1)*p*nsdim+j] + w[i*p*nsdim+j]);
                }
            }



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


        Vec lssSolver::APinvx(Vec x)
        {
            Vec y;
            VecDuplicate(x,&y);
            VecAssemblyBegin(y);
            VecAssemblyEnd(y);


            Vec temp1;
            VecCreate(PETSC_COMM_WORLD, &temp1);
            VecSetSizes(temp1,PETSC_DECIDE,(nn-1)*(p+1));
            VecSetFromOptions(temp1);
            VecAssemblyBegin(temp1);
            VecAssemblyEnd(temp1);

            Vec temp2;
            VecCreate(PETSC_COMM_WORLD, &temp2);
            VecDuplicate(temp1, &temp2);
            VecAssemblyBegin(temp2);
            VecAssemblyEnd(temp2);

            Vec temp3;
            VecCreate(PETSC_COMM_WORLD, &temp3);
            VecSetSizes(temp3, PETSC_DECIDE, (nn-1)*(p*f+c));
            VecSetFromOptions(temp3);
            VecAssemblyBegin(temp3);
            VecAssemblyEnd(temp3);
        
            MatMultTranspose(E, x, temp1);
            MatMult(wEinv, temp1, temp2);
            MatMult(E, temp2, temp3);
            MatMult(Pinv,temp3, y);

            VecDestroy(&temp1);
            VecDestroy(&temp2);
            VecDestroy(&temp3);

            VecCreate(PETSC_COMM_WORLD, &temp1);
            VecSetSizes(temp1,PETSC_DECIDE,(nn-1)*(p*f+c)+f+c);
            VecSetFromOptions(temp1);
            VecAssemblyBegin(temp1);
            VecAssemblyEnd(temp1);

            VecCreate(PETSC_COMM_WORLD, &temp2);
            VecDuplicate(temp1, &temp2);
            VecAssemblyBegin(temp2);
            VecAssemblyEnd(temp2);

            VecCreate(PETSC_COMM_WORLD, &temp3);
            VecSetSizes(temp3, PETSC_DECIDE, (nn-1)*(p*f+c));
            VecSetFromOptions(temp3);
            VecAssemblyBegin(temp3);
            VecAssemblyEnd(temp3);


            Vec Ax;
            VecCreate(PETSC_COMM_WORLD, &Ax);
            VecSetSizes(Ax, PETSC_DECIDE, (nn-1)*(p*f+c));
            VecSetFromOptions(Ax);
            VecAssemblyBegin(Ax);
            VecAssemblyEnd(Ax);

            MatMultTranspose(B, x, temp1);
            MatMult(wBinv, temp1, temp2);
            MatMult(B, temp2, temp3);
            MatMult(Pinv, temp3, Ax);


            VecDestroy(&temp1);
            VecDestroy(&temp2);
            VecDestroy(&temp3);

            VecAXPY(Ax, 1, y);

            VecDestroy(&y);
            return Ax;
        }

        void lssSolver::CGSolve(Vec b, Vec* w)
        {

            double Tol = 1.e-7;
            int maxiter = 2500;
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
                        std::cout << "exceeded maxiter" << std::endl;
                    }

                    VecDestroy(&r);
                    VecDestroy(&p);
                    VecDestroy(&Ap);

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

        double* lssSolver::lss(int maxIter,double atol,double rtol)
        {

            Vec b,w,v,eta,temp4;

            VecCreate(PETSC_COMM_WORLD, &b);
            VecSetSizes(b, PETSC_DECIDE, (nn-1)*(p*f+c));
            VecSetFromOptions(b);

            double* funcvalf = F(uMidf, 1, ntdim, nsdim, wf);
            double* funcvalc = F(uMidc, 1, nn, nsdim, wc);

            for (int i=0; i< nn-1; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    VecSetValue(b, i*(p*f+c)+j, dudtc[i*c+j] - funcvalc[i*nsdim+mapcoarse[j]], INSERT_VALUES);
		}
                for (int j = 0; j < p; j++)
                {
                    for (int k = 0; k < f; k++)
                    {
                        VecSetValue(b, i*(p*f+c)+c+j*f+k, dudtf[(i*p+j)*f+k] - funcvalf[(i*p+j)*nsdim+mapfine[k]], INSERT_VALUES);
		   }
                }
            }
            

            VecAssemblyBegin(b);
            VecAssemblyEnd(b);

            PetscReal norm_b0;
            VecNorm(b,NORM_2, &norm_b0);

            std::cout << "Norm of b is " << norm_b0 << std::endl;        


            for (int iNewton=0; iNewton<maxIter; iNewton++)
            {

                MatCreate(PETSC_COMM_WORLD, &B);
                MatSetType(B, MATSEQAIJ);
                MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, (nn-1)*(p*f+c),(nn-1)*(p*f+c)+f+c);
                MatSetFromOptions(B);
                MatSetUp(B);
                MatSeqAIJSetPreallocation(B,2*nsdim,NULL);


		double* val = 0;
		int* ia = 0;
		int* ja = 0;
		int nnz;

                for (int i = 0; i < nn-1; i++)
                {
                    dfdu(uMidc,wc,nsdim,i,val,ia,ja,nnz);
                    for( int k = 0; k < nnz; k++)
                    {
                        if ((coarse[ia[k]]!=-1)&&(coarse[ja[k]]!=-1))
                        {
                            MatSetValue(B, i*(p*f+c)+coarse[ia[k]], i*(p*f+c)+coarse[ja[k]], val[k], INSERT_VALUES);
                            MatSetValue(B, i*(p*f+c)+coarse[ia[k]], (i+1)*(p*f+c)+coarse[ja[k]], val[k], INSERT_VALUES);
                        }
                        else { 
                            if ((coarse[ia[k]]!=-1)&&(fine[ja[k]]!=-1))
                            {
                                MatSetValue(B, i*(p*f+c)+coarse[ia[k]], i*(p*f+c)+c+fine[ja[k]], val[k], INSERT_VALUES);
                                MatSetValue(B, i*(p*f+c)+coarse[ia[k]], (i+1)*(p*f+c)+c+fine[ja[k]], val[k], INSERT_VALUES);
                            }
                        }

                    }

                    dfdu(uMidf,wf,nsdim,i*p,val,ia,ja,nnz);
                    for (int k =0; k<nnz; k++)
                    {
                        if ((fine[ia[k]]!=-1)&&(fine[ja[k]]!=-1))
                        {
                            MatSetValue(B, i*(p*f+c)+c+fine[ia[k]], i*(p*f+c)+c+fine[ja[k]], val[k], INSERT_VALUES);
                            MatSetValue(B, i*(p*f+c)+c+fine[ia[k]], i*(p*f+c)+c+f+fine[ja[k]], val[k], INSERT_VALUES);
                        }
                        
                        else { 
                            if ((fine[ia[k]]!=-1)&&(coarse[ja[k]]!=-1))
                            {
                                MatSetValue(B, i*(p*f+c)+c+fine[ia[k]], i*(p*f+c)+coarse[ja[k]], (2-dt[i*p]/ds[i])*val[k], INSERT_VALUES);
                                MatSetValue(B, i*(p*f+c)+c+fine[ia[k]], (i+1)*(p*f+c)+coarse[ja[k]], dt[i*p]/ds[i]*val[k], INSERT_VALUES);
                            }
                        }
                    }
                    double sub=0;
                    for (int j = 1; j < p; j++)
                    {
                        sub += dt[i*p+j-1];
                        dfdu(uMidf,wf,nsdim,i*p+j,val,ia,ja,nnz);    
                        for (int k = 0; k < nnz; k++)
                        {
                            if ((fine[ia[k]]!=-1)&&(fine[ja[k]]!=-1))
                            {
                                MatSetValue(B, i*(p*f+c)+c+j*f+fine[ia[k]], i*(p*f+c)+c+j*f+fine[ja[k]], val[k], INSERT_VALUES);
                                if (j != p-1)
                                {
                                    MatSetValue(B, i*(p*f+c)+c+j*f+fine[ia[k]], i*(p*f+c)+c+(j+1)*f+fine[ja[k]], val[k], INSERT_VALUES);
                                }
                                else
                                {
                                    MatSetValue(B, i*(p*f+c)+c+j*f+fine[ia[k]], (i+1)*(p*f+c)+c+fine[ja[k]], val[k], INSERT_VALUES);
                                }
                            }

                            else if ((fine[ia[k]]!=-1)&&(coarse[ja[k]]!=-1))
                            {
                                MatSetValue(B, i*(p*f+c)+c+j*f+fine[ia[k]], i*(p*f+c)+coarse[ja[k]], (2-(2*sub+dt[i*p+j])/ds[i])*val[k], INSERT_VALUES);
                                MatSetValue(B, i*(p*f+c)+c+j*f+fine[ia[k]], (i+1)*(p*f+c)+coarse[ja[k]], (2*sub+dt[i*p+j])/ds[i]*val[k], INSERT_VALUES);
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
                        MatSetValue(B, i*(p*f+c)+j, i*(p*f+c)+j, -1/ds[i], ADD_VALUES);
                        MatSetValue(B, i*(p*f+c)+j, (i+1)*(p*f+c)+j, 1/ds[i], ADD_VALUES);
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            MatSetValue(B, i*(p*f+c)+c+j*f+k, i*(p*f+c)+c+j*f+k, -1/dt[i*p+j], ADD_VALUES);
                            if (j != p-1)
                            {
                                MatSetValue(B, i*(p*f+c)+c+j*f+k, i*(p*f+c)+c+(j+1)*f+k, 1/dt[i*p+j], ADD_VALUES);
                            }
                            else {
                                MatSetValue(B, i*(p*f+c)+c+j*f+k, (i+1)*(p*f+c)+c+k, 1/dt[i*p+j], ADD_VALUES);
                            }
                        }
                    }
                }
        
                MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);


                MatCreate(PETSC_COMM_WORLD, &E);
                MatSetType(E, MATSEQAIJ);
                MatSetSizes(E, PETSC_DECIDE, PETSC_DECIDE, (nn-1)*(p*f+c), (nn-1)*(p+1));
                MatSetFromOptions(E);
                MatSetUp(E);
                MatSeqAIJSetPreallocation(E,1,NULL);
                for (int i = 0; i < (nn-1); i++)
                {
                    for (int j = 0; j< c; j++)
                    {
                        MatSetValue(E, i*(p*f+c)+j, i*(p+1), dudtc[i*c+j], INSERT_VALUES);
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            MatSetValue(E, i*(p*f+c)+c+j*f+k, i*(p+1)+1+j, dudtf[i*p*f+j*f+k], INSERT_VALUES);
                        }
                    }
                }

                MatAssemblyBegin(E, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(E, MAT_FINAL_ASSEMBLY);



                MatCreate(PETSC_COMM_WORLD, &wEinv);
                MatSetType(wEinv, MATSEQAIJ);
                MatSetSizes(wEinv, PETSC_DECIDE, PETSC_DECIDE, (nn-1)*(p+1), (nn-1)*(p+1));
                MatSetFromOptions(wEinv);
                MatSetUp(wEinv);
                MatSeqAIJSetPreallocation(wEinv,1,NULL);



                for (int i = 0; i< (nn-1); i++)
                {
                    MatSetValue(wEinv, i*(p+1), i*(p+1), (t[ntdim-1]-t[0])/(ds[i]*pow(alpha,2)),INSERT_VALUES);
                    for (int j = 0; j<p; j++)
                    {
                        MatSetValue(wEinv, i*(p+1)+1+j, i*(p+1)+1+j, (t[ntdim-1]-t[0])/(dt[i*p+j]*pow(alpha,2)), INSERT_VALUES);
                    }
                }
                    
                MatAssemblyBegin(wEinv, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(wEinv, MAT_FINAL_ASSEMBLY);
                

                MatCreate(PETSC_COMM_WORLD, &wBinv);
                MatSetType(wBinv, MATSEQAIJ);
                MatSetSizes(wBinv, PETSC_DECIDE, PETSC_DECIDE, (nn-1)*(p*f+c)+f+c, (nn-1)*(p*f+c)+f+c);
                MatSetFromOptions(wBinv);
                MatSetUp(wBinv);
                MatSeqAIJSetPreallocation(wBinv,1,NULL);
                for (int i = 0; i < nn-1; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        if (i != 0)
                        {    
                            MatSetValue(wBinv, i*(p*f+c)+j,i*(p*f+c)+j, (t[ntdim-1]-t[0])/((ds[i]+ds[i-1])/2), INSERT_VALUES);
                        }
                        else {
                            MatSetValue(wBinv, i*(p*f+c)+j,i*(p*f+c)+j, (t[ntdim-1]-t[0])/(ds[0]/2), INSERT_VALUES);
                        }
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k ++)
                        {
                            if ((i+j)!=0)
                            {
                                MatSetValue(wBinv, i*(p*f+c)+c+j*f+k, i*(p*f+c)+c+j*f+k, (t[ntdim-1]-t[0])/((dt[i*p+j]+dt[i*p+j-1])/2), INSERT_VALUES);
                            }
                            else
                            {
                                MatSetValue(wBinv, i*(p*f+c)+c+j*f+k, i*(p*f+c)+c+j*f+k, (t[ntdim-1]-t[0])/(dt[0]/2), INSERT_VALUES);
                            }
                        }
                    }
                }
                for (int j = 0; j < c; j++)
                {
                    MatSetValue(wBinv, (nn-1)*(p*f+c)+j, (nn-1)*(p*f+c)+j, (t[ntdim-1]-t[0])/(ds[nn-2]/2), INSERT_VALUES);
                }
                for (int j = 0; j < f; j++)
                {
                    MatSetValue(wBinv, (nn-1)*(p*f+c)+c+j, (nn-1)*(p*f+c)+c+j, (t[ntdim-1]-t[0])/(dt[ntdim-2]/2), INSERT_VALUES);
                }

                MatAssemblyBegin(wBinv, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(wBinv, MAT_FINAL_ASSEMBLY);


                VecDuplicate(b,&w);
                VecSet(w, 0.);
                VecAssemblyBegin(w);
                VecAssemblyEnd(w);
            
                Vec b1;

                VecDuplicate(w, &b1);
                VecAssemblyBegin(b1);
                VecAssemblyEnd(b1);

                MatCreate(PETSC_COMM_WORLD, &Pinv);
                MatSetType(Pinv, MATSEQAIJ);
                MatSetSizes(Pinv, PETSC_DECIDE,PETSC_DECIDE, (nn-1)*(p*f+c), (nn-1)*(p*f+c));
                MatSetFromOptions(Pinv);
                MatSetUp(Pinv);
                MatSeqAIJSetPreallocation(Pinv,1,NULL);


                for (int i = 0; i < (nn-1)*(p*f+c); i++)
                {
                    PetscInt ncols;
                    const PetscInt *cols;
                    const PetscScalar *vals;
                    MatGetRow(B,i,&ncols,&cols,&vals);

                    double sum = 0;
                    for (int j = 0; j <ncols; j++)
                    {
                        int blah = cols[j];
                        PetscScalar wBinvval;
                        MatGetValues(wBinv,1,&blah,1,&blah,&wBinvval);
                        sum += vals[j]*vals[j]*wBinvval;
                    }


                    MatRestoreRow(B,i,&ncols,&cols,&vals);

                    MatGetRow(E,i,&ncols,&cols,&vals);


                    for (int j = 0; j <ncols; j++)
                    {
                        int blah = cols[j];
                        PetscScalar wEinvval;
                        MatGetValues(wEinv,1,&blah,1,&blah,&wEinvval);
                        sum += vals[j]*vals[j]*wEinvval;
                    }

                    MatRestoreRow(E, i, &ncols, &cols, &vals);

                    MatSetValue(Pinv, i, i, 1/sum, INSERT_VALUES);
                }

                MatAssemblyBegin(Pinv,MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(Pinv,MAT_FINAL_ASSEMBLY);


                MatMult(Pinv,b,b1);

                CGSolve(b1,&w);


                VecDestroy(&b1);
                VecDestroy(&b);
                MatDestroy(&Pinv); 

                VecCreate(PETSC_COMM_WORLD, &temp4);
                VecSetSizes(temp4,PETSC_DECIDE,(nn-1)*(p*f+c)+f+c);
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
                VecSetSizes(temp4, PETSC_DECIDE, (nn-1)*(p+1));
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
                        int loc = i*(p*f+c)+j;
                        VecGetValues(v, 1 , &loc, &y);
                        u[i*p*nsdim+mapcoarse[j]] = u[i*p*nsdim+mapcoarse[j]]-y;
                    }
                    for (int j = 0; j <p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            int loc = i*(p*f+c)+c+j*f+k;
                            VecGetValues(v, 1, &loc, &y);
                            u[(i*p+j)*nsdim + mapfine[k]] = u[(i*p+j)*nsdim + mapfine[k]] - y;
                        }
                    }
                }
                for (int j = 0; j < c; j++)
                {
                    int loc = (nn-1)*(p*f+c)+j;
                    VecGetValues(v,1,&loc,&y);
                    u[(nn-1)*p*nsdim+mapcoarse[j]] = u[(nn-1)*p*nsdim+mapcoarse[j]] - y;
                }
                for (int j = 0; j < f; j++)
                {
                    int loc = (nn-1)*(p*f+c)+c+j;
                    VecGetValues(v,1,&loc,&y);
                    u[(nn-1)*p*nsdim+mapfine[j]] = u[(nn-1)*p*nsdim+mapfine[j]] - y;
                }

                    
                for (int i = 0; i < nn-1; i++)    
                {
                    for (int j =0; j < p; j++)
                    {
                        int loc1 = i*(p+1)+1+j;
                        int loc2 = i*(p+1);
                        VecGetValues(eta,1,&loc1, &y);
                        dt[i*p+j] *= exp(y);
                        VecGetValues(eta,1,&loc2, &y);
                        dt[i*p+j] *= exp(y*dt[i*p+j]/ds[i]);
                    }
                }


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
                            sub += dt[i*p+j-1];
                            u[(i*p+j)*nsdim+mapcoarse[k]] = (1-sub/ds[i])*u[i*p*nsdim+mapcoarse[k]] + (sub/ds[i])*u[(i+1)*p*nsdim+mapcoarse[k]];
                        }
                    }
                }

                for (int i = 0; i<ntdim - 1; i++)
                {
                    for (int j = 0; j< nsdim; j++)
                    { 
                        uMidf[i*nsdim+j] = 0.5*(u[(i+1)*nsdim+j] + u[i*nsdim+j]);
                    }
                }
    
                for (int i = 0; i < nn-1; i++)
                {
                    for (int j = 0; j < nsdim; j++)
                    {
                        uMidc[i*nsdim+j] = 0.5*(u[(i+1)*p*nsdim+j] + u[i*p*nsdim+j]);
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
                VecSetSizes(b, PETSC_DECIDE, (nn-1)*(p*f+c));
                VecSetFromOptions(b);

                funcvalf = F(uMidf, 1, ntdim, nsdim,wf);
                funcvalc = F(uMidc, 1, nn, nsdim,wc);


                for (int i=0; i< nn-1; i++)
                {
                    for (int j = 0; j < c; j++)
                    {
                        VecSetValue(b, i*(p*f+c)+j, dudtc[i*c+j] - funcvalc[i*nsdim+mapcoarse[j]], INSERT_VALUES);
                    }
                    for (int j = 0; j < p; j++)
                    {
                        for (int k = 0; k < f; k++)
                        {
                            VecSetValue(b, i*(p*f+c)+c+j*f+k, dudtf[(i*p+j)*f+k] - funcvalf[(i*p+j)*nsdim+mapfine[k]], INSERT_VALUES);
                        }
                    }    
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
