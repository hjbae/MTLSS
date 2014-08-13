#include "Lss2.cpp"
#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include "eulervortex.hpp"

int main(int argc, char** args)
{



    int ntdim = 31;
    double* t = (double*)malloc(ntdim*sizeof(double));
    for (int i = 0; i<ntdim; i++)
    {
        t[i] = 0.1*i;
    }


    int nsdim = 2244;
    double* u0 = (double*)malloc(ntdim*nsdim*sizeof(double));    


    std::string line;             
    std::ifstream myfile ("vortex_51.txt");

    double a;
    int linenum = -1;

    if(myfile.is_open())
    {

        while(! myfile.eof())
        {
            getline(myfile,line);
            linenum ++;
            if (linenum >= nsdim)
            { break; }
            std::istringstream iss(line);
            for (int i = 0; i<ntdim; i++)
            {
                iss >> a;
                u0[i*nsdim+linenum] = a;
            }
            
        }
    }

    double* s;
    
    double alpha = 10;


    lssSolver A(t,u0,s,ntdim,nsdim,alpha);


    MPI_Init(&argc,&args);
    PetscInitialize(&argc,&args,(char*)0,NULL);
    double* uf = A.lss();


    std::ofstream outfile ("uf_31.txt");

    if(outfile.is_open())
    {
        for(int i = 0 ; i < ntdim; i++)
        {
            for (int j = 0; j < nsdim; j++)
            {
                outfile << uf[i*nsdim+j] << ", ";
            }
            outfile <<  std::endl;
        }
    }


    PetscFinalize();
    MPI_Finalize();

    return 0;
}
