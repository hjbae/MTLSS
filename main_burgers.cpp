#include "LSSm.hpp"
#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

int main(int argc, char** args)
{


    int ntdim = 61;
    double* t = new double[ntdim];
    for (int i = 0; i<ntdim; i++)
    {
        t[i] = 0.05*i;
    }


    int nsdim = 2319;
    double* u0 = new double[ntdim*nsdim];    


    std::string line;             
    std::ifstream myfile ("burgers_81.txt");

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

    double* k;
    
    double alpha = 10;
    int c = 294;
    int f = nsdim - c;

    int* mapcoarse = new int[c];
    int* mapfine = new int[f];

    int pos = 0;
    for (int i = 0; i< 7*21; i++)
    {
        mapcoarse[i] = pos;
        pos ++;
    }
    for (int i = 0; i < 81*25; i ++)
    {
        mapfine[i] = pos;
        pos ++;
    }
    for (int i = 0; i < 7*21; i++)
    {
        mapcoarse[7*21+i] = pos;
        pos ++;
    }

    lssSolver A(t,u0,k,4,c,f,mapcoarse,mapfine,ntdim,nsdim,alpha);
    
    PetscInitialize(&argc,&args,(char*)0,NULL);
    double* uf = A.lss();


    std::ofstream outfile ("uf_burger_61_025.txt");

    if(outfile.is_open())
    {
        for(int i = 0 ; i < ntdim; i++)
        {
            for (int j = 0; j < nsdim; j++)
            {
                outfile << uf[i*nsdim+j] << " ";
            }
            outfile <<  std::endl;
        }
    }

    return 0;
}
