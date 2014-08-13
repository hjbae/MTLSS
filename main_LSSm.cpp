#include "LSSm.hpp"
#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

int main(int argc, char** args)
{


    int ntdim = 81;
    double* t = new double[ntdim];
    for (int i = 0; i<ntdim; i++)
    {
        t[i] = 0.1*i;
    }


    int nsdim = 2244;
    double* u0 = new double[ntdim*nsdim];    


    std::string line;             
    std::ifstream myfile ("vortex_101.txt");

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
    int c = 1568;
    int f = nsdim - c;

    int* mapcoarse = new int[c];
    int* mapfine = new int[f];

    int pos = 0;
    int posf = 0;
    for (int ii = 0; ii < 4; ii++)
    {
        for (int i = 0; i< 12*21; i++)
        {
            mapcoarse[pos] = ii*nsdim/4+i;
            pos ++;
        }
        for (int i = 0; i <7; i++)
        {
            for (int j = 0; j < 12; j ++)
            {
                mapcoarse[pos] = ii*nsdim/4+12*21+i*40+j;
                pos ++;
            }
            for (int j = 12; j < 25; j++)
            {
                mapfine[posf] = ii*nsdim/4+12*21+i*40+j;
                posf ++;
            }
            for (int j = 25; j <27; j++)
            {
                mapcoarse[pos] = ii*nsdim/4+12*21+i*40+j;
                pos ++;
            }
            if (i != 6)
            {
                for (int j = 27; j < 40; j ++)
                {
                    mapfine[posf] = ii*nsdim/4+12*21+i*40+j;
                    posf ++;
                }
            }
        }
        for (int i = 21*12+40*6+27; i < 21*12+40*6+27+21*2; i++)
        {
            mapcoarse[pos] = ii*nsdim/4+i;
            pos ++;
        }
    }

    lssSolver A(t,u0,k,2,c,f,mapcoarse,mapfine,ntdim,nsdim,alpha);

    PetscInitialize(&argc,&args,(char*)0,NULL);
    double* uf = A.lss(20);


    std::ofstream outfile ("uf_81.txt");

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
