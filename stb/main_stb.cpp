#include "LSSm_st.hpp"
#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

int main(int argc, char** args)
{


    int ntdim = 1001;
    double* t = new double[ntdim];
    for (int i = 0; i<ntdim; i++)
    {
        t[i] = 0.01*i;
    }


    int nsdim = 553;
    double* u0 = new double[ntdim*nsdim];    

    std::string line;             
    std::ifstream myfile ("stb_u.txt");

    double a;
    int linenum = -1;

    if(myfile.is_open())
    {

        while(! myfile.eof())
        {
            getline(myfile,line);
            linenum ++;
            if (linenum >= ntdim)
            { break; }
            std::istringstream iss(line);
            for (int i = 0; i<nsdim; i++)
            {
                iss >> a;
                u0[linenum*nsdim+i] = a;
            }
            
        }
    }



    double* w = new double[ntdim*nsdim];

    std::string line2;
    std::ifstream myfile2 ("stb_w.txt");

    linenum = -1;

    if(myfile2.is_open())
    {

        while(! myfile2.eof())
        {
            getline(myfile2,line2);
            linenum ++;
            if (linenum >= ntdim)
            { break; }
            std::istringstream iss(line2);
            for (int i = 0; i<nsdim; i++)
            {
                iss >> a;
                w[linenum*nsdim+i] = a;
            }

        }
    }



    double alpha = 10;
    int c = 471;
    int f = nsdim-c;

    int* mapcoarse = new int[c];
    int* mapfine = new int[f];

    for (int i = 0; i < 471; i++)
    {
        mapcoarse[i] = i+41;
    }

    for (int i = 0; i < 41; i++)
    {
        mapfine[i] = i;
    }
    for (int i = 0; i < 41; i++)
    {
        mapfine[i+41] = i + 512;
    }

    lssSolver A(t,u0,w,2,c,f,mapcoarse,mapfine,ntdim,nsdim,alpha);
    
    PetscInitialize(&argc,&args,(char*)0,NULL);
    double* uf = A.lss();


    std::ofstream outfile ("uf_stb_25.txt");

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
