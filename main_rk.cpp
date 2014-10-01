#include "RK.hpp"
#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

int main(int argc, char** args)
{


    int ntdim = 180;
    double* dt = new double[ntdim];

    std::string line;
    std::ifstream mytime ("dt.txt");

    double a;
    int linenum = -1;
    if(mytime.is_open())
    {
      while( linenum < ntdim)
      {
        getline(mytime,line);
        linenum ++;
        std::istringstream iss(line);
        iss >> a;
        dt[linenum] = a;
      }
    }
   

    int nsdim = 2319;
    double* u0 = new double[nsdim];    

    linenum = -1;
    std::ifstream myfile ("uf_burger_61_025.txt");

    if(myfile.is_open())
    {

        while(linenum < 0)
        {
            getline(myfile,line);
            linenum ++;
            std::istringstream iss(line);
            for (int i = 0; i<nsdim; i++)
            {
                iss >> a;
                u0[i] = a;
            }
            
        }
    }

    double* uf = RK(u0,dt,ntdim,nsdim);


    std::ofstream outfile ("uf_rk.txt");
    outfile.precision(15);
    if(outfile.is_open())
    {
            for (int j = 0; j < nsdim; j++)
            {
                outfile << uf[j] << " ";
            }
            outfile <<  std::endl;
       
    }

    return 0;
}
