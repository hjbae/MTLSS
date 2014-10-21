#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <mpi.h>

int main(int argc, char** args)
{

    int ntdim = 41;
    int nsdim = 2319;

    MPI_Init(&argc,&args);

    double* u0 = new double[ntdim*nsdim];
    std::string line;             
    std::ifstream myfile ("u0_burger_81.txt");

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



    int localsize[2] = {ntdim,nsdim};
    int globalsize[2] = {ntdim,nsdim};
    int start[2] = {0,0}; 
    MPI_File Myfile;
    MPI_Status status;
    MPI_Datatype filetype;

    MPI_Type_create_subarray(2, globalsize,localsize, start, MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(MPI_COMM_WORLD,"uf_burger_61_025.txt",MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&Myfile);
    MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
    MPI_File_write_all(Myfile,u0,nsdim*ntdim,MPI_DOUBLE,&status);
    MPI_File_close(&Myfile);

    MPI_Finalize();

    return 0;
}
