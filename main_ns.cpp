#include "LSSm2.hpp"
#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <mpi.h>

int main(int argc, char** args)
{

    //Parameters////////////////////////////////////////////////////////////////


    // initial file
    char* input_file = "uf_burger_61_025.txt";
    // ntdim should be a (power of 2) + 1
    int ntdim = 65;
    // initial dt
    double dt = 3./80.;

    // spacial dimensions
    int nsdim_x = ;
    int nsdim_y = ;
    int nsdim_z = ;
    int nsdim = nsdim_x * nsdim_y * nsdim_z;
    // this is the buffer size for parallization in space
    int bufsize = 1;
    // p is the variable timestep factor difference. This should be a power of 2.
    int p = 2;


    //////////////////////////////////////////////////////////////////////////////

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


    ///////////////////////////////////////////////////////////////////////////////

    int size, rank, ssize, srank, tsize, trank, locntdim, startrow, locnsdim, startcol;

    MPI_Init(&argc,&args);
    // get size and rank of MPI communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // get the size of spatial and temporal parallelization
    tsize = (ntdim-1)/p;
    if (size < tsize)
    {
        tsize = size;
        ssize = 1;
    }
    else{
        ssize = size / tsize;
    }

    trank = rank / ssize;

    if (tsize == 1)
    { locntdim = ntdim;}
    else {
        if (trank != tsize-1)
            locntdim = ntdim/tsize+1;
        else
            locntdim = ntdim - (ntdim/tsize)*(tsize-1);
    }

    double* t = new double[locntdim];
    for (int i = 0; i < locntdim; i++)
    {
        t[i] = (trank*(locntdim-1)+i)*dt;
    }

    // Divide MPI Group
    int ranks[ssize];
    for (int i = 0; i < ssize; i++)
    {
        ranks[i] = trank*ssize + i;
    }

    MPI_Group orig_group, new_group;
    MPI_Comm new_comm;

    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    MPI_Group_incl(orig_group, ssize, ranks, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
    MPI_Group_rank(new_group, &srank);


    NS = nssolver();




    // Read in initial condition
    double* u0 = new double[locntdim*locnsdim];
    int localsize[2] = {locntdim,locnsdim};
    startrow = trank * (locntdim-1);
    startcol = srank * (locnsdim-1);
    int globalsize[2] = {ntdim,nsdim};
    int start[2] = {startrow,startcol}; 
    MPI_File Myfile;
    MPI_Status status;
    MPI_Datatype memtype;
    MPI_Datatype filetype;

    MPI_Type_create_subarray(2, globalsize,localsize, start, MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(MPI_COMM_WORLD,input_file,MPI_MODE_RDONLY,MPI_INFO_NULL,&Myfile);
    MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
    MPI_File_read_all(Myfile,u0,locnsdim*locntdim,MPI_DOUBLE,&status);
    MPI_File_close(&Myfile);

    
    PetscInitialize(&argc,&args,(char*)0,NULL);
    lssSolver A(t,u0,p,c,c,0,f,f,0,mapcoarse,mapfine,locntdim,ntdim,trank,tsize,locnsdim,nsdim,srank,ssize,alpha,new_comm);
    
    double* uf = A.lss();


    MPI_File_open(MPI_COMM_WORLD,"uf_burger_61_025.txt",MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&Myfile);
    MPI_File_set_view(Myfile,0,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
    MPI_File_write_all(Myfile,uf,locnsdim*locntdim,MPI_DOUBLE,&status);
    MPI_File_close(&Myfile);


    return 0;
}
