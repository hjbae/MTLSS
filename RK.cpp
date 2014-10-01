#include "burgers2D.hpp"
#include "stdlib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

double* RK(double* u0,double* dt, int ntdim, int nsdim)
{
  double* u = new double[nsdim];
  for (int j = 0; j < nsdim; j++)
  {
    u[j] = u0[j];
  }
  for (int i = 0; i < ntdim; i++)
  {
    double* k1 = F(u,1,2,nsdim);
    double* u_int = new double[nsdim];
    for (int j = 0; j < nsdim; j++)
    {
      u_int[j] = u[j] + 0.5 * dt[i] * k1[j];
    }
    double* k2 = F(u_int,1,2,nsdim);
    for (int j = 0; j< nsdim; j++)
    {
      u_int[j] = u[j] + 0.5* dt[i] * k2[j];
    }
    double* k3 = F(u_int,1,2,nsdim);
    for (int j = 0; j < nsdim; j++)
    {
      u_int[j] = u[j] + dt[i] * k3[j];
    }
    double* k4 = F(u_int,1,2,nsdim);
    for (int j = 0; j < nsdim; j++)
    {
      u[j] = u[j] + dt[i]*(1./6.*k1[j] + 1./3.*k2[j] + 1./3.*k3[j] + 1./6.*k4[j]);
    }
    delete[] u_int;
  }
  return u;
}
