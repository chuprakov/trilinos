#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_TEUCHOS)
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "Teuchos_BLAS.hpp"

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Creating an instance of the BLAS class for double-precision kernels looks like:
  Teuchos::BLAS<int, double> blas;

  // This instance provides the access to all the BLAS kernels listed in Figure \ref{blas_kernels}:
  const int n = 10;
  double alpha = 2.0;
  double x[ n ];
  for ( int i=0; i<n; i++ ) { x[i] = i; }
  blas.SCAL( n, alpha, x, 1 );
  int max_idx = blas.IAMAX( n, x, 1 );
  cout<< "The index of the maximum magnitude entry of x[] is the "
      <<  max_idx <<"-th and x[ " << max_idx-1 << " ] = "<< x[max_idx-1] 
      << endl;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
       "--enable-teuchos");

  return 0;
}
#endif

