#include "TSF.h"
#include "MatlabReader.h"
#include <iostream.h>
#include <fstream.h>

// Find the infinity norm of a vector
double normInf(TSFVector x)
{
  double norm = x[0];
  for(int i=1; i<x.space().dim(); i++)
    {
      if(abs(x[i]) > norm) norm = x[i];
    }
  return norm;
}

// Find an eigenvalue using the power method
int main(int argc, void** argv)
{
  TSFMPI::init(&argc, &argv);

  ofstream fout("powerMethodTest.csv");

  try
    {
      const char matrixFile [] = "matrices/twobytwo.m";
      const double tol = 1.0e-10;

      TSFVectorType vt = new DenseSerialVectorType();

      // Read A from Matrix Market file
      TSFMatrixReader reader = new MatlabReader(matrixFile);
      TSFLinearOperator A = reader.read(vt);

      // Create vectors
      TSFVector x = A.domain().createMember();
      x[0] = 0;
      x[1] = 1;
      TSFVector y = x.copy();
      
      TSFReal yn0;
      TSFReal yn1 = normInf(y);

      do {

	yn0 = yn1;
	TSFVector y = A * x;

	// yn1 = y.normInf();
	yn1 = normInf(y);

	x.scalarMult(1/yn1, y);
      } while(abs(yn1-yn0) > tol + abs(yn1)*tol);
      
      // Print results
      cerr << "x: " << endl << x << endl;
      cerr << "e: " << yn1 << endl;

      if(abs(yn1-yn0) <= tol + abs(yn1)*tol)
	{
	  fout << " PASSED" << endl;
	}
      else
	{
	  fout << " FAILED error>tol" << endl;
	}
    }
  catch(exception& e)
    {
      cerr << e.what() << endl;
      fout << " EXCEPTION" << endl;
    }

  fout.close();

  TSFMPI::finalize();
}


