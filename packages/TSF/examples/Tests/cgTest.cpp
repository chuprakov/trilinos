#include "MatlabReader.h"
#include "TSF.h"
#include <iostream>
#include <fstream>

// Michael Boldt
// June 5, 2002

// Solve Ax = b using CG with TSF
int main(int argc, void** argv)
{
  TSFMPI::init(&argc, &argv);

  ofstream fout("cgTest.csv");

  try
    {
      const char matrixFile [] = "matrices/tenbyten.m";
      const double tol = 1.0e-10;

      TSFVectorType vt = new DenseSerialVectorType();

      // Read A from Matrix file
      cout << "here\n";
      TSFMatrixReader reader = new MatlabReader(matrixFile);
      cout << "here\n";
      TSFLinearOperator A = reader.read(vt);
      cout << "here\n";

      // Create a known solution to check computed solution
      TSFVector soln = A.domain().createMember();
      soln.randomize();
      TSFVector b = A * soln;

      // Create and initialize vectors for CG
      TSFVector x = A.domain().createMember();
      x.zero();
      TSFVector r = b.copy();
      TSFVector p = r.copy();

      // Do CG
      int i = 0;
      while (r.norm2() > tol)
	{
	  
	  cerr << "iteration: " << ++i << endl;
	  
	  TSFReal r2 = r * r;
	  TSFReal pAp = p*(A*p);
	  
	  if(TSFUtils::chop(pAp) == 0)
	    {
	      TSFError::raise("failure in CG");
	    }
	  TSFReal alpha = r2 / pAp;

	  x = x + alpha * p;
	  r = r - alpha * A * p;
	  TSFReal beta = (r * r) / r2;
	  p = r + beta * p;
	};
      
      // Print results
      cerr << "RHS: " << endl << b << endl;
      cerr << "Solution: " << endl << soln << endl;
      cerr << "Result: " << endl << x << endl;

      TSFReal error = (x-soln).norm2();
      cerr << "error norm = " << error << endl;

      if(error < tol)
	{
	  fout << " PASSED" << endl;
	}
      else
	{
	  fout << " FAILED error>tolerance" << endl;
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
