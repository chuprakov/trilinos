#include "MatlabReader.h"
#include "TSF.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace TSF;
using namespace std;

int main(int argc, void** argv)
{
  TSFMPI::init(&argc, &argv);

  try
    {
      int n = 5;
      TSFVectorType epetra = new PetraVectorType();
      TSFVectorSpace space = epetra.createSpace(n);
      TSFLinearOperator A = epetra.createMatrix(space, space);
      TSFMatrixView M = A;
      
      vector<int> ncols(n);
      vector<vector<int> > colIndices(n);
      vector<const int*> colIndicesPtrs(n);
      vector<vector<double> > Avals(n);

      for (int row=0; row<n; row++)
        {
          if (row==0)
            {
              ncols[row] = 2;
              colIndices[row].resize(2);
              colIndices[row][0] = row;
              colIndices[row][1] = row+1;
              Avals[row].resize(2);
              Avals[row][0] = 2.0;
              Avals[row][1] = -1.0;
            }
          else if (row==(n-1))
            {
              ncols[row] = 2;
              colIndices[row].resize(2);
              colIndices[row][0] = row;
              colIndices[row][1] = row-1;
              Avals[row].resize(2);
              Avals[row][0] = 2.0;
              Avals[row][1] = -1.0;
            }
          else
            {
              ncols[row] = 3;
              colIndices[row].resize(3);
              Avals[row].resize(3);
              colIndices[row][0] = row-1;
              colIndices[row][1] = row;
              colIndices[row][2] = row+1;
              Avals[row][0] = -1.0;
              Avals[row][1] = 2.0;
              Avals[row][2] = -1.0;
            }
          colIndicesPtrs[row] = &(colIndices[row][0]);
        }

      M.setGraph(n, &(ncols[0]), &(colIndicesPtrs[0]));

      M.freezeStructure();

      M.zero();

      for (int row=0; row<n; row++)
        {
          M.addToRow(row, ncols[row], &(colIndices[row][0]), &(Avals[row][0]));
        }

      M.freezeValues();
      cerr << "A = " << endl << A << endl;

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
      double convTol = 1.0e-12;
      int i = 0;
      while (r.norm2() > convTol)
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

      double tol = 1.0e-10;
      if(error < tol)
        {
          cerr << " PASSED" << endl;
        }
      else
        {
          cerr << " FAILED error>tolerance" << endl;
        }
    }
  catch(exception& e)
    {
      cerr << e.what() << endl;
      cerr << " EXCEPTION" << endl;
    }
  

  TSFMPI::finalize();
}
