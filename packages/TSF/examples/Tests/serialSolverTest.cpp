#include "TSF.h"
#include <iostream>
#include <fstream>

using namespace TSF;


class PoissonTest1D : public BVP1D
{
public:
  PoissonTest1D(int n, TSFMatrixOperator* matrix)
    : BVP1D(n, -1.0, 1.0, matrix) {;}
  virtual ~PoissonTest1D(){;}
  
  TSFReal gamma(TSFReal x) const {return -6.0*x;}
  
  TSFReal solution(TSFReal x) const {return x*(1.0-x*x);}
};


int main(int argc, void** argv)
{
  TSFMPI::init(&argc, &argv);
  
  ofstream fout("serialSolverTest.csv");

  try
    {
      char matrixFile [] = "matrices/fourbyfour.mtx";
      int n = 4;

      TSFVectorType petra = new PetraVectorType();
      TSFVectorSpace space = petra.createSpace(n);

      MatrixMarketReader reader(matrixFile);
      TSFLinearOperator A = reader.read(petra);

      TSFVector soln = space.createMember();
      soln.randomize();
      TSFVector b = A * soln;

      TSFLinearSolver solver = new AZTECSolver();
      TSFLinearOperator aInv = A.inverse(solver);
      TSFVector x = aInv * b;

      cerr << "A: " << endl;
      A.getMatrix()->print(cout);
      cerr << "RHS: " << endl << b << endl;
      cerr << "Solution: " << endl << soln << endl;
      cerr << "Result: " << endl << x << endl;

      TSFReal error = (x-soln).norm2();
      cerr << "error norm = " << error << endl;

      TSFReal tol = 1.0e-10;
      if(error < tol)
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






