#include "TSF.h"
#include "BICGSTABSolver.h"
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

  ofstream fout("petraPoissonTest.csv");
  
  try
    {
      int n = 3;
      TSFVectorType petra = new PetraVectorType();
      TSFVectorSpace space = petra.createSpace(n);
      TSFMatrixOperator* mat = petra.createMatrix(space, space);
      
      TSFLinearProblem prob = new PoissonTest1D(space.dim(), mat);
      
      TSFLinearOperator A = prob.getOperator();
      
      TSFVector b = prob.getRHS();
      cerr << "RHS: " << endl << b << endl;
      cerr << "A: " << endl << A << endl;
      
      TSFVector soln = prob.getKnownSolution();
      
      TSFLinearSolver solver = new AZTECSolver();
      // TSFPreconditionerFactory prec = new ILUKPreconditionerFactory(1);
      // TSFLinearSolver solver = new BICGSTABSolver(prec, 1.0e-10, 100);
      
      TSFLinearOperator aInv = A.inverse(solver);

      TSFVector x = aInv*b;

      TSFLinearOperator At = A.adjoint();
			
      cerr << "RHS: " << endl << b << endl;
      cerr << "Solution: " << endl << soln << endl;
      cerr << "Result: " << endl << x << endl;
      
      TSFReal error = (x-soln).norm2();
      cerr << "error norm = " << error << endl;

      double tol = 1.0e-10;
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






