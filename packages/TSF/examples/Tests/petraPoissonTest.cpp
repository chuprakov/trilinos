#include "BICGSTABSolver.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearSolver.h"

#include "PetraMatrix.h"
#include "LAPACKGeneralMatrix.h"
#include "TSFLinearProblem.h"
#include "TSFMPI.h"
#include "TSFOut.h"
#include "BVP1D.h"
#include "ILUKPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "TSFVectorType.h"
#include "PetraVectorType.h"
#include "DenseSerialVectorType.h"
#include "AZTECSolver.h"
#include "TSFDeferredLinearCombination.h"

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
  TSFMPI::init(argc, argv);
  
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
      
			TSFPreconditionerFactory prec = new ILUKPreconditionerFactory(1);
			TSFLinearSolver solver = new BICGSTABSolver(prec, 1.0e-10, 100);
      TSFLinearOperator At = A.adjoint();
			
			

			TSFLinearOperator aInv = At.inverse(solver);
      
      TSFVector x;
      aInv.apply(b, x);
      
      cerr << "RHS: " << endl << b << endl;
      cerr << "Solution: " << endl << soln << endl;
      cerr << "Result: " << endl << x << endl;
      
      TSFReal error = (x-soln).norm2();
      cerr << "error norm = " << error << endl;
    }
  catch(exception& e)
    {
      cerr << e.what() << endl;
    }
  TSFMPI::finalize();
}






