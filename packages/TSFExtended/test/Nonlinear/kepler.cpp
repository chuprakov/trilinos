//@HEADER
// ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "Teuchos_MPISession.hpp"
#include "TSFVector.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFNonlinearOperator.hpp"
#include "Teuchos_Time.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFBICGSTABSolver.hpp"
#include <cmath>

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;


namespace TSFExtended
{
  class Kepler : public NonlinearOperatorBase<double>
  {
  public:
    /** */
    Kepler(const double& e, const double& m, const VectorType<double>& type)
      : NonlinearOperatorBase<double>(type.createSpace(1, 1, &(tuple(0)[0])), 
                                      type.createSpace(1, 1, &(tuple(0)[0]))),
        e_(e), m_(m), type_(type)
    {;}

    /** */
    void setM(double m) {m_ = m;}

    /** */
    void apply(const Vector<double>& in, Vector<double>& out) const 
    {
      const double& xIn = in.getElement(0);
      double xOut = xIn + e_*::sin(xIn) - m_;
      out.setElement(0, xOut);
    }
    /** */
    RefCountPtr<TSFCore::LinearOp<double> > jacobian(const Vector<double>& in) const 
    {
      /* create a new operator */
      LinearOperator<double> J = type_.createMatrix(domain(), range());
      /* get a "view" of a loadable matrix underneath the operator */
      RefCountPtr<LoadableMatrix<double> > matview = J.matrix();

      /* compute the derivative of the residual */
      const double& xIn = in.getElement(0);
      double jVal = 1 + e_*::cos(xIn);

      /* insert the derivative into the (0,0) element of the matrix */
      Array<int> colIndices = tuple(0);
      Array<double> colValues = tuple(jVal);
      matview->setRowValues(0, colIndices.size(), &(colIndices[0]),
                            &(colValues[0]));
      matview->freezeValues();
    
      return J.ptr();
    }

    /* */
    GET_RCP(NonlinearOperatorBase<double>);

  private:
    /* eccentricity of orbit */
    double e_;
    /* mean anomaly */
    double m_;
    /* vector type */
    VectorType<double> type_;
  };
}


int main(int argc, void *argv[]) 
{
  try
    {
      int verbosity = 2;

      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      /* a number of occasional interest */
      const double pi = 4.0*atan(1.0);

      /* eccentricity */
      double e = 0.1;

      Kepler* kepler = new Kepler(e, 0.0, type);
      NonlinearOperator<double> F  = kepler;
      
      
      Vector<double> x0 = F.domain().createMember();
      x0.setElement(0, 0.0);

      Vector<double> resid =  F.range().createMember();


      ParameterList solverParams;

      solverParams.set(LinearSolverBase<double>::verbosityParam(), 2);
      solverParams.set(IterativeSolver<double>::maxitersParam(), 100);
      solverParams.set(IterativeSolver<double>::tolParam(), 1.0e-14);

      LinearSolver<double> solver = new BICGSTABSolver<double>(solverParams);
      
      /* simple Newton solver with no line search. This should work for 
       * small values of the eccentricity */
      for (double m=0.0; m<=2.0*pi; m+=pi/10.0)
        {
          kepler->setM(m);
          int maxiters = 20;
          int iters = 0;
          double tol = 1.0e-12;
          Vector<double> step;
          while(iters < maxiters)
            {
              F.apply(x0, resid);
              LinearOperator<double> J = F.jacobian(x0);
              solver.solve(J, resid, step);
              cerr << "iter=" << iters << " " << " |step|=" << step.norm2() 
                   << " |F|=" << resid.norm2() << endl;
              x0 = x0 - step;
              if (step.norm2() < tol) break;
            }
          if (iters >= maxiters)
            {
              cerr << "Newton's method failed to converge!" << endl;
            }
          else
            {
              cerr << "Converged! Yippee!" << endl;
              cerr << "m= " << m << " x=" << x0.getElement(0) << endl;
            }
        }
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();

}
