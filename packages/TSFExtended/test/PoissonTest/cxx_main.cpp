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
#include "TSFLinearOperator.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFLinearSolverBuilder.hpp"

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;


int main(int argc, void *argv[]) 
{
  try
    {
      int verbosity = 1;

      MPISession::init(&argc, &argv);

      MPIComm::world().synchronize();

      VectorType<double> type = new EpetraVectorType();

      /* create the range space  */
      int nLocalRows = 50;

      int rank = MPISession::getRank();
      int nProc = MPISession::getNProc();

      int spaceDimension = nProc*nLocalRows;
      int lowRow = nLocalRows*rank;

      std::vector<int> localRows(nLocalRows);
      for (int i=0; i<nLocalRows; i++)
        {
          localRows[i] = lowRow + i;
        }

      VectorSpace<double> space = type.createSpace(spaceDimension, 
                                                   nLocalRows, 
                                                   &(localRows[0]));

      /* create an empty operator, sized to map domain to range */
      LinearOperator<double> A = type.createMatrix(space, space);
      
      RefCountPtr<LoadableMatrix<double> > mat = A.matrix();

      /* fill in with the Laplacian operator */
      for (int i=0; i<nLocalRows; i++)
        {
          Array<int> colIndices;
          Array<double> colVals;
          if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
            {
              colIndices = tuple(localRows[i]);
              colVals = tuple(1.0);
            }
          else
            {
              colIndices = tuple(localRows[i]-1, localRows[i], localRows[i]+1);
              colVals = tuple(-1.0, 2.0, -1.0);
            }
          mat->setRowValues(localRows[i], colIndices.size(), 
                            &(colIndices[0]), &(colVals[0]));
        }

      mat->freezeValues();

      //      cerr << A << endl;

      Vector<double> x = space.createMember();
      Vector<double> y;
      Vector<double> ans = space.createMember();

      for (int i=0; i<nLocalRows; i++)
        {
          x.setElement(localRows[i], localRows[i]*localRows[i]);

          if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
            {
              ans.setElement(localRows[i], localRows[i]*localRows[i]);
            }
          else
            {
              ans.setElement(localRows[i], -2.0);
            }
        }

      MPIComm::world().synchronize();

      cerr << "x=" << x << endl;

      MPIComm::world().synchronize();

      A.apply(x, y);
      cerr << "y=" << y << endl;

      cerr << endl << "error 1-norm = " << (y-ans).norm1() << endl;
      cerr << endl << "error 2-norm = " << (y-ans).norm2() << endl;
      cerr << endl << "error inf-norm = " << (y-ans).normInf() << endl;

      ParameterList params;
      ParameterList solverParams;
      solverParams.set("Type", "TSF");
      solverParams.set("Method", "BICGSTAB");
      solverParams.set("Max Iterations", 100);
      solverParams.set("Tolerance", 1.0e-12);
      solverParams.set("Precond", "ILUK");
      solverParams.set("Graph Fill", 1);
      solverParams.set("Verbosity", 4);

      params.set("Linear Solver", solverParams);


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(params);

      SolverState<double> state = solver.solve(A, y, ans);
      
      cerr << state << endl;

      cerr << endl << "solve error 2-norm = " << (x-ans).norm2() << endl;
      
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}

