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
#include "TSFEpetraVectorSpace.hpp"
//#include "TSFCoreEpetraVectorSpace.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
//#include "TSFLinearSolver.hpp"
#include "TSFBICGSTABSolver.hpp"
#include "TSFProductVectorSpace.hpp"
#include "TSFTransposeOperator.hpp"
//#include "TSFInverseOperator.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFEpetraMatrix.hpp"
//#include "TSFCoreLinearOp.hpp"
#include "TSFIdentityOperator.hpp"
#include "TSFZeroOperator.hpp"
#include "TSFBlockOperator.hpp"
#include "TSFDiagonalOperator.hpp"
#include "TSFScaledOperator.hpp"
#include "TSFSumOperator.hpp"
#include "TSFComposedOperator.hpp"


//using namespace Teuchos;
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

      /* create two vectors  */
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

      int nLocalRows2 = 25;

      //int rank = MPISession::getRank();
      //int nProc2 = MPISession::getNProc();

      int spaceDimension2 = nProc*nLocalRows2;
      int lowRow2 = nLocalRows2*rank;

      std::vector<int> localRows2(nLocalRows);
      for (int i=0; i<nLocalRows2; i++)
        {
          localRows2[i] = lowRow2 + i;
        }


      VectorSpace<double> space = type.createSpace(spaceDimension, 
                                                   nLocalRows, 
                                                   &(localRows[0]));
      cerr << "Created space" << endl << space.describe() << endl;
      VectorSpace<double> space2 = 
	type.createSpace(spaceDimension2,nLocalRows2, 
			 &(localRows2[0]));


      /* Create operator */


      LinearOperator<double> A = type.createMatrix(space, space2);
      cerr << "Created Epetra Operator using new\n";
      cerr << A.describe() << endl;
      



      LinearOperator<double> AZZ = new ZeroOperator<double>(space, space2);
      cerr << "Created ZeroOperator using new\n";
      cerr << AZZ.describe() << endl;
      
      

      LinearOperator<double> I = new IdentityOperator<double>(space);
      cerr << "Created IdentityOperator" << endl; 
      cerr << I.describe() << endl;
      
      
      LinearOperator<double> AT = A.transpose();
      LinearOperator<double> ATT = new TransposeOperator<double>(AT);

      cerr << "Created TransposeOperator both ways" << endl;
      cerr << AT.describe() << endl;

      

      /* create vectors and initialize */
      Vector<double> vec = space.createMember();
      Vector<double> vec2 = space2.createMember();
      cerr << "Created two vectors\n";
      cerr << vec.describe() << endl;
      

      /* trivial to debug*/
      vec.zero();
      vec2.zero();

      VectorSpace<double> pvs = 
	new ProductVectorSpace<double>(space, space2);
      cerr << "Created a product vector space\n";
      cerr << pvs.describe() << endl;

      cerr << "Creating ProductVector" << endl;
      Vector<double> pv = pvs.createMember();
      cerr << pv.describe() << endl;
      
      pv.setBlock(0, vec);
      pv.setBlock(1, vec2);
      cerr << "Set the blocks of the pv\n";
      cerr << pv.describe() << endl;


      cerr << "Setting up block Operator" << endl;
      LinearOperator<double> B = new BlockOperator<double>(pvs, pvs);
      cerr << "B = " << B.describe() << endl;

      cerr << "Getting nBlockRows = " << B.numBlockRows() << endl;


      cerr << "Setting up the blocks" << endl;
      B.setBlock(0, 0, I);

      
      LinearOperator<double> I2 = new IdentityOperator<double>(space2);
      B.setBlock(1, 1, I2);
     
      B.finalize(true);
      cerr << "B set up and finalized" << endl;
      cerr << B.describe() << endl;


      vec.setToConstant(1.0);
      vec2.setToConstant(1.0);
      /* Let's do a matrix-vect mult  */
      Vector<double> r = B*pv;
      cerr << "Did a mv mult; norm of r = " << pow(r.norm2(), 2) << endl;

      cerr << endl << "Trying a getRow " << endl;
      Teuchos::Array<int> indices;
      Teuchos::Array<double> values;
      cerr << "calling" << endl;
      B.getRow(55, indices, values);
      for (int i = 0; i < indices.size(); i++)
	{
	  cerr << "   incex = " << indices[i] << " value = " << values[i] 
	       << endl;
	}

      cerr << endl << "Setting up a diagonal matrix and calling getRow" 
	   << endl;
      LinearOperator<double> D = new DiagonalOperator<double>(vec);
      cerr << "   " << vec.describe() << endl;

      D.getRow(5, indices, values);
      for (int i = 0; i < indices.size(); i++)
	{
	  cerr << "   incex = " << indices[i] << " value = " << values[i] 
	       << endl;
	}
      
      cerr << endl << "Setting up a scaled matrix and calling getRow" 
	   << endl;
      LinearOperator<double> Sc = new ScaledOperator<double>(D, 3.14);

      Sc.getRow(5, indices, values);
      for (int i = 0; i < indices.size(); i++)
	{
	  cerr << "   incex = " << indices[i] << " value = " << values[i] 
	       << endl;
	}
      

      /* Trying SumOperator  */
      cerr << endl << "Creating a SumOperator and calling getRow" << endl;
      LinearOperator<double> SumOp = new SumOperator<double>(Sc, D);
      SumOp.getRow(5, indices, values);
      for (int i = 0; i < indices.size(); i++)
	{
	  cerr << "   incex = " << indices[i] << " value = " << values[i] 
	       << endl;
	}


      /* Trying ComposedOperator  */
      cerr << endl << "Creating a ComposedOperator and calling getRow" << endl;
      LinearOperator<double> ComOp = new ComposedOperator<double>(Sc, D);
      ComOp.getRow(5, indices, values);
      for (int i = 0; i < indices.size(); i++)
	{
	  cerr << "   incex = " << indices[i] << " value = " << values[i] 
	       << endl;
	}



      /* Trying to form a matrix  */
      cerr << endl << "Forming a matrix from a row accessible matrix" 
	   << endl;
      LinearOperator<double> formSumOp = SumOp.form(type);
      cerr << "    Formed it: testing with getRow" << endl;
      formSumOp.getRow(5, indices, values);
      for (int i = 0; i < indices.size(); i++)
	{
	  cerr << "   incex = " << indices[i] << " value = " << values[i] 
	       << endl;
	}



      

      /* Try Block operators  */

//       LinearOperator<double> Empty = new BlockOperator<double>();
//       cerr << "set up empty block Operator" << endl << "Empty = " << 
// 	Empty.describe() << endl << endl;


//       // set some blocks
//       cerr << "NumBlocks = " << Empty.numBlockRows() << endl;
//       cerr << "Setting a block \n";
//       Empty.setBlock(0, 0, I);
      

      



//       /* create an empty operator, sized to map domain to range */
//       LinearOperator<double> A = type.createMatrix(space, space);
      
//       RefCountPtr<LoadableMatrix<double> > mat = A.matrix();

//       /* fill in with the Laplacian operator */
//       for (int i=0; i<nLocalRows; i++)
//         {
//           Array<int> colIndices;
//           Array<double> colVals;
//           if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
//             {
//               colIndices = tuple(localRows[i]);
//               colVals = tuple(1.0);
//             }
//           else
//             {
//               colIndices = tuple(localRows[i]-1, localRows[i], localRows[i]+1);
//               colVals = tuple(-1.0, 2.0, -1.0);
//             }
//           mat->setRowValues(localRows[i], colIndices.size(), 
//                             &(colIndices[0]), &(colVals[0]));
//         }

//       mat->freezeValues();

//       //      cerr << A << endl;

//       Vector<double> x = space.createMember();
//       Vector<double> y;
//       Vector<double> ans = space.createMember();

//       for (int i=0; i<nLocalRows; i++)
//         {
//           x.setElement(localRows[i], localRows[i]*localRows[i]);

//           if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
//             {
//               ans.setElement(localRows[i], localRows[i]*localRows[i]);
//             }
//           else
//             {
//               ans.setElement(localRows[i], -2.0);
//             }
//         }

//       MPIComm::world().synchronize();

//       cerr << "x=" << x << endl;

//       MPIComm::world().synchronize();

//       A.apply(x, y);
//       cerr << "y=" << y << endl;

//       cerr << endl << "error 1-norm = " << (y-ans).norm1() << endl;
//       cerr << endl << "error 2-norm = " << (y-ans).norm2() << endl;
//       cerr << endl << "error inf-norm = " << (y-ans).normInf() << endl;

      cerr << "end of try\n";
      
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  cerr << "About to finalize MPISession\n";
  MPISession::finalize();
  exit(0);
}

