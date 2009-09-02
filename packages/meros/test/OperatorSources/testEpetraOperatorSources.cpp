//@HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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

// Test Meros operator sources given Epetra matrices.


#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_LinearCombinationImpl.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "Meros_LSCOperatorSource.h"
#include "Meros_PCDOperatorSource.h"
#include "Meros_AztecSolveStrategy.hpp"
#include "Meros_InverseOperator.hpp"
#include "Meros_ZeroOperator.hpp"
#include "Meros_IdentityOperator.hpp"
#include "Meros_LinearSolver.hpp"


#include "AztecOO.h"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"

using namespace Teuchos;
using namespace EpetraExt;
using namespace Thyra;
using namespace Meros;


int main(int argc, char *argv[]) 
{
  bool success = false;
  
  GlobalMPISession mpiSession(&argc, &argv);
  typedef Teuchos::ScalarTraits<double> ST;

#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif


  // Get stream that can print to just root or all streams!
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  try
    {
      CommandLineProcessor  clp;
      clp.throwExceptions(false);
      clp.addOutputSetupOptions(true);
      bool verbose = true;
      clp.setOption( "verbose", "quiet", &verbose, 
                     "Determines if any output is printed or not." );

      
      CommandLineProcessor::EParseCommandLineReturn parse_return 
        = clp.parse(argc,argv);
      if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) 
	return parse_return;

      if (!verbose) out = rcp(new FancyOStream(rcp(new oblackholestream())));

      // Make necessary Epetra maps.
      // Need a velocity space map and a pressure space map.
      const Epetra_Map* velocityMap = new Epetra_Map(578, 0, Comm); 
      const Epetra_Map* pressureMap = new Epetra_Map(192, 0, Comm); 
 

      // Read matrix and vector blocks into Epetra_Crs matrices
      Epetra_CrsMatrix* FMatrix(0);
      char * filename = "../../../../../packages/meros/example/data/q1/Aq1.mm";
      MatrixMarketFileToCrsMatrix(filename,
				  *velocityMap, *velocityMap,
				  *velocityMap, *velocityMap,
				  FMatrix);

      // cerr << "If you get Epetra ERROR -1 here, the file was not found"
      //  << "Need to make symbolic link to data directory" 
      //  << endl;

      Epetra_CrsMatrix* BtMatrix(0);
      filename = "../../../../../packages/meros/example/data/q1/Btq1.mm";
      MatrixMarketFileToCrsMatrix(filename,
				  *velocityMap, *pressureMap,
				  *velocityMap, *pressureMap,
				  BtMatrix);
      Epetra_CrsMatrix* BMatrix(0);
      filename = "../../../../../packages/meros/example/data/q1/Bq1.mm";
      MatrixMarketFileToCrsMatrix(filename,
				  *pressureMap, *velocityMap,
				  *pressureMap, *velocityMap,
				  BMatrix);


      // Make an empty C matrix for the 2,2 block
      Epetra_CrsMatrix* CMatrix = new Epetra_CrsMatrix(View, *pressureMap, 1);
      CMatrix->Scale(0.0);
      CMatrix->FillComplete();



      // Wrap Epetra operators into Thyra operators so we can test
      // against them. To do this, we first wrap as Thyra core
      // operators, then convert to the handle layer LinearOperators.
      
      RCP<LinearOpBase<double> >
        tmpF = nonconstEpetraLinearOp(rcp(FMatrix,false));
      const LinearOperator<double> F = tmpF;

      RCP<LinearOpBase<double> >
        tmpBt = nonconstEpetraLinearOp(rcp(BtMatrix,false));
      const LinearOperator<double> Bt = tmpBt;

      RCP<LinearOpBase<double> >
        tmpB = nonconstEpetraLinearOp(rcp(BMatrix,false));
      const LinearOperator<double> B = tmpB;

      RCP<LinearOpBase<double> >
        tmpC = nonconstEpetraLinearOp(rcp(CMatrix,false));
      const LinearOperator<double> C = tmpC;

      const LinearOperator<double> S = block2x2(F, Bt, B, C);

      // Make a random vector for testing.
      Vector<double> svec = S.domain().createMember();
      Thyra::randomize(-1.0, 1.0, svec.ptr().get());



//       // Set up a PCD operator source with these three operators
//       RCP<const PCDOperatorSource> pcdOpSrcRcp 
// 	= rcp(new PCDOperatorSource(FMatrix, BtMatrix, BMatrix, CMatrix));

//       // Get the blocks back from the operator source
//       ConstLinearOperator<double> blockOpOut = pcdOpScrRcp->getSaddleOp();
//       ConstLinearOperator<double> Fout = blockOpOut.getBlock(0,0);
//       ConstLinearOperator<double> Btout = blockOpOut.getBlock(0,1);
//       ConstLinearOperator<double> Bout = blockOpOut.getBlock(1,0);
//       ConstLinearOperator<double> Cout = blockOpOut.getBlock(1,1);

//       RCP<const LinearOpBase<double> > tmpFpOut = pcdOpSrcRcp->getFp();
//       ConstLinearOperator<double> FpOut = tmpFpOut;
//       RCP<const LinearOpBase<double> > tmpApOut = pcdOpSrcRcp->getAp();
//       ConstLinearOperator<double> ApOut = tmpApOut;

//       Vector<double> rangevec = S.range().createMember();
//       Thyra::randomize(-1.0, 1.0, rangevec.ptr().get());

//       double err_pcd1 = norm2(S*svec - blockOpOut*svec);
//       cerr << "norm2(S*svec - blockOpOut*svec) = "
// 	   << err_pcd1 << endl;
//       bool ok_pcd1 = abs(err_pcd1) < 1.0e-10;

//       double err_pcd2 = norm2(Fp*fpvec - FpOut*fpvec);
//       cerr << "norm2(Fp*fpvec - FpOut*fpvec) = "
// 	   << err_pcd2 << endl;
//       bool ok_pcd2 = abs(err_pcd2) < 1.0e-10;

//       double err_pcd3 = norm2(Ap*apvec - ApOut*apvec);
//       cerr << "norm2(Ap*apvec - ApOut*apvec) = "
// 	   << err_pcd3 << endl;
//       bool ok_pcd3 = abs(err_pcd3) < 1.0e-10;

//       ConstLinearOperator<double> newblockOp = block2x2(Fout,Btout,Bout,Cout);
//       double err_pcd4 = norm2(S*svec - newblockOp*svec);
//       cerr << "norm2(S*svec - newblockOp*svec) = "
// 	   << err_pcd4 << endl;
//       bool ok_pcd4 = abs(err_pcd4) < 1.0e-10;

//       cout << "first pcd test status: " << ok_pcd1 << endl;
//       cout << "second pcd test status: " << ok_pcd2 << endl;
//       cout << "third pcd test status: " << ok_pcd3 << endl;
//       cout << "forth pcd test status: " << ok_pcd4 << endl;



      /* ------------- test LSC Operator Source --------------*/
      // Tests of LSCOperatorSource
      RCP<const LSCOperatorSource> lscOpSrcRcp 
	= rcp(new LSCOperatorSource(FMatrix, BtMatrix, BMatrix, CMatrix));

      // Get the blocks back from the operator source
      ConstLinearOperator<double> blockOpOut2 = lscOpSrcRcp->getSaddleOp();
      ConstLinearOperator<double> Fout2 = blockOpOut2.getBlock(0,0);
      ConstLinearOperator<double> Btout2 = blockOpOut2.getBlock(0,1);
      ConstLinearOperator<double> Bout2 = blockOpOut2.getBlock(1,0);
      ConstLinearOperator<double> Cout2 = blockOpOut2.getBlock(1,1);

      double err_lsc1 = norm2(S*svec - blockOpOut2*svec);
      cerr << "norm2(S*svec - blockOpOut2*svec) = "
	   << err_lsc1 << endl;
      bool ok_lsc1 = abs(err_lsc1) < 1.0e-10;

      ConstLinearOperator<double> newblockOp2 
	= block2x2(Fout2,Btout2,Bout2,Cout2);
      double err_lsc2 = norm2(S*svec - newblockOp2*svec);
      cerr << "norm2(S*svec - newblockOp*svec) = "
	   << err_lsc2 << endl;
      bool ok_lsc2 = abs(err_lsc2) < 1.0e-10;

      cout << "first lsc test status: " << ok_lsc1 << endl;
      cout << "second lsc test status: " << ok_lsc2 << endl;


      //  success = ok_pcd1 && ok_pcd2 && ok_pcd3 && ok_pcd4 
      //     && ok_lsc1 && ok_lsc2;
      
      success = ok_lsc1 && ok_lsc2;
      
    }

  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,out.get()?*out:std::cerr,success)

    if (success)
      {
        *out << "all tests PASSED!" << std::endl;
	return 0;
      }
    else
      {
	*out << "at least one test FAILED!" << std::endl;
        return 1;
      }
}




