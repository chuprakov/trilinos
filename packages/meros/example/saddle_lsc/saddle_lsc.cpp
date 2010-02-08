// @HEADER
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
// @HEADER

// saddle_lsc.cpp

// Example program that reads in Epetra matrices from files, creates a
// Meros Least Squares Commutator (LSC) preconditioner and does a
// solve.

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"

#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_Vector.hpp"
#include "Thyra_VectorSpace.hpp"
#include "Thyra_LinearOperator.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"

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

#include "AztecOO.h"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"

#include "Meros_ConfigDefs.h"
#include "Meros_Version.h"
#include "Meros_LSCPreconditionerFactory.h"
#include "Meros_LSCOperatorSource.h"
#include "Meros_AztecSolveStrategy.hpp"
#include "Meros_InverseOperator.hpp"
#include "Meros_ZeroOperator.hpp"
#include "Meros_IdentityOperator.hpp"
#include "Meros_LinearSolver.hpp"

using namespace Teuchos;
using namespace EpetraExt;
using namespace Thyra;
using namespace Meros;


int main(int argc, char *argv[]) 
{
  GlobalMPISession mpiSession(&argc, &argv);

  // DEBUG 0 = no extra tests or printing
  // DEBUG 1 = print out basic diagnostics as we go.
  //           prints outer saddle system iterations, but not inner solves
  // DEBUG > 1 = test usage of some of the operators before proceeding
  //           prints inner and outer iterations
  int DEBUG = 1;

  // Get stream that can print to just root or all streams!
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  //  Epetra_Comm* Comm;
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  const int myRank = Comm.MyPID();
  const int numProcs = Comm.NumProc();
#else
  Epetra_SerialComm Comm;
  const int myRank = 0;
  const int numProcs = 1;
#endif

  if(DEBUG > 0)
    {
      if (myRank == 0)
	{
	  cout << "Proc " << myRank << ": " 
	       << "Number of processors = " 
	       << numProcs << endl;
	  cout << "Proc " << myRank << ": " 
	       << Meros::Meros_Version() << endl;
	}
    }
  
  try
    {
      /* ------ Read in epetra matrices and rhs ------------------- */

      // Using Q1 matrix example; see Meros/examples/data/q1
      // (2,2 block (C block) is zero in this example)
  
      // Make necessary Epetra maps.
      // Need a velocity space map and a pressure space map.
      const Epetra_Map* velocityMap = new Epetra_Map(578, 0, Comm); 
      const Epetra_Map* pressureMap = new Epetra_Map(192, 0, Comm); 
 

      // Read matrix and vector blocks into Epetra_Crs matrices
      Epetra_CrsMatrix* FMatrix(0);
      char * filename = "../../../../../packages/meros/example/data/q1/Aq1.mm";
      MatrixMarketFileToCrsMatrix(filename,
                                  *velocityMap,
                                  *velocityMap, *velocityMap,
                                  FMatrix);

      // cerr << "If you get Epetra ERROR -1 here, the file was not found"
      //  << "Need to make symbolic link to data directory" 
      //  << endl;

      Epetra_CrsMatrix* BtMatrix(0);
      filename = "../../../../../packages/meros/example/data/q1/Btq1.mm";
      MatrixMarketFileToCrsMatrix(filename,
                                  *velocityMap, 
                                  *velocityMap, *pressureMap,
                                  BtMatrix);

      Epetra_CrsMatrix* BMatrix(0);
      filename = "../../../../../packages/meros/example/data/q1/Bq1.mm";
      MatrixMarketFileToCrsMatrix(filename,
                                  *pressureMap,
                                  *pressureMap, *velocityMap,
                                  BMatrix);


      Epetra_Vector* rhsq1_vel(0);
      filename = "../../../../../packages/meros/example/data/q1/rhsq1_vel.mm";
      MatrixMarketFileToVector(filename,
                               *velocityMap, 
                               rhsq1_vel);

      Epetra_Vector* rhsq1_press(0);
      filename = "../../../../../packages/meros/example/data/q1/rhsq1_press.mm";
      MatrixMarketFileToVector(filename,
                               *pressureMap, 
                               rhsq1_press);

      FMatrix->FillComplete();
      BtMatrix->FillComplete(*pressureMap, *velocityMap);
      BMatrix->FillComplete(*velocityMap, *pressureMap);


      // Wrap Epetra operators into Thyra operators. To do this, we
      // first wrap as Thyra core operators, then convert to the
      // handle layer LinearOperators.
      
      RCP<LinearOpBase<double> >
        tmpF = nonconstEpetraLinearOp(rcp(FMatrix,false));
      const LinearOperator<double> F = tmpF;

      RCP<LinearOpBase<double> >
        tmpBt = nonconstEpetraLinearOp(rcp(BtMatrix,false));
      const LinearOperator<double> Bt = tmpBt;

      RCP<LinearOpBase<double> >
        tmpB = nonconstEpetraLinearOp(rcp(BMatrix,false));
      const LinearOperator<double> B = tmpB;


      // Wrap Epetra vectors into Thyra vectors similarly by first
      // wrapping to the Thyra core vector layer, then converting to
      // the handle layer

      RCP<const Thyra::VectorSpaceBase<double> > epetra_vs_press
        = Thyra::create_VectorSpace(rcp(pressureMap,false));
      RCP<const Thyra::VectorSpaceBase<double> > epetra_vs_vel
        = Thyra::create_VectorSpace(rcp(velocityMap,false));


      RCP<VectorBase<double> > rhs1
        = create_Vector(rcp(rhsq1_press, false), epetra_vs_press);
      RCP<VectorBase<double> > rhs2
        = create_Vector(rcp(rhsq1_vel, false), epetra_vs_vel);

      // Convert the vectors to handled vectors
      RCP<VectorBase<double> > tmp1 = rhs1;
      const Vector<double> rhs_press = tmp1;

      RCP<VectorBase<double> > tmp2 = rhs2;
      const Vector<double> rhs_vel = tmp2;


      if(DEBUG > 1){
        // Test the wrapped operators
        cerr << "F matrix: description  " << F.description() << endl;
        cerr << "F domain: " << F.domain().dim() 
             << ",  F range: " << F.range().dim() 
             << endl;
        cerr << "B matrix: description  " << B.description() << endl;
        cerr << "B domain: " << B.domain().dim() 
             << ",  B range: " << B.range().dim() 
             << endl;
        cerr << "Bt matrix: description  " << Bt.description() << endl;
        cerr << "Bt domain: " << Bt.domain().dim() 
             << ",  Bt range: " << Bt.range().dim() 
             << endl;

        Vector<double> testvec1 = Bt * rhs_press;
        cerr << "Bt * rhs_press = " << norm2(testvec1) <<endl;

	cerr << "domain size of B " << B.domain().dim() << endl;
	cerr << "range size of B " << B.range().dim() << endl;
	cerr << "size of testvec1 " << dim(testvec1) << endl;

        Vector<double> testvec2 = B * testvec1;
        cerr << "B * Bt * rhs_press = " << norm2(testvec2) <<endl;

        Vector<double> testvec3 = F * rhs_vel;
        cerr << "F * rhs_vel = " << norm2(testvec3) <<endl;

        Vector<double> testvec4 = B * F * Bt * rhs_press;
        cerr << "BFBt * rhs_press = " << norm2(testvec4) <<endl;

      }



      // Get Thyra velocity and pressure spaces from the operators.
      VectorSpace<double> velocitySpace = F.domain();
      VectorSpace<double> pressureSpace = Bt.domain();

      if(DEBUG > 0){
        cerr << "P" << myRank 
             << ": vel space dim = " << velocitySpace.dim() << endl;
        cerr << "P" << myRank 
             << ": press space dim = " << pressureSpace.dim() << endl;
      }

      // Make a zero operator on the small (pressure) space since the
      // 2,2 block (C) is zero in this example.
      // RCP<Thyra::LinearOpBase<double> > tmpZ = 
      // rcp(new DefaultZeroLinearOp<double>(tmpBt->domain(), 
      //                                   tmpBt->domain()));
      // const LinearOperator<double> Z = tmpZ;
      const LinearOperator<double> Z;

      // Build the block saddle operator with F, Bt, B, and Z.  This
      // operator will be blockOp = [F Bt; B zero] (in Matlab
      // notation).
      LinearOperator<double> blockOp = block2x2(F, Bt, B, Z);

      if(DEBUG > 1) 
        {
          // Test getting subblock components out of the block operator.
          LinearOperator<double> testF1 = blockOp.getBlock(0,0);
          LinearOperator<double> testBt1 = blockOp.getBlock(0,1);
          LinearOperator<double> testB1 = blockOp.getBlock(1,0);
          cerr << "Checking blocks after extracting them out the block op"
               << endl;
          cerr << "testF domain: " << testF1.domain().dim() << endl;
          cerr << "testF range: " << testF1.range().dim() << endl;
          cerr << "testBt domain: " << testBt1.domain().dim() << endl;
          cerr << "testBt range: " << testBt1.range().dim() << endl;
          cerr << "testB domain: " << testB1.domain().dim() << endl;
          cerr << "testB range: " << testB1.range().dim() << endl;

          Vector<double> testvec1 = testBt1 * rhs_press;
          cerr << "Bt * rhs_press = " << norm2(testvec1) <<endl;
	  
          Vector<double> testvec2 = testB1 * testvec1;
          cerr << "B * Bt * rhs_press = " << norm2(testvec2) <<endl;
	  
          Vector<double> testvec3 = testF1 * rhs_vel;
          cerr << "F * rhs_vel = " << norm2(testvec3) <<endl;
	  
          Vector<double> testvec4 = testB1 * testF1 * testBt1 * rhs_press;
          cerr << "BFBt * rhs_press = " << norm2(testvec4) <<endl;
        }

      // Get the domain and range product spaces from the block operator. 
      VectorSpace<double> domain = blockOp.domain();
      VectorSpace<double> range = blockOp.range();

      // Build a product vector for the rhs and set the velocity and
      // pressure components of the rhs with the vectors we previously
      // read in from files and wrapped in Thyra.
      Vector<double> rhs = range.createMember();
      rhs.setBlock(0,rhs_vel);
      rhs.setBlock(1,rhs_press);

      // Build a solution vector and initialize it to zero.
      Vector<double> solnblockvec = domain.createMember();
      zeroOut(solnblockvec);


      // We now have Thyra block versions of the saddle point system,
      // rhs, and solution vector. Next we build the Meros LSC
      // preconditioner.


      /* -------- Build the Meros preconditioner factory ---------*/

      // Build a Least Squares Commutator (LSC) block preconditioner
      // with Meros
      // 
      // | inv(F) 0 | | I  -Bt | | I        |
      // | 0      I | |     I  | |   -inv(X)|
      // 
      // where inv(X) = inv(B*Bt) * B * F * Bt * inv(B*Bt)
      // (velocity mass matrix is the identity in this example)
      //
      // We'll do this in 4 steps:
      // 1) Build an AztecOO ParameterList for inv(F) solve
      // 2) Build an AztecOO ParameterList for inv(B*Bt) solve
      //    The Schur complement approximation inverse requires solves
      //    on the composed operator B*Bt.
      // 3) Make an LSCOperatorSource with blockOp  (and Qu if needed)
      // 4) Build the LSC block preconditioner factory 


      // 1) Build an AztecOO ParameterList for inv(F) solve
      //    This one corresponds to (unpreconditioned) GMRES

      RCP<ParameterList>
        aztecFParams = rcp(new ParameterList("aztecOOFSolverFactory"));

      RCP<LinearOpWithSolveFactoryBase<double> > aztecFLowsFactory;

      if(DEBUG> 1)
        {
          // Print out valid parameters and the existing default params.
          aztecFLowsFactory = rcp(new AztecOOLinearOpWithSolveFactory());
          cerr << "\naztecFLowsFactory.getValidParameters():\n" << endl;
          aztecFLowsFactory->getValidParameters()->print(cerr, 0, true, false);
          cerr << "\nPrinting initial parameters. " << endl;
          aztecFLowsFactory->setParameterList(aztecFParams);
          aztecFLowsFactory->getParameterList()->print(cerr, 0, true, false);
        }


      // forward solve settings
      aztecFParams->sublist("Forward Solve").set("Max Iterations", 100);
      aztecFParams->sublist("Forward Solve").set("Tolerance", 10e-8);
      // aztecOO solver settings
      aztecFParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Aztec Solver", "GMRES");
      aztecFParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Aztec Preconditioner", "none");
      aztecFParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Size of Krylov Subspace", 100);


      if(DEBUG > 1)
        {
          // turn on AztecOO output
          aztecFParams->sublist("Forward Solve")
            .sublist("AztecOO Settings").set("Output Frequency", 10);
        }

      if(DEBUG > 1)
        {
          // Print out the parameters we just set
          aztecFLowsFactory->setParameterList(aztecFParams);
          aztecFLowsFactory->getParameterList()->print(cerr, 0, true, false);
        }


      // 2) Build an AztecOO ParameterList for inv(Ap) solve
      //    This one corresponds to unpreconditioned CG.

      RCP<ParameterList>
        aztecBBtParams = rcp(new ParameterList("aztecOOBBtSolverFactory"));

      // forward solve settings
      aztecBBtParams->sublist("Forward Solve").set("Max Iterations", 100);
      aztecBBtParams->sublist("Forward Solve").set("Tolerance", 10e-8);
      // aztecOO solver settings
      aztecBBtParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Aztec Solver", "CG");
      aztecBBtParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Aztec Preconditioner", "none");


      if(DEBUG > 1)
        {
          // turn on AztecOO output
          aztecBBtParams->sublist("Forward Solve")
            .sublist("AztecOO Settings").set("Output Frequency", 10);
        }


      if(DEBUG > 1)
        {
          // Print out the parameters we just set
          RCP<LinearOpWithSolveFactoryBase<double> > 
            aztecBBtLowsFactory = rcp(new AztecOOLinearOpWithSolveFactory());
          aztecBBtLowsFactory->setParameterList(aztecBBtParams);
          aztecBBtLowsFactory->getParameterList()->print(cerr, 0, true, false);
        }


      // 3) Make an LSCOperatorSource that contains  blockOp
      //    The velocity mass matrix Qu is the identity in this example.
      RCP<const LinearOpSourceBase<double> > myLSCopSrcRcp 
        = rcp(new LSCOperatorSource(blockOp));


      if(DEBUG > 1) 
        {
          // Test getting subblock components out of the operator source
          RCP<const LSCOperatorSource> lscOpSrcPtr 
            = rcp_dynamic_cast<const LSCOperatorSource>(myLSCopSrcRcp);  
	  
          // Retrieve operators from the LSC operator source
          RCP<const LinearOpBase<double> > tmpBlockOp 
            = lscOpSrcPtr->getOp();
          ConstLinearOperator<double> blockOp = tmpBlockOp;
          ConstLinearOperator<double> testF2 = blockOp.getBlock(0,0);
          ConstLinearOperator<double> testBt2 = blockOp.getBlock(0,1);
          ConstLinearOperator<double> testB2 = blockOp.getBlock(1,0);	 
	  
          cerr << "Checking blocks after extracting them out the block op"
               << endl;
          cerr << "testF domain: " << testF2.domain().dim() << endl;
          cerr << "testF range: " << testF2.range().dim() << endl;
          cerr << "testBt domain: " << testBt2.domain().dim() << endl;
          cerr << "testBt range: " << testBt2.range().dim() << endl;
          cerr << "testB domain: " << testB2.domain().dim() << endl;
          cerr << "testB range: " << testB2.range().dim() << endl;

          Vector<double> testvec11 = testBt2 * rhs_press;
          cerr << "Bt * rhs_press = " << norm2(testvec11) <<endl;
	  
          Vector<double> testvec22 = testB2 * testvec11;
          cerr << "B * Bt * rhs_press = " << norm2(testvec22) <<endl;
	  
          Vector<double> testvec33 = testF2 * rhs_vel;
          cerr << "F * rhs_vel = " << norm2(testvec33) <<endl;
	  
          Vector<double> testvec44 = testB2 * testF2 * testBt2 * rhs_press;
          cerr << "BFBt * rhs_press = " << norm2(testvec44) <<endl;
        }


      // 4) Build the LSC block preconditioner factory.
      RCP<PreconditionerFactoryBase<double> > merosPrecFac
        = rcp(
	      new LSCPreconditionerFactory(
		rcp(new Thyra::AztecOOLinearOpWithSolveFactory(aztecFParams)),
		rcp(new Thyra::AztecOOLinearOpWithSolveFactory(aztecBBtParams))
		)
	      );
      
      RCP<PreconditionerBase<double> > Prcp 
        = merosPrecFac->createPrec();
      
      merosPrecFac->initializePrec(myLSCopSrcRcp, &*Prcp);


      // Checking that isCompatible and uninitializePrec at least
      // compile and throw the intended exceptions for now.
      // merosPrecFac->uninitializePrec(&*Prcp, &myLSCopSrcRcp);
      // bool passed;
      // passed = merosPrecFac->isCompatible(*(myLSCopSrcRcp.get()));


      /* --- Now build a solver factory for outer saddle point problem --- */

      // Set up parameter list and AztecOO solver
      RCP<ParameterList> aztecSaddleParams 
        = rcp(new ParameterList("aztecOOSaddleSolverFactory"));
      
      RCP<LinearOpWithSolveFactoryBase<double> >
        aztecSaddleLowsFactory = rcp(new AztecOOLinearOpWithSolveFactory());
      
      double saddleTol = 1.0e-6;

      // forward solve settings
      aztecSaddleParams->sublist("Forward Solve").set("Max Iterations", 500);
      aztecSaddleParams->sublist("Forward Solve").set("Tolerance", saddleTol);
      // aztecOO solver settings
      aztecSaddleParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Aztec Solver", "GMRES");
      aztecSaddleParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Aztec Preconditioner", "none");
      aztecSaddleParams->sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Size of Krylov Subspace", 500);

      if(DEBUG > 0)
        {
          // turn on AztecOO output
          aztecSaddleParams->sublist("Forward Solve")
            .sublist("AztecOO Settings").set("Output Frequency", 1);
        }

      aztecSaddleLowsFactory->setParameterList(aztecSaddleParams);

      if(DEBUG > 0)
        {
          // Print out the parameters we've set.
          aztecSaddleLowsFactory->getParameterList()->print(cerr, 0, 
                                                            true, false);
        }



      // Set up the preconditioned inverse object and do the solve!
      RCP<LinearOpWithSolveBase<double> > rcpAztecSaddle 
        = aztecSaddleLowsFactory->createOp();

      // LinearOperator<double> epetraBlockOp = makeEpetraOperator(blockOp);

      // initializePreconditionedOp<double>(*aztecSaddleLowsFactory, 
      //			 epetraBlockOp.ptr(), 
      //				 Prcp,
      //				 &*rcpAztecSaddle );

      //       initializePreconditionedOp<double>(*aztecSaddleLowsFactory, 
      // 					 blockOp.ptr(), 
      // 					 Prcp,
      // 					 &*rcpAztecSaddle );
      
      //       RCP<LinearOpBase<double> > tmpSaddleInv 
      // 	= rcp(new DefaultInverseLinearOp<double>(rcpAztecSaddle));
      
      //       LinearOperator<double> saddleInv = tmpSaddleInv;
      //       saddleInv.description();

      RCP<const LinearOpBase<double> > tmpPinv 
        = Prcp->getRightPrecOp();
      ConstLinearOperator<double> Pinv = tmpPinv;

      LinearSolveStrategy<double> azSaddle 
        = new AztecSolveStrategy(*(aztecSaddleParams.get()));

      ConstLinearOperator<double> saddleInv 
	= new InverseOperator<double>(blockOp * Pinv, azSaddle);

      // Do the solve!
      solnblockvec = saddleInv * rhs;

      // Check our results.
      Vector<double> residvec = blockOp * Pinv * solnblockvec - rhs;

      cerr << "norm of resid " << norm2(residvec) << endl;
      
      double normResvec = norm2(residvec);

      if(normResvec < 10.0*saddleTol)
	{
	  cerr << "Example PASSED!" << endl;
	  return 0;
	}
      else
	{
	  cerr << "Example FAILED!" << endl;
	  return 1;
	}



    } // end of try block


  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }

  MPISession::finalize();

} // end of main()

