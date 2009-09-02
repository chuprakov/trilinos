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
using namespace Thyra;
using namespace Meros;

template <class Scalar>
LinearOperator<Scalar> makeRandomDenseOperator(int nc, const VectorSpace<Scalar>& rowSp);


LinearOperator<double> makeOp();


int main(int argc, char *argv[]) 
{
  bool success = false;
  
  GlobalMPISession mpiSession(&argc, &argv);
  typedef Teuchos::ScalarTraits<double> ST;
  
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

      // Make a test block operator
      LinearOperator<double> S = makeOp();

      // Grabbing subblocks as tests for Fp and Ap operators
      LinearOperator<double> Ap = S.getBlock(0,0);
      LinearOperator<double> Fp = S.getBlock(1,1);

      Vector<double> svec = S.domain().createMember();
      Thyra::randomize(-1.0, 1.0, svec.ptr().get());
      Vector<double> fpvec = Fp.domain().createMember();
      Thyra::randomize(-1.0, 1.0, fpvec.ptr().get());
      Vector<double> apvec = Ap.domain().createMember();
      Thyra::randomize(-1.0, 1.0, apvec.ptr().get());

      // Set up a PCD operator source with these three operators
      RCP<const PCDOperatorSource> pcdOpSrcRcp 
	= rcp(new PCDOperatorSource(S, Fp, Ap));

      // Get the blocks back from the operator source
      ConstLinearOperator<double> blockOpOut = pcdOpSrcRcp->getSaddleOp();
      ConstLinearOperator<double> Fout = blockOpOut.getBlock(0,0);
      ConstLinearOperator<double> Btout = blockOpOut.getBlock(0,1);
      ConstLinearOperator<double> Bout = blockOpOut.getBlock(1,0);
      ConstLinearOperator<double> Cout = blockOpOut.getBlock(1,1);

      ConstLinearOperator<double> FpOut = pcdOpSrcRcp->getFp();
      ConstLinearOperator<double> ApOut = pcdOpSrcRcp->getAp();

      Vector<double> rangevec = S.range().createMember();
      Thyra::randomize(-1.0, 1.0, rangevec.ptr().get());

      double err_pcd1 = norm2(S*svec - blockOpOut*svec);
      cerr << "norm2(S*svec - blockOpOut*svec) = "
	   << err_pcd1 << endl;
      bool ok_pcd1 = abs(err_pcd1) < 1.0e-10;

      double err_pcd2 = norm2(Fp*fpvec - FpOut*fpvec);
      cerr << "norm2(Fp*fpvec - FpOut*fpvec) = "
	   << err_pcd2 << endl;
      bool ok_pcd2 = abs(err_pcd2) < 1.0e-10;

      double err_pcd3 = norm2(Ap*apvec - ApOut*apvec);
      cerr << "norm2(Ap*apvec - ApOut*apvec) = "
	   << err_pcd3 << endl;
      bool ok_pcd3 = abs(err_pcd3) < 1.0e-10;

      ConstLinearOperator<double> newblockOp = block2x2(Fout,Btout,Bout,Cout);
      double err_pcd4 = norm2(S*svec - newblockOp*svec);
      cerr << "norm2(S*svec - newblockOp*svec) = "
	   << err_pcd4 << endl;
      bool ok_pcd4 = abs(err_pcd4) < 1.0e-10;

      cout << "first pcd test status: " << ok_pcd1 << endl;
      cout << "second pcd test status: " << ok_pcd2 << endl;
      cout << "third pcd test status: " << ok_pcd3 << endl;
      cout << "forth pcd test status: " << ok_pcd4 << endl;


      // Tests of LSCOperatorSource
      RCP<const LSCOperatorSource> lscOpSrcRcp 
	= rcp(new LSCOperatorSource(S));

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


      success = ok_pcd1 && ok_pcd2 && ok_pcd3 && ok_pcd4 
	&& ok_lsc1 && ok_lsc2;
      
      
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




template <class Scalar> inline
LinearOperator<Scalar> makeRandomDenseOperator(int nc, const VectorSpace<Scalar>& rowSp)
{
  typedef typename Teuchos::ScalarTraits<Scalar> ST;
  RCP<Thyra::MultiVectorBase<Scalar> > mv = rowSp.createMembers(nc);
  Thyra::randomize(-ST::one(), ST::one(), &*mv);
  RCP<Thyra::LinearOpBase<Scalar> > rtn = mv;
  return rtn;
}



LinearOperator<double> makeOp()
{
  Array<int> rangeSpaceSizes = tuple(2, 2);
  Array<int> domainSpaceSizes = tuple(2, 2);
  
  Array<VectorSpace<double> > rangeBlocks(rangeSpaceSizes.size());
  
  Array<Array<LinearOperator<double> > > blocks(rangeSpaceSizes.size());
  
  for (unsigned int br=0; br<rangeSpaceSizes.size(); br++)
    {
      int n = rangeSpaceSizes[br];
      rangeBlocks[br] 
        = new DefaultSpmdVectorSpace<double>(DefaultComm<Index>::getComm(),n,-1);
      
      blocks[br].resize(domainSpaceSizes.size());
      
      for (unsigned int bc=0; bc<domainSpaceSizes.size(); bc++)
        {
          blocks[br][bc] = makeRandomDenseOperator<double>(domainSpaceSizes[bc],
                                                           rangeBlocks[br]);
        }
    }
  
  
  LinearOperator<double> A = block2x2(blocks[0][0], blocks[0][1],
                                      blocks[1][0], blocks[1][1]);
  
  return A;
}

