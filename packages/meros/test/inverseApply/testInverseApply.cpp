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

      LinearOperator<double> A = makeOp();
      LinearOperator<double> A11 = A.getBlock(0,0);

      LinearOperator<double> I = new IdentityOperator<double>(A11.range());

      Vector<double> rangevec = A11.range().createMember();
      Thyra::randomize(-1.0, 1.0, rangevec.ptr().get());
      Vector<double> domainvec = A11.domain().createMember();
      Thyra::randomize(-1.0, 1.0, domainvec.ptr().get());

      double err1 = norm2(domainvec - I*domainvec);
      cerr << "norm(domainvec - I * domainvec) = " 
	   << err1 << endl;

      bool ok1 = abs(err1) < 1.0e-10;


      ParameterList azParams = AztecSolveStrategy::getParameters();
      azParams.sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Aztec Preconditioner", "none");
      azParams.sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Output Frequency", 1);

      LinearSolveStrategy<double> az = new AztecSolveStrategy(azParams);

      // LinearOperator<double> Iinv 
      // = new InverseOperator<double>(makeEpetraOperator(I), az);
      LinearOperator<double> Iinv 
	= new InverseOperator<double>(I, az);

      double err2 = norm2(domainvec - Iinv*domainvec);
      cerr << "norm(domainvec - Iinv*domainvec) = " 
	   << err2 << endl;

      bool ok2 = abs(err2) < 1.0e-10;


      LinearOperator<double> Zero 
	= new ZeroOperator<double>(A11.range(), A11.domain());
      
      LinearOperator<double> blockI = block2x2(I, Zero, Zero, I);
      Vector<double> blockvec = blockI.range().createMember();
      Thyra::randomize(-1.0, 1.0, blockvec.ptr().get());

      double err3 = norm2(blockvec - blockI*blockvec);
      cerr << "norm2(blockvec - blockI*blockvec) = "
	   << err3 << endl;

      bool ok3 = abs(err3) < 1.0e-10;


      // LinearOperator<double> blockIinv 
      // = new InverseOperator<double>(makeEpetraOperator(blockI), az);
      LinearOperator<double> blockIinv 
	= new InverseOperator<double>(blockI, az);

      Vector<double> tmpvec = blockIinv * blockvec;
      double err4 = norm2(blockvec - tmpvec);

      cerr << "norm2(blockvec - blockIinv*blockvec) = "
	   << err4 << endl;

      bool ok4 = abs(err4) < 1.0e-10;


      cout << "first test status: " << ok1 << endl;
      cout << "second test status: " << ok2 << endl;
      cout << "third test status: " << ok3 << endl;
      cout << "forth test status: " << ok4 << endl;

      success = ok1 && ok2 && ok3 && ok4;
      
      
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
      rangeBlocks[br] = defaultSpmdVectorSpace<double>(
        DefaultComm<Index>::getComm(),n,-1);
      
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

