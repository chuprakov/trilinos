// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

#include "TSFCoreSerialVectorSpaceStd.hpp"
#include "TSFCoreSerialMultiVectorStd.hpp"
#include "TSFCoreLinearOpScalarProd.hpp"
#include "TSFCoreEuclideanScalarProd.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "TSFCoreDiagonalLinearOp.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreLinearOpTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

///
/** Main test driver function for scalar products
 */
template <class Scalar>
bool run_scalar_product_tests(
	const int                                                     n
	,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &tol
	,const bool                                                   dumpAll
	,std::ostream                                                 *out
	)
{

	using TSFCore::relErr;
	typedef TSFCore::LinOpNonPersisting<Scalar> LONP;
	typedef Teuchos::ScalarTraits<Scalar> ST;
	typedef typename ST::magnitudeType    ScalarMag;
	using Teuchos::RefCountPtr;
	using Teuchos::rcp;

	if(out) *out << "\n*** Entering run_scalar_product_tests<"<<ST::name()<<">(...) ...\n";

	bool success = true, result;

	RefCountPtr<TSFCore::ScalarProdVectorSpaceBase<Scalar> >
		domain = rcp(new TSFCore::SerialVectorSpaceStd<Scalar>(n/2)),
		range  = rcp(new TSFCore::SerialVectorSpaceStd<Scalar>(n));
	RefCountPtr<TSFCore::SerialMultiVectorStd<Scalar> >
		op_coeff = rcp(new TSFCore::SerialMultiVectorStd<Scalar>(range,domain)),
		op       = rcp(new TSFCore::SerialMultiVectorStd<Scalar>(range,domain));
	TSFCore::seed_randomize<Scalar>(0);
	TSFCore::randomize( Scalar(Scalar(-1)*ST::one()), Scalar(Scalar(+1)*ST::one()), &*op_coeff );
	if(out && dumpAll) *out << "\nop_coeff =\n" << *op_coeff;
	RefCountPtr<TSFCore::Vector<Scalar> >
		domainScalarProdDiag = domain->createMember(),
		rangeScalarProdDiag  = range->createMember();
	TSFCore::randomize( Scalar(Scalar(+1)*ST::one()), Scalar(Scalar(+2)*ST::one()), &*domainScalarProdDiag );
	TSFCore::randomize( Scalar(Scalar(+1)*ST::one()), Scalar(Scalar(+2)*ST::one()), &*rangeScalarProdDiag );

	const ScalarMag warning_tol = ScalarMag(1e-2)*tol, error_tol = tol;
	TSFCore::LinearOpTester<Scalar> linearOpTester(warning_tol,error_tol);
	
	if(out) *out << "\nTesting LinearOp with Euclidean domain and range scalar products ...\n";
	TSFCore::assign( &*op, *op_coeff );
	if(out && dumpAll) *out << "\nop =\n" << LONP(*op,TSFCore::NOTRANS);
	if(out && dumpAll) *out << "\nop' =\n" << LONP(*op,TSFCore::CONJTRANS);
	result = linearOpTester.check(*op,out);
	if(!result) success = false;
	
	if(out) *out << "\nTesting LinearOp with non-Euclidean domain and Euclidean range scalar products ...\n";
	range->setScalarProd(rcp(new TSFCore::EuclideanScalarProd<Scalar>()));
	domain->setScalarProd(
		rcp(
			new TSFCore::LinearOpScalarProd<Scalar>(
				TSFCore::LinOpPersisting<Scalar>(rcp(new TSFCore::DiagonalLinearOp<Scalar>(domainScalarProdDiag)))
				)
			)
		);
	op->initialize(range,domain);
	TSFCore::assign( &*op, *op_coeff );
	if(out && dumpAll) *out << "\nop =\n" << LONP(*op,TSFCore::NOTRANS);
	if(out && dumpAll) *out << "\nop' =\n" << LONP(*op,TSFCore::CONJTRANS);
	result = linearOpTester.check(*op,out);
	if(!result) success = false;
	
	if(out) *out << "\nTesting LinearOp with Euclidean domain and non-Euclidean range scalar products ...\n";
	range->setScalarProd(
		rcp(
			new TSFCore::LinearOpScalarProd<Scalar>(
				TSFCore::LinOpPersisting<Scalar>(rcp(new TSFCore::DiagonalLinearOp<Scalar>(rangeScalarProdDiag)))
				)
			)
		);
	domain->setScalarProd(rcp(new TSFCore::EuclideanScalarProd<Scalar>()));
	op->initialize(range,domain);
	TSFCore::assign( &*op, *op_coeff );
	if(out && dumpAll) *out << "\nop =\n" << LONP(*op,TSFCore::NOTRANS);
	if(out && dumpAll) *out << "\nop' =\n" << LONP(*op,TSFCore::CONJTRANS);
	result = linearOpTester.check(*op,out);
	if(!result) success = false;
	
	if(out) *out << "\nTesting LinearOp with non-Euclidean domain and non-Euclidean range scalar products ...\n";
	range->setScalarProd(
		rcp(
			new TSFCore::LinearOpScalarProd<Scalar>(
				TSFCore::LinOpPersisting<Scalar>(rcp(new TSFCore::DiagonalLinearOp<Scalar>(rangeScalarProdDiag)))
				)
			)
		);
	domain->setScalarProd(
		rcp(
			new TSFCore::LinearOpScalarProd<Scalar>(
				TSFCore::LinOpPersisting<Scalar>(rcp(new TSFCore::DiagonalLinearOp<Scalar>(domainScalarProdDiag)))
				)
			)
		);
	op->initialize(range,domain);
	TSFCore::assign( &*op, *op_coeff );
	if(out && dumpAll) *out << "\nop =\n" << LONP(*op,TSFCore::NOTRANS);
	if(out && dumpAll) *out << "\nop' =\n" << LONP(*op,TSFCore::CONJTRANS);
	result = linearOpTester.check(*op,out);
	if(!result) success = false;

  if(out) *out << "\n*** Leaving run_scalar_product_tests<"<<ST::name()<<">(...) ...\n";

	return success;

} // end run_scalar_product_tests() [Doxygen looks for this!]

int main( int argc, char* argv[] ) {

	using Teuchos::CommandLineProcessor;

	bool success = true;
	bool verbose = true;

	std::ostream &out = std::cout;

	try {

		//
		// Read options from commandline
		//

		int n         = 4;
		bool dumpAll  = false;

		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		clp.setOption( "n", &n, "Number of elements in each constituent vector." );
		clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

		//
		// Run the tests
		//

		if( !run_scalar_product_tests<float>(n,float(1e-6),dumpAll,verbose?&out:NULL) ) success = false;
		//if( !run_scalar_product_tests<double>(n,double(1e-14),dumpAll,verbose?&out:NULL) ) success = false;
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
		//if( !run_scalar_product_tests<std::complex<float> >(n,float(1e-6),dumpAll,verbose?&out:NULL) ) success = false;
		//if( !run_scalar_product_tests<std::complex<double> >(n,double(1e-14),dumpAll,verbose?&out:NULL) ) success = false;
#endif
#ifdef HAVE_TEUCHOS_GNU_MP
		//if( !run_scalar_product_tests<mpf_class>(n,mpf_class(1e-14),dumpAll,verbose?&out:NULL) ) success = false;
		// Above commented out code will not compile because its ScalarTraits specializatioin does not support eps()
#endif

	} // end try
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
		success = false;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception!\n";
		success = false;
	}

	if(verbose) {
		if(success)
			out << "\nAll of the tests seem to have run sucessfully!\n";
		else
			out << "\nOh no! at least one of the test failed!\n";	
	}
	
	return success ? 0 : 1;

} // end main() [Doxygen looks for this!]
