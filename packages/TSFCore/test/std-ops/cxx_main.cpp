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

// ///////////////////////////////
// cxx_main.cpp

#include "TSFCoreSerialVectorSpaceStd.hpp"
#include "TSFCoreProductVectorSpace.hpp"
#include "TSFCoreVectorStdOpsTester.hpp"
#include "TSFCoreTestingTools.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_arrayArg.hpp"
#ifdef RTOp_USE_MPI
#  include "TSFCoreSimpleMPIVectorSpace.hpp"
#endif

namespace TSFCore {

///
/** Main test driver that runs tests on all standard operators
 */
template <class Scalar>
bool run_tests(
	const int                                                       n
	,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType    max_rel_err
	,const bool                                                     dumpAll
	,std::ostream                                                   *out
	)
{

	VectorStdOpsTester<Scalar> vectorStdOpsTester(max_rel_err);

	typedef Teuchos::ScalarTraits<Scalar> ST;

	if(out) *out << "\n*** Entering run_tests<"<<ST::name()<<">(...) ...\n";

	bool success = true;

	if(out) *out << "\nCreating a serial vector space svp with n="<<n<<" vector elements ...\n";
	const SerialVectorSpaceStd<Scalar>  svp(n);

	if(!vectorStdOpsTester.checkStdOps(svp,out,dumpAll)) success = false;

	const int numBlocks = 2;

	if(out) *out << "\nCreating a product space ps with numBlocks="<<numBlocks<<" and n="<<n<<"vector elements per block ...\n";

	std::vector<Teuchos::RefCountPtr<const TSFCore::VectorSpace<Scalar> > >
		vecSpaces(numBlocks);
	Teuchos::RefCountPtr<const TSFCore::VectorSpace<Scalar> >
		spaceBlock = Teuchos::rcp(new TSFCore::SerialVectorSpaceStd<Scalar>(n));
	for( int i = 0; i < numBlocks; ++i )
		vecSpaces[i] = spaceBlock;

	TSFCore::ProductVectorSpace<Scalar> ps(numBlocks,&vecSpaces[0]);

	if(!vectorStdOpsTester.checkStdOps(ps,out,dumpAll)) success = false;

	return success;

}

} // namespace TSFCore

int main( int argc, char* argv[] ) {

	using Teuchos::CommandLineProcessor;

	bool success = true;
	bool verbose = true;
	bool dumpAll = false;

	std::ostream &out = std::cout;

	try {

		//
		// Read options from the commandline
		//

		int     local_dim         = 4;
		double  max_rel_err       = 1e-13;
		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
		clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );
		clp.setOption( "local-dim", &local_dim, "Number of vector elements per process." );
		clp.setOption( "max-rel-err", &max_rel_err, "Maximum relative error for tests." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

		//
		// Run the tests
		//

		if( !TSFCore::run_tests<float>(local_dim,float(max_rel_err),dumpAll,verbose?&out:NULL) ) success = false;
		if( !TSFCore::run_tests<double>(local_dim,double(max_rel_err),dumpAll,verbose?&out:NULL) ) success = false;
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
		if( !TSFCore::run_tests<std::complex<float> >(local_dim,float(max_rel_err),dumpAll,verbose?&out:NULL) ) success = false;
		if( !TSFCore::run_tests<std::complex<double> >(local_dim,double(max_rel_err),dumpAll,verbose?&out:NULL) ) success = false;
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

}
