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

#include "TSFCoreSerialVectorSpace.hpp"
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
bool run_tests( const int n, const Scalar max_rel_err, const bool dumpAll, std::ostream* out )
{

	using Teuchos::arrayArg;

	typedef Teuchos::ScalarTraits<Scalar> ST;

	if(out) *out << "\n*** Entering run_tests<"<<ST::name()<<">) ...\n";

	bool success = true;
	//bool result;
	//Scalar sresult1, sresult2;

	if(out) *out << "\nCreating a serial vector space svp with n="<<n<<" vector elements ...\n";
	const SerialVectorSpace<Scalar>  svp(n);
	if(out) *out << "\nsvp.dim() = " << svp.dim() << std::endl;

	if(out) *out << "\nCreating serial vectors v1, v2, v3 and z ...\n";
	Teuchos::RefCountPtr<Vector<Scalar> >
		v1 = svp.createMember(),
		v2 = svp.createMember(),
		v3 = svp.createMember(),
		z  = svp.createMember();

	if(out) *out << "\nassign(&*v1,-2.0);\n";
	assign(&*v1,Scalar(-2.0));
	if(out) *out << "\nassign(&*v2,-3.0);\n";
	assign(&*v2,Scalar(-3.0));
	if(out) *out << "\nassign(&*v3,-4.0);\n";
	assign(&*v3,Scalar(-4.0));
	
	if(out) *out << "\nabs(&*z,*v1);\n";
	abs(&*z,*v1);
	if(!testRelErr("sum(*z)",sum(*z),"2.0*svp.dim()",Scalar(2.0)*Scalar(svp.dim()),"max_rel_err",max_rel_err,out)) success=false;
	
	if(out) *out << "\nreciprocal(&*z,*v1);\n";
	reciprocal(&*z,*v1);
	if(!testRelErr("sum(*z)",sum(*z),"-0.5*svp.dim()",Scalar(-0.5)*Scalar(svp.dim()),"max_rel_err",max_rel_err,out)) success=false;

	if(out) *out << "\nlinear_combination(2,{0.5,0.25},{&*v1,&*v2},0.0,&*z);\n";
	linear_combination(2,arrayArg<Scalar>(0.5,0.25)(),arrayArg<const Vector<Scalar>*>(&*v1,&*v2)(),Scalar(0.0),&*z);
	if(!testRelErr("sum(*z)",sum(*z),"(-0.5*2.0-0.25*3.0)*svp.dim()",Scalar(-0.5*2.0-0.25*3.0)*Scalar(svp.dim()),"max_rel_err",max_rel_err,out)) success=false;

	if(out) *out << "\nassign(&*z,2.0);\n";
	assign(&*z,Scalar(2.0));
	if(out) *out << "\nlinear_combination(3,{0.5,0.25,0.125},{&*v1,&*v2,&*v2},0.5,&*z);\n";
	linear_combination(3,arrayArg<Scalar>(0.5,0.25,0.125)(),arrayArg<const Vector<Scalar>*>(&*v1,&*v2,&*v3)(),Scalar(0.5),&*z);
	if(!testRelErr(
			 "sum(*z)",sum(*z)
			 ,"(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*svp.dim()",Scalar(0.5*2.0-0.5*2.0-0.25*3.0-0.125*4.0)*Scalar(svp.dim())
			 ,"max_rel_err",max_rel_err,out
			 )
		) success=false;

	if(out) *out << "\nassign(&*z,2.0);\n";
	assign(&*z,Scalar(2.0));
	if(!testRelErr(
			 "norm_2(*z,*v2)",norm_2(*z,*v2)
			 ,"sqrt(2.0*3.0*3.0*svp.dim())",ST::squareroot(Scalar(2.0*3.0*3.0)*Scalar(svp.dim()))
			 ,"max_rel_err",max_rel_err,out
			 )
		) success=false;

	// ToDo: Add tests for *all* standard operators!

  if(out) *out
		<< "\n*** Leaving run_tests<"<<ST::name()<<">) ...\n";

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
		if( !TSFCore::run_tests<std::complex<float> >(local_dim,std::complex<float>(max_rel_err),dumpAll,verbose?&out:NULL) ) success = false;
		if( !TSFCore::run_tests<std::complex<double> >(local_dim,std::complex<double>(max_rel_err),dumpAll,verbose?&out:NULL) ) success = false;
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
