// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
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

#include "RTOpUnaryFuncPtr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

//
// Here are some example unary functions for use with RTOpUnaryFuncPtr class
//

namespace MyDumbPack {

void sin( const double x[], int x_dim, double out[] )
{
	for( int i = 0; i < x_dim; ++i ) out[i] = ::sin(x[i]);
}

void sqrt( const double x[], int x_dim, double out[] )
{
	for( int i = 0; i < x_dim; ++i ) out[i] = ::sqrt(x[i]);
}

} // namespace MyDumbPack

int main( int argc, char* argv[] ) {

	using Teuchos::CommandLineProcessor;

	bool verbose = true;

	try {

		//
		// Read options from the commandline
		//

		CommandLineProcessor  clp(false); // Don't throw exceptions
		clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

		//
		// Create a couple of vectors to use
		//
		
		const int n = 10;

		double x_store[n], y_store[n];

		std::fill_n( x_store, n, 1.0 );

		RTOpPack::SubVectorT<double>         sub_vecs[1];
		RTOpPack::MutableSubVectorT<double>  targ_sub_vecs[1];

		sub_vecs[0].initialize( 0, n, x_store, 1 );
		targ_sub_vecs[0].initialize( 0, n, y_store, 1 );

		//
		// Create a RTOpUnaryFuncPtr and use it
		//

		RTOpPack::RTOpUnaryFuncPtr<double>  unary_rtop;

		unary_rtop.initialize( MyDumbPack::sin, "sin" );

		unary_rtop.apply_op( 1, sub_vecs, 1, targ_sub_vecs, NULL );

		// ToDo: Test the output!

		unary_rtop.initialize( MyDumbPack::sqrt, "sqrt" );

		unary_rtop.apply_op( 1, sub_vecs, 1, targ_sub_vecs, NULL );

		// ToDo: Test the output!

		if(verbose)
			std::cout << "\nAll of the tests seem to have run sucessfully!\n";

	} // end try
	catch( const std::exception &excpt ) {
		if(verbose)
			std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
		return -1;
	}
	catch( ... ) {
		if(verbose)
			std::cerr << "*** Caught an unknown exception!\n";
		return -1;
	}

	return 0;

}
