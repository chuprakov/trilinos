// //////////////////////////////////////////////////////////////////
// ExampleNLPBandedMain.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <iostream>

#include "ExampleNLPBanded.h"
#include "ReducedSpaceSQPPack/Configurations/rSQPppSolver.h"

int main( int argc, char* argv[] )
{
	namespace rcp   = MemMngPack;
	namespace rsqp  = ReducedSpaceSQPPack;
	namespace nlpip = NLPInterfacePack;
	using rsqp::rSQPppSolver;
	using nlpip::NLP;
	using nlpip::ExampleNLPBanded;
	typedef nlpip::size_type  size_type;
	typedef nlpip::value_type value_type;

	try {

		// Get options!
		
		size_type    nD = 2;
		size_type    nI = 1;
		size_type    bw = 1;
		value_type   xo = 0.1;
		bool         nlp_selects_basis = true;
		value_type   xl = -NLP::infinite_bound();
		value_type   xu = +NLP::infinite_bound();
		size_type    mU = 0;
		size_type    mI = 0;
		value_type   hl = -NLP::infinite_bound();
		value_type   hu = +NLP::infinite_bound();
		value_type   diag_scal = 1.0;
		value_type   diag_vary = 1.0;
		bool         sym_basis = false;
		
		// Read from the arguments
		if(argc > 1)
			nD = ::atoi(argv[1]);
		if(argc > 2)
			nI = ::atoi(argv[2]);
		if(argc > 3)
			bw = ::atoi(argv[3]);
		if(argc > 4)
			xo = ::atof(argv[4]);
		if(argc > 5)
			nlp_selects_basis = (::atoi(argv[5]) != 0 );
		if(argc > 6)
			diag_scal  = ::atof(argv[6]);
		if(argc > 7)
			diag_vary  = ::atof(argv[7]);
		if(argc > 8)
			sym_basis = (::atoi(argv[8]) != 0 );
		
		// ToDo: readin more the arguments from argv[] when options are supported
		
		ExampleNLPBanded
			nlp(nD,nI,bw,mU,mI,xo,xl,xu,hl,hu,nlp_selects_basis,diag_scal,diag_vary,sym_basis);

		rSQPppSolver  solver;

		solver.set_nlp( rcp::rcp(&nlp,false) );

		const rSQPppSolver::ESolutionStatus
			solution_status = solver.solve_nlp();
		
		return solution_status;

	}
	catch(const std::exception& excpt) {
		std::cerr << "\nCaught a std::exception " << excpt.what() << std::endl;
	}
	catch(...) {
		std::cerr << "\nCaught an unknown exception\n";
	}

	return rSQPppSolver::SOLVE_RETURN_EXCEPTION;
}
