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

int main( int argc, char argv[] )
{
	namespace rcp   = ReferenceCountingPack;
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
		size_type    mU = 0;
		size_type    mI = 0;
		value_type   xl = -NLP::infinite_bound();
		value_type   xu = +NLP::infinite_bound();
		value_type   hl = -NLP::infinite_bound();
		value_type   hu = +NLP::infinite_bound();
		
		// ToDo: readin the arguments from argv[]
		
		ExampleNLPBanded
			nlp(nD,nI,bw,mU,mI,xl,xu,hl,hu);

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
