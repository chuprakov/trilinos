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

#include "NLPInterfacePack_ExampleNLPBanded.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

int main( int argc, char* argv[] )
{
	namespace mp  = MoochoPack;
	namespace nlpip = NLPInterfacePack;
	using mp::MoochoSolver;
	using nlpip::NLP;
	using nlpip::ExampleNLPBanded;
	typedef nlpip::size_type  size_type;
	typedef nlpip::value_type value_type;
	using Teuchos::CommandLineProcessor;

	try {
	
		//
		// Get options from the command line
		//
		
		int      nD = 1;
		int      nI = 1;
		int      bw = 1;
		int      mI = 0;
		double   xo = 1;
		bool     nlp_selects_basis = true;
		double   xDl = -NLP::infinite_bound();
		double   xDu = +NLP::infinite_bound(); 
		double   xIl = -NLP::infinite_bound(); 
		double   xIu = +NLP::infinite_bound();
		/*double   xDl = -1000.;
		double   xDu = +1000.;
		double   xIl = -1000.;
		double   xIu = +1000.;*/

		int      mU = 0;
		double   hl = -NLP::infinite_bound();
		double   hu = +NLP::infinite_bound();
		double   diag_scal = 1.0;
		double   diag_vary = 1.0;
		bool     sym_basis = false;
		double   f_offset  = 0.0;
		double   co        = 0.0;
		bool     ignore_constraints = false;
		
		CommandLineProcessor  command_line_processor(false); // don't throw exceptions

		command_line_processor.setOption( "nD",  &nD, "Number of dependent variables" );
		command_line_processor.setOption( "nI",  &nI, "Number of independent variables" );
		command_line_processor.setOption( "bw",  &bw, "Band width of the basis matrix" );
		command_line_processor.setOption( "mI",  &mI, "Number of general inequality constriants" );
		command_line_processor.setOption( "xo",  &xo, "Initial guess for x" );
		command_line_processor.setOption( "xDl", &xDl, "Lower bounds on xD" );
		command_line_processor.setOption( "xDu", &xDu, "Upper bounds on xD" );
		command_line_processor.setOption( "xIl", &xIl, "Lower bounds on xI" );
		command_line_processor.setOption( "xIu", &xIu, "Upper bounds on xI" );
//		command_line_processor.setOption( "mU",  &mU,  "Number of dependent equality constriants" );
		command_line_processor.setOption( "hl", &hl, "Lower bounds on general inequalities" );
		command_line_processor.setOption( "hu", &hu, "Upper bounds on general inequalities" );
		command_line_processor.setOption( "diag-scal", &diag_scal, "Scaling of the basis diagonal" );
		command_line_processor.setOption( "diag-vary", &diag_vary, "Variation of the basis diagonal scaling" );
		command_line_processor.setOption(
			"nlp-selects-basis", "no-nlp-selects-basis", &nlp_selects_basis
			,"Determine if the NLP will select basis" );
		command_line_processor.setOption(
			"sym-basis", "unsym-basis", &sym_basis
			,"Determine if the basis is symmetric" );
		command_line_processor.setOption( "f_offset", &f_offset, "Constant offset for objective function" );
		command_line_processor.setOption( "co", &co, "Constant term in general equalities" );
		command_line_processor.setOption(
			"ignore-constraints", "no-ignore-constraints", &ignore_constraints
			,"Determine if constraints are ignored or not" );
	
		CommandLineProcessor::EParseCommandLineReturn
			parse_return = command_line_processor.parse(argc,argv,&std::cerr);

		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
			return parse_return;
    
    std::cout << "\nPrinting commandline options used (options used shown as (default: \?\?\?) ...\n\n";
    command_line_processor.printHelpMessage(argv[0],std::cout);
		
		//
		// Create the NLP
		//

		ExampleNLPBanded
			nlp(nD,nI,bw,mU,mI,xo,xDl,xDu,xIl,xIu,hl,hu
				,nlp_selects_basis,diag_scal,diag_vary
				,sym_basis,f_offset,co,ignore_constraints
				);

		MoochoSolver  solver;

		solver.set_nlp( Teuchos::rcp(&nlp,false) );

		const MoochoSolver::ESolutionStatus
			solution_status = solver.solve_nlp();
		
		return solution_status;

	}
	catch(const std::exception& excpt) {
		std::cerr << "\nCaught a std::exception " << excpt.what() << std::endl;
	}
	catch(...) {
		std::cerr << "\nCaught an unknown exception\n";
	}

	return MoochoSolver::SOLVE_RETURN_EXCEPTION;
}
