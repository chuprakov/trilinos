// /////////////////////////////////////////////////////////////////
// ExampleNLPFirstOrderRun.cpp
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

#ifdef _INTEL_CXX
// disable Intel C++ 5.0 warnings about debugger limitations
#pragma warning(disable : 985)
#endif

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "NLPInterfacePack_ExampleNLPFirstOrderRun.hpp"
#include "NLPInterfacePack_ExampleNLPFirstOrder.hpp"
#include "NLPInterfacePack_ExampleBasisSystem.hpp"
#include "MoochoPack_NLPAlgoConfigMamaJama.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"

MoochoPack::MoochoSolver::ESolutionStatus
NLPInterfacePack::ExampleNLPFirstOrderRun(
	const VectorSpace&   vec_space
	,value_type          xo
	,bool                has_bounds
	,bool                dep_bounded
	,std::ostream*       console_out
	,std::ostream*       error_out
	,bool                throw_solve_exception
	,std::ostream*       algo_out
	,std::ostream*       summary_out
	,std::ostream*       journal_out
	)
{
	using std::endl;
	using std::setw;
	namespace rcp = MemMngPack;
	using Teuchos::RefCountPtr;
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::OptionsFromStream;
	namespace rsqp = MoochoPack;
	using rsqp::MoochoSolver;
	using rsqp::NLPAlgoConfigMamaJama;

	MoochoSolver::ESolutionStatus
		solve_return = MoochoSolver::SOLVE_RETURN_EXCEPTION;

	int prec = 8;

	if(console_out)
		*console_out
			<< std::setprecision(prec)
			<< std::scientific
			<< "*************************************************\n"
			<< "*** Running Tests on ExampleNLPFirstOrder ***\n"
			<< "*************************************************\n"
			<< "\nUsing a vector space of type \'" << typeid(vec_space).name() << "\'"
			<< "\nwith a dimension of vec_space.dim() = " << vec_space.dim()
			<< std::endl;

	// Create the nlp
	ExampleNLPFirstOrder
		nlp(VectorSpace::space_ptr_t(&vec_space,false),xo,has_bounds,dep_bounded);

	// Create the solver object and set it up
	MoochoSolver solver;
	solver.set_nlp(Teuchos::rcp(&nlp,false));                  // Set nlp
	// set up outputting
	solver.set_error_handling(
		throw_solve_exception
		,Teuchos::rcp(error_out,false)
		);
	solver.set_console_out(Teuchos::rcp(console_out,false));
	solver.set_summary_out(Teuchos::rcp(summary_out,false));
	solver.set_journal_out(Teuchos::rcp(journal_out,false));
	solver.set_algo_out(   Teuchos::rcp(algo_out,false)   );

	// Run MOOCHO using the MamaJama configuration
	solve_return = solver.solve_nlp();

	return solve_return;
}
