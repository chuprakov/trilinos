// //////////////////////////////////////////////////////////
// ExampleNLPFirstOrderInfoMain.cpp
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
//

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "ExampleNLPFirstOrderInfoRun.hpp"
#include "ExampleVectorLib/src/MPIDenseVector.hpp"
#include "SparseLinAlgPack/src/VectorSpaceSerial.hpp"
#include "OptionsFromStream.hpp"
#include "WorkspacePack.hpp"
#include "oblackholestream.hpp"
#include "CommandLineProcessor.hpp"

int main(int argc, char* argv[] ) {

	using std::endl;
	using std::setw;
	namespace rcp = MemMngPack;
	using rcp::ref_count_ptr;
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::OptionsFromStream;
	typedef AbstractLinAlgPack::size_type size_type;
	typedef AbstractLinAlgPack::value_type value_type;
	namespace NLPIP = NLPInterfacePack;
	namespace rsqp = ReducedSpaceSQPPack;
	using rsqp::rSQPppSolver;

	using AbstractLinAlgPack::VectorSpace;
	using AbstractLinAlgPack::Vector;
	using AbstractLinAlgPack::VectorMutable;

	using CommandLineProcessorPack::CommandLineProcessor;

	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

/*
	// Print out the input arguments
	printf("argc = %d\n",argc);
	{for( int i = 0; i < argc; ++i) {
		printf("argv[%d] = %s\n",i,argv[i]);
	}}
*/
	// Get an idea of what processors we have.
	MPI_Init(&argc,&argv);
	int num_proc, proc_rank;
	MPI_Comm_size( MPI_COMM_WORLD, &num_proc );
	MPI_Comm_rank( MPI_COMM_WORLD, &proc_rank );
/*
	// Print out the input arguments
	printf("\nproc_rank = %d\n",proc_rank);
	printf("argc = %d\n",argc);
	{for( int i = 0; i < argc; ++i) {
		printf("argv[%d] = %s\n",i,argv[i]);
	}}
*/

	// Define program return values
	const int
		PROG_SUCCESS				=  0,
		PROG_NLP_TEST_ERR			= -1,
		PROG_EXCEPTION				= -2,
		PROG_MAX_ITER_EXEEDED		= -3,
		PROG_MAX_TIME_EXEEDED		= -4;

	int prog_return = PROG_SUCCESS;

	// Get the size of the problem to solve
	size_type n = 4;
	// Get the starting point
	value_type xo = 0.1;
	// Determine if the NLP has bounds or not.
	bool has_bounds = false;
	// Make the dependent or independent variables bounded.
	bool dep_bounded = true;
	// Serial or parallel?
	bool in_parallel = false;
	
	CommandLineProcessor  command_line_processor;
	
	command_line_processor.set_option( "n",  &n,   "Global number of dependent (and independent) variables" );
	command_line_processor.set_option( "xo", &xo,  "Initial guess of the solution" );
	command_line_processor.set_option(
		"has-bounds", "no-has-bounds", &has_bounds
		,"Determine if the NLP has bounds or not" );
	command_line_processor.set_option(
		"dep-bounded", "indep-bounded", &dep_bounded
		,"Determine if the dependent or independent variables are bounded" );
	command_line_processor.set_option(
		"in-parallel", "in-serial", &in_parallel
		,"Determine if computations are performed in parallel or not" );
	
	CommandLineProcessor::EParseCommandLineReturn
		parse_return = command_line_processor.parse_command_line(argc,argv,&std::cerr);
	
	if( parse_return != CommandLineProcessor::PARSE_SUCCESSFULL )
		return parse_return;

	// Set the output stream
	std::ostream &out  = std::cout;
	std::ostream &eout = std::cerr;
	IOStreamHelperPack::oblackholestream  blackhole;

	try {
	
	// Create the vector space object to use.
	VectorSpace::space_ptr_t    vec_space;

	if(in_parallel) {
		//
		// Use parallel vectors!
		//
		// Determine the mapping of elements to processors for MPIDenseVectorSpace
		RTOp_index_type local_dim = n/num_proc; // assume n > num_proc
		RTOp_index_type *ind_map  = new RTOp_index_type[num_proc];
		RTOp_index_type i_u = local_dim;
		for( int p = 0; p < num_proc; ++p, i_u += local_dim )
			ind_map[p] = i_u;
		ind_map[num_proc-1] = n; // Make sure we don't go past n
		local_dim = ( proc_rank > 0
					  ? ind_map[proc_rank]-ind_map[proc_rank-1]
					  : ind_map[0] );
		vec_space = rcp::rcp(new MPIDenseVectorSpace(MPI_COMM_WORLD,ind_map,true,1,n));
	}
	else {
		//
		// Use serial vectors
		//
		vec_space = rcp::rcp(new SparseLinAlgPack::VectorSpaceSerial(n));
	}

	// Create and test the NLP using this vector space object
	const rSQPppSolver::ESolutionStatus
		solve_return = NLPIP::ExampleNLPFirstOrderInfoRun(
			*vec_space, xo, has_bounds, dep_bounded
			,proc_rank == 0 ? &out  : &blackhole  // console_out
			,proc_rank == 0 ? &eout : &blackhole  // error_out
			,proc_rank == 0 ? false : true        // throw_solve_exception
			,proc_rank == 0 ? NULL  : &blackhole  // algo_out
			,proc_rank == 0 ? NULL  : &blackhole  // summary_out
			,proc_rank == 0 ? NULL  : &blackhole  // journal_out
			);

	switch(solve_return) {
		case rSQPppSolver::SOLVE_RETURN_SOLVED:
			prog_return = PROG_SUCCESS;
			break;
		case rSQPppSolver::SOLVE_RETURN_MAX_ITER:
			prog_return = PROG_MAX_ITER_EXEEDED;
			break;
		case rSQPppSolver::SOLVE_RETURN_MAX_RUN_TIME:
			prog_return = PROG_MAX_TIME_EXEEDED;
			break;
		case rSQPppSolver::SOLVE_RETURN_NLP_TEST_FAILED:
			prog_return = PROG_NLP_TEST_ERR;
			break;
		case rSQPppSolver::SOLVE_RETURN_EXCEPTION:
			prog_return = PROG_EXCEPTION;
			break;
		default:
			assert(0);
	}

	}	// end try
	catch(const std::exception& excpt) {
		eout << "\nCaught a std::exception on process " << proc_rank<< ": " << excpt.what() << endl;
		prog_return = PROG_EXCEPTION;
	}
	catch(...) {
		eout << "\nCaught an unknown exception on process " << proc_rank<< "\n";
		prog_return = PROG_EXCEPTION;
	}

 	MPI_Finalize();

	return prog_return;
}
