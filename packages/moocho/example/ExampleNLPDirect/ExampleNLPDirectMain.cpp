// //////////////////////////////////////////////////////////
// ExampleNLPFirstOrderDirectMain.cpp
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

#include "ExampleNLPFirstOrderDirectRun.h"
#include "ExampleVectorLib/include/MPIDenseVector.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "OptionsFromStream.h"
#include "WorkspacePack.h"

int main(int argc, char* argv[] ) {

	using std::endl;
	using std::setw;
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::OptionsFromStream;
	typedef AbstractLinAlgPack::size_type size_type;
	typedef AbstractLinAlgPack::value_type value_type;
	namespace NLPIP = NLPInterfacePack;

	using AbstractLinAlgPack::VectorSpace;
	using AbstractLinAlgPack::VectorWithOp;
	using AbstractLinAlgPack::VectorWithOpMutable;

	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	int err = 0;

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
	bool has_bounds = true;
	// Make the dependent or independent variables bounded.
	bool dep_bounded = true;

	// Read from the arguments
	if(argc > 1)
		n = ::atoi(argv[1]);
	if(argc > 2)
		xo = ::atof(argv[2]);
	if(argc > 3)
		has_bounds = (::atoi(argv[3]) == 1);
	if(argc > 4)
		dep_bounded = (::atoi(argv[4]) == 1);

	// Set the output stream
	char console_out_name[20];
	::sprintf( console_out_name, "console.%d.out", proc_rank );
	std::ofstream console_out(console_out_name);
	std::ostream // Only send output to the root process!
		&out  = ( proc_rank == 0 ? std::cout : console_out );
	std::ostream // Only send output to the root process!
		&eout = ( proc_rank == 0 ? std::cerr : console_out );

	try {
	
	int w = 15;
	int prec = 8;

	out
		<< std::setprecision(prec)
		<< std::scientific
		<< "********************************************************\n"
		<< "*** Running Tests on Example First Order Dirrect NLP ***\n"
		<< "********************************************************\n";

	// Create the vector space object to use.
	VectorSpace::space_ptr_t    vec_space;

	// Determine the mapping of elements to processors for MPIDenseVectorSpace
	RTOp_index_type local_dim = n/num_proc; // assume n > num_proc
	wsp::Workspace<RTOp_index_type>  ind_map(wss,num_proc);
	{
		RTOp_index_type i_u = local_dim;
		{for( int p = 0; p < num_proc; ++p, ++i_u ) {
			ind_map[p] = i_u;
		}}
		ind_map[num_proc-1] = n;
		local_dim = ( proc_rank > 0
					  ? ind_map[proc_rank]-ind_map[proc_rank-1]
					  : ind_map[0] );
	}
	
	vec_space = rcp::rcp_implicit_cast<const VectorSpace>(
		ref_count_ptr<MPIDenseVectorSpace>(
			new MPIDenseVectorSpace(MPI_COMM_WORLD,&ind_map[0],1,n)
			)
		);

	// Create and test the NLP using this vector space object
	const bool
		result = NLPIP::ExampleNLPFirstOrderDirectRun(
			*vec_space, xo, has_bounds, dep_bounded
			,&out,&eout
			);
	if(!result)
		prog_return = PROG_NLP_TEST_ERR;

	}	// end try
	catch(const std::exception& excpt) {
		out << "\nCaught a std::exception: " << excpt.what() << endl;
		prog_return = PROG_EXCEPTION;
	}
	catch(...) {
		out << "\nCaught an unknown exception\n";
		prog_return = PROG_EXCEPTION;
	}

	if(prog_return == PROG_SUCCESS) {
		eout   << "Congradulations!  NLP class seems to check out!\n";
		out    << "\nCongradulations!  NLP class seems to check out!\n";
	}
	else {
		eout   << "Oh No!  Something did not checkout!\n";
		out    << "\nOh No!  Something did not checkout!\n";
	}

 	MPI_Finalize();

	return prog_return;
}
