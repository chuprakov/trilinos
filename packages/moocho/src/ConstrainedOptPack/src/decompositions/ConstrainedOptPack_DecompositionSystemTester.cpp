// ///////////////////////////////////////////////////////////
// DecompositionSystemTester.cpp
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

#include "ConstrainedOptimizationPack/include/DecompositionSystemTester.h"

namespace ConstrainedOptimizationPack {

DecompositionSystemTester::DecompositionSystemTester(
	EPrintTestLevel  print_tests
	,bool            dump_all
	,bool            throw_exception
	,size_type       num_random_tests
	,value_type      mult_warning_tol
	,value_type      mult_error_tol
	,value_type      solve_warning_tol
	,value_type      solve_error_tol
	)
	:print_tests_(print_tests)
	,dump_all_(dump_all)
	,throw_exception_(throw_exception)
	,num_random_tests_(num_random_tests)
	,mult_warning_tol_(mult_warning_tol)
	,mult_error_tol_(mult_error_tol)
	,solve_warning_tol_(solve_warning_tol)
	,solve_error_tol_(solve_error_tol)
{}
 
bool DecompositionSystemTestertest_decomp_system(
	const DecompositionSystem       &decomp_sys
	,const MatrixWithOp             &Gc
	,const MatrixWithOp             *Gh
	,const MatrixWithOp             *Z
	,const MatrixWithOp             *Y
	,const MatrixWithOpNonsingular  *R
	,const MatrixWithOp             *Uz
	,const MatrixWithOp             *Uy
	,const MatrixWithOp             *Vz
	,const MatrixWithOp             *Vy
	,std::ostream                   *out
	)
{
	assert(0); // ToDo: Implement!
	return false;
}

} // end namespace ConstrainedOptimizationPack
