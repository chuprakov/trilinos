// ///////////////////////////////////////////////////////////
// BasisSystemTester.cpp
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

#include <ostream>

#include "AbstractLinAlgPack/include/BasisSystemTester.h"
#include "AbstractLinAlgPack/include/BasisSystem.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"

namespace AbstractLinAlgPack {

BasisSystemTester::BasisSystemTester(
	EPrintTestLevel  print_tests
	,bool            dump_all
	,bool            throw_exception
	,size_type       num_random_tests
	,value_type      warning_tol
	,value_type      error_tol
	)
{}

bool BasisSystemTester::test_basis_system(
	const BasisSystem               &bs
	,const MatrixWithOp             *Gc
	,const MatrixWithOp             *Gh
	,const MatrixWithOpNonsingular  *C
	,const MatrixWithOp             *D
	,bool                           print_all_warnings
	,std::ostream                   *out
	)
{
	namespace rcp = ReferenceCountingPack;
	
	// ToDo: Check the preconditions
	
	bool success = false, result;
	
	// Print the input?
	if( out && print_tests() != PRINT_NONE ) {
		if( print_tests() >= PRINT_BASIC )
			*out << "\n*************************************************"
				 << "\n*** BasisSystemTester::test_basis_system(...) ***"
				 << "\n*************************************************\n";
		if(dump_all()) {
			if(Gc)
				*out << "\nGc =\n" << *Gc;
			if(Gh)
				*out << "\nGh =\n" << *Gh;
			if(C)
				*out << "\nC =\n"  << *C;
			if(D)
				*out << "\nD =\n"  << *D;
		}
	}

	//
	// Check the dimensions of everything
	//

	if( out && print_tests() >= PRINT_MORE ) {
		*out
			<< "\nbs.var_dep()        = ["<<bs.var_dep().lbound()<<","<<bs.var_dep().ubound()<<"]"
			<< "\nbs.var_indep( )     = ["<<bs.var_indep().lbound()<<","<<bs.var_indep().ubound()<<"]"
			<< "\nbs.equ_decomp()     = ["<<bs.equ_decomp().lbound()<<","<<bs.equ_decomp().ubound()<<"]"
			<< "\nbs.equ_undecomp()   = ["<<bs.equ_undecomp().lbound()<<","<<bs.equ_undecomp().ubound()<<"]"
			<< "\nbs.inequ_decomp()   = ["<<bs.inequ_decomp().lbound()<<","<<bs.inequ_decomp().ubound()<<"]"
			<< "\nbs.inequ_undecomp() = ["<<bs.inequ_undecomp().lbound()<<","<<bs.inequ_undecomp().ubound()<<"]"
			<< std::endl;
	}

	if( out && print_tests() >= PRINT_BASIC )
		*out << "\ncheck: bs.var_dep().size() != bs.equ_decomp().size() + bs.inequ_decomp().size() : ";
	result = bs.var_dep().size() == bs.equ_decomp().size() + bs.inequ_decomp().size();
	if( out && print_tests() >= PRINT_BASIC )
		*out << std::boolalpha << result << std::endl;

	if(Gc) {
		if( out && print_tests() >= PRINT_BASIC )
			*out << "\ncheck: bs.var_dep().size() + bs.var_indep().size() == Gc->rows() : ";
		result = bs.var_dep().size() + bs.var_indep().size() == Gc->rows();
		if( out && print_tests() >= PRINT_BASIC )
			*out << std::boolalpha << result << std::endl;
	}	

	if(Gh) {
		if( out && print_tests() >= PRINT_BASIC )
			*out << "\ncheck: bs.var_dep().size() + bs.var_indep().size() == Gh->rows() : ";
		result = bs.var_dep().size() + bs.var_indep().size() == Gh->rows();
		if( out && print_tests() >= PRINT_BASIC )
			*out << std::boolalpha << result << std::endl;
	}	

	if(Gc) {
		if( out && print_tests() >= PRINT_BASIC )
			*out << "\ncheck: bs.equ_decomp().size() + bs.equ_undecomp().size() == Gc->cols() : ";
		result = bs.equ_decomp().size() + bs.equ_undecomp().size() == Gc->cols();
		if( out && print_tests() >= PRINT_BASIC )
			*out << std::boolalpha << result << std::endl;
	}	

	if(Gh) {
		if( out && print_tests() >= PRINT_BASIC )
			*out << "\ncheck: bs.inequ_decomp().size() + bs.inequ_undecomp().size() == Gh->cols() : ";
		result = bs.inequ_decomp().size() + bs.inequ_undecomp().size() == Gh->cols();
		if( out && print_tests() >= PRINT_BASIC )
			*out << std::boolalpha << result << std::endl;
	}	

	// Create the N matrix
	rcp::ref_count_ptr<AbstractLinAlgPack::MatrixCompositeStd>
		N = new AbstractLinAlgPack::MatrixCompositeStd(bs.var_dep().size(),bs.var_indep().size());
	if( (Gc || Gh) && C ) {
		if( bs.equ_decomp().size() )
			N->add_matrix( 0, 0, 1.0, bs.equ_decomp(), Gc, NULL, BLAS_Cpp::trans, bs.var_indep() );
		if( Gh && bs.inequ_decomp().size() )
			N->add_matrix( bs.equ_decomp().size(), 0, 1.0, bs.inequ_decomp(), Gh, NULL, BLAS_Cpp::trans, bs.var_indep() );
		N->finish_construction(
			Gc->space_cols().sub_space(bs.var_indep())->clone()
			,Gc->space_rows().sub_space(bs.equ_decomp())->clone()
			);
		if( out && dump_all() )
			*out << "\nN =\n" << *N;
	}
	
	//
	// Perform the numerical tests
	//
	
	if(out)
		*out << "Warning, implementation is not finished!\n";
	return true;

	assert(0); // ToDo: Implement!
	return false;
}

} // end namespace AbstractLinAlgPack
