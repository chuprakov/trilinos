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

#include <math.h>

#include <ostream>

#include "ConstrainedOptimizationPack/include/DecompositionSystemTester.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystem.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingularTester.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"

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
 
bool DecompositionSystemTester::test_decomp_system(
	const DecompositionSystem       &ds
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
	namespace rcp = MemMngPack;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using AbstractLinAlgPack::sum;
	using AbstractLinAlgPack::dot;
	using AbstractLinAlgPack::Vp_StV;
	using AbstractLinAlgPack::Vp_StMtV;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using AbstractLinAlgPack::random_vector;
	using LinAlgOpPack::V_StMtV;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_VpV;
	using LinAlgOpPack::Vp_V;

	bool success = true, result, lresult, llresult;
	const value_type
		rand_y_l  = -1.0,
		rand_y_u  = 1.0,
		small_num = ::pow(std::numeric_limits<value_type>::epsilon(),0.25),
		alpha     = 2.0,
		beta      = 3.0;

	// Print the input?
	if( out && print_tests() != PRINT_NONE ) {
		if( print_tests() >= PRINT_BASIC )
			*out << "\n**********************************************************"
				 << "\n*** DecompositionSystemTester::test_decomp_system(...) ***"
				 << "\n**********************************************************\n";
	}

	const size_type
		n = ds.n(),
		m = ds.m(),
		r = ds.r();
	const Range1D
		con_decomp       = ds.con_decomp(),
		con_undecomp     = ds.con_undecomp();

	// print dimensions, ranges
	if( out && print_tests() >= PRINT_MORE ) {
		*out
			<< "\nds.n()                  = " << n
			<< "\nds.m()                  = " << m
			<< "\nds.r()                  = " << r
			<< "\nds.con_decomp()         = ["<<con_decomp.lbound()<<","<<con_decomp.ubound()<<"]"
			<< "\nds.con_undecomp()       = ["<<con_undecomp.lbound()<<","<<con_undecomp.ubound()<<"]"
			<< "\nds.space_range()->dim() = " << ds.space_range()->dim()
			<< "\nds.space_null()->dim()  = " << ds.space_null()->dim()
			<< std::endl;
	}

	// validate input matrices
	THROW_EXCEPTION(
	    Z==NULL&&Y==NULL&&R==NULL&&Uz==NULL&&Uy==NULL&&Vz==NULL&&Vy==NULL
		, std::invalid_argument
		,"DecompositionSystemTester::test_decomp_system(...) : Error, "
		"at least one of Z, Y, R, Uz, Uy, Vz or Vy can not be NULL NULL!" );
	THROW_EXCEPTION(
		m == r && Uz != NULL, std::invalid_argument
		,"DecompositionSystemTester::test_decomp_system(...) : Error, "
		"Uz must be NULL if m==r is NULL!" );
	THROW_EXCEPTION(
		m == r && Uy != NULL, std::invalid_argument
		,"DecompositionSystemTester::test_decomp_system(...) : Error, "
		"Uy must be NULL if m==r is NULL!" );
	THROW_EXCEPTION(
	    Gh == NULL && Vz != NULL, std::invalid_argument
		,"DecompositionSystemTester::test_decomp_system(...) : Error, "
		"Vz must be NULL if Gh is NULL!" );
	THROW_EXCEPTION(
	    Gh == NULL && Vy != NULL, std::invalid_argument
		,"DecompositionSystemTester::test_decomp_system(...) : Error, "
		"Vy must be NULL if Gh is NULL!" );

	// Print the input?
	if( out && print_tests() != PRINT_NONE ) {
		if(dump_all()) {
			*out << "\nGc =\n"       << Gc;
			if(Gh)
				*out << "\nGh =\n"   << *Gh;
			if(Z)
				*out << "\nZ =\n"    << *Z;
			if(Y)
				*out << "\nY =\n"    << *Y;
			if(R)
				*out << "\nR =\n"    << *R;
			if(Uz)
				*out << "\nUz =\n"   << *Uz;
			if(Uy)
				*out << "\nUy =\n"   << *Uy;
			if(Vz)
				*out << "\nVz =\n"   << *Vz;
			if(Vy)
				*out << "\nVy =\n"   << *Vy;
		}
	}

	//
	// Check the dimensions of everything
	//

	if( out && print_tests() >= PRINT_BASIC )
		*out << "\n1) Check the partitioning ranges and vector space dimensions ...";
	lresult = true;

	if( out && print_tests() >= PRINT_MORE )
		*out << "\n\n1.a) check: con_decomp.size() + con_undecomp.size() == ds.m() : ";
	result = con_decomp.size() + con_undecomp.size() == ds.m();
	if(out && print_tests() >= PRINT_MORE)
		*out << ( result ? "passed" : "failed" );
	if(!result) lresult = false;

	if( out && print_tests() >= PRINT_MORE )
		*out << "\n\n1.b) check: con_decomp.size() == ds.r() : ";
	result = con_decomp.size() == ds.r();
	if(out && print_tests() >= PRINT_MORE)
		*out << ( result ? "passed" : "failed" );
	if(!result) lresult = false;

	if( out && print_tests() >= PRINT_MORE )
		*out << "\n\n1.c) check: ds.space_range()->dim() == ds.r() : ";
	result = ds.space_range()->dim() == ds.r();
	if(out && print_tests() >= PRINT_MORE)
		*out << ( result ? "passed" : "failed" );
	if(!result) lresult = false;

	if( out && print_tests() >= PRINT_MORE )
		*out << "\n\n1.d) check: ds.space_null()->dim() == ds.n() - ds.r() : ";
	result = ds.space_null()->dim() == ds.n() - ds.r();
	if(out && print_tests() >= PRINT_MORE)
		*out << ( result ? "passed" : "failed" );
	if(!result) lresult = false;

	if(out && print_tests() >= PRINT_MORE)
		*out << std::endl;

	if(!lresult) success = false;
	if( out && print_tests() == PRINT_BASIC )
		*out << " : " << ( lresult ? "passed" : "failed" );

	//
	// Perform the tests
	//

	if(out && print_tests() >= PRINT_BASIC)
		*out
			<< "\n2) Check the compatibility of the vector spaces for Gc, Gh Z, Y, R, Uz, Uy, Vz and Vy  ...";
	lresult = true;
	
	if(Z) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< "\n2.a) Check consistency of the vector spaces for:"
				<< "\n    Z.space_cols() == Gc.space_cols() and Z.space_rows() == ds.space_null()";
		llresult = true;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.a.1) Z->space_cols().is_compatible(Gc.space_cols()) == true : ";
		result = Z->space_cols().is_compatible(Gc.space_cols());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.a.2) Z->space_cols().is_compatible(*ds.space_null()) == true : ";
		result = Z->space_rows().is_compatible(*ds.space_null());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );
	}

	if(Y) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< "\n2.b) Check consistency of the vector spaces for:"
				<< "\n    Y.space_cols() == Gc.space_cols() and Y.space_rows() == ds.space_range()";
		llresult = true;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.b.1) Y->space_cols().is_compatible(Gc.space_cols()) == true : ";
		result = Y->space_cols().is_compatible(Gc.space_cols());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.b.2) Y->space_cols().is_compatible(*ds.space_range()) == true : ";
		result = Y->space_rows().is_compatible(*ds.space_range());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );
	}

	if(R) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< "\n2.c) Check consistency of the vector spaces for:"
				<< "\n    R.space_cols() == Gc.space_cols()(con_decomp) and R.space_rows() == ds.space_range()";
		llresult = true;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.c.1) R->space_cols().is_compatible(*Gc.space_cols().sub_space(con_decomp)) == true : ";
		result = R->space_cols().is_compatible(*Gc.space_cols().sub_space(con_decomp));	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.c.2) R->space_cols().is_compatible(*ds.space_range()) == true : ";
		result = R->space_rows().is_compatible(*ds.space_range());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );
	}

	if(Uz) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< "\n2.d) Check consistency of the vector spaces for:"
				<< "\n    Uz.space_cols() == Gc.space_cols()(con_undecomp) and Uz.space_rows() == ds.space_null()";
		llresult = true;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.d.1) Uz->space_cols().is_compatible(*Gc.space_cols().sub_space(con_undecomp)) == true : ";
		result = Uz->space_cols().is_compatible(*Gc.space_cols().sub_space(con_undecomp));	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.d.2) Uz->space_cols().is_compatible(*ds.space_null()) == true : ";
		result = Uz->space_rows().is_compatible(*ds.space_null());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );
	}

	if(Uy) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< "\n2.e) Check consistency of the vector spaces for:"
				<< "\n    Uy.space_cols() == Gc.space_cols()(con_undecomp) and Uy.space_rows() == ds.space_range()";
		llresult = true;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.e.1) Uy->space_cols().is_compatible(*Gc.space_cols().sub_space(con_undecomp)) == true : ";
		result = Uy->space_cols().is_compatible(*Gc.space_cols().sub_space(con_undecomp));	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.e.2) Uy->space_cols().is_compatible(*ds.space_range()) == true : ";
		result = Uy->space_rows().is_compatible(*ds.space_range());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );
	}

	if(Vz) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< "\n2.f) Check consistency of the vector spaces for:"
				<< "\n    Vz.space_cols() == Gh.space_cols() and Vz.space_rows() == ds.space_null()";
		llresult = true;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.f.1) Vz->space_cols().is_compatible(Gh->space_cols()) == true : ";
		result = Vz->space_cols().is_compatible(Gh->space_cols());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.f.2) Vz->space_cols().is_compatible(*ds.space_null()) == true : ";
		result = Vz->space_rows().is_compatible(*ds.space_null());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );
	}

	if(Vy) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< "\n2.g) Check consistency of the vector spaces for:"
				<< "\n    Vy.space_cols() == Gh.space_cols() and Vy.space_rows() == ds.space_range()";
		llresult = true;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.g.1) Vy->space_cols().is_compatible(Gh->space_cols()) == true : ";
		result = Vy->space_cols().is_compatible(Gh->space_cols());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(out && print_tests() >= PRINT_ALL)
			*out << "\n\n2.g.2) Vy->space_cols().is_compatible(*ds.space_range()) == true : ";
		result = Vy->space_rows().is_compatible(*ds.space_range());	
		if(out && print_tests() >= PRINT_ALL)
			*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );
	}

	if(!lresult) success = false;
	if( out && print_tests() == PRINT_BASIC )
		*out << " : " << ( lresult ? "passed" : "failed" );

	if(out && print_tests() >= PRINT_BASIC)
		*out
			<< "\n3) Check the compatibility of the matrices Gc, Gh Z, Y, R, Uz, Uy, Vz and Vy numerically ...";
	
	if(Z) {

		if(out && print_tests() >= PRINT_MORE)
			*out
				<< std::endl
				<< "\n3.a) Check consistency of:"
				<< "\n     op ( alpha*[ Gc(:,con_decomp)'   ]"
				<< "\n                [ Gc(:,con_undecomp)' ] * beta*Z ) * v"
				<< "\n         \\_____________________________________/"
				<< "\n                         A"
				<< "\n    ==  op( alpha*beta*[ 0  ]"
				<< "\n                       [ Uz ] ) * v"
				<< "\n            \\_______________/"
				<< "\n                     B"
				<< "\nfor random vectors v ...";

		VectorSpace::vec_mut_ptr_t
			v_c       = Gc.space_rows().create_member(),
			v_c_tmp   = v_c->space().create_member(),
			v_x       = Gc.space_cols().create_member(),
			v_x_tmp   = v_x->space().create_member(),
			v_z       = ds.space_null()->create_member(),
			v_z_tmp   = v_z->space().create_member();
		
		if(out && print_tests() >= PRINT_MORE)
			*out << "\n\n3.a.1) Testing non-transposed A*v == B*v ...";
		if(out && print_tests() > PRINT_MORE)
			*out << std::endl;
		llresult = true;
		{for( int k = 1; k <= num_random_tests(); ++k ) {
			random_vector( rand_y_l, rand_y_u, v_z.get() );
			if(out && print_tests() >= PRINT_ALL) {
				*out
					<< "\n3.a.1."<<k<<") random vector " << k << " ( ||v_z||_1 / n = " << (v_z->norm_1() / v_z->dim()) << " )\n";
				if(dump_all() && print_tests() >= PRINT_ALL)
					*out << "\nv_z =\n" << *v_z;
			}
			V_StMtV( v_x.get(), beta, *Z, no_trans, *v_z );
			V_StMtV( v_c.get(), alpha, Gc, trans, *v_x );
			*v_c_tmp->sub_view(con_decomp) = 0.0;
			if(con_undecomp.size()) {
				if(Uz)
					V_StMtV( v_c_tmp->sub_view(con_undecomp).get(), alpha*beta, *Uz, no_trans, *v_z );
				else
					*v_c_tmp->sub_view(con_undecomp).get() = *v_c->sub_view(con_undecomp);
			}
			const value_type
				sum_Bv  = sum(*v_c_tmp), // should be zero if con_undecomp.size() == 0 so scale by 1.0
				sum_Av  = sum(*v_c);
			assert_print_nan_inf(sum_Bv, "sum(B*v_z)",true,out);
			assert_print_nan_inf(sum_Av, "sum(A*v_z)",true,out);
			const value_type
				calc_err = ::fabs( ( sum_Av - sum_Bv )
								   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + (con_undecomp.size() ? small_num : 1.0) ) );
			if(out && print_tests() >= PRINT_ALL)
				*out
					<< "\nrel_err(sum(A*v_z),sum(B*v_z)) = "
					<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
					<< calc_err << std::endl;
			if( calc_err >= mult_warning_tol() ) {
				if(out && print_tests() >= PRINT_ALL)
					*out
						<< std::endl
						<< ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
						<< ", rel_err(sum(A*v_z),sum(B*v_z)) = "
						<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
						<< calc_err
						<< " exceeded "
						<< ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
						<< " = "
						<< ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
						<< std::endl;
				if(calc_err >= mult_error_tol()) {
					if(dump_all() && print_tests() >= PRINT_ALL) {
						*out << "\nalpha = " << alpha << std::endl;
						*out << "\nbeta  = " << beta  << std::endl;
						*out << "\nv_z =\n"           << *v_z;
						*out << "\nbeta*Z*v_z =\n"    << *v_x;
						*out << "\nA*v_z =\n"         << *v_c;
						*out << "\nB*v_z =\n"         << *v_c_tmp;
					}
					llresult = false;
				}
			}
		}}
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" )
				 << std::endl;
		
		if(out && print_tests() >= PRINT_MORE)
			*out << "\n\n3.a.2) Testing transposed A'*v == B'*v ...";
		if(out && print_tests() > PRINT_MORE)
			*out << std::endl;
		llresult = true;
		{for( int k = 1; k <= num_random_tests(); ++k ) {
			random_vector( rand_y_l, rand_y_u, v_c.get() );
			if(out && print_tests() >= PRINT_ALL) {
				*out
					<< "\n3.a.2."<<k<<") random vector " << k << " ( ||v_c||_1 / n = " << (v_c->norm_1() / v_c->dim()) << " )\n";
				if(dump_all() && print_tests() >= PRINT_ALL)
					*out << "\nv_c =\n" << *v_c;
			}
			V_StMtV( v_x.get(), alpha, Gc, no_trans, *v_c );
			V_StMtV( v_z.get(), beta,  *Z, trans,    *v_x );
			*v_z_tmp = 0.0;
			if(con_undecomp.size()) {
				if(Uz)
					V_StMtV( v_z_tmp.get(), alpha*beta, *Uz, trans, *v_c->sub_view(con_undecomp) );
				else
					*v_z_tmp = *v_z;
			}
			const value_type
				sum_Bv  = sum(*v_z_tmp), // should be zero so scale by 1.0
				sum_Av  = sum(*v_z);
			assert_print_nan_inf(sum_Bv, "sum(B'*v_c)",true,out);
			assert_print_nan_inf(sum_Av, "sum(A'*v_c)",true,out);
			const value_type
				calc_err = ::fabs( ( sum_Av - sum_Bv )
								   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + (con_undecomp.size() ? small_num : 1.0) ) );
			if(out && print_tests() >= PRINT_ALL)
				*out
					<< "\nrel_err(sum(A'*v_c),sum(B'*v_c)) = "
					<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
					<< calc_err << std::endl;
			if( calc_err >= mult_warning_tol() ) {
				if(out && print_tests() >= PRINT_ALL)
					*out
						<< std::endl
						<< ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
						<< ", rel_err(sum(A'*v_c),sum(B'*v_c)) = "
						<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
						<< calc_err
						<< " exceeded "
						<< ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
						<< " = "
						<< ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
						<< std::endl;
				if(calc_err >= mult_error_tol()) {
					if(dump_all() && print_tests() >= PRINT_ALL) {
						*out << "\nalpha = " << alpha << std::endl;
						*out << "\nbeta  = " << beta  << std::endl;
						*out << "\nv_c =\n"           << *v_c;
						*out << "\nalpha*Gc*v_c =\n"  << *v_x;
						*out << "\nA'*v_c =\n"        << *v_z;
						*out << "\nB'*v_c =\n"        << *v_z_tmp;
					}
					llresult = false;
				}
			}
		}}
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" )
				 << std::endl;

	}
	else {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< std::endl
				<< "\n3.a) Warning! Z ==NULL; Z, and Uz are not checked numerically ...\n";
	}

	if(Y) {

		if(out && print_tests() >= PRINT_MORE)
			*out
				<< std::endl
				<< "\n3.b) Check consistency of:"
				<< "\n     op ( alpha*[ Gc(:,con_decomp)'   ]"
				<< "\n                [ Gc(:,con_undecomp)' ] * beta*Y ) * v"
				<< "\n         \\_____________________________________/"
				<< "\n                         A"
				<< "\n    ==  op( alpha*beta*[ R  ]"
				<< "\n                       [ Uy ] ) * v"
				<< "\n            \\_______________/"
				<< "\n                     B"
				<< "\nfor random vectors v ...";

		VectorSpace::vec_mut_ptr_t
			v_c       = Gc.space_rows().create_member(),
			v_c_tmp   = v_c->space().create_member(),
			v_x       = Gc.space_cols().create_member(),
			v_x_tmp   = v_x->space().create_member(),
			v_y       = ds.space_range()->create_member(),
			v_y_tmp   = v_y->space().create_member();
		
		if(out && print_tests() >= PRINT_MORE)
			*out << "\n\n3.b.1) Testing non-transposed A*v == B*v ...";
		if(out && print_tests() > PRINT_MORE)
			*out << std::endl;
		llresult = true;
		{for( int k = 1; k <= num_random_tests(); ++k ) {
			random_vector( rand_y_l, rand_y_u, v_y.get() );
			if(out && print_tests() >= PRINT_ALL) {
				*out
					<< "\n3.b.1."<<k<<") random vector " << k << " ( ||v_y||_1 / n = " << (v_y->norm_1() / v_y->dim()) << " )\n";
				if(dump_all() && print_tests() >= PRINT_ALL)
					*out << "\nv_y =\n" << *v_y;
			}
			V_StMtV( v_x.get(), beta, *Y, no_trans, *v_y );
			V_StMtV( v_c.get(), alpha, Gc, trans, *v_x );
			V_StMtV( v_c_tmp->sub_view(con_decomp).get(), alpha*beta, *R, no_trans, *v_y );
			if(con_undecomp.size()) {
				if(Uy)
					V_StMtV( v_c_tmp->sub_view(con_undecomp).get(), alpha*beta, *Uy, no_trans, *v_y );
				else
					*v_c_tmp->sub_view(con_undecomp) = *v_c->sub_view(con_undecomp);
			}
			const value_type
				sum_Bv  = sum(*v_c_tmp),
				sum_Av  = sum(*v_c);
			assert_print_nan_inf(sum_Bv, "sum(B*v_y)",true,out);
			assert_print_nan_inf(sum_Av, "sum(A*v_y)",true,out);
			const value_type
				calc_err = ::fabs( ( sum_Av - sum_Bv )
								   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + small_num ) );
			if(out && print_tests() >= PRINT_ALL)
				*out
					<< "\nrel_err(sum(A*v_y),sum(B*v_y)) = "
					<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
					<< calc_err << std::endl;
			if( calc_err >= mult_warning_tol() ) {
				if(out && print_tests() >= PRINT_ALL)
					*out
						<< std::endl
						<< ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
						<< ", rel_err(sum(A*v_y),sum(B*v_y)) = "
						<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
						<< calc_err
						<< " exceeded "
						<< ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
						<< " = "
						<< ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
						<< std::endl;
				if(calc_err >= mult_error_tol()) {
					if(dump_all() && print_tests() >= PRINT_ALL) {
						*out << "\nalpha = " << alpha << std::endl;
						*out << "\nbeta  = " << beta  << std::endl;
						*out << "\nv_y =\n"           << *v_y;
						*out << "\nbeta*Y*v_y =\n"    << *v_x;
						*out << "\nA*v_y =\n"         << *v_c;
						*out << "\nB*v_y =\n"         << *v_c_tmp;
					}
					llresult = false;
				}
			}
		}}
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" )
				 << std::endl;
		
		if(out && print_tests() >= PRINT_MORE)
			*out << "\n\n3.b.2) Testing transposed A'*v == B'*v ...";
		if(out && print_tests() > PRINT_MORE)
			*out << std::endl;
		llresult = true;
		{for( int k = 1; k <= num_random_tests(); ++k ) {
			random_vector( rand_y_l, rand_y_u, v_c.get() );
			if(out && print_tests() >= PRINT_ALL) {
				*out
					<< "\n3.a.2."<<k<<") random vector " << k << " ( ||v_c||_1 / n = " << (v_c->norm_1() / v_c->dim()) << " )\n";
				if(dump_all() && print_tests() >= PRINT_ALL)
					*out << "\nv_c =\n" << *v_c;
			}
			V_StMtV( v_x.get(), alpha, Gc, no_trans, *v_c );
			V_StMtV( v_y.get(), beta,  *Y, trans,    *v_x );
			V_StMtV( v_y_tmp.get(), alpha*beta, *R, trans, *v_c->sub_view(con_decomp) );
			if(con_undecomp.size()) {
				if(Uy)
					Vp_StMtV( v_y_tmp.get(), alpha*beta, *Uy, trans, *v_c->sub_view(con_undecomp) );
				else
					Vp_V( v_y_tmp.get(), *v_y );
			}
			const value_type
				sum_Bv  = sum(*v_y_tmp), // should be zero so scale by 1.0
				sum_Av  = sum(*v_y);
			assert_print_nan_inf(sum_Bv, "sum(B'*v_c)",true,out);
			assert_print_nan_inf(sum_Av, "sum(A'*v_c)",true,out);
			const value_type
				calc_err = ::fabs( ( sum_Av - sum_Bv )
								   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + small_num ) );
			if(out && print_tests() >= PRINT_ALL)
				*out
					<< "\nrel_err(sum(A'*v_c),sum(B'*v_c)) = "
					<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
					<< calc_err << std::endl;
			if( calc_err >= mult_warning_tol() ) {
				if(out && print_tests() >= PRINT_ALL)
					*out
						<< std::endl
						<< ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
						<< ", rel_err(sum(A'*v_c),sum(B'*v_c)) = "
						<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
						<< calc_err
						<< " exceeded "
						<< ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
						<< " = "
						<< ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
						<< std::endl;
				if(calc_err >= mult_error_tol()) {
					if(dump_all() && print_tests() >= PRINT_ALL) {
						*out << "\nalpha = " << alpha << std::endl;
						*out << "\nbeta  = " << beta  << std::endl;
						*out << "\nv_c =\n"           << *v_c;
						*out << "\nalpha*Gc*v_c =\n"  << *v_x;
						*out << "\nA'*v_c =\n"        << *v_y;
						*out << "\nB'*v_c =\n"        << *v_y_tmp;
					}
					llresult = false;
				}
			}
		}}
		if(!llresult) lresult = false;
		if( out && print_tests() == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" )
				 << std::endl;

	}
	else {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< std::endl
				<< "\n3.b) Warning! Y ==NULL; Y, R and Uy are not checked numerically ...\n";
	}

	assert(Vz == NULL && Vy == NULL); // ToDo: 3.c) Check Vz and Vy

	if(R) {
		if(out && print_tests() >= PRINT_MORE)
			*out
				<< std::endl
				<< "\n3.b) Check consistency of: op(op(inv(R))*op(R)) == I ...\n";
		typedef MatrixWithOpNonsingularTester  MWONST_t;
		MWONST_t::EPrintTestLevel
			olevel;
		switch(print_tests()) {
			case PRINT_NONE:
			case PRINT_BASIC:
				olevel = MWONST_t::PRINT_NONE;
				break;
			case PRINT_MORE:
				olevel = MWONST_t::PRINT_MORE;
				break;
			case PRINT_ALL:
				olevel = MWONST_t::PRINT_ALL;
				break;
			default:
				assert(0); // Should not get here
		}
		MWONST_t
			R_tester(
				MWONST_t::TEST_LEVEL_2_BLAS
				,olevel
				,dump_all()
				,throw_exception()
				,num_random_tests()
				,solve_warning_tol()
				,solve_error_tol()
				);
		lresult = R_tester.test_matrix(*R,"R",out);
	}

	if(!lresult) success = false;
	if( out && print_tests() == PRINT_BASIC )
		*out << " : " << ( lresult ? "passed" : "failed" );
	
	if( out && print_tests() != PRINT_NONE ) {
		if(success)
			*out << "\nCongradulations! The DecompositionSystem object and its associated matrix objects seem to check out!\n";
		else
			*out << "\nOops! At last one of the tests did not check out!\n";
		if( print_tests() >= PRINT_BASIC )
			*out << "\nEnd DecompositionSystemTester::test_decomp_system(...)\n";
	}

	return success;
}

} // end namespace ConstrainedOptimizationPack
