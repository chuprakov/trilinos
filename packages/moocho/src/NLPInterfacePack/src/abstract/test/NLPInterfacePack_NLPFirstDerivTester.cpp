// ///////////////////////////////////////////////////////////
// NLPFirstDerivativesTester.cpp

#include <assert.h>
#include <iomanip>
#include <sstream>

#include "NLPFirstDerivativesTester.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "NLPInterfacePack/include/CalcFiniteDiffFirstDerivativeProduct.h"
#include "NLPInterfacePack/include/CalcFiniteDiffFirstDerivatives.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/test/CompareDenseVectors.h"
#include "LinAlgPack/include/random_vector.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/MatVecCompare.h"
#include "LinAlgPack/include/assert_print_nan_inf.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StMtV;
}

namespace NLPInterfacePack {
namespace TestingPack {

NLPFirstDerivativesTester::NLPFirstDerivativesTester(
		  ETestingMethod		fd_testing_method
		, size_type				num_fd_directions
		, value_type			warning_tol
		, value_type			error_tol
		, const comp_Gc_ptr_t&	comp_Gc		)
	:
		  fd_testing_method_(fd_testing_method)
		, num_fd_directions_(num_fd_directions)
		, warning_tol_(warning_tol)
		, error_tol_(error_tol)
		, comp_Gc_(comp_Gc)
{}	

bool NLPFirstDerivativesTester::finite_diff_check(
	  NLP						*nlp
	, const VectorSlice			&xo
	, const SpVectorSlice		*xl
	, const SpVectorSlice		*xu
	, const value_type			&max_var_bounds_viol
	, const MatrixWithOp		*Gc
	, const VectorSlice			*Gf
	, bool						print_all_warnings
	, std::ostream				*out
	) const
{
	// ///////////////////////////////////
	// Validate the input

	// Check the input vectors
	assert_print_nan_inf(xo, "xo",true,out);
	if(Gf)
		assert_print_nan_inf(*Gf, "Gf",true,out); 

	bool success = true;

	try {

	switch(fd_testing_method_) {
		case FD_COMPUTE_ALL:
			return fd_check_all(nlp,xo,xl,xu,max_var_bounds_viol,Gc,Gf
				,print_all_warnings,out);
		case FD_DIRECTIONAL:
			return fd_directional_check(nlp,xo,xl,xu,max_var_bounds_viol,Gc,Gf
				,print_all_warnings,out);
		default:
			assert(0);
	}

	} // end try
	catch( const LinAlgPack::NaNInfException& except ) {
		if(out)
			*out
				<< "Error, found a NaN or Inf.  Stoping tests\n";
		success = false;
	}
	
	return success;	// will never be executed
}

// private

bool NLPFirstDerivativesTester::fd_check_all(
	  NLP						*nlp
	, const VectorSlice			&xo
	, const SpVectorSlice		*xl
	, const SpVectorSlice		*xu
	, const value_type			&max_var_bounds_viol
	, const MatrixWithOp		*Gc
	, const VectorSlice			*Gf
	, bool						print_all_warnings
	, std::ostream				*out
	) const
{
	using std::setw;
	using std::endl;
	using std::right;

	using LinAlgPack::Vp_StV;
	using LinAlgPack::sqrt_eps;
	using LinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_StV;

	using SparseLinAlgPack::TestingPack::CompareDenseVectors;

	using NLPInterfacePack::CalcFiniteDiffFirstDerivatives;

	using TestingHelperPack::update_success;

	bool success = true, result;

	const size_type
		n = nlp->n(),
		m = nlp->m();
	
	// //////////////////////////////////////
	// Compute the derivatives and test them.

	if(out)
		*out
			<< "\nComputing derivatives of objective f(x) and constraints c(x) "
				"by finite differences ...\n";

	Vector FDGf;
	if(Gf)
		FDGf.resize(n);
	GenMatrix FDGc;
	if(Gc)
		FDGc.resize(n,m);

	CalcFiniteDiffFirstDerivatives fd_deriv_computer;
	fd_deriv_computer.calc_deriv(xo,xl,xu,max_var_bounds_viol,Range1D(),nlp
		, Gf ? &FDGf() : NULL
		, Gc ? &FDGc() : NULL	,BLAS_Cpp::no_trans
		,out
		);

	// /////////////////////////////////////////
	// Compare results

	if(Gf) {
		if(out)
			*out
				<< "\nComparing derivatives of objective f(x)\n"
				<< "where u(i) = finite_d(f(x))/d(x(i)), v(i) = d(f(x))/d(x(i)) ...\n";
		CompareDenseVectors comp_Gf;
		result = comp_Gf.comp( FDGf, *Gf, warning_tol(), error_tol()
			, print_all_warnings, out );
		update_success( result, &success );
	}

	if(Gc) {
		if(out)
			*out
				<< "\nComparing derivatives of constraints c(x)\n"
				<< "where D(i,j) = finite_d(c(j))/d(x(i)), M(i,j) = d(c(j))/d(x(i)) ...\n";
		result = comp_Gc().comp( FDGc, *Gc, BLAS_Cpp::no_trans
			, CompareDenseSparseMatrices::FULL_MATRIX
			, warning_tol(), error_tol(), print_all_warnings, out );
		update_success( result, &success );
	}

	return success;
}

bool NLPFirstDerivativesTester::fd_directional_check(
	  NLP						*nlp
	, const VectorSlice			&xo
	, const SpVectorSlice		*xl
	, const SpVectorSlice		*xu
	, const value_type			&max_var_bounds_viol
	, const MatrixWithOp		*Gc
	, const VectorSlice			*Gf
	, bool						print_all_warnings
	, std::ostream				*out
	) const
{
	using std::setw;
	using std::endl;
	using std::right;

	using LinAlgPack::Vp_StV;
	using LinAlgPack::sqrt_eps;
	using LinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_MtV;

	using SparseLinAlgPack::TestingPack::CompareDenseVectors;

	using NLPInterfacePack::CalcFiniteDiffFirstDerivativeProduct;

	using TestingHelperPack::update_success;

	bool success = true, result;

	const size_type
		n = nlp->n(),
		m = nlp->m();

	CalcFiniteDiffFirstDerivativeProduct
		fd_deriv_prod;
	CompareDenseVectors	
		comp_v;

	const value_type
		rand_y_l = -1.0, rand_y_u = 1.0,
		small_num = ::sqrt(std::numeric_limits<value_type>::epsilon());

	if(out)
		*out
			<< "\nComparing products Gf'*y and/or Gc'*y with finite difference values "
				" FDGf'*y and/or FDGc'*y for random y's ...\n";
	Vector y(n);
	value_type max_Gf_warning_viol = 0.0;
	int num_Gf_warning_viol = 0;
	Vector Gc_prod, FDGc_prod;
	if(Gc) {
		Gc_prod.resize(m);
		FDGc_prod.resize(m);
	}
	for( int direc_i = 1; direc_i <= num_fd_directions(); ++direc_i ) {
		random_vector( rand_y_l, rand_y_u, &y() );
		// Compute exact??? values
		value_type
			Gf_y = Gf ? dot( *Gf, y() ) : 0.0;
		if(Gc)
			V_MtV( &Gc_prod(), *Gc, BLAS_Cpp::trans, y() );
		// Compute finite difference values
		value_type
			FDGf_y;
		fd_deriv_prod.calc_deriv_product(
			xo,xl,xu,max_var_bounds_viol
			,y(),nlp
			, Gf ? &FDGf_y : NULL
			, Gc ? &FDGc_prod() : NULL
			,out
			);
		// Compare the quantities
		assert_print_nan_inf(FDGf_y, "FDGf'*y",true,out);
		const value_type
			Gf_err = ::fabs( ( Gf_y - FDGf_y )/( Gf_y + FDGf_y + small_num ) );
		if( Gf_err >= warning_tol() ) {
			max_Gf_warning_viol = std::_MAX( max_Gf_warning_viol, Gf_err );
			++num_Gf_warning_viol;
		}
		if( Gf_err >= error_tol() ) {
			if(out)
				*out
					<< "\nError, rel_err(Gf'*y,FDGf'*y) = "
					<< "rel_err(" << Gf_y << "," << FDGf_y << ") = "
					<< Gf_err << endl
					<< "exceeded Gf_error_tol = " << error_tol() << endl
					<< "Stoping the tests!\n";
			return false;
		}
		if(out)
			*out
				<< "\nChecking the computed FDGc'*y\n"
				<< "where u(i) = (FDGc'*y)(i), v(i) = (Gc'*y)(i) ...\n";
		bool result = comp_v.comp( FDGc_prod(), Gc_prod()
						, warning_tol(), error_tol()
						, print_all_warnings, out );
		update_success( result, &success );
		if(!result) return false;

	}
	if(out && num_Gf_warning_viol)
		*out
			<< "\nFor Gf, there were " << num_Gf_warning_viol << " warning tolerance "
			<< "violations out of num_fd_directions = " << num_fd_directions()
			<< " computations of FDGf'*y\n"
			<< "and the maximum violation was " << max_Gf_warning_viol
			<< " > Gf_waring_tol = " << warning_tol() << endl;

	return true;
}

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack

