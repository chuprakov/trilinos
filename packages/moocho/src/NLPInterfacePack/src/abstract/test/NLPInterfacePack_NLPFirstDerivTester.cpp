// ///////////////////////////////////////////////////////////
// NLPFirstDerivativesTester.cpp

#include <assert.h>
#include <iomanip>
#include <sstream>

#include "NLPFirstDerivativesTester.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "NLPInterfacePack/include/CalcFiniteDiffFirstDerivatives.h"
#include "SparseLinAlgPack/test/CompareDenseVectors.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/MatVecCompare.h"
#include "LinAlgPack/include/assert_print_nan_inf.h"

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
	  NLPFirstOrderInfo*					nlp
	, const MatrixWithOp&					Gc
	, const VectorSlice&					Gf
	, const VectorSlice&					xo
	, bool									print_all_warnings
	, std::ostream*							out
	) const
{
	switch(fd_testing_method_) {
		case FD_COMPUTE_ALL:
			return fd_check_all(nlp,Gc,Gf,xo,print_all_warnings,out);
		case FD_DIRECTIONAL:
			return fd_directional_check(nlp,Gc,Gf,xo,print_all_warnings,out);
		default:
			assert(0);
	}
	return false;	// will never be executed
}

// private

bool NLPFirstDerivativesTester::fd_check_all(
	  NLPFirstOrderInfo*					nlp
	, const MatrixWithOp&					Gc
	, const VectorSlice&					Gf
	, const VectorSlice&					xo
	, bool									print_all_warnings
	, std::ostream*							out
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
	
	Vector FDGf(0.0,n);
	GenMatrix FDGc(0.0,n,m);

	// Check the input vectors
	assert_print_nan_inf(xo, "xo",true,out); 
	assert_print_nan_inf(Gf, "Gf",true,out); 

	// Compute derivatives by finite differences
	if(out) {
		*out
			<< "\nComputing derivatives of objective f(x) and constraints c(x) "
				"by finite differences ...\n";
	}
	CalcFiniteDiffFirstDerivatives fd_deriv_computer;
	fd_deriv_computer.calc_deriv(xo,nlp,&FDGf(),&FDGc(),out);

	// Compare results
	// 
	if(out) {
		*out
			<< "\nComparing derivatives of objective f(x)\n"
			<< "where u(i) = finite_d(f(x))/d(x(i)), v(i) = d(f(x))/d(x(i)) ...\n";
	}
	CompareDenseVectors comp_Gf;
	result = comp_Gf.comp( FDGf, Gf, warning_tol(), error_tol()
		, print_all_warnings, out );
	update_success( result, &success );

	if(out) {
		*out
			<< "\nComparing derivatives of constraints c(x)\n"
			<< "where D(i,j) = finite_d(c(j))/d(x(i)), M(i,j) = d(c(j))/d(x(i)) ...\n";
	}
	result = comp_Gc().comp( FDGc, Gc, BLAS_Cpp::no_trans
		, CompareDenseSparseMatrices::FULL_MATRIX
		, warning_tol(), error_tol(), print_all_warnings, out );
	update_success( result, &success );

	return success;
}

bool NLPFirstDerivativesTester::fd_directional_check(
	  NLPFirstOrderInfo*					nlp
	, const MatrixWithOp&					Gc
	, const VectorSlice&					Gf
	, const VectorSlice&					xo
	, bool									print_all_warnings
	, std::ostream*							out
	) const
{
	// ToDo: Finish this!
	assert(0);
}

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack

