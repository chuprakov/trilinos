// ///////////////////////////////////////////////////////////////////////////////////////
// VariableBoundsTester.cpp

#include "ConstrainedOptimizationPack/include/VariableBoundsTester.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "SparseLinAlgPack/test/CompareDenseVectors.h"
#include "SparseLinAlgPack/include/sparse_bounds_diff.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/LinAlgOpPack.h"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Vp_StV;
	using SparseLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptimizationPack {

// public

VariableBoundsTester::VariableBoundsTester(
		  value_type	warning_tol
		, value_type	error_tol
		)
	:
		warning_tol_(warning_tol)
		,error_tol_(error_tol)
{}

bool VariableBoundsTester::check_in_bounds(
	  std::ostream* out, bool print_all_warnings, bool print_vectors
	, const SpVectorSlice& xL, const char xL_name[]
	, const SpVectorSlice& xU, const char xU_name[]
	, const VectorSlice& x, const char x_name[]
	)
{
	using std::endl;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::upper;
	using BLAS_Cpp::lower;
	using LinAlgPack::norm_inf;
	using LinAlgPack::Vt_S;
	using SparseLinAlgPack::imp_sparse_bnd_diff;

	value_type scale = 0.0;
	Vector
		u;	// hold the result to pass to comparison function
	const value_type
		x_norm_inf = norm_inf(x);

	SparseLinAlgPack::TestingPack::CompareDenseVectors comp_v;

	if(out)
		*out
			<< "\n*** Checking that variables are in bounds\n";

	///////////////////////////////////
	// xL - x <= 0
	if(out)
		*out
			<< "\nChecking "<<xL_name<<" - "<<x_name<<" <= 0 ...\n";
	u.resize(x.size());
	imp_sparse_bnd_diff( +1, xL, lower, x, &u() );
	if(out && print_vectors)
		*out
			<<xL_name<<" - "<<x_name<<" =\n" << u();
	Vt_S( &u(), 1.0/(1.0+x_norm_inf) );
	if(out) {
		*out
			<< "Comparing: u - v <= 0\n"
			<< "u = ("<<xL_name<<" - "<<x_name<<") | / (1 + ||"<<x_name<<"||inf ), v = 0 ...\n";
	}
	if(!comp_v.comp_less( u(), 0.0, warning_tol(), error_tol()
		, print_all_warnings, out )) return false;

	///////////////////////////////////
	// x - xU <= 0
	if(out)
		*out
			<< "\nChecking "<<x_name<<" - "<<xU_name<<" <= 0 ...\n";
	u.resize(x.size());
	imp_sparse_bnd_diff( -1, xU, upper, x, &u() );
	if(out && print_vectors)
		*out
			<<x_name<<" - "<<xU_name<<" =\n" << u();
	Vt_S( &u(), 1.0/(1.0+x_norm_inf) );
	if(out) {
		*out
			<< "Comparing: u - v <= 0\n"
			<< "u = ("<<x_name<<" - "<<xU_name<<") | / (1 + ||"<<x_name<<"||inf ), v = 0 ...\n";
	}
	if(!comp_v.comp_less( u(), 0.0, warning_tol(), error_tol()
		, print_all_warnings, out )) return false;

	if(out)
		*out
			<< "\n*** Congradulations, the variables are within the bounds by the specified"
			<< "\n*** error tolerances!\n";

	return true;	// If we get here then success

}

}	// end namespace ConstrainedOptimizationPack
