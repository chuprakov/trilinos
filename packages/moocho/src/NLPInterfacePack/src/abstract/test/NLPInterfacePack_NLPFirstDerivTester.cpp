// ///////////////////////////////////////////////////////////
// NLPFirstDerivativesTester.cpp

#include <assert.h>
#include <iomanip>
#include <sstream>

#include "NLPFirstDerivativesTester.h"
#include "../include/NLPFirstOrderInfo.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/MatVecCompare.h"
#include "LinAlgPack/include/assert_print_nan_inf.h"

namespace NLPInterfacePack {
namespace TestingPack {

NLPFirstDerivativesTester::NLPFirstDerivativesTester(
		  value_type			warning_tol
		, value_type			error_tol
		, const comp_Gc_ptr_t&	comp_Gc		)
	:
		  warning_tol_(warning_tol)
		, error_tol_(error_tol)
		, comp_Gc_(comp_Gc)
{}	

bool NLPFirstDerivativesTester::finite_diff_check(
	  NLPFirstOrderInfo*					nlp
	, const MatrixWithOp&					Gc
	, const VectorSlice&					Gf
	, const VectorSlice&					xo
	, std::ostream*							out		) const
{
	using std::setw;
	using std::endl;
	using std::right;

	using LinAlgPack::Vp_StV;
	using LinAlgPack::sqrt_eps;
	using LinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_StV;

	using TestingHelperPack::update_success;

	bool success = true, result;

	// Gc = [ Gc1, Gc2, ..., Gcm ]
	//
	//		[	dc1/dx1		dc2/dx1	...		dcm/dx1	]
	//		[	dc1/dx2		dc2/dx2	...		dcm/dx2	]
	// Gc= 	[	.			.		...		.		]
	//		[	dc1/dxn		dc2/dxn	...		dcm/dxn	]
	//
	//		[	dc/dx1	]
	//		[	dc/dx2	]
	// Gc =	[	.		]
	//		[	dc/dxn	]
	//
	// Here we will compute dc/dxi using a divided difference formula:
	//
	//		dc/dxi	= ( c(xo + ei*h) - c(xo - ei*h) ) / 2*h
	//		

	const value_type small_num = ::sqrt(std::numeric_limits<value_type>::epsilon());

	value_type h = sqrt_eps;

	size_type	n = nlp->n(),
				m = nlp->m();
	
	GenMatrix FDGc(0.0,n,m);
	value_type f;
	Vector x = xo, c(m), FDGf(0.0,n);

	// Check the input vectors
	assert_print_nan_inf(xo, "xo",true,out); 
	assert_print_nan_inf(xo, "Gf",true,out); 

	// Remember what was set

	value_type		*f_saved;
	Vector			*c_saved;

	f_saved = nlp->get_f();
	c_saved = nlp->get_c();

	// Set the quanitities used to compute with

	nlp->set_f( &f );
	nlp->set_c( &c );

	// Output header
	bool printed_header = false;
	const int int_w = 10, dbl_w = 25, dbl_p = 15;
	const char int_ul[] = "--------", dbl_ul[] = "-----------------------"; 

	int p_saved;
	if(out) {
		p_saved = out->precision();
		*out << std::setprecision(dbl_p);
	}

	// For each x(i) compute the finite difference approximations to
	// Gf and Gc
	if(out)
		*out << "Checking derivatives of objective function f(x)...\n";
	for( size_type i = 1; i <= n; ++i ) {

		x(i) = x(i) + h;

		nlp->calc_c( x );
		assert_print_nan_inf(c(), "c(xo+h*e(i))",true,out); 
		V_StV( &FDGc.row(i), 1/(2*h), c() );

		nlp->calc_f( x, false );
		FDGf(i) = f;

		x(i) = x(i) - 2 * h;

		nlp->calc_c( x );
		assert_print_nan_inf(c(), "c(xo-h*e(i))",true,out); 
		Vp_StV( &FDGc.row(i), -1/(2*h), c() );

		nlp->calc_f( x, false );
		FDGf(i) = (FDGf(i) - f) / (2 * h);

		x(i) = x(i) + h;

		const value_type
			ele_err = ::fabs( (FDGf(i) - Gf(i)) / (FDGf(i) + Gf(i) + small_num ) );

		if( ele_err >= warning_tol() || ele_err >= error_tol() ) {
			if(out) {
				if(!printed_header) {
					*out
						<< "Warning, the following relative error between "
						<< "the elements is greater than warning_tol = "
						<< warning_tol() << endl << endl
						<< right << setw(int_w) << "i"
						<< right << setw(dbl_w) << "d(f(x))/d(x(i))"
						<< right << setw(dbl_w) << "finite_d(f(x))/d(x(i))"
						<< right << setw(dbl_w) << "relative error\n"
						<< right << setw(int_w) << int_ul
						<< right << setw(dbl_w) << dbl_ul
						<< right << setw(dbl_w) << dbl_ul
						<< right << setw(dbl_w) << dbl_ul << endl;
					printed_header = true;
				}
				*out
					<< right << setw(int_w) << i
					<< right << setw(dbl_w) << Gf(i)
					<< right << setw(dbl_w) << FDGf(i)
					<< right << setw(dbl_w) << ele_err << endl;
			}
			if( ele_err >= error_tol() ) {
				if(out)
					*out
						<< "*** Error, the last error is larger than error_tol = "
						<< error_tol() << endl; 
				return false;	// no more testing!
			}
		}
	}
	if(out) {
		if(printed_header) *out << std::endl;
		*out
			<< "Checking derivatives of constraints c(x)\n"
			<< "where D(i,j) = finite_d(c(j))/d(x(i)), M(i,j) = d(c(j))/d(x(i)) ...\n";
	}
	result = comp_Gc().comp( FDGc, Gc, CompareDenseSparseMatrices::FULL_MATRIX
		, warning_tol(), error_tol(), out );
	update_success( result, &success );

	return success;
}

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack

