// ///////////////////////////////////////////////////////////
// NLPFirstDerivativesTester.cpp

#include "NLPFirstDerivativesTester.h"
#include "../include/NLPFirstOrderInfo.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/MatVecCompare.h"

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
	using LinAlgPack::Vp_StV;
	using LinAlgPack::sqrt_eps;

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

	// Remember what was set

	value_type		*f_saved;
	Vector			*c_saved;

	f_saved = nlp->get_f();
	c_saved = nlp->get_c();

	// Set the quanitities used to compute with

	nlp->set_f( &f );
	nlp->set_c( &c );

	// For each x(i) compute the finite difference approximations to
	// Gf and Gc
	if(out)
		*out << "Checking derivatives of objective function f(x)...\n";
	for( size_type i = 1; i <= n; ++i ) {

		x(i) = x(i) + h;

		nlp->calc_c( x );
		V_StV( &FDGc.row(i), 1/(2*h), c() );

		nlp->calc_f( x, false );
		FDGf(i) = f;

		x(i) = x(i) - 2 * h;

		nlp->calc_c( x );
		Vp_StV( &FDGc.row(i), -1/(2*h), c() );

		nlp->calc_f( x, false );
		FDGf(i) = (FDGf(i) - f) / (2 * h);

		x(i) = x(i) + h;

		const value_type
			err = ::fabs( (FDGf(i) - Gf(i)) / (FDGf(i) + Gf(i) + small_num ) );

		if( err >= warning_tol() || err >= error_tol() ) {
			if(out) {
				std::ostringstream omsg;
				if( err >= error_tol() ) {
					omsg
						<< "*** Error, ";
				}
				else {	// ele_err >= warning_tol
					omsg
						<< "*** Warning, ";
				}
				omsg	<< "d(f(x))/d(x(" << i << ")) = " << Gf(i)
						<< " not close enough to finite_diff(f(x))/d(x(" << i << ")) = "
						<< FDGf(i) << "\n";
				*out << omsg.str();
			}
			if( err >= error_tol() )
				return false;	// no more testing!
		}
	}

	if(out)
		*out << "Checking derivatives of constraints c(x)...\n";
	result = comp_Gc().comp( FDGc, Gc, CompareDenseSparseMatrices::FULL_MATRIX
		, warning_tol(), error_tol(), out );
	if( !result && out )
		*out << "where:\n"
				"    M(i,j) = d(c(j))/d(x(i)), D(i,j) = finite_d(c(j))/d(x(i)) ...\n";
	update_success( result, &success );

	return success;
}

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack

