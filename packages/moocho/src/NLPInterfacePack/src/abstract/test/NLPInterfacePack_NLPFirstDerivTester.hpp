// ///////////////////////////////////////////////////////////
// NLPFirstDerivativesTester.h

#ifndef NLP_FIRST_DERIVATIVES_TESTER_H
#define NLP_FIRST_DERIVATIVES_TESTER_H

#include <iosfwd>

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"
#include "SparseLinAlgPack/test/CompareDenseSparseMatrices.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace NLPInterfacePack {
namespace TestingPack {

///
/** Concrete class that tests the derivatives using finite differences.
  *
  * This class computes the derivatives using two sided finite differences.
  * These test can be very expensive in storage space and run time.
  *
  * The amount of storage space needed is approximatly n*m and f(x) and
  * c(x) will be computed n times more.
  *
  * The client can set a set of tolerances to measure if the anylitical
  * values of Gf and Gc are close enough to the finite difference
  * values.  Let the function h(x) be f(x) or any cj(x), j = 1...m.  Let
  * gh(i) = d(h(x))/d(x(i)) and fdh(i) = finite_diff(h(x))/d(x(i)).  Then
  * let's define the relative error between the anylitic value and the
  * finite difference value to be:
  *
  *   err(i) = | (gh(i) - fdh(i)) /  (gh(i) + fdh(i) + sqrt(epsilon)) |
  *
  * The above error takes into account the relative sizes of the elements and also
  * allows one or both of the elements to be zero without ending up with 0/0
  * or something like 1e-16 not comparing with zero.
  *
  * All errors err(i) >= warning_tol are reported to *out if out != NULL.
  * The first error err(i) >= error_tol that is found is reported is reported
  * to *out if out != NULL and immediatly finite_diff_check(...) returns false.
  * If all errors err(i) < error_tol then finite_diff_check(...) will return true.
  *
  * Given these two tolerances the client can do many things:
  *
  * 1) Print out all the comparisons that are not equal by setting warning_tol
  *    <= epsilon and error_tol = very_large_number.
  *
  * 2) Print out all suspect comparisons by setting epsilon < warning_tol < 1
  *    and error_tol = very_large_number.
  *
  * 3) Just validate that matrices are approximatly equal and report the first
  *    discrepency if not by setting epsilon < error_tol < 1 and warning_tol
  *    >= error_tol.
  *
  * 4) Print out any suspect comparisons by setting epsilon < warning_tol < 1
  *    but also quit if the error is too large by setting error_tol > warning_tol.
  *
  */
class NLPFirstDerivativesTester {
public:

	///
	typedef SparseLinAlgPack::TestingPack::CompareDenseSparseMatrices
		CompareDenseSparseMatrices;

	/// «std comp» members for comparision object compatible with Gc
	STANDARD_COMPOSITION_MEMBERS( CompareDenseSparseMatrices, comp_Gc )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol )

	/// Constructor
	NLPFirstDerivativesTester(
			  value_type			warning_tol = 1e-6
			, value_type			error_tol	= 1e-3
			, const comp_Gc_ptr_t&	comp_Gc		= new CompareDenseSparseMatrices
		);

	///
	/** This function takes an NLP object and its computed derivatives
	  * and function values and validates
	  * the functions and the derivatives by evaluating them
	  * about the given point #x#.  If all the checks as described in the
	  * intro checkout then this function will return true, otherwise it
	  * will return false.
	  *
	  * @param	nlp		NLP object used to compute and test derivatives for.
	  *	@param	Gc		A matrix object for the Gc computed at x
	  *	@param	Gf		Gradient of f(x) computed at x
	  *	@param	x		Point at which the derivatives are computed at.
	  *	@param	out		If != null then some summary information is printed to it
	  *					and if a derivative does not match up then it prints which
	  *					derivative failed.  If #out == 0# then no output is printed.
	  *
	  * @return Returns #true# if all the derivatives check out, and false
	  *	otherwise.
	  */
	bool finite_diff_check(
		  NLPFirstOrderInfo*					nlp
		, const MatrixWithOp&					Gc
		, const VectorSlice&					Gf
		, const VectorSlice&					x
		, std::ostream*							out		) const;
};

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_DERIVATIVES_TESTER_H