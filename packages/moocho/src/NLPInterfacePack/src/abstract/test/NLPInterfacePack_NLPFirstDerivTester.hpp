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
  * There are two options from testing the derivatives by finite differences.
  * 
  * The first option (fd_testing_method==FD_COMPUTE_ALL) is to compute all of
  * them as dense vectors and matrices. This option can be very expensive
  * in runtime and storage costs.
  * The amount of storage space needed is O(n*m) and f(x) and
  * c(x) will be computed O(n) times.
  * 
  * The other option (fd_testing_method==FD_DIRECTIONAL) computes products
  * of the form g'*v and compares them to the finite difference computed
  * value g_fd'*v.  This method only costs O(n) storage and two function
  * evaluations per direction.  The directions v are computed randomly between
  * [1,10] so that they are well scaled and should give good results.
  * The option num_fd_directions determines how many directions are
  * used.
  *
  * This class computes the derivatives using two sided finite differences
  * which is more accurate than one sided differences.
  *
  * The client can set a set of tolerances to measure if the anylitical
  * values of Gf and Gc are close enough to the finite difference
  * values.  Let the function h(x) be f(x) or any cj(x), j = 1...m.  Let
  * gh(i) = d(h(x))/d(x(i)) and fdh(i) = finite_diff(h(x))/d(x(i)).  Then
  * let's define the relative error between the anylitic value and the
  * finite difference value to be:
  *
  *   err(i) = |(gh(i) - fdh(i))| /  (||gh||inf + ||fdh||inf + (epsilon)^(1/4))
  *
  * The above error takes into account the relative sizes of the elements and also
  * allows one or both of the elements to be zero without ending up with 0/0
  * or something like 1e-16 not comparing with zero.
  *
  * All errors err(i) >= warning_tol are reported to *out if out != NULL
  * and print_all_warnings==true.  Otherwise, if out != NULL,
  * only the number of elements and the maxinum violation of the
  * warning tolerance will be printed.
  * The first error err(i) >= error_tol that is found is reported is reported
  * to *out if out != NULL and immediatly finite_diff_check(...) returns false.
  * If all errors err(i) < error_tol then finite_diff_check(...) will return true.
  *
  * Given these two tolerances the client can do many things:
  *
  * 1) Print out all the comparisons that are not equal by setting warning_tol
  *    == 0.0 and error_tol = very_large_number.
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
  * There is one minor hitch to this testing.  For many NLPs, there is a
  * strict region of x where f(x) or c(x) are not defined.  In order to
  * help ensure that we stay out of these regions, variable bounds can be
  * included and a scalar max_var_bounds_viol so that the testing software
  * will never evaluation f(x) or c(x) outside the region:
  * 
  * xl - max_var_bounds_viol <= x <= xu + max_var_bounds_viol
  * 
  * This is an important agreement made with the user.
  */
class NLPFirstDerivativesTester {
public:

	///
	enum ETestingMethod { FD_COMPUTE_ALL, FD_DIRECTIONAL };

	///
	typedef SparseLinAlgPack::TestingPack::CompareDenseSparseMatrices
		CompareDenseSparseMatrices;

	/// «std comp» members for comparision object compatible with Gc
	STANDARD_COMPOSITION_MEMBERS( CompareDenseSparseMatrices, comp_Gc )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ETestingMethod, fd_testing_method )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_fd_directions )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol )

	/// Constructor
	NLPFirstDerivativesTester(
			  ETestingMethod		fd_testing_method	= FD_DIRECTIONAL
			, size_type				num_fd_directions	= 3
			, value_type			warning_tol 		= 1e-6
			, value_type			error_tol			= 1e-1
			, const comp_Gc_ptr_t&	comp_Gc				= new CompareDenseSparseMatrices
		);

	///
	/** This function takes an NLP object and its computed derivatives
	  * and function values and validates
	  * the functions and the derivatives by evaluating them
	  * about the given point #x#.  If all the checks as described in the
	  * intro checkout then this function will return true, otherwise it
	  * will return false.
	  *
	  * @param	nlp		[I]	NLP object used to compute and test derivatives for.
	  *	@param	xo		[I]	Point at which the derivatives are computed at.
	  *	@param	xl		[I] If != NULL then this is the lower variable bounds.
	  *	@param	xu		[I] If != NULL then this is the upper variable bounds.
	  *						If xl != NULL then xu != NULL must also be true
	  *						and visa-versa or a std::invalid_arguement exceptions
	  *						will be thrown.
	  *	@param	max_var_bounds_viol
	  *					[I] If the bounds are set then this is the maximum
	  *						violation in the bounds allowed. when computing
	  *						points.  If xl==NULL and xu==NULL then this
	  *						number is not important.
	  *	@param	Gc		[I] A matrix object for the Gc computed at xo.
	  *						If Gc==NULL then this is not tested for.
	  *	@param	Gf		[I] Gradient of f(x) computed at xo.
	  *						If Gf==NULL then this is not tested for.
	  *	@param	print_all_warnings
	  *					[I] If true then all errors greater than warning_tol
	  *						will be printed if out!=NULL
	  *	@param	out		If != null then some summary information is printed to it
	  *					and if a derivative does not match up then it prints which
	  *					derivative failed.  If #out == 0# then no output is printed.
	  *
	  * @return Returns #true# if all the derivatives check out, and false
	  *	otherwise.
	  */
	bool finite_diff_check(
		  NLP						*nlp
		, const VectorSlice			&xo
		, const SpVectorSlice		*xl
		, const SpVectorSlice		*xu
		, const value_type			&max_var_bounds_viol
		, const MatrixWithOp		*Gc
		, const VectorSlice			*Gf
		, bool						print_all_warnings
		, std::ostream				*out
		) const;

private:

	///
	bool fd_check_all(
		  NLP						*nlp
		, const VectorSlice			&xo
		, const SpVectorSlice		*xl
		, const SpVectorSlice		*xu
		, const value_type			&max_var_bounds_viol
		, const MatrixWithOp		*Gc
		, const VectorSlice			*Gf
		, bool						print_all_warnings
		, std::ostream				*out
		) const;

	///
	bool fd_directional_check(
		  NLP						*nlp
		, const VectorSlice			&xo
		, const SpVectorSlice		*xl
		, const SpVectorSlice		*xu
		, const value_type			&max_var_bounds_viol
		, const MatrixWithOp		*Gc
		, const VectorSlice			*Gf
		, bool						print_all_warnings
		, std::ostream				*out
		) const;
};

}	// end namesapce TestingPack
}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_DERIVATIVES_TESTER_H
