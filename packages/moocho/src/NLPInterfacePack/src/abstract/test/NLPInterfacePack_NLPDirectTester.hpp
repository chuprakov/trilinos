// ///////////////////////////////////////////////////////////
// NLPFirstOrderDirectTester.h
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

#ifndef NLP_FIRST_ORDER_DIRECT_TESTER_H
#define NLP_FIRST_ORDER_DIRECT_TESTER_H

#include <iosfwd>

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"
#include "StandardMemberCompositionMacros.h"
#include "StandardCompositionMacros.h"

namespace NLPInterfacePack {

///
/** Concrete class that tests the computed values of the
  * NLPrSQPTailoredApproach interface using finite differences.
  *
  * There are two options for testing the derivatives by finite differences.
  * Each option can be picked independently for the computations with
  * the objective f(x) and its gradient Gf and for the constraints
  * c(x) and its Jacobian Gc' = [ C, N ].  The tests involving the objective
  * and the constraints will be discussed separatly.
  * 
  * For testing the gradient of the objective function, two options
  * are available.  The first option (Gf_testing_method==FD_COMPUTE_ALL)
  * is to compute FDGf by brute force which requries 2*n evaluations of
  * f(x) using central differences.  Given FDGf the following comparison
  * is then make:
  * 
  *  (1)    FDGf \approx Gf
  * 
  * The other option (Gf_testing_method==FD_DIRECTIONAL) is to compute
  * random dot products of the form DFGf'*y where y is a randomly generated
  * vector.  Using central differences DFGf'*y can be computed using
  * two evaluations of f(x) per random y.  The number of random
  * y's used is determined by the option num_fd_directions.  So the
  * number of evaluations of f(x) for this option is 2*num_fd_directions.
  * 
  * The test for the quantity py = -inv(C)*c is shown below:
  * 
  *  (2)  - FDC * (-inv(C)*c) \approx c
  *               \_________/
  *                   py
  * 
  * Computing -FDC * py requires only two evaluations of
  * c(x) using central differences.  There is no other option needed
  * for test.
  * 
  * Lastly, we have to test D = -inv(C)*N.  The first option
  * (Gc_testing_method==FD_COMPUTE_ALL) is to directly compute
  * N using central differences (2*(n-m) evaluations of c(x)) as
  * FDN and then perform the comparison:
  * 
  *  (3)  - FDC * (-inv(C)*N) \approx FDN
  *               \_________/
  *                    D
  * 
  * The matrix -FDC * D can be computed using 2*(n-m) evaluations
  * with c(x) using central differences.
  * Therefore, the total number of evaluations with c(x) for
  * comparing (3) is 4*(n-m).  If n-m is not too large then this
  * is definitly the preferred method to use.
  * 
  * The other option for testing D = -inv(C)*N is to compute
  * directional derivatives using finite differences.  In this approach
  * for the random vector y, we can compute:
  * 
  *  (4)  - FDC * (-inv(C)*N) * y \approx FDN * y
  *               \_________/
  *                    D
  * 
  * Using central differences, (4) can be computed with 4 evaluations
  * of c(x).  The number of random y's used is determined by the option
  * num_fd_directions.  So the number of evaluations of c(x) for this
  * option is 4*num_fd_directions.
  * 
  * The client can pick a set of tolerances to measure if the
  * values of the above comparisons are close enough to the finite difference
  * values.  Let's define the relative error between the computed value and the
  * finite difference value to be:
  *
  *   err(i) = | (h(i) - fdh(i)) | /  (||h||inf + ||fdh||inf + sqrt(epsilon))
  *
  * The above error takes into account the relative sizes of the elements and also
  * allows one or both of the elements to be zero without ending up with 0/0
  * or something like 1e-16 not comparing with zero.
  *
  * All errors err(i) >= warning_tol are reported to *out if out != NULL.
  * The first error err(i) >= error_tol that is found is reported
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
  * The tolerances Gf_warning_tol and Gf_error_tol are applied to the tests for
  * Gf shown in (1) for instance.  The tolerances Gc_warning_tol and Gc_error_tol
  * are used for the comparisions (2), (3) and (4).
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
class NLPFirstOrderDirectTester {
public:

	///
	enum ETestingMethod { FD_COMPUTE_ALL, FD_DIRECTIONAL };

	///
	typedef const MatrixSpace<MatrixWithOpMutable>  mat_space_t;

	///
	STANDARD_COMPOSITION_MEMBERS( mat_space_t, mat_space )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ETestingMethod, Gf_testing_method )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( ETestingMethod, Gc_testing_method )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gf_warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gf_error_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gc_warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gc_error_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gh_warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, Gh_error_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_fd_directions )

	/// Constructor
	NLPFirstOrderDirectTester(
		const mat_space_ptr_t&  mat_space           = NULL
		,ETestingMethod         Gf_testing_method   = FD_DIRECTIONAL
		,ETestingMethod         Gc_testing_method   = FD_DIRECTIONAL
		,value_type             Gf_warning_tol      = 1e-6
		,value_type             Gf_error_tol        = 1e-1
		,value_type             Gc_warning_tol      = 1e-6
		,value_type             Gc_error_tol        = 1e-1
		,value_type             Gh_warning_tol      = 1e-6
		,value_type             Gh_error_tol        = 1e-1
		,size_type              num_fd_directions   = 3
		);

	///
	/** This function takes an NLP object and its computed derivatives
	 * and function values and validates
	 * the functions and the derivatives by evaluating them
	 * about the given point #xo#.
	 * 
	 * If all the checks as described in the
	 * intro checkout then this function will return true, otherwise it
	 * will return false.
	 * 
	 * If the finite difference steps are limited by relaxed variable
	 * bounds then a warning message is printed and the derivatives
	 * computed could be very inaccurate.
	 *
	 * @param  nlp     [in] NLP object used to compute and test derivatives for.
	 * @param  xo      [in] Point at which the derivatives are computed at.
	 * @param  xl      [in] If != NULL then this is the lower variable bounds.
	 * @param  xu      [in] If != NULL then this is the upper variable bounds.
	 *	                If xl != NULL then xu != NULL must also be true
	 *                 and visa-versa or a std::invalid_arguement exceptions
	 *                 will be thrown.
	 * @param  max_var_bounds_viol
	 *                 [in] If the bounds are set then this is the maximum
	 *                 violation in the bounds allowed. when computing
	 *                 points.  If xl==NULL and xu==NULL then this
	 *                 number is ignored.
	 * @param  c       [in] Value of c(x) computed at xo.
	 *                 If NULL, then none of the tests involving it will
	 *                 be performed.
	 * @param  h       [in] Value of h(x) computed at xo.
	 *                 If NULL, then none of the tests involving it will
	 *                 be performed.  Should be NULL if nlp->mI() == 0.
	 * @param  Gf      [in] Gradient of f(x) computed at xo.
	 *                 If NULL, then none of the tests involving it will
	 *                 be performed.
	 * @param  py      [in] Newton step #py = -inv(C) * c(con_decomp)
	 *                 If NULL, then none of the tests involving it will
	 *                 be performed.
	 * @param  rGf     [in] Reduced gradient of the objective function
	 *                 #rGf = Gf(var_indep) - D' * Gf(var_dep).  If NULL,
	 *                 then none of the tests involving it will be performed.
	 * @param  GcU     [in]  Auxiliary jacobian matrix #del(c(con_undecomp),x)#.
	 *                 If NULL, htne none of the tests involving it will be performed.
	 * @param  Gh      [in] Auxiliary jacobian matrix #del(h,x)#.  If NULL, then none
	 *                 of the tests involving it will be performed.
	 * @param  D       [in] Direct sensitivity matrix #D = -inv(C)*N#.  If NULL,
	 *                 none of the tests involving it will be performed.
	 * @param  V       [in] #V = F + E * D#, which is the an auxiliary sensitivity matrix.
	 *                 If NULL, then none of the tests involving it will be performed.
	 * @param  P       [in]  #P = GhI' + GhD'* D#, which is the an auxiliary sensitivity matrix.
	 *                 If NULL, then none of the tests involving it will be performed.
	 * @param  print_all_warnings
	 *                 [in] If true then all errors greater than warning_tol
	 *                 will be printed if out!=NULL
	 * @param  out     [in/out] If != null then some summary information is printed to it
	 *                 and if a derivative does not match up then it prints which
	 *                 derivative failed.  If #out == 0# then no output is printed.
	 *
	 * @return Returns #true# if all the derivatives comparisons are
	 * within the error tolerances or returns false
	 *	otherwise.  This function will return false if any NaN or Inf values
	 *	where encountered.
	 */
	bool finite_diff_check(
		NLPFirstOrderDirect            *nlp
		,const VectorWithOp            &xo
		,const VectorWithOp            *xl
		,const VectorWithOp            *xu
		,const value_type              &max_var_bounds_viol
		,VectorWithOpMutable           *c
		,VectorWithOpMutable           *h
		,VectorWithOpMutable           *Gf
		,VectorWithOpMutable           *py
		,VectorWithOpMutable           *rGf
		,MatrixWithOp                  *GcU
		,MatrixWithOp                  *Gh
		,MatrixWithOp                  *D
		,MatrixWithOp                  *V
		,MatrixWithOp                  *P
		,bool                          print_all_warnings
		,std::ostream                  *out
		) const;

};	// end class NLPFirstOrderDirectTester

}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_ORDER_DIRECT_TESTER_H
