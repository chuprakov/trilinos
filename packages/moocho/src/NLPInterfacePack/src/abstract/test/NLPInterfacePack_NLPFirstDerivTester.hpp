// ///////////////////////////////////////////////////////////
// NLPFirstDerivativesTester.h
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

#ifndef NLP_FIRST_DERIVATIVES_TESTER_H
#define NLP_FIRST_DERIVATIVES_TESTER_H

#include <iosfwd>

#include "NLPInterfacePack/include/NLPInterfacePackTypes.h"
//#include "SparseLinAlgPack/test/CompareDenseSparseMatrices.h"
#include "StandardCompositionMacros.h"
#include "StandardMemberCompositionMacros.h"

namespace NLPInterfacePack {

///
/** Concrete class that tests the derivatives using finite differences.
 *
 * There are two options for testing the derivatives by finite differences.
 * 
 * The first option (<tt>fd_testing_method==FD_COMPUTE_ALL</tt>) is to compute all of
 * them as dense vectors and matrices. This option can be very expensive in runtime
 * and storage costs.  The amount of storage space needed is <tt>O(n*m)</tt> and
 * \c f(x) and \c c(x) will be computed <tt>O(n)</tt> times.
 * 
 * The other option (<tt>fd_testing_method==FD_DIRECTIONAL</tt>) computes products
 * of the form <tt>g'*v</tt> and compares them to the finite difference computed
 * value <tt>g_fd'*v</tt>.  This method only costs <tt>O(n)</tt> storage and two function
 * evaluations per direction (assuming central differences are used.  The directions \c v
 * are computed randomly between [1,10] so that they are well scaled and should give good results.
 * The option <tt>num_fd_directions()</tt> determines how many random directions are used.
 *
 * This class computes the derivatives using two sided (central )finite differences
 * which is more accurate than one sided differences.
 *
 * The client can set the tolerances used to measure if the anylitical
 * values of \c Gf and \c Gc are close enough to the finite difference
 * values.  Let the function \a h(x) be \c f(x) or any \c c<sub>j</sub>(x),
 * <tt>j = 1...m</tt>.  Let <tt>gh(i) = d(h(x))/d(x(i))</tt> and
 * <tt>fdh(i) = finite_diff(h(x))/d(x(i))</tt>.  Then let's define the
 * relative error between the anylitic value and the finite difference value to be:
 \verbatim

     err(i) = |(gh(i) - fdh(i))| /  (||gh||inf + ||fdh||inf + (epsilon)^(1/4))
 \endverbatim
 * The above error takes into account the relative sizes of the elements and also
 * allows one or both of the elements to be zero without ending up with <tt>0/0</tt>
 * or something like <tt>1e-16</tt> not comparing with zero.
 *
 * All errors <tt>err(i) >= warning_tol</tt> are reported to <tt>*out</tt> if
 * <tt>out != NULL</tt> and <tt>print_all_warnings==true</tt>.  Otherwise, if
 * <tt>out != NULL</tt>, only the number of elements and the maxinum violation of the
 * warning tolerance will be printed.  The first error <tt>err(i) >= error_tol</tt>
 * that is found is reported is reported to <tt>*out</tt> if <tt>out != NULL</tt> and
 * immediatly \c finite_diff_check() returns \c false.  If all errors
 * <tt>err(i) < error_tol</tt> then \c finite_diff_check() will return \c true.
 *
 * Given these two tolerances the client can do many things:
 * <ol>
 * <li> Print out all the comparisons that are not equal by setting warning_tol
 *    == 0.0 and error_tol = very_large_number.
 *
 * <li> Print out all suspect comparisons by setting epsilon < warning_tol < 1
 *    and error_tol = very_large_number.
 *
 * <li> Just validate that matrices are approximatly equal and report the first
 *    discrepency if not by setting epsilon < error_tol < 1 and warning_tol
 *    >= error_tol.
 *
 * <li> Print out any suspect comparisons by setting epsilon < warning_tol < 1
 *    but also quit if the error is too large by setting error_tol > warning_tol.
 * </ol>
 * There is one minor hitch to this testing.  For many NLPs, there is a
 * strict region of \a x where \a f(x) or \a c(x) are not defined.  In order to
 * help ensure that we stay out of these regions, variable bounds can be
 * included and a scalar \c max_var_bounds_viol so that the testing software
 * will never evaluate \a f(x) or \a c(x) outside the region:
 \verbatim
   
     xl - max_var_bounds_viol <= x <= xu + max_var_bounds_viol
 \endverbatim
 * This is an important agreement made with the user.
 */
class NLPFirstDerivativesTester {
public:

	///
	enum ETestingMethod { FD_COMPUTE_ALL, FD_DIRECTIONAL };

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
		);

	///
	/** This function takes an NLP object and its computed derivatives
	 * and function values and validates
	 * the functions and the derivatives by evaluating them
	 * about the given point <tt>x</tt>.  If all the checks as described in the
	 * intro checkout then this function will return true, otherwise it
	 * will return false.
	 *
	 * @param  nlp      [in] NLP object used to compute and test derivatives for.
	 * @param  xo       [in] Point at which the derivatives are computed at.
	 * @param  xl       [in] If != NULL then this is the lower variable bounds.
	 * @param  xu       [in] If != NULL then this is the upper variable bounds.
	 *                  If xl != NULL then xu != NULL must also be true
	 *                  and visa-versa or a std::invalid_arguement exceptions
	 *                  will be thrown.
	 * @param  max_var_bounds_viol
	 *                  [in] If the bounds are set then this is the maximum
	 *                  violation in the bounds allowed. when computing
	 *                  points.  If xl==NULL and xu==NULL then this
	 *                  number is not important.
	 * @param  Gc       [in] A matrix object for the Gc computed at xo.
	 *                  If Gc==NULL then this is not tested for.
	 * @param  Gh       [in] A matrix object for the Gh computed at xo.
	 *                  If Gh==NULL then this is not tested for.
	 * @param  Gf       [in] Gradient of f(x) computed at xo.
	 *                  If Gf==NULL then this is not tested for.
	 * @param  print_all_warnings
	 *                  [in] If true then all errors greater than warning_tol
	 *                  will be printed if out!=NULL
	 * @param  out      [in/out] If != null then some summary information is printed to it
	 *                  and if a derivative does not match up then it prints which
	 *                  derivative failed.  If <tt>out == 0</tt> then no output is printed.
	 *
	 * @return Returns <tt>true</tt> if all the derivatives check out, and false
	 * otherwise.
	 */
	bool finite_diff_check(
		NLP                     *nlp
		,const VectorWithOp     &xo
		,const VectorWithOp     *xl
		,const VectorWithOp     *xu
		,const value_type       &max_var_bounds_viol
		,const MatrixWithOp     *Gc
		,const MatrixWithOp     *Gh
		,const VectorWithOp     *Gf
		,bool                   print_all_warnings
		,std::ostream           *out
		) const;

private:

	///
	bool fd_check_all(
		NLP                     *nlp
		,const VectorWithOp     &xo
		,const VectorWithOp     *xl
		,const VectorWithOp     *xu
		,const value_type       &max_var_bounds_viol
		,const MatrixWithOp     *Gc
		,const MatrixWithOp     *Gh
		,const VectorWithOp     *Gf
		,bool                   print_all_warnings
		,std::ostream           *out
		) const;

	///
	bool fd_directional_check(
		NLP                     *nlp
		,const VectorWithOp     &xo
		,const VectorWithOp     *xl
		,const VectorWithOp     *xu
		,const value_type       &max_var_bounds_viol
		,const MatrixWithOp     *Gc
		,const MatrixWithOp     *Gh
		,const VectorWithOp     *Gf
		,bool                   print_all_warnings
		,std::ostream           *out
		) const;
};

}	// end namespace NLPInterfacePack

#endif	// NLP_FIRST_DERIVATIVES_TESTER_H
