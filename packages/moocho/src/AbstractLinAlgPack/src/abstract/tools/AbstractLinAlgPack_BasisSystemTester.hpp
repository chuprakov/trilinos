// ///////////////////////////////////////////////////////////
// BasisSystemTester.h
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

#ifndef BASIS_SYSTEM_TESTER_H
#define BASIS_SYSTEM_TESTER_H

#include <iosfwd>

#include "AbstractLinAlgPackTypes.h"
#include "StandardMemberCompositionMacros.h"

namespace AbstractLinAlgPack {

///
/** Testing class for \c BasisSystem interface.
 *
 * This testing class is basically a unit tester for \c BasisSystem.  The method \c test_basis_system()
 * runs many different tests to validate the interface and the objects allocated with the interface.
 * The method \c test_basis_system() should only be called after
 * <tt>basis_sys\ref BasisSystem::update_basis ".update_basis(...)"</tt> is called on the <tt>BasisSystem</tt>
 * object <tt>basis_sys</tt>.  The output basis matrix \a C and/or direct sensitivity matrix \a D are passed through
 * a series of tests using the testing classes <tt>MatrixWithOpNonsingularTester</tt> and <tt>MatrixWithOpTester</tt>
 * respectively.  The compatibility of the matrices \c Gc, \c Gh, \c C and/or \c D are also checked in a series of
 * tests.  If the method \c test_basis_system() returns true, then the client can feel fairly confident that the
 * basis matrix object is functioning properly.
 *
 * The tests performed by this testing class are designed to allow some validation for even the larges systems
 * and will produce various levels of output so as to be usefull in debugging.
 *
 * ToDo:  Finish documentation!
 */
class BasisSystemTester {
public:

	///
	enum EPrintTestLevel { PRINT_NONE=0, PRINT_BASIC=1, PRINT_MORE=2, PRINT_ALL=3 };
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EPrintTestLevel, print_tests )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_random_tests )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol )
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol )

	///	Constructor (default options)
	BasisSystemTester(
		EPrintTestLevel  print_tests      = PRINT_NONE
		,bool            dump_all         = false
		,bool            throw_exception  = true
		,size_type       num_random_tests = 1
		,value_type      warning_tol      = 1e-14
		,value_type      error_tol        = 1e-8
		);
 
	///
	/** Test a \c BasisSystem object after <tt>BasisSystem::update_basis()</tt> is called.
	 *
	 * @param  print_all_warnings
	 *              [in] Determines if warnings for all of the comparison tests are printed or not.
	 *              Warning! may cause as much as <i>O(</i><tt>bs->var_dep().size())<tt><i>)</i> output.
	 * @param  out  [in/out] If <tt>out != NULL</tt> any and all output will be sent here.  If
	 *              <tt>out == NULL</tt> then no output will be produced.
	 *
	 * @return Returns \c true if all of the tests checked out and no unexpected exceptions were
	 * thrown.
	 *
	 * The behavior of this method depends on a set of options and the input arguments.
	 * <ul>
	 * <li> <b><tt>throw_exception(bool)</tt></b>:
	 * If <tt>throw_exception()</tt> == true</tt>, then if any of the objects within
	 * this function throw exceptions, these exceptions will be be thrown clean
	 * out of this function for the caller to handle.  If <tt>throw_exception()</tt> == false</tt>,
	 * then if any object throws an exception, the exception is caught and this this function will
	 * return <tt>false</tt>.  In any case an error message will be printed
	 * to <tt>*out</tt> (if <tt>out != NULL</tt) before leaving the function (by \c return or \c throw).
	 * <li> <b><tt>dump_all(bool)</tt></b>:
	 * If <tt>dump_all() == true</tt> then all of the computed quantities will but dumped to \c out.
	 * Note that this is a useful option for initial debugging of small systems but not a good idea for
	 * larger systems as it will result in an excessive amount of output.
	 * <li> ToDo: Add rest of options!
	 * </ul>
	 */
	bool test_basis_system(
		const BasisSystem               &basis_sys
		,const MatrixWithOp             *Gc
		,const MatrixWithOp             *Gh
		,const MatrixWithOpNonsingular  *C
		,const MatrixWithOp             *D
		,bool                           print_all_warnings
		,std::ostream                   *out
		);

}; // end class BasisSystemTester

} // end namespace AbstractLinAlgPack

#endif // BASIS_SYSTEM_TESTER_H
