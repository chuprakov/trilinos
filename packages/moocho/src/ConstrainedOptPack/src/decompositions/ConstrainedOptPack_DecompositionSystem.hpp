// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystem.h
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

#ifndef DECOMPOSITION_SYSTEM_H
#define DECOMPOSITION_SYSTEM_H

#include <stdexcept>

#include "ConstrainedOptimizationPackTypes.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"

namespace ConstrainedOptimizationPack {

///
/** This class abstracts a decomposition choice for the range space \a Y,
 * and null space \a Z, matrices for a linearly independent set of columns of \a Gc.
 *
 * <tt>Gc = [ Gc(:,con_decomp),  Gc(:,con_undecomp) ]</tt>
 *
 * where \c Gc is <tt>n x m</tt>, \c Gc(:,con_decomp) is <tt>n x r</tt> and
 * \c Gc(:,con_undecomp) is <tt>n x (m - r)</tt>.
 *
 * Note that the columns in <tt>Gc(:,con_undecomp)</tt> may be linearly dependent with
 * the columns in <tt>Gc(:,con_undecomp)</tt> but they may just be undecomposed
 * linearly independent equality constraints.
 *
 * The decomposition formed by subclasses must have the properties:
 \verbatim
	 Z s.t. Gc(:,con_deomp)' * Z = 0
	 Y s.t. [Z  Y] is nonsingular
	 R = Gc(:,con_decomp)' * Y is nonsingular
	 Uz = Gc(:,con_undecomp)' * Z
	 Uy = Gc(:,con_undecomp)' * Y
	 Vz = Gh' * Z
	 Vy = Gh' * Y
 \endverbatim
 *
 * The matrix factory objects returned by ??? are ment to have a lifetime that is
 * independent of \c this.
 *
 * The decomposition matrices \c Z, \c Y, \c R, \c Uz, \c Uy, \c Vz and \c Vy which
 * are updated in <tt>this->update_decomp()</tt> must be completely independent from
 * \c this and from each other and \c Gc and \c Gh that they based on.  For example,
 * Once \c update_decomp() is called, \c this, \c Gc and \c Gh can be destroyed and
 * the behaviors of the decomposition matrices must not be altered.  In this respect
 * the <tt>%DecompositionSystem</tt> interface is really nothing more than a "Strategy"
 * interface (with some state data of course) for computing range/null decompositions.
 * This gives the client great flexibility in how the decomposition matrices are used. 
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystem {
public:

	/** @name Public types */
	//@{

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractFactoryPack::AbstractFactory<MatrixWithOpNonsingular> >    mat_nonsing_fcty_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractFactoryPack::AbstractFactory<MatrixWithOp> >               mat_fcty_ptr_t;
	///
	class InvalidMatrixType : public std::logic_error
	{public: InvalidMatrixType(const std::string& what_arg) : std::logic_error(what_arg) {}};
	///
	class TestFailed : public std::runtime_error
	{public: TestFailed(const std::string& what_arg) : std::runtime_error(what_arg) {}};
	/// Enumeration for the amount of output to create from <tt>update_decomp()</tt>.
	enum EOutputLevel {
		PRINT_NONE          = 0,
		PRINT_BASIC_INFO    = 1,
		PRINT_MORE_INFO     = 2,
		PRINT_VECTORS       = 3,
		PRINT_EVERY_THING   = 4
		};
	/// Enumeration for if to run internal tests or not.
	enum ERunTests { RUN_TESTS, NO_TESTS };
	///
	enum EMatRelations { MATRICES_INDEP_IMPS, MATRICES_ALLOW_DEP_IMPS };

	//@}

	/** @name Dimensionality of the decomposition */
	//@{

	///
	/** Return the number of rows in \c Gc and \c Gh.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>n > m</tt>
	 * </ul>
	 *
	 * The default implementation returns
	 * <tt>this->space_range()->dim() + this->space_null()->dim()</tt>.
	 */
	virtual size_type n() const;

	///
	/** Return the number of columns in \c Gc.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>m > 0</tt>
	 * </ul>
	 */
	virtual size_type m() const = 0;

	///
	/** Returns the rank of \c Gc(:,con_decomp()).
	 *
	 * Postconditions:<ul>
	 * <li> <tt>r =< m</tt>
	 * </ul>
	 *
	 * The default implementation returns
	 * <tt>this->space_range()->dim()</tt>.
	 */
	virtual size_type r() const;

	///
	/** Returns the range of the decomposed equalities.
	 *
	 * The default implementation returns <tt>Range1D(1,this->r())</tt>.
	 */
	virtual Range1D con_decomp() const;

	///
	/** Returns the range of the undecomposed equalities.
	 *
	 * The default implementation returns <tt>Range1D(this->r()+1,this->m())</tt>
	 * or <tt>Range1D::Invalid</tt> if <tt>this->r() == this->m()<tt>
	 */
	virtual Range1D con_undecomp() const;

	//@}

	/** @name Range and null vector spaces */
	//@{

	///
	/** Return a \c VectorSpace object for the range space.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>return->dim() == this->r()</tt>
	 * </ul>
	 */
	virtual const VectorSpace::space_ptr_t space_range() const = 0;

	///
	/** Return a \c VectorSpace object for the range space.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>return->dim() == this->n() - this->r()</tt>
	 * </ul>
	 */
	virtual const VectorSpace::space_ptr_t space_null() const = 0;

	//@}

	/** @name Matrix factories */
	//@{

	///
	/** Return a matrix factory object for <tt>Z</tt>
	 */
	virtual const mat_fcty_ptr_t factory_Z() const = 0;

	///
	/** Return a matrix factory object for <tt>Y</tt>
	 */
	virtual const mat_fcty_ptr_t factory_Y() const = 0;

	///
	/** Return a matrix factory object for <tt>R</tt>.
	 */
	virtual const mat_nonsing_fcty_ptr_t factory_R() const = 0;
	
	///
	/** Return a matrix factory object for <tt>Uz</tt>
	 */
	virtual const mat_fcty_ptr_t factory_Uz() const = 0;

	///
	/** Return a matrix factory object for <tt>Uy</tt>
	 */
	virtual const mat_fcty_ptr_t factory_Uy() const = 0;

	///
	/** Return a matrix factory object for <tt>Vz</tt>
	 */
	virtual const mat_fcty_ptr_t factory_Vz() const = 0;

	///
	/** Return a matrix factory object for <tt>Vy</tt>
	 */
	virtual const mat_fcty_ptr_t factory_Vy() const = 0;

	//@}

	/** @name Update range/null decomposition */
	//@{

	///
	/** Creates the range/null decomposition for <tt>Gc(:,con_decomp)'</tt>.
	 *
	 * The decomposition is based on the linearly independent columns \c Gc(:,con_decomp)
	 * of \c Gc
	 *
	 * <tt>Gc = [ Gc(:,con_decomp),  Gc(:,con_undecomp) ]</tt>
	 *
	 * Specifically this operation finds the matrices:
	 \verbatim
	 Z s.t. Gc(:,con_deomp)' * Z = 0
	 Y s.t. [Z  Y] is nonsingular
	 R = Gc(:,con_decomp)' * Y is nonsingular
	 Uz = Gc(:,con_undecomp)' * Z
	 Uy = Gc(:,con_undecomp)' * Y
	 Vz = Gh' * Z
	 Vy = Gh' * Y
	 \endverbatim
	 * If there is some problem creating the decomposition then exceptions
	 * with the base class \c std::exception may be thrown.  The meaning
	 * of these exceptions are more associated with the subclasses
	 * that implement this operation.
	 *
	 * The concrete types for <tt>Gc</tt>, <tt>Gh</tt>, <tt>Z</tt>, <tt>Y</tt>,
	 * <tt>Uz</tt>, <tt>Uy</tt>, <tt>Vz</tt> and <tt>Vy</tt> must be compatable with
	 * the concrete implementation of \c this or an <tt>InvalidMatrixType</tt> exeption
	 * will be thrown..
	 *
	 * Preconditions:<ul>
	 * <li> <tt>Gc.rows() == this->n()</tt> (throw \c std::invalid_argument)
	 * <li> <tt>Gc.cols() == this->m()</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>this->m() == this->r()</tt>] <tt>Uz == NULL</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>this->m() == this->r()</tt>] <tt>Uy == NULL</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>Gh == NULL</tt>] <tt>Gh->space_cols().is_compatible(Gc.space_cols()) == true</tt> (throw \c ???)
	 * <li> [<tt>Gh == NULL</tt>] <tt>Vz == NULL</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>Gh == NULL</tt>] <tt>Vy == NULL</tt> (throw \c std::invalid_argument)
	 * <li> <tt>Z!=NULL || Y!=NULL || R!=NULL || Uz!=NULL || Uy!=NULL || Vz!=NULL | Vy!=NULL</tt>
	 *      (throw \c std::invalid_argument)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>Gc(:,con_decomp())' * Z = 0</tt>
	 * <li> <tt>[ Y  Z ]</tt> nonsingular
	 * <li> [<tt>Z != NULL</tt>] <tt>Z.space_cols().is_compatible(Gc.space_cols()) == true)</tt>
	 * <li> [<tt>Z != NULL</tt>] <tt>Z.cols() == this->n() - this->r()</tt>
	 * <li> [<tt>Y != NULL</tt>] <tt>Y.space_cols().is_compatible(Gc.space_cols()) == true)</tt>
	 * <li> [<tt>Y != NULL</tt>] <tt>Y.cols() == this->r()</tt>
	 * <li> [<tt>R != NULL</tt>] <tt>R->space_cols().is_compatible(*Gc.space_cols()->sub_space(con_decomp())) == true</tt>
	 * <li> [<tt>R != NULL</tt>] <tt>R->space_rows().is_compatible(Y->space_rows()) == true</tt>
	 * <li> [<tt>Uz != NULL</tt>] <tt>Uz.space_cols().is_compatible(*Gc.space_rows()->sub_space(con_undecomp())) == true</tt>
	 * <li> [<tt>Uz != NULL</tt>] <tt>Uz.space_rows().is_compatible(Z.space_rows()) == true</tt>
	 * <li> [<tt>Uy != NULL</tt>] <tt>Uy.space_cols().is_compatible(*Gc.space_rows()->sub_space(con_undecomp())) == true</tt>
	 * <li> [<tt>Uy != NULL</tt>] <tt>Uy.space_rows().is_compatible(Y.space_rows()) == true</tt>
	 * <li> [<tt>Vz != NULL</tt>] <tt>Vz.space_cols().is_compatible(*Gh->space_rows()) == true</tt>
	 * <li> [<tt>Vz != NULL</tt>] <tt>Vz.space_rows().is_compatible(Z.space_rows())</tt>
	 * <li> [<tt>Vy != NULL</tt>] <tt>Vy.space_cols().is_compatible(*Gh->space_rows()) == true</tt>
	 * <li> [<tt>Vy != NULL</tt>] <tt>Vy.space_rows().is_compatible(Y.space_rows())</tt>
	 * <li> The behaviors of all of the participating matrices must not be altered by changes to the other matrices.
	 * </ul>
	 *
	 * @param  out [out] If <tt>out!=NULL</tt> then output is printed to this stream
	 *             depending on the value of \c olevel.
	 * @param olevel
	 *             [in] Determines the amount of output to print to \c *out.
	 *             The exact type of output is determined by the implementing
	 *             subclass but here is the sugguested behavior:<ul>
	 *             <li> \c PRINT_NONE : Don't print anything (same as <tt>out==NULL</tt>).
	 *             <li> \c PRINT_BASIC_INFO : Only print basic information about
	 *                  how the decomposition is formed.
	 *                  Amount of output = \c O(1).
	 *             <li> \c PRINT_VECTORS : Prints out important vectors computed
	 *                  durring the computations (usually only durring testing).
	 *                  This level is only useful for debugging.
	 *                  Amount of output = \c O(n).
	 *             <li> \c PRINT_EVERY_THING : Print out nearly every important
	 *                  quantity that is computed (except for the output matrices
	 *                  themselves, clients can do that) while the matrices are being
	 *                  formed or tests are being conducted.  This level is only
	 *                  useful for debugging.
	 *                  Amount of output = <tt>O(m*n)</tt>
	 *             </ul>
	 * @param test_what
	 *             [in] Determines if internal validation tests are performed.
	 *             The post conditions for the output matrices are not checked
	 *             internally.  This is something that client can (and should)
	 *             do independently (see \c DecompositionSystemTester).  Values:<ul>
	 *             <li> \c RUN_TESTS : As many validation/consistency tests
	 *                  are performed internally as possible.  If a test
	 *                  fails then a \c TestFailed execption will be thrown.
	 *                  The subclasses determine what the tests are and
	 *                  what failing a test means.
	 *             <li> \c NO_TEST : No tests are performed internally.  This is
	 *                  to allow the fastest possible execution.
	 *             </ul>
	 *             If a test fails, then a \c TestFailed exception will be thrown with
	 *             a helpful error message.
	 * @param  Gc  [in] The matrix for which the range/null decomposition is defined.
	 * @param  Gh  [in] An auxlillary matrix that will have the range/null decompositon applied to.
	 *             It is allowed for <tt>Gh == NULL</tt>.
	 * @param  Z   [out] On output represents the <tt>n x (n-r)</tt> null space	matrix such that
	 *             <tt>Gc(:,con_decomp) * Z == 0</tt>.  This matrix object must have been created
	 *             by <tt>this->factory_Z()->create()</tt>.
	 * @param  Y   [out] On output represents the <tt>n x r</tt> range space matrix	such that
	 *             <tt>[ Y  Z ]</tt> is nonsingular.  This matrix object must have been created
	 *             by <tt>this->factory_Y()->create()</tt>.
	 * @param  R   [out] On output represents the nonsingular <tt>r x r</tt> matrix <tt>Gc(:,con_decomp) * Y</tt>.
	 *             This matrix object must have been created by <tt>this->factory_R()->create()</tt>.
	 * @param  Uz  [in/out] If <tt>Uz != NULL</tt> (<tt>this->m() > this->r()</tt> only) then on output
	 *             <tt>*Uz</tt> represents the <tt>(m-r) x (n-r)</tt> matrix <tt>Gc(:,con_undecomp) * Z</tt>.
	 *             If <tt>this->m() == this->r()</tt> then <tt>Uz == NULL</tt> must be true.
	 *             If <tt>Uz!=NULL</tt>, then this matrix object must have been created by
	 *             <tt>this->factory_Uz()->create()</tt>.
	 * @param  Uy  [in/out] If <tt>Uy != NULL</tt> (<tt>this->m() > this->r()</tt> only) then on output
	 *             <tt>*Uy</tt> represents the <tt>(m-r) x r</tt> matrix <tt>Gc(:,con_undecomp) * Y</tt>.
	 *             If <tt>this->m() == this->r()</tt> then <tt>Uy == NULL</tt> must be true.
	 *             If <tt>Uy!=NULL</tt>, then this matrix object must have been created by
	 *             <tt>this->factory_Uy()->create()</tt>.
	 * @param  Vz  [in/out] If <tt>Vz != NULL</tt> (<tt>Gh != NULL</tt> only) then on output
	 *             <tt>*Vz</tt> represents the <tt>mI x (n-r)</tt> matrix <tt>Gh * Z</tt>.
	 *             If <tt>Gh == NULL</tt> then <tt>Vz == NULL</tt> must be true.
	 *             If <tt>Vz!=NULL</tt>, then this matrix object must have been created by
	 *             <tt>this->factory_Vz()->create()</tt>.
	 * @param  Vy  [in/out] If <tt>Vy != NULL</tt> (<tt>Gh != NULL</tt> only) then on output
	 *             <tt>*Vy</tt> represents the <tt>mI x r</tt> matrix <tt>Gh * Y</tt>.
	 *             If <tt>Gh == NULL</tt> then <tt>Vy == NULL</tt> must be true.
	 *             If <tt>Vy!=NULL</tt>, then this matrix object must have been created by
	 *             <tt>this->factory_Vy()->create()</tt>.
	 *
	 * Note that this method requires that all of the output matrix objects \c Z, \c Y, \c R,
	 * \c Uz, \c Uy, \c Vz and \c Vy must be independent of \c this and of each other.
	 * For example, the behavior of \c R must be be altered if \c this is destroyed or if
	 * \c Z is modified of destroyed.  This requirment constrains the implementations somewhat
	 * but makes things much easier for the client and gives the client much more power.
	 */
	virtual void update_decomp(
		std::ostream              *out
		,EOutputLevel             olevel
		,ERunTests                test_what
		,const MatrixWithOp       &Gc
		,const MatrixWithOp       *Gh
		,MatrixWithOp             *Z
		,MatrixWithOp             *Y
		,MatrixWithOpNonsingular  *R
		,MatrixWithOp             *Uz
		,MatrixWithOp             *Uy
		,MatrixWithOp             *Vz
		,MatrixWithOp             *Vy
		,EMatRelations            mat_rel = MATRICES_INDEP_IMPS
		) const = 0;
	
	///
	/** Print the sub-algorithm by which the decomposition is formed
	 *
	 */
	virtual void print_update_decomp(
		std::ostream& out, const std::string& leading_str ) const = 0;

	//@}
	
};	// end class DecompositionSystem

}	// end namespace ConstrainedOptimizationPack

#endif // DECOMPOSITION_SYSTEM_H
