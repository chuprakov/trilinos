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

namespace ConstrainedOptimizationPack {

///
/** This class abstracts a decomposition choice for the
 * range space \a Y, and null space \a Z, matrices for a linearly
 * independent set of columns of \a Gc.
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
 * <ul>
 * <li> <tt>Gc(:,con_decomp)' * Z = 0</tt> (null space property)
 * <li> <tt>[Y , Z]</tt> is nonsingular
 * </ul>
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystem {
public:

	/** @name Public types */
	//@{

	///
	class InvalidMatrixType : public std::logic_error
	{public: InvalidMatrixType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	/** @name Dimensionality of the decomposition */
	//@{

	///
	/** Return the number of rows in \c Gc and \c Gh.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>n > m</tt>
	 * </ul>
	 */
	virtual size_type n() const = 0;

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
	 */
	virtual size_type r() const = 0;

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
	 *
	 * \c Z s.t. <tt>Gc(:,con_deomp)' * Z = 0</tt><br>
	 * \c Y s.t. <tt>[Z  Y]</tt> nonsingular<br>
	 * <tt>R = Gc(:,con_decomp)' * Y</tt> nonsingular<br>
	 * <tt>Uz = Gc(:,con_undecomp)' * Z</tt><br>
	 * <tt>Uy = Gc(:,con_undecomp)' * Y</tt><br>
	 * <tt>Vz = Gh' * Z</tt><br>
	 * <tt>Vy = Gh' * Y</tt><br>
	 *
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
	 * <li> <tt>Z != NULL</tt> (throw \c std::invalid_argument)
	 * <li> <tt>Y != NULL</tt> (throw \c std::invalid_argument)
	 * <li> <tt>R != NULL</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>this->m() == this->r()</tt>] <tt>Uz == NULL</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>this->m() == this->r()</tt>] <tt>Uy == NULL</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>Gh == NULL</tt>] <tt>Gh->space_cols().is_compatible(Gc.space_cols()) == true</tt> (throw \c ???)
	 * <li> [<tt>Gh == NULL</tt>] <tt>Vz == NULL</tt> (throw \c std::invalid_argument)
	 * <li> [<tt>Gh == NULL</tt>] <tt>Vy == NULL</tt> (throw \c std::invalid_argument)
	 * </ul> 
	 *
	 * Postconditions:<ul>
	 * <li> <tt>Gc(:,con_decomp())' * Z = 0</tt>
	 * <li> <tt>[ Y  Z ]</tt> nonsingular
	 * <li> <tt>Z.space_cols().is_compatible(Gc.space_cols()) == true)</tt>
	 * <li> <tt>Z.cols() == this->n() - this->r()</tt>
	 * <li> <tt>Y.space_cols().is_compatible(Gc.space_cols()) == true)</tt>
	 * <li> <tt>Y.cols() == this->r()</tt>
	 * <li> <tt>R->space_cols().is_compatible(*Gc.space_cols()->sub_space(con_decomp())) == true</tt>
	 * <li> <tt>R->space_rows().is_compatible(Y->space_rows()) == true</tt>
	 * <li> [<tt>Uz != NULL</tt>] <tt>Uz.space_cols().is_compatible(*Gc.space_rows()->sub_space(con_undecomp())) == true</tt>
	 * <li> [<tt>Uz != NULL</tt>] <tt>Uz.space_rows().is_compatible(Z.space_rows()) == true</tt>
	 * <li> [<tt>Uy != NULL</tt>] <tt>Uy.space_cols().is_compatible(*Gc.space_rows()->sub_space(con_undecomp())) == true</tt>
	 * <li> [<tt>Uy != NULL</tt>] <tt>Uy.space_rows().is_compatible(Y.space_rows()) == true</tt>
	 * <li> [<tt>Vz != NULL</tt>] <tt>Vz.space_cols().is_compatible(*Gh->space_rows()) == true</tt>
	 * <li> [<tt>Vz != NULL</tt>] <tt>Vz.space_rows().is_compatible(Z.space_rows())</tt>
	 * <li> [<tt>Vy != NULL</tt>] <tt>Vy.space_cols().is_compatible(*Gh->space_rows()) == true</tt>
	 * <li> [<tt>Vy != NULL</tt>] <tt>Vy.space_rows().is_compatible(Y.space_rows())</tt>
	 * </ul>
	 *
	 * @param  Gc  [in] The matrix for which the range/null decomposition is defined.
	 * @param  Gh  [in] An auxlillary matrix that will have the range/null decompositon applied to.
	 *             It is allowed for <tt>Gh == NULL</tt>.
	 * @param  Z   [out] On output represents the <tt>n x (n-r)</tt> null space	matrix such that
	 *             <tt>Gc(:,con_decomp) * Z == 0</tt>.
	 * @param  Y   [out] On output represents the <tt>n x r</tt> range space matrix	such that
	 *             <tt>[ Y  Z ]</tt> is nonsingular.
	 * @param  R   [out] On output represents the nonsingular <tt>r x r</tt> matrix <tt>Gc(:,con_decomp) * Y</tt>.
	 * @param  Uz  [in/out] If <tt>Uz != NULL</tt> (<tt>this->m() > this->r()</tt> only) then on output
	 *             <tt>*Uz</tt> represents the <tt>(m-r) x (n-r)</tt> matrix <tt>Gc(:,con_undecomp) * Z</tt>.
	 *             If <tt>this->m() == this->r()</tt> then <tt>Uz == NULL</tt> must be true.
	 * @param  Uy  [in/out] If <tt>Uy != NULL</tt> (<tt>this->m() > this->r()</tt> only) then on output
	 *             <tt>*Uy</tt> represents the <tt>(m-r) x r</tt> matrix <tt>Gc(:,con_undecomp) * Y</tt>.
	 *             If <tt>this->m() == this->r()</tt> then <tt>Uy == NULL</tt> must be true.
	 * @param  Vz  [in/out] If <tt>Vz != NULL</tt> (<tt>Gh != NULL</tt> only) then on output
	 *             <tt>*Vz</tt> represents the <tt>mI x (n-r)</tt> matrix <tt>Gh * Z</tt>.
	 *             If <tt>Gh == NULL</tt> then <tt>Vz == NULL</tt> must be true.
	 * @param  Vy  [in/out] If <tt>Vy != NULL</tt> (<tt>Gh != NULL</tt> only) then on output
	 *             <tt>*Vy</tt> represents the <tt>mI x r</tt> matrix <tt>Gh * Y</tt>.
	 *             If <tt>Gh == NULL</tt> then <tt>Vy == NULL</tt> must be true.
	 */
	virtual void update_decomp(
		const MatrixWithOp        &Gc
		,const MatrixWithOp       *Gh
		,MatrixWithOp             *Z
		,MatrixWithOp             *Y
		,MatrixWithOpNonsingular  *R
		,MatrixWithOp             *Uz
		,MatrixWithOp             *Uy
		,MatrixWithOp             *Vz
		,MatrixWithOp             *Vy
		) const = 0;

	//@}
	
};	// end class DecompositionSystem

}	// end namespace ConstrainedOptimizationPack

#endif // DECOMPOSITION_SYSTEM_H
