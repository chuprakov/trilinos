// ////////////////////////////////////////////////////////////////////
// DecompositionSystemCoordinate.hpp
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

#ifndef DECOMPOSITION_SYSTEM_COORDINATE_H
#define DECOMPOSITION_SYSTEM_COORDINATE_H

#include "DecompositionSystemVarReductImp.hpp"
#include "StandardCompositionMacros.hpp"

namespace ConstrainedOptimizationPack {

///
/** Coordinate variable reduction subclass.
 *
 * This is the interface for the coordinate variable reduction decomposition
 * where:
 \verbatim

  Y = [ I ]   (class MatrixIdentConcatStd with MatrixZero)
      [ 0 ]

  R = Gc(:,con_decomp)'*Y = [ C   N ] * [ I ] = C
                                        [ 0 ]

  Uy = Gc(:,con_undecomp)'*Y = [ E  F ] * [ I ] = E
                                          [ 0 ]

  Vy = Gh'*Y = [ Gh(var_dep,:)'   Gh(var_indep,:)' ] * [ I ] = Gh(var_dep,:)'
                                                       [ 0 ]
 \endverbatim
 * The solution of the
 *
 * For now the copy constructor and the assignment operator are not defined.
 */
class DecompositionSystemCoordinate : public DecompositionSystemVarReductImp {
public:

	/** @name Constructors / initializers */
	//@{

	///
	DecompositionSystemCoordinate(
		const VectorSpace::space_ptr_t     &space_x          = MemMngPack::null
		,const VectorSpace::space_ptr_t    &space_c          = MemMngPack::null
		,const VectorSpace::space_ptr_t    &space_h          = MemMngPack::null
		,const basis_sys_ptr_t             &basis_sys        = MemMngPack::null
		,const basis_sys_tester_ptr_t      &basis_sys_tester = MemMngPack::null
		,EExplicitImplicit                 D_imp             = MAT_IMP_AUTO
		,EExplicitImplicit                 Uz_imp            = MAT_IMP_AUTO
		,EExplicitImplicit                 Vz_imp            = MAT_IMP_AUTO
		);

	//@}

	/** @name Overridden from DecompositionSystem */
	//@{

	///
	const mat_fcty_ptr_t factory_Y() const;
	///
	const mat_nonsing_fcty_ptr_t factory_R() const;
	///
	const mat_fcty_ptr_t factory_Uy() const;
	///
	const mat_fcty_ptr_t factory_Vy() const;

	//@}

protected:

	/** @name Overridden from DecompositionSystemVarReductImp */
	//@{

	///
	mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t	uninitialize_matrices(
		std::ostream                                           *out
		,EOutputLevel                                          olevel
		,MatrixWithOp                                          *Y
		,MatrixWithOpNonsingular                               *R
		,MatrixWithOp                                          *Uy
		,MatrixWithOp                                          *Vy
		) const;
	///
	void initialize_matrices(
		std::ostream                                           *out
		,EOutputLevel                                          olevel
		,const mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t &C
		,const mat_fcty_ptr_t::element_type::obj_ptr_t         &D
		,MatrixWithOp                                          *Y
		,MatrixWithOpNonsingular                               *R
		,MatrixWithOp                                          *Uy
		,MatrixWithOp                                          *Vy
		,EMatRelations                                         mat_rel
		) const;
	///
	void print_update_matrices(
		std::ostream& out, const std::string& leading_str ) const;

	//@}

};	// end class DecompositionSystemCoordinate

}	// end namespace ConstrainedOptimizationPack

#endif	// DECOMPOSITION_SYSTEM_COORDINATE_H
