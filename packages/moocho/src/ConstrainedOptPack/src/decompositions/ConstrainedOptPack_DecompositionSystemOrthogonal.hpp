// ////////////////////////////////////////////////////////////////////
// DecompositionSystemOrthogonal.h
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

#ifndef DECOMPOSITION_SYSTEM_ORTHOGONAL_H
#define DECOMPOSITION_SYSTEM_ORTHOGONAL_H

#include "DecompositionSystemVarReduct.h"
#include "StandardCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Orthogonal variable reduction subclass.
 *
 * This is the interface for the coordinate variable reduction decomposition
 * where:
 \verbatim

  Y = [      I     ] = [  I  ]  (class MatrixIdentConcatStd)
      [ N'*inv(C') ]   [ -D' ]

  R = Gc(:,con_decomp)'*Y = [ C   N ] * [     I      ] = (C + N*N'*inv(C'))
                                        [ N'*inv(C') ]
    = C*(I + inv(C)*N*N'*inv(C'))
    = C*(I + D*D')

  Uy = Gc(:,con_undecomp)'*Y = [ E  F ] * [  I  ]
                                          [ -D' ]
      = E - F*D'

  Vy = Gh'*Y = [ Gh(var_dep,:)'   Gh(var_indep,:)' ] * [  I  ]
                                                       [ -D' ]
     =  = Gh(var_dep,:)' -  = Gh(var_indep,:)'*D'
 \endverbatim
 * The solution of the
 *
 * For now the copy constructor and the assignment operator are not defined.
 */
class DecompositionSystemOrthogonal : public DecompositionSystemVarReduct {
public:

	/** @name Constructors / initializers */
	//@{

	/// Set the strategy object for forming <tt>S = I + D'*D</tt>.
	STANDARD_CONST_COMPOSITION_MEMBERS( VarReductOrthog_Strategy, var_reduct_orthog_strategy )

	///
	DecompositionSystemOrthogonal(
		const VectorSpace::space_ptr_t           &space_x                    = MemMngPack::null
		,const VectorSpace::space_ptr_t          &space_c                    = MemMngPack::null
		,const VectorSpace::space_ptr_t          &space_h                    = MemMngPack::null
		,const basis_sys_ptr_t                   &basis_sys                  = MemMngPack::null
		,const basis_sys_tester_ptr_t            &basis_sys_tester           = MemMngPack::null
		,const var_reduct_orthog_strategy_ptr_t  &var_reduct_orthog_strategy = MemMngPack::null
		,EExplicitImplicit                       D_imp                       = MAT_IMP_EXPLICIT
		,EExplicitImplicit                       Uz_imp                      = MAT_IMP_EXPLICIT
		,EExplicitImplicit                       Vz_imp                      = MAT_IMP_EXPLICIT
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

	/** @name Overridden from DecompositionSystemVarReduct */
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

private:
	
	// ////////////////////////
	// Private types

	typedef MemMngPack::ref_count_ptr<MatrixSymWithOpNonsingular>  S_ptr_t;

	// ////////////////////////
	// Private data members

	mutable S_ptr_t        S_ptr_;
	// We just hold on to this between calls of uninitialize_matrices() and initialize_matrices()

	// ////////////////////////
	// Private member functions

	// Not defined and not to be called
	DecompositionSystemOrthogonal(const DecompositionSystemOrthogonal&);
	DecompositionSystemOrthogonal& operator=(const DecompositionSystemOrthogonal&);

};	// end class DecompositionSystemOrthogonal

}	// end namespace ConstrainedOptimizationPack

#endif	// DECOMPOSITION_SYSTEM_ORTHOGONAL_H
