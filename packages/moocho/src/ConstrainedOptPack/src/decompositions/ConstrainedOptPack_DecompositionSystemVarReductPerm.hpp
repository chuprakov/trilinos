// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystemVarReductPerm.h
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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H

#include <stdexcept>

#include "DecompositionSystemVarReduct.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"

namespace ConstrainedOptimizationPack {

///
/** Specialization interface of \c DecompositonSystem that allows basis permutations.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemVarReductPerm : public DecompositionSystemVarReduct {
public:

	/** @name Public types */
	//@{

	///
	typedef MemMngPack::ref_count_ptr<
		const MemMngPack::AbstractFactory<Permutation> >         perm_fcty_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	///
	DecompositionSystemVarReductPerm(
		EExplicitImplicit     D_imp    = MAT_IMP_AUTO
		,EExplicitImplicit    Uz_imp   = MAT_IMP_AUTO
		,EExplicitImplicit    Vz_imp   = MAT_IMP_AUTO
		)
		:DecompositionSystemVarReduct(D_imp,Uz_imp,Vz_imp)
	{}

	//@}

	/** @name Permutation factories */
	//@{

	///
	virtual const perm_fcty_ptr_t   factory_P_var() const = 0;
	///
	virtual const perm_fcty_ptr_t   factory_P_equ() const = 0;

	//@}

	/** @name Setting or selecting a new decomposition */
	//@{

	///
	/** Query to see if a current basis is already selected.
	 */
	virtual bool has_basis() const = 0;

	///
	/** Client selects the basis for <tt>Gc(:,con_decomp)'</tt>.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void set_decomp(
		std::ostream              *out
		,EOutputLevel             olevel
		,ERunTests                test_what
		,const Permutation        &P_var
		,const Range1D            &var_dep
		,const Permutation        *P_equ
		,const Range1D            *equ_decomp
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
		) = 0;
	
	///
	/** Client asks decompostion system object to select the basis for <tt>Gc(:,con_decomp)'</tt>.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void select_decomp(
		std::ostream              *out
		,EOutputLevel             olevel
		,ERunTests                test_what
		,const VectorWithOp       *nu
		,MatrixWithOp             *Gc
		,MatrixWithOp             *Gh
		,Permutation              *P_var
		,Range1D                  *var_dep
		,Permutation              *P_equ
		,Range1D                  *equ_decomp
		,MatrixWithOp             *Z
		,MatrixWithOp             *Y
		,MatrixWithOpNonsingular  *R
		,MatrixWithOp             *Uz
		,MatrixWithOp             *Uy
		,MatrixWithOp             *Vz
		,MatrixWithOp             *Vy
		,EMatRelations            mat_rel = MATRICES_INDEP_IMPS
		) = 0;

	//@}
	
};	// end class DecompositionSystemVarReductPerm

}	// end namespace ConstrainedOptimizationPack

#endif // DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_H
