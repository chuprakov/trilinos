// ////////////////////////////////////////////////////////////
// BasisSystemPerm.h
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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H

#include "BasisSystem.h"

namespace AbstractLinAlgPack {

///
/** Interface for setting and selecting a basis from the Jacobian
 * from a set of equation.
 *
 * ToDo: Finish documentation!
 */
class BasisSystemPerm : public BasisSystem {
public:

	/** @name Public types */
	//@{

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractFactoryPack::AbstractFactory<Permutation> >         perm_fcty_ptr_t;

	//@}

	/** @name Basis selection / manipulation */
	//@{

	///
	/** Factor a basis selected by the client.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void set_basis(
		const Permutation          &P_var
		,const Range1D             &var_dep
		,const Permutation         *P_equ
		,const Range1D             *equ_decomp
		,const Permutation         *P_inequ
		,const Range1D             *inequ_decomp
		,const MatrixWithOp        *Gc
		,const MatrixWithOp        *Gh
		,MatrixWithOpNonsingular   *C
		,MatrixWithOp              *D
		,MatrixWithOp              *GcUP
		,MatrixWithOp              *GhUP
		,EMatRelations              mat_rel = MATRICES_INDEP_IMPS
		) const = 0;

	///
	/** Select a basis.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void select_basis(
		const VectorWithOp          &nu
		,const VectorWithOp         &lambdaI
		,MatrixWithOp               *Gc
		,MatrixWithOp               *Gh
		,Permutation                *P_var
		,Range1D                    *var_dep
		,Permutation                *P_equ
		,Range1D                    *equ_decomp
		,Permutation                *P_inequ
		,Range1D                    *inequ_decomp
		,MatrixWithOpNonsingular    *C
		,MatrixWithOp               *D
		,MatrixWithOp               *GcUP
		,MatrixWithOp               *GhUP
		,EMatRelations              mat_rel = MATRICES_INDEP_IMPS
		) const = 0;
	
	//@}

}; // end class BasisSystem

}; // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H
