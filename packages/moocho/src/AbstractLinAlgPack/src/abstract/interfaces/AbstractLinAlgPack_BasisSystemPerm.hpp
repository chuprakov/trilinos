// ////////////////////////////////////////////////////////////
// BasisSystemPerm.hpp
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

#include "BasisSystem.hpp"

namespace AbstractLinAlgPack {

///
/** Interface for setting and selecting a basis from the Jacobian
 * from a set of equations.
 *
 * ToDo: Finish documentation!
 */
class BasisSystemPerm : public BasisSystem {
public:

	/** @name Public types */
	//@{

	///
	typedef Teuchos::RefCountPtr<
		const Teuchos::AbstractFactory<Permutation> >   perm_fcty_ptr_t;

	//@}


	///
	/** Required constructor (calls <tt>initialize()</tt>).
	 */
	BasisSystemPerm(
		const mat_sym_fcty_ptr_t             &factory_transDtD
		,const mat_sym_nonsing_fcty_ptr_t    &factory_S
		)
		:BasisSystem(factory_transDtD,factory_S)
			{}

	/** @name Permutation factories */
	//@{

	///
	virtual const perm_fcty_ptr_t   factory_P_var() const = 0;
	///
	virtual const perm_fcty_ptr_t   factory_P_equ() const = 0;

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
		,const MatrixOp            &Gc
		,MatrixOpNonsing           *C
		,MatrixOp                  *D
		,MatrixOp                  *GcUP
		,EMatRelations             mat_rel = MATRICES_INDEP_IMPS
		,std::ostream              *out    = NULL
		) = 0;

	///
	/** Select a basis.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void select_basis(
		const Vector               *nu
		,MatrixOp                  *Gc
		,Permutation               *P_var
		,Range1D                   *var_dep
		,Permutation               *P_equ
		,Range1D                   *equ_decomp
		,MatrixOpNonsing           *C
		,MatrixOp                  *D
		,MatrixOp                  *GcUP
		,EMatRelations             mat_rel = MATRICES_INDEP_IMPS
		,std::ostream              *out    = NULL
		) = 0;
	
	//@}

private:
	// not defined and not to be called
	BasisSystemPerm();

}; // end class BasisSystemPerm

}  // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_PERM_SYSTEM_H
