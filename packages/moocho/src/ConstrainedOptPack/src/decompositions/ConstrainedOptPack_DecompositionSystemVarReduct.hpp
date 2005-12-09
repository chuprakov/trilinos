// /////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_DecompositionSystemVarReduct.hpp
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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_H

#include "ConstrainedOptPack_DecompositionSystem.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"

namespace ConstrainedOptPack {

///
/** Specialization of \c DecompositionSystem for variable reduction decompositions.
 *
 * This interface abstracts a variable reduction decomposition where:
 *
 \verbatim
  
  Gc' = [ C  N ] 
        [ E  F ]

  Z   = [ D ]
        [ I ]

  Uz  = F + E * D

      where:
           C = Gc(var_dep,con_decomp)'     [nonsingular]
           N = Gc(var_indep,con_decomp)'
           E = Gc(var_dep,con_undecomp)'
           F = Gc(var_indep,con_undecomp)'
           D = -inv(C) * N
 \endverbatim
 *
 * This interface simply allows clients to determine how \c D and \c Uz
 * are implemented (implicitly or explicity).
 */
class DecompositionSystemVarReduct : public DecompositionSystem {
public:

	/** @name Public types */
	//@{

	///
	enum EExplicitImplicit {
		MAT_IMP_EXPLICIT
		,MAT_IMP_IMPLICIT
		,MAT_IMP_AUTO
	};

	//@}

	/** @name Matrix representations */
	//@{

	/// Set whether to use explicit or implicit <tt>D = -inv(C)*N</tt> matrix.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,D_imp)
	/// Set whether to use explicit or implicit <tt>Uz = F + E * D</tt> matrix.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,Uz_imp)

		// ToDo: The above could be implemented as pure virtual funtions if needed later!

	//@}
		
	/** @name Constructors / initializers */
	//@{

	///
	DecompositionSystemVarReduct(
		EExplicitImplicit     D_imp    = MAT_IMP_AUTO
		,EExplicitImplicit    Uz_imp   = MAT_IMP_AUTO
		)
		:D_imp_(D_imp), Uz_imp_(Uz_imp)
	{}

	//@}

	/** @name Variable partitions. */
	//@{

	///
	virtual Range1D var_indep() const = 0;
	///
	virtual Range1D var_dep() const = 0;

	//@}

private:

	// not defined and not to be called!
	DecompositionSystemVarReduct(const DecompositionSystemVarReduct&);
	DecompositionSystemVarReduct& operator=(const DecompositionSystemVarReduct&);

};	// end class DecompositionSystemVarReduct

}	// end namespace ConstrainedOptPack

#endif	// DECOMPOSITION_SYSTEM_VAR_REDUCT_H
