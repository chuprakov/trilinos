// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystemVarReduct.h
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

#include "DecompositionSystem.h"
#include "AbstractLinAlgPack/include/BasisSystem.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "StandardCompositionMacros.h"
#include "StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

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

  Vz  = Gh(var_indep,:)' + Gh(var_dep,:)'*D

      where:
           C = Gc(var_dep,con_decomp)'     [nonsingular]
           N = Gc(var_indep,con_decomp)'
           E = Gc(var_dep,con_undecomp)'
           F = Gc(var_indep,con_undecomp)'
           D = -inv(C) * N
 \endverbatim
 *
 * Above, \a C is a <tt>r x r</tt> nonsingular matrix.  Subclasses define
 * how \c Y is defined which in turn determines how \c R, \c Uy and \c Vy are
 * defined.
 *
 * The implementation of this subclass is completly determined by an aggregate
 * <tt>BasisSytem</tt> object.  Since the <tt>BasisSystem</tt> interface does
 * not allow any permutations of the basis to be performed ???.
 *
 * <b>Subclass implementors notes:</b>
 *
 * It is up to subclasses to override \c factory_R(), \c factory_Uy() and
 * \c factory_Vy() in order to define the types for \c R, \c Uy and \c Vy
 * respectively.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemVarReduct : public DecompositionSystem {
public:

	/** @name Public types */
	//@{

	///
	typedef DecompositionSystem	                                          inherited;
	///
	typedef BasisSystem::NumericallySingular                              NumericallySingular;
	///
	typedef ReferenceCountingPack::ref_count_ptr<const BasisSystem>       basis_sys_ptr_t;
	///
	enum EExplicitImplicit {
		MAT_IMP_EXPLICIT
		,MAT_IMP_IMPLICIT
		,MAT_IMP_AUTO
	};

	//@}	

	/** @name Constructors / initializers */
	//@{

	/// Set whether to use explicit or implicit <tt>D = -inv(C)*N</tt> matrix.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,D_imp)
	/// Set whether to use explicit or implicit <tt>Uz = F + E * D</tt> matrix.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,Uz_imp)
	/// Set whether to use explicit or implicit <tt>Vz = Gh(var_indep,:)' + Gh(var_dep,:)'*D</tt> matrix.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(EExplicitImplicit,Vz_imp)

	///
	/** Construct a variable reduction decomposition.
	 *
	 * Calls <tt>this->initialize()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> ???
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->factory_Z().get() != NULL</tt>
	 * </ul>
	 */
	DecompositionSystemVarReduct(
		const VectorSpace::space_ptr_t     &space_x   = NULL
		,const VectorSpace::space_ptr_t    &space_c   = NULL
		,const VectorSpace::space_ptr_t    &space_h   = NULL
		,const basis_sys_ptr_t             &basis_sys = NULL
		,EExplicitImplicit                 D_imp      = MAT_IMP_AUTO
		,EExplicitImplicit                 Uz_imp     = MAT_IMP_AUTO
		,EExplicitImplicit                 Vz_imp     = MAT_IMP_AUTO
		);

	//@}

	///
	/** Initialize.
	 *
	 * @param  space_x
	 *             [in] Vector space for variables \c x.
	 * @param  space_c
	 *             [in] Vector space for general equalities \c c.
	 * @param  space_h
	 *             [in] Vector space for general inequalities \c h.
	 * @param  basis_sys
	 *             [in] The <tt>BasisSystem</tt> object that defines the
	 *             variable reduction that the decomposition is based on.
	 *             This object must be fully initialized before this
	 *             method is called.  The object <tt>*basis_sys</tt> must
	 *             not be altered while \c this is still using it.  It is
	 *             allowed for <tt>basis_sys.get() == NULL</tt> in which
	 *             case \c this will not be fully initialized.
	 *
	 * Preconditions:<ul>
	 * <li> ToDo: Spell these out!
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> ToDo: Spell these out!
	 * </ul>
	 */
	void initialize(
		const VectorSpace::space_ptr_t     &space_x
		,const VectorSpace::space_ptr_t    &space_c
		,const VectorSpace::space_ptr_t    &space_h
		,const basis_sys_ptr_t             &basis_sys
		);

	/** @name Access */
	//@{
	
	///
	const VectorSpace::space_ptr_t& space_x() const;
	///
	const VectorSpace::space_ptr_t& space_c() const;
	///
	const VectorSpace::space_ptr_t& space_h() const;
	///
	const basis_sys_ptr_t& basis_sys() const;

	//@}

	/** @name Overridden from DecompositionSystem */
	//@{

	///
	size_type n() const;
	///
	size_type m() const;
	///
	size_type r() const;
	///
	const mat_fcty_ptr_t factory_Z() const;
	///
	const mat_fcty_ptr_t factory_Uz() const;
	///
	const mat_fcty_ptr_t factory_Vz() const;
	///
	/**
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->space_x().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> <tt>this->space_c().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li> [<tt>Gh != NULL</tt>] <tt>this->space_h().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
	 * </ul>
	 */
	void update_decomp(
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
		) const;

	//@}

protected:

	// ToDo: Add pure virtual method to uninitialize Y, R, Uy and Vy and return C if references!
	// ToDo: Add pure virtual method to initialize Y, R, Uy and Vy given C and D!

private:

	// //////////////////////////////////
	// Private data members
	
#ifdef DOXYGEN_COMPILE
	AbstractLinAlgPack::BasisSystem       *basis_sys;
	VectorSpace                           *space_x;
	VectorSpace                           *space_c;
	VectorSpace                           *space_h;
#else
	basis_sys_ptr_t                       basis_sys_;
	VectorSpace::space_ptr_t              space_x_;
	VectorSpace::space_ptr_t              space_c_;
	VectorSpace::space_ptr_t              space_h_;
#endif
	// //////////////////////////////////
	// Private member functions

	// not defined and not to be called!
	DecompositionSystemVarReduct(const DecompositionSystemVarReduct&);
	DecompositionSystemVarReduct& operator=(const DecompositionSystemVarReduct&);

};	// end class DecompositionSystemVarReduct

// //////////////////////////////////////////
// Inline members

inline
const VectorSpace::space_ptr_t&
DecompositionSystemVarReduct::space_x() const
{
	return space_x_;
}

inline
const VectorSpace::space_ptr_t&
DecompositionSystemVarReduct::space_c() const
{
	return space_c_;
}

inline
const VectorSpace::space_ptr_t&
DecompositionSystemVarReduct::space_h() const
{
	return space_h_;
}

inline
const DecompositionSystemVarReduct::basis_sys_ptr_t&
DecompositionSystemVarReduct::basis_sys() const
{
	return basis_sys_;
}

}	// end namespace ConstrainedOptimizationPack

#endif	// DECOMPOSITION_SYSTEM_VAR_REDUCT_H
