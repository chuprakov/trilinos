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
#include "AbstractLinAlgPack/include/BasisSystemTester.h"
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
	typedef BasisSystem::SingularBasis                                    SingularBasis;
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

	/// Set the BasisSystem tester object
	STANDARD_COMPOSITION_MEMBERS( BasisSystemTester, basis_sys_tester )
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
		const VectorSpace::space_ptr_t     &space_x
		,const VectorSpace::space_ptr_t    &space_c
		,const VectorSpace::space_ptr_t    &space_h
		,const basis_sys_ptr_t             &basis_sys
		,const basis_sys_tester_ptr_t      &basis_sys_tester
		,EExplicitImplicit                 D_imp
		,EExplicitImplicit                 Uz_imp
		,EExplicitImplicit                 Vz_imp
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

	/// Returns <tt>this->space_x()->dim()</tt>.
	size_type n() const;
	/// Returns <tt>this->space_c()->dim()</tt>.
	size_type m() const;
	/// Returns <tt>this->basis_sys()->equ_decomp().size()</tt>.
	size_type r() const;
	/// Returns <tt>this->space_x()->sub_space(var_dep)</tt>
	const VectorSpace::space_ptr_t space_range() const;
	/// Returns <tt>this->space_x()->sub_space(var_indep)</tt>
	const VectorSpace::space_ptr_t space_null() const;
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
	 * <li> [<tt>test_what == TEST</tt>] <tt>this->basis_sys_tester().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
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
		,EMatRelations            mat_rel
		) const;
	///
	void print_update_decomp(
		std::ostream& out, const std::string& leading_str ) const;

	//@}

protected:

	///
	/** Overridden by subclasses to uninitialized Y, R, Uy and Vy and return C if referenced.
	 *
	 * Note that the returned smart pointer to \c C may have <tt>return.has_ownership() == false</tt>
	 * in which case this will not be a shared resource with any other object (at least not
	 * on this level).
	 *
	 * ToDo: Finish documentatation!
`	 */
	virtual	mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t	uninitialize_matrices(
		std::ostream                             *out
		,EOutputLevel                            olevel
		,MatrixWithOp                            *Y
		,MatrixWithOpNonsingular                 *R
		,MatrixWithOp                            *Uy
		,MatrixWithOp                            *Vy
		) const = 0;

	///
	/** Overridden by subclasses to initialize Y, R, Uy and Vy given C and D.
	 *
 	 * If C_ptr.has_ownership() == false, then the subclass implementation of this
	 * method will use clone_mwons() to clone it so that the output matrices are
	 * independent.
	 *
	 * ToDo: Finish documentatation!
`	 */
	virtual void initialize_matrices(
		std::ostream                                           *out
		,EOutputLevel                                          olevel
		,const mat_nonsing_fcty_ptr_t::element_type::obj_ptr_t &C_ptr
		,const mat_fcty_ptr_t::element_type::obj_ptr_t         &D_ptr
		,MatrixWithOp                                          *Y
		,MatrixWithOpNonsingular                               *R
		,MatrixWithOp                                          *Uy
		,MatrixWithOp                                          *Vy
		,EMatRelations                                         mat_rel
		) const = 0;

	///
	/** Print the sub-algorithm by which the matrices Y, R, Uy and Uy are updated.
	 *
	 * ToDo: Finish documentatation!
	 */
	virtual void print_update_matrices(
		std::ostream& out, const std::string& leading_str ) const = 0;

private:

	// //////////////////////////////////
	// Private data members
	
#ifdef DOXYGEN_COMPILE
	AbstractLinAlgPack::BasisSystem       *basis_sys;
	VectorSpace                           *space_x;
	VectorSpace                           *space_c;
	VectorSpace                           *space_h;
	VectorSpace                           *space_range;
	VectorSpace                           *space_null;
#else
	basis_sys_ptr_t                       basis_sys_;
	VectorSpace::space_ptr_t              space_x_;
	VectorSpace::space_ptr_t              space_c_;
	VectorSpace::space_ptr_t              space_h_;
	VectorSpace::space_ptr_t              space_range_;
	VectorSpace::space_ptr_t              space_null_;
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
