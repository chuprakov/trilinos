// ////////////////////////////////////////////////////////////
// BasisSystem.h
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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_H

#include "AbstractLinAlgPackTypes.h"
#include "ref_count_ptr.h"

namespace AbstractLinAlgPack {

///
/** Interface for the creation and maintainance of a basis matrix
 * for a decomposition of linearlized constriants.
 *
 * <b>Overview:</b>
 *
 * This interface is designed to take the Jacobian for a sub-set of
 * equality and inequality constraints and to create a basis matrix.
 * Assume we have the folloing linealrized equality and inequality constraints:
 *
 \f[ \nabla c^T d + c = 0 \f]
 \f[ h^l \le \nabla h^T d + h \le h^u \f]
 *
 * The C++ identifiers given to \f$ \nabla c \f$ and \f$ \nabla h \f$ are <tt>Gc</tt> and <tt>Gh</tt>
 * respectively.
 *
 * In this basis interface we will assume that <tt>d</tt>, <tt>c</tt> and <tt>h</tt> are
 * sorted such that we define the following sets (given the partitioning matrices
 * \f$ Q^{x} = \left[\begin{array}{cc} Q^{xD} & Q^{xI} \end{array}\right] \f$,
 * \f$ Q^{c} = \left[\begin{array}{cc} Q^{cD} & Q^{cU} \end{array}\right] \f$,
 * \f$ Q^{h} = \left[\begin{array}{cc} Q^{hD} & Q^{hU} \end{array}\right] \f$
 * ):
 * <ul>
 * <li> <tt>d(var_dep)</tt> (\f$ d^D = (Q^{xD})^T d\f$) : Dependent (i.e. basis) variables.
 * <li> <tt>d(var_indep)</tt> (\f$ d^I = (Q^{xI})^T d\f$) : Independent (i.e. nonbasic) variables.
 * <li> <tt>c(equ_decomp)</tt> (\f$ c^D = (Q^{cD})^T c\f$) : Decomposed equality constriants.
 * <li> <tt>c(equ_undecomp)</tt> (\f$ c^U = (Q^{cU})^T c\f$): Undecomposed equality constriants.
 * <li> <tt>h(inequ_decomp)</tt> (\f$ h^D = (Q^{hD})^T h\f$) : Decomposed inequality constriants.
 * <li> <tt>h(inequ_undecomp)</tt> (\f$ h^U = (Q^{hU})^T h\f$) : Undecomposed inequality constriants.
 * </ul>
 * Given these partitionings we can define a basis matrix \a C for the
 * following Jacobian sub-matrices (in mathematical and Matlab-like notation):
 \f[
    C = \left[\begin{array}{c}
            (Q^{cD})^T \nabla c^T Q^{xD} \\
            (Q^{cD})^T \nabla c^T Q^{xD}
	    \end{array}\right]
 \f]
 \verbatim

 C = [ Gc(var_dep,equ_decomp)'   ]
     [ Gh(var_dep,inequ_decomp)' ]
 \endverbatim
 * We can also define a nonbasis matrix \a N for the decomposed constraints as:
 \f[
    N = \left[\begin{array}{c}
            (Q^{cD})^T \nabla c^T Q^{xI} \\
            (Q^{cD})^T \nabla c^T Q^{xI}
	    \end{array}\right]
 \f]
 \verbatim

 N = [ Gc(var_indep,equ_decomp)'   ]
     [ Gh(var_indep,inequ_decomp)' ]
 \endverbatim
 * Given the definitions of \a C and \a N above, we can define the following
 * matrix <tt>D</tt>:
 \f[
   D = - C^{-1} N
 \f]
 \verbatim

  D = -inv(C)*N
 \endverbatim
 * Given this matrix \a D, we can define some other projected sensistivity matrices:
 * <ul>
 * <li> <tt>GcUP = Gc(var_indep,equ_undecomp)   + D'*Gc(var_dep,equ_undecomp)</tt>
 * <li> <tt>GhUP = Gh(var_indep,inequ_undecomp) + D'*Gh(var_dep,inequ_undecomp)</tt>
 * </ul>
 *
 * This interface allows a client to create the basis matrix <tt>C</tt> and optionally
 * the direct sensitivity matrix <tt>D = -inv(C)*N</tt> and the auxiliary projected
 * sensistivity matrices <tt>GcUP</tt> and <tt>GhUP</tt> (shown above).  These matrix
 * objects are independent from \c this \c BasisSystem object or from other \a C, \a D, \c GcUP
 * or \c GhUP objects.  Therefore, a <tt>BasisSystem</tt> object can be thought of
 * as an "Abstract Factory" for basis matrices and auxillary matrices.
 *
 * Note that the purpose of this interface is to abstract client code away from the
 * details of how the matrices \c Gc and \c Gh are represented and implemented and how
 * the basis matrix \a C is formed and implemented.  The complexity of these matrices could
 * vary from simple dense serial matrices all the way up massively parallel matrices using
 * iterative solvers for \a C running on computers with thousands of nodes.
 *
 * <b>Client usage:</b>
 *
 * The matrix objects for the basis matrix \a C and the direct sensitivity matrix
 * \a D are created by the client using the \c MatrixSpace objects returned from
 * \c space_C() and \c space_D().  These methods return smart pointers to these
 * <tt>%MatrixSpace</tt> objects and these objects are ment to have a lifetime that extends
 * up to and beyond the lifetime of the <tt>%BasisSystem</tt> object that created them.
 * Note that the matrix objects returned by this matrix space objects are not to be
 * considered usable until they have passed through \c update_basis().
 *
 * The ranges of the dependent and independent variables, decomposed and undecomposed
 * equality constriants and decomposed and undecomposed inequality constriants are returned
 * by the methods \c var_dep(), \c var_indep(), \c equ_decomp(), \c equ_undecomp(),
 * \c inequ_decomp() and \c inequ_undecomp() respectively.  There are a few obvious assertions
 * for the values that these ranges can take on.  Assuming that \c Gc and \c Gh are non-null
 * when passed to \c update_basis(), the following assertions apply:
 *
 * <A NAME="ranges_assertions"></A>
 * Assertions:<ul>
 * <li> <tt>var_dep().size() == equ_decomp().size() + inequ_decomp().size()</tt>
 * <li> <tt>var_dep().size() + var_indep().size() == Gc.rows() == Gh.rows()</tt>
 * <li> <tt>equ_decomp().size() + equ_undecomp().size() == Gc.cols()</tt>
 * <li> <tt>inequ_decomp().size() + inequ_undecomp().size() == Gh.cols()</tt>
 * </ul>
 *
 * Note that the client should not rely on \c var_dep(), \c var_indep(), \c equ_decomp(),
 * \c equ_undecomp(), \c inequ_decomp() or \c inequ_undecomp() until after the first
 * call to \c update_basis().  This allows a <tt>%BasisSystem</tt> object to adjust itself
 * to accommodate the input matrices \c Gc and \c Gh.  However, in general, a <tt>%BasisSystem</tt>
 * will be ready to deal with only specific \c Gc and/or \c Gh matrices when it is
 * given to a client. The return value <ttt>var_dep() == Range1D::Invalid</tt> is used to
 * flag that \c this <tt>%BasisSystem</tt> object is not initialized for a specific
 * set of matrices yet.  This interface is completely worthless unless \c var_dep() returns
 * some valid range (i.e. a basis matrix exists).
 *
 * The method \c update_basis() is used by the client to update the basis matrix \a C and
 * perhaps the direct sensitivity matrix \a D and it's auxillary projected sensistivity
 * matrices \a V and \a P.  Strictly speaking, it would be possible to
 * form the matrix \a D externally through the <tt>MatrixNonsingular</tt> interface using
 * the returned \a C and an \a N matrix object, but this may not take advantage of any
 * special application specific tricks that
 * can be used to form \a D.  Note that this interface does not return a nonbasis matrix
 * object for \a N.  However, this matrix object will be needed for an implicit \a D matrix
 * object that the client will undoubtably want to create.  Creating such a matrix object
 * is simple given the <tt>MatrixCompositeStd</tt> subclass.  The following code example
 * shows how to create a matrix object for \a N (given the matrices \c Gc and \c Gh input to
 * <tt>bs.update_basis(Gc,Gh,...)</tt> and \c bs):
 \code

 ReferenceCountingPack::ref_count_ptr<const MatrixWithOp>
 create_N(
     const AbstractLinAlgPack::MatrixWithOp*   Gc
     ,const AbstractLinAlgPack::MatrixWithOp*  Gh
     ,const AbstractLinAlgPack::BasisSystem&   bs
     )
 {
     namespace rcp = ReferenceCountingPack;
     rcp::ref_count_ptr<AbstractLinAlgPack::MatrixCompositeStd>
         N = new AbstractLinAlgPack::MatrixCompositeStd(bs.var_dep().size(),bs.var_indep().size());
	 if( Gc && bs.equ_decomp().size() )
         N->add_matrix( 0, 0, 1.0, bs.equ_decomp(), Gc, NULL, BLAS_Cpp::trans, bs.var_indep() );
	 if( Gh && bs.inequ_decomp().size() )
         N->add_matrix( bs.equ_decomp().size(), 0, 1.0, bs.inequ_decomp(), Gh, NULL, BLAS_Cpp::trans, bs.var_indep() );
     N->finish_construction(
         Gc->space_rows().sub_space(bs.equ_decomp())->clone()
         ,Gc->space_cols().sub_space(bs.var_indep())->clone()
         );
     return N;
 }
 \endcode
 * Note that the above nonbasis matrix object \a N returned from the above function depends the matrix objects
 * \c Gc and \c Gh not being modified while \a N is in use.  To make \c N independnet of \c Gc and \c Gh we would
 * of had to clone then (which is not part of the <tt>MatrixWithOp</tt> interface) and may have resulted in a 
 * large allocation of memory. 
 * Given the nonbasis matrix object for \a N returned by the above function, this matrix object could be used
 * to form an explicit \a D matrix object (but perhaps not very efficiently) or be used to implicitly implement
 * matrix vector products with \a D as:
 \verbatim

 op(D)*v = op(-inv(C)*N)*v = -inv(C)*(N*v) or -N'*(inv(C')*v)
 \endverbatim
 *
 * <b>Subclass developer's notes:</b>
 *
 * The default implementation (of the methods that have default implementations) assume
 * that there are no undecomposed equality constriants and no general inequality constriants at all.
 * 
 * ToDo: Finish documentation!
 *
 */
class BasisSystem {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractLinAlgPack::MatrixSpace<MatrixWithOpNonsingular> >    mat_nonsing_space_ptr_t;
	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractLinAlgPack::MatrixSpace<MatrixWithOp> >               mat_space_ptr_t;

	///
	virtual ~BasisSystem() {}

	/** @name MatrixSpace objects */
	//@{

	///
	/** Return a matrix space object for basis <tt>C = [ Gc(var_dep,equ_decomp)';  Gh(var_dep,inequ_decomp)' ]</tt>.
	 */
	virtual const mat_nonsing_space_ptr_t& space_C() const = 0;
	
	///
	/** Return a matrix space object for sensitivity matrix <tt>D = -inv(C)*N</tt>.
	 */
	virtual const mat_space_ptr_t& space_D() const = 0;

	///
	/** Return a matrix space object for auxiliary sensitivity matrix <tt>GcUP = Gc(var_indep,equ_undecomp) + D'*Gc(var_dep,equ_undecomp)</tt>.
	 */
	virtual const mat_space_ptr_t& space_GcUP() const;

	///
	/** Return a matrix space object for auxiliary sensitivity matrix <tt>GhUP = Gh(var_indep,inequ_undecomp) + D'*Gh(var_dep,inequ_undecomp)</tt>.
	 */
	virtual const mat_space_ptr_t& space_GhUP() const;

	//@}

	/** @name Return the ranges for variable and constraint partitioning */
	//@{

	///
	/** Range of dependent (basic) variables.
	 *
	 * If there are no dependent variables then <tt>return.size() == 0</tt>.
	 * This would be a strange case where there was no basis matrix in which
	 * case this whole interface would be worthless.  Therefore, to be useful
	 * <tt>return.size() > 0</tt> must be true.
	 */
	virtual Range1D var_dep() const = 0;
	///
	/** Range of independnet (nonbasic) variables.
	 *
	 * It is possible that the basis matrix may take up all of the degrees of
	 * freedom with <tt>var_dep().size() == Gc->rows()</tt>.  In this case, there
	 * is no nonbasis matrix \a N and no direct sensitivity matrix \a D.
	 * In this case <tt>return.size() == 0</tt>.  In the more general case
	 * however, <tt>return.size() > 0</tt>.
	 */
	virtual Range1D var_indep() const = 0;
	///
	/** Range of decomposed general equality constraints.
	 *
	 * If there are no decomposed general equality constriants then
	 * <tt>return.size() == 0</tt>.  Otherwise, <tt>return.size() > 0</tt>.
	 *
	 * The default implementation return <tt>Range1D(1,this->var_dep().size())</tt>
	 */
	virtual Range1D equ_decomp() const;
	///
	/** Range of undecomposed general equality constriants.
	 *
	 * If there are no undecomposed equality constriants then
	 * <tt>return.size() == 0</tt>.  Otherwise, <tt>return.size() > 0</tt>.
	 *
	 * The default implementation return <tt>Range1D::Invalid</tt>
	 */
	virtual Range1D equ_undecomp() const;
	///
	/** Range for decomposed general inequality constraints.
	 *
	 * If there are no decomposed general inequality constriants then
	 * <tt>return.size() == 0</tt>.  Otherwise, <tt>return.size() > 0</tt>.
	 *
	 * The default implementation return <tt>Range1D::Invalid</tt>
	 */
	virtual Range1D inequ_decomp() const;
	///
	/** Range for undecomposed general inequality constraints.
	 *
	 * If there are no undecomposed general inequality constriants then
	 * <tt>return.size() == 0</tt>.  Otherwise, <tt>return.size() > 0</tt>.
	 *
	 * The default implementation return <tt>Range1D::Invalid</tt>
	 */
	virtual Range1D inequ_undecomp() const;

	//@}

	/** @name Update matrices */
	//@{

	///
	/** Update a basis and posssibly the direct sensitivity matrix for a 
	 * set of Jacobian matrices.
	 *
	 * @param  Gc    [in] Jacobian of the equality constriants.
	 * @param  Gh    [in] Jacobian of the inequality constraints.
	 * @param  C     [out] Basis matrix.  If <tt>C == NULL</tt> on input, then this
	 *               quantity is not updated.  If <tt>C != NULL</tt> then this must
	 *               have been created by <tt>this->space_C()->create_member()</tt>.
	 *               This basis matrix object must be independent of the input
	 *               matrices \c Gc and/or \c Gh.  Therefore, it must be legal to
	 *               destroy \c Gc and/or \c Gh without affecting the behavior of
	 *               the basis matrix object \c C.
	 * @param  D     [out] Direct sensitivity matrix <tt>D = -inv(C)*N</tt>.  If
	 *               <tt>D == NULL</tt> on input then this quantity is not updated.
	 *               If <tt>D != NULL</tt> then this must have been created by
	 *               <tt>this->space_D()->create_member()</tt>.  This matrix object
	 *               is meaningless if <tt>this->var_indep() == Range1D::Invalid</tt>
	 *               on return.
	 *               This matrix object must be independent of the input
	 *               matrices \c Gc and/or \c Gh.  Therefore, it must be legal to
	 *               destroy \c Gc and/or \c Gh without affecting the behavior of
	 *               the direct sensitivity matrix object \c D.
	 * @param  GcUP  [out] Auxiliary sensistivity matrix
	 *               <tt>GcUP = Gc(var_indep,equ_undecomp) + D'*Gc(var_dep,equ_undecomp)</tt>.
	 *               If <tt>GcUP == NULL</tt> on input then this quantity is not updated.
	 *               If <tt>GcUP != NULL</tt> then this must have been created by
	 *               <tt>this->space_GcUP()->create_member()</tt>.  This matrix object
	 *               is meaningless if <tt>this->var_indep() == Range1D::Invalid</tt>
	 *               on return.
	 *               This matrix object must be independent of the input
	 *               matrices \c Gc and/or \c Gh.  Therefore, it must be legal to
	 *               destroy \c Gc and/or \c Gh without affecting the behavior of
	 *               the direct sensitivity matrix object \c D.
	 * @param  GhUP  [out] Auxiliary sensistivity matrix
	 *               <tt>GhUP = Gh(var_indep,inequ_undecomp) + D'*Gh(var_dep,inequ_undecomp)</tt>.
	 *               If <tt>GhUP == NULL</tt> on input then this quantity is not updated.
	 *               If <tt>GhUP != NULL</tt> then this must have been created by
	 *               <tt>this->space_GhUP()->create_member()</tt>.  This matrix object
	 *               is meaningless if <tt>this->var_indep() == Range1D::Invalid</tt>
	 *               on return.
	 *               This matrix object must be independent of the input
	 *               matrices \c Gc and/or \c Gh.  Therefore, it must be legal to
	 *               destroy \c Gc and/or \c Gh without affecting the behavior of
	 *               the direct sensitivity matrix object \c D.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>Gc != NULL || Gh != NULL</tt>
	 * <li> [<tt>Gc != NULL && Gh != NULL</tt>]
	 *      <tt>Gc->space_cols().is_compatible(Gh->space_cols()) == true</tt>
	 * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space_cols().sub_space(var_dep()).get() != NULL</tt>
	 * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space_cols().sub_space(var_indep()).get() != NULL</tt>
	 * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space_rows().sub_space(equ_decomp()).get() != NULL</tt>
	 * <li> [<tt>Gc != NULL && equ_decomp().size() < Gc.cols()</tt>]
	 *      <tt>Gc.space_rows().sub_space(equ_undecomp()).get() != NULL</tt>
	 * <li> [<tt>Gh != NULL</tt>]  <tt>Gh->space_rows().sub_space(inequ_decomp()).get() != NULL</tt>
	 * <li> [<tt>Gh != NULL && inequ_decomp().size() < Gh->cols()</tt>]
	 *      <tt>Gh.space_rows().sub_space(inequ_undecomp()).get() != NULL</tt>
	 * <li> [<tt>Gc == NULL</tt>] <tt>GcUP == NULL</tt>
	 * <li> [<tt>Gh == NULL</tt>] <tt>GhUP == NULL</tt>
	 * <li> <tt>C != NULL || D != NULL || GcUP != NULL || GhUP != NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->var_dep() != Range1D::Invalid && !this->var_dep().full_range()</tt>
	 * <li> [<tt>C != NULL && Gc != NULL && equ_decomp().size() > 0</tt>]
	 *      <tt>C->space_cols().sub_space(equ_decomp())->is_compatible(Gc->space_rows().sub_space(equ_decomp()))
	 *      && C->space_rows().is_compatible(Gc->space_cols().sub_space(var_dep()))</tt>
	 * <li> [<tt>C != NULL && Gh != NULL && inequ_decomp().size() > 0</tt>]
	 *      <tt>C->space_cols().sub_space(equ_decomp().size()+inequ_decomp())->is_compatible(Gh->space_rows().sub_space(inequ_decomp()))
	 *      && C->space_rows().is_compatible(Gh->space_cols().sub_space(var_dep()))</tt>
	 * <li> [<tt>D != NULL && Gc != NULL && var_indep().size() > 0 && equ_decomp().size() > 0</tt>]
	 *      <tt>D->space_cols().sub_space(equ_decomp())->is_compatible(Gc->space_rows().sub_space(equ_decomp()))
	 *      && D->space_rows().is_compatible(Gc->space_cols().sub_space(var_indep()))</tt>
	 * <li> [<tt>D != NULL && Gh != NULL && var_indep().size() > 0 && inequ_decomp().size() > 0</tt>]
	 *      <tt>D->space_cols().sub_space(equ_decomp().size()+inequ_decomp())->is_compatible(Gh->space_rows().sub_space(inequ_decomp()))
	 *      && D->space_rows().is_compatible(Gh->space_cols().sub_space(var_indep()))</tt>
	 * <li> [<tt>GcUP != NULL && var_indep().size() > 0 && equ_undecomp().size() > 0</tt>]
	 *      <tt>GcUP->space_cols()->is_compatible(Gc->space_cols().sub_space(var_indep()))
	 *      && GcUP->space_rows()->is_compatible(Gc->space_rows().sub_space(equ_undecomp()))</tt>
	 * <li> [<tt>GhUP != NULL && var_indep().size() > 0 && inequ_undecomp().size() > 0</tt>]
	 *      <tt>GhUP->space_cols()->is_compatible(Gh->space_cols().sub_space(var_indep()))
	 *      && GhUP->space_rows()->is_compatible(Gh->space_rows().sub_space(inequ_undecomp()))</tt>
	 * </ul>
	 *
	 * This method is declared non-const because it may change what is returned by the partitioning functions such as
	 * \c var_dep().
	 */
	virtual void update_basis(
		const MatrixWithOp*         Gc
		,const MatrixWithOp*        Gh
		,MatrixWithOpNonsingular*   C
		,MatrixWithOp*              D
		,MatrixWithOp*              GcUP
		,MatrixWithOp*              GhUP
		) = 0;

	//@}

private:

	static const mat_space_ptr_t   mat_space_ptr_null_;

}; // end class BasisSystem

}; // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_SYSTEM_H
