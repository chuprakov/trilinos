// ////////////////////////////////////////////////////////////////////////////////////
// rSQPState.h
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

#ifndef RSQP_STATE_H
#define RSQP_STATE_H

#include "ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/include/IterQuantityAccess.h"
#include "GeneralIterationPack/include/AlgorithmState.h"
#include "GeneralIterationPack/include/cast_iq.h"
#include "GeneralIterationPack/include/IterQuantityAccessContiguous.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
//#include "ConstrainedOptimizationPack/include/DecompositionSystem.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"
//#include "LinAlgPack/include/IVector.h"
#include "StandardCompositionMacros.h"
#include "StandardMemberCompositionMacros.h"
#include "Range1D.h"

namespace ReducedSpaceSQPPack {

/** \defgroup rSQPState_IQ_names_grp rSQPState iteration quantities names.
 *
 * ToDo: Finish documentation!
 */
//@{

// Iteration Info
extern const std::string num_basis_name;
// NLP Problem Info 
extern const std::string x_name;
extern const std::string f_name;
extern const std::string Gf_name;
extern const std::string HL_name;
extern const std::string c_name;
extern const std::string h_name;
extern const std::string Gc_name;
extern const std::string Gh_name;
// Constraint Gradient Null Space / Range Space Decomposition Info
extern const std::string Y_name;
extern const std::string Z_name;
extern const std::string R_name;
extern const std::string Uy_name;
extern const std::string Uz_name;
extern const std::string Vy_name;
extern const std::string Vz_name;
// Search Direction Info
extern const std::string py_name;
extern const std::string Ypy_name;
extern const std::string pz_name;
extern const std::string Zpz_name;
extern const std::string d_name;
// Reduced QP Subproblem Info
extern const std::string rGf_name;
extern const std::string rHL_name;
extern const std::string w_name;
extern const std::string zeta_name;
extern const std::string qp_grad_name;
extern const std::string eta_name;
// Global Convergence Info
extern const std::string alpha_name;
extern const std::string merit_func_nlp_name;
extern const std::string mu_name;
extern const std::string phi_name;
// KKT Info
extern const std::string opt_kkt_err_name;
extern const std::string feas_kkt_err_name;
extern const std::string GL_name;
extern const std::string rGL_name;
extern const std::string lambda_name;
extern const std::string lambdaI_name;
extern const std::string nu_name;

//@}

/** \defgroup rSQPState_Macros_grp Macros for adding IterQuantity objects to rSQPState.
 *
 * These macros make it easy to add and remove iteration quantities of any type
 * to the state class.  These macros can even be used by subclasses of rSQPState
 * (only any <tt>GeneralIterationPack::AlgorithmState</tt> subclass) to add
 * iteration quantities.  Since scalars and vectors are so pervasive they have
 * there own special macros that take care of actually instantiating the
 * iteration quantities.
 */
//@{

///
/** Add class declarations for an arbitrary iteration quantity
 */
#define RSQP_STATE_IQ_DECL(TYPE,NAME)                                     \
	virtual IterQuantityAccess<TYPE>&       NAME();                       \
	virtual const IterQuantityAccess<TYPE>& NAME() const;                 \
private:                                                                  \
	iq_id_encap NAME ## _iq_id_;                                          \
public:

///
/** Add class declarations for an index (i.e. index_type) iteration quantity
 */
#define RSQP_STATE_INDEX_IQ_DECL(NAME)                                    \
    RSQP_STATE_IQ_DECL(index_type,NAME)                                   \

///
/** Add class declarations for a scalar (i.e. value_type) iteration quantity
 */
#define RSQP_STATE_SCALAR_IQ_DECL(NAME)                                   \
    RSQP_STATE_IQ_DECL(value_type,NAME)                                   \

///
/** Add class declarations for a VectorWithOpMutable iteration quantity.
 */
#define RSQP_STATE_VECTOR_IQ_DECL(NAME)                                   \
    RSQP_STATE_IQ_DECL(VectorWithOpMutable,NAME)                          \

///
/** Add class definitions for an arbitrary iteration quantity.
 *
 * This implementation can not instantate the underlying iteration quantity
 * so it is the responsibility of some client to do this.  This implementation
 * just initializes the iq_id for the iteration quantity on the fly and
 * then casts the iteration quantity. 
 */
#define RSQP_STATE_IQ_DEF(CLASS,TYPE,NAME,NAME_STR)                       \
IterQuantityAccess<TYPE>&                                                 \
CLASS::NAME()                                                             \
{                                                                         \
	update_iq_id( NAME_STR, &NAME ## _iq_id_ );                           \
	return GeneralIterationPack::cast_iq<TYPE>(                           \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<TYPE>&                                           \
CLASS::NAME() const                                                       \
{                                                                         \
	return const_cast<rSQPState*>(this)->NAME();                          \
}

///
/** Add class definitions for a index_type iteration quantity.
 *
 * This implementation will instantiate the IterQuantity object on the fly.
 */
#define RSQP_STATE_INDEX_IQ_DEF(CLASS,NAME,NAME_STR)                      \
IterQuantityAccess<index_type>&                                           \
CLASS::NAME()                                                             \
{                                                                         \
	update_index_type_iq_id( NAME_STR, &NAME ## _iq_id_ );                \
	return GeneralIterationPack::cast_iq<index_type>(                     \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<index_type>&                                     \
CLASS::NAME() const                                                       \
{                                                                         \
	return const_cast<rSQPState*>(this)->NAME();                          \
}

///
/** Add class definitions for a value_type iteration quantity.
 *
 * This implementation will instantiate the IterQuantity object on the fly.
 */
#define RSQP_STATE_SCALAR_IQ_DEF(CLASS,NAME,NAME_STR)                     \
IterQuantityAccess<value_type>&                                           \
CLASS::NAME()                                                             \
{                                                                         \
	update_value_type_iq_id( NAME_STR, &NAME ## _iq_id_ );                \
	return GeneralIterationPack::cast_iq<value_type>(                     \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<value_type>&                                     \
CLASS::NAME() const                                                       \
{                                                                         \
	return const_cast<rSQPState*>(this)->NAME();                          \
}

///
/** Add class definitions for a VectorWithOpMutable iteration quantity.
 *
 * This implementation will instantiate the IterQuantity object on the fly
 * given the VectorSpace (VEC_SPC).  Note that this VEC_SPC can be any
 * code that will returns a smart pointer to a AbstractFactory<VectorWithOpMutable>
 * object from within the class body.  It is best if VEC_SPC is some function
 * that is called on *this for the maximum safety and to avoid strage
 * behavior.
 */
#define RSQP_STATE_VECTOR_IQ_DEF(CLASS,NAME,NAME_STR,VEC_SPC)             \
IterQuantityAccess<VectorWithOpMutable>&                                  \
CLASS::NAME()                                                             \
{                                                                         \
    update_vector_iq_id( NAME_STR, VEC_SPC, &NAME ## _iq_id_ );           \
	return GeneralIterationPack::cast_iq<VectorWithOpMutable>(            \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<VectorWithOpMutable>&                            \
CLASS::NAME() const                                                       \
{                                                                         \
	return const_cast<rSQPState*>(this)->NAME();                          \
}

//@}

///
/** Reduced space SQP state encapsulation interface.
 *
 * This in an interface to a set of data specific to a reduced space SQP algorithms.
 * The iteration quantites are abstracted within <tt>IterQuantityAccess<></tt> objects.
 * A set of boilerplate macros are used to add the necessary declarations and
 * implemetations of these iteration quantity access functions.  As shown by
 * these macros the access methods are declared virtual so that subclasses can
 * override these methods.  Otherwise, much of these could have been declared
 * inline.
 *
 * The implementation defined in this class uses <tt>IterQuantityAccessContiguous<></tt>
 * for iteration quantities of type <tt>index_type</tt>, <tt>value_type</tt> and
 * <tt>VectorWithOpMutable</tt> with a default of one storage location.  The default
 * implementation is able to create the <tt>VectorWithOpMutable</tt> iteration quantities
 * by using <tt>VectorSpace</tt> objects that the client sets \c this up with.
 *
 * For all other types of iteration quantities (i.e. <tt>MatrixWithOp</tt> etc.) the
 * client is responsible for setting the iteration quantity object of type
 * <tt>IterQuantityAccess<></tt>.  The client can also change the type of class
 * used for any iteration quantity by simply calling
 * <tt>AlgorithmState::set_iter_quant(...)</tt>.
 *
 * The number of storage locations for any iteration quantity of type
 * <tt>IterQuantityAccessContiguous<></tt> can be changed by fetching
 * the iteration quantity using the access methods defined here and then
 * using <tt>dynamic_cast<></tt> and calling the
 * <tt>IterQuantityAccessContiguous<>::resize(...)</tt> method.
 *
 * Note that the underlying <tt>AlgorithmState</tt> object will not know
 * about the iteration quantity objects with default implementations
 * until the access functions have been called at least once.
 *
 * ToDo: Finish documentation.
 */
class rSQPState
	: public GeneralIterationPack::AlgorithmState // doxygen needs full path
{
public:

	/** @name Public Types */
	//@{

	/// Thrown if an iteration quantity is of an invalid type.
	class InvalidType : public std::logic_error
	{public: InvalidType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

protected:

	// /////////////////////////////
	// Protected types.

	///
	struct iq_id_encap {
		iq_id_encap() : iq_id(DOES_NOT_EXIST) {}
		iq_id_type iq_id;
	};

public:

	/** @name Constructors/initializers */
	//@{

	/// <<std comp>> members for the VectorSpace of x
	STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpace, space_x )
	/// <<std comp>> members for the VectorSpace of c
	STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpace, space_c )
	/// <<std comp>> members for the VectorSpace of h
	STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpace, space_h )
	/// <<std comp>> members for the VectorSpace of py
	STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpace, space_py )
	/// <<std comp>> members for the VectorSpace of pz
	STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpace, space_pz )

	///
	/** Construct
	 *
	 * Initializes num_basis() == 0
	 */
	rSQPState(
		const space_x_ptr_t&   space_x   = NULL
		,const space_x_ptr_t&  space_c   = NULL
		,const space_x_ptr_t&  space_h   = NULL
		,const space_x_ptr_t&  space_py  = NULL
		,const space_x_ptr_t&  space_pz  = NULL
		);

	///
	virtual ~rSQPState() {}

	//@}

	/** @name Iteration Info */
	//@{

    /// num_basis: Counts basis changes durring the algorithm
	RSQP_STATE_INDEX_IQ_DECL(num_basis)
	
	//@}

	/** @name NLP Problem Info */
	//@{

	/// x:  The current NLP point
	RSQP_STATE_VECTOR_IQ_DECL(x)
	/// f:  Objective function value
	RSQP_STATE_SCALAR_IQ_DECL(f)
	/// Gf:  Gradient of the objective function sorted according to current basis selection ( n x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(Gf)
	/// HL:  Hessian of the Lagrangian ( n x n 
	RSQP_STATE_IQ_DECL(MatrixSymWithOp,HL)
	/// c:  Vector of general nonlinear equality constraints ( m x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(c)
	/// h:  Vector of general nonlinear inequality constraints ( mI x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(h)
	/// Gc:  Gradient of equality constraints ('c') matrix ( n x m )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Gc)
	/// Gh:  Gradient of inequality constraints ('h') matrix ( n x mI )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Gh)

	//@}

	/** @name Constraint Gradient Null Space / Range Space Decomposition Info */
	//@{

	/// Y:  Range space matrix for Gc ([Y  Z] is non-singular) ( n x r )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Y)
	/// Z:  Null space matrix for Gc(con_decomp)' (Gc(con_decomp)' * Z) ( n x (n-r) )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Z)
	/// R:  Represents the nonsingular matrix Gc(con_decomp)' * Y ( r x r )
	RSQP_STATE_IQ_DECL(MatrixWithOpNonsingular,R)
	/// Uy:  Represents Gc(con_undecomp)' * Y ( (m-r) x r )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Uy)
	/// Uz:  Represents Gc(con_undecomp)' * Z ( (m-r) x (m-r) )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Uz)
	/// Vy:  Represents Gh' * Y ( mI x r )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Vy)
	/// Vz:  Represents Gh' * Z ( mI x (m-r) )
	RSQP_STATE_IQ_DECL(MatrixWithOp,Vz)

	//@}

	/** @name Search Direction Info */
	//@{

	/// py:  Range space (dependent) QP solution component ( space_py, m x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(py)
	/// Ypy:  Range space (dependent) contribution to search direction (Ypy = Y * py) ( n x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(Ypy)
	/// pz:  Null space (independent) QP solution component ( (n-m) x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(pz)
	/// Zpz:  Null space (independent) contribution to the search direction (Zpz = Z * pz) ( n x 1)
	RSQP_STATE_VECTOR_IQ_DECL(Zpz)
	/// d:  Search direction (d = Zpz + Ypy) ( n x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(d)

	//@}

	/** @name QP Subproblem Info */
	//@{

	/// rGf:  Reduced gradient of the objective function ( (n-r) x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(rGf)
	/// rHL:  Reduced Hessian of the Lagrangian function ( (n-r) x (n-r) )
	RSQP_STATE_IQ_DECL(MatrixSymWithOp,rHL)
	/// w:  QP gradient crossterm correction (Z' * HL * Y * py) ( (n-r) x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(w)
	/// zeta:  QP crossterm dampening parameter [0, 1]
	RSQP_STATE_SCALAR_IQ_DECL(zeta)
	/// qp_grad:  QP gradient (qp_grad = rGf + zeta * w) ( (n-m) x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(qp_grad)
	/// eta:  QP relaxation parameter [0, 1]
	RSQP_STATE_SCALAR_IQ_DECL(eta)

	//@}

	/** @name Global Convergence Info */
	//@{

	/// alpha:  Line seach parameter
	RSQP_STATE_SCALAR_IQ_DECL(alpha)
	/// merit_func_nlp: Primary merit function for the NLP
	RSQP_STATE_IQ_DECL(MeritFuncNLP,merit_func_nlp)
	/// mu:  Merit function penalty parameter
	RSQP_STATE_SCALAR_IQ_DECL(mu)
	/// phi:  Merit function value
	RSQP_STATE_SCALAR_IQ_DECL(phi)

	//@}

	/** @name KKT Info */
	//@{

	/// Scaled KKT error for optimality ||rGL||
	RSQP_STATE_SCALAR_IQ_DECL(opt_kkt_err)
	/// Scaled KKT error for feasibility ||c|| and ||hl <= h <= hu||
	RSQP_STATE_SCALAR_IQ_DECL(feas_kkt_err)
	/// GL:  Gradient of the Lagrangian ( n x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(GL)
	/// rGL:  Reduced gradient of the Lagrangian ( (n-m) x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(rGL)
	/// lambda:  Lagrange multipliers for the equality constraints 'c' ( m x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(lambda)
	/// lambdaI:  Lagrange multipliers for the ineequality constraints 'h' ( mI x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(lambdaI)
	/// nu:  Difference between Lagrange multipiers for the upper and lower bounds ( n x 1 )
	RSQP_STATE_VECTOR_IQ_DECL(nu)

	//@}

	/** @name Basis Pivot Info */
	//@{

	//@}

protected:

	// /////////////////////////////
	// Protected member functions

	// These implementations are used to avoid code blot and help in debugging
	// (can't debug macros very well).

	///
	void update_iq_id(
		const std::string&                iq_name
		,iq_id_encap*                     iq_id
		) const;
	///
	void update_index_type_iq_id(
		const std::string&                iq_name
		,iq_id_encap*                     iq_id
		);
	///
	void update_value_type_iq_id(
		const std::string&                iq_name
		,iq_id_encap*                     iq_id
		);
	///
	void update_vector_iq_id(
		const std::string&                iq_name
		,const VectorSpace::space_ptr_t&  vec_space
		,iq_id_encap*                     iq_id
		);

private:

	// ////////////////////////////
	// Private member functions.

	// not defined and not to be called
	rSQPState(const rSQPState&);
	rSQPState& operator=(const rSQPState&);

};	// end class rSQPState

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_STATE_H
