// //////////////////////////////////////////////////////////////////////////////////////
// RTOpCppDecl.hpp
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

#ifndef RTOPPACK_RTOP_CPP_DECL_HPP
#define RTOPPACK_RTOP_CPP_DECL_HPP

#include <assert.h>

#include "RTOpPackTypes.hpp"
#include "RTOp.h"

namespace RTOpPack {

/** \defgroup RTOpCpp_grp Templated interfaces for generalized vector
 * reduction/transformation operators in C++.
 */
//@{

extern "C" {
///
typedef int (*op_create_func_t)( const RTOp_obj_type_vtbl_t* dummy1, const void* dummy2, void** op );
///
typedef int (*op_free_func_t)( const RTOp_obj_type_vtbl_t* dummy1, const void* dummy2, void** op );
} // extern "C"

///
/** C++ interface to vector reduction/transformation operators {abstract}.
 *
 * This class provides a superior interface to reduction/transformation operators
 * to C++ clients than is provided by the access functions for the C struct
 * <tt>RTOp_RTOp</tt>.  In addition, this class can have direct C++ subclasses.  This
 * allows new reduction/transformation operator classes to be developed much easier
 * by those who know C++.  However, this interface is not entirely type safe because
 * of the opaque <tt>RTOp_ReductTarget</tt> objects still included in the interface.
 * Unfortunately, this can not be avoided however because of the the need for compatibility
 * with C <tt>RTOp_RTOp</tt> operator object classes and C clients.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class RTOpT {
public:

	///
	virtual ~RTOpT() {}
	///
	/** Return a function pointer to a function that will create an
	 * RTOp object.
	 *
	 * The operator object pointed to by <tt>op</tt> returned from
	 * <tt>this->get_op_create_func()(NULL,&op)</tt> can be deleted
	 * only by calling <tt>this->get_op_free_func()(NULL,NULL,&op)</tt>.
	 * Do not attempt to call <tt>delete</tt> on the pointer returned
	 * from this function.  The prototype for this function may
	 * seem strange, but this is required so that it will be
	 * compatible with C implementations (see RTOp_obj_type_vtbl_t).
	 *
	 * Factory objects will be considered equivalent if the
	 * pointers returned from this function are equal.  In fact,
	 * comparison of the return value from this function is used
	 * as the test to see if operator factories are equivalent.
	 */
	virtual op_create_func_t get_op_create_func() const;
	///
	/** Return a function pointer to a function that will free an
	 * RTOp object.
	 *
	 * The operator object pointed to by <tt>op</tt> returned from
	 * <tt>this->get_op_create_func()(NULL,NULL,&op)</tt> must be deleted
	 * by calling <tt>this->get_op_free_func()(NULL,NULL,&op)</tt>.  Again, the
	 * prototype for this function may seem strange but it is required
	 * to be this way for C compatibility.
	 */
	virtual op_free_func_t get_op_free_func() const;
	///
	/** Return the number of entries of each type of basic data type in the externalized state
	 * for the operator object.
	 *
	 * The default implementation returns zeros for <tt>*num_values</tt>, <tt>*num_indexes</tt> and <tt>*num_chars</tt>
	 * (i.e. the default reduction/transformation operator has no state).
	 */
	virtual void get_op_type_num_entries(
		int*  num_values
		,int* num_indexes
		,int* num_chars
		) const;
	///
	/** Extract the state of the object in a portable format.
	 *
	 * This method allows the state of a reduction/transformation operator to be
	 * externalized so that it can be passed over a heterogeneous network of computers.
	 *
	 * The default implemenation does nothing (i.e. the default reduction/transformation
	 * operator has no state).
	 *
	 * @param  num_values
	 *              [in] Value returned from <tt>this->\ref RTOp::get_op_type_num_entries "get_op_type_num_entries(...)"</tt>.
	 * @param  value_data
	 *              [out] Array (size <tt>num_values</tt>) of floating point data.
	 * @param  num_indexes
	 *              [in] Value returned from <tt>this->get_op_type_num_entries(...)</tt>.
	 * @param  index_data
	 *              [out] Array (size <tt>num_indexes</tt>) of integral data.
	 * @param  num_chars
	 *              [in] Value returned from <tt>this->get_op_type_num_entries(...)</tt>.
	 * @param  char_data
	 *              [out] Array (size <tt>num_chars</tt>) of character data.
	 */
	virtual void extract_op_state(
		int               num_values
		,Scalar           value_data[]
		,int              num_indexes
		,RTOp_index_type  index_data[]
		,int              num_chars
		,RTOp_char_type   char_data[]
		) const;
	///
	/** Load the state of the object from a portable format.
	 *
	 * This method allows the state of the operator object to be set given an the externalized
	 * state as extracted using <tt>\ref extract_op_state "extract_op_state(...)"</tt>
	 * called on a compatible operator object (possibly on a different heterogeneous computer).
	 *
	 * The default implementation does nothing (i.e. the default reduction/transformation
	 * operator has no state).
	 *
	 * @param  num_values
	 *              [in] Value returned from <tt>this->\ref get_op_type_num_entries "get_op_type_num_entries(...)"</tt>.
	 * @param  value_data
	 *              [out] Array (size <tt>num_values</tt>) of floating point data.
	 * @param  num_indexes
	 *              [in] Value returned from <tt>this->get_op_type_num_entries(...)</tt>.
	 * @param  index_data
	 *              [out] Array (size <tt>num_indexes</tt>) of integral data.
	 * @param  num_chars
	 *              [in] Value returned from <tt>this->get_op_type_num_entries(...)</tt>.
	 * @param  char_data
	 *              [out] Array (size <tt>num_chars</tt>) of character data.
	 */
	virtual void load_op_state(
		int                       num_values
		,const Scalar             value_data[]
		,int                      num_indexes
		,const RTOp_index_type    index_data[]
		,int                      num_chars
		,const RTOp_char_type     char_data[]
		);
	///
	/** Get the number of entries of each basic data type in the externalized state for
	 * a reduction object.
	 *
	 * Note that a specific reduction object is not passed in as an argument.  This is because
	 * the structure of a reduction object is completely determined by its associated operator object and
	 * this structure can not change as a result of a reduction operation (this is needed to simplyfy
	 * difficult networking code).
	 *
	 * The default implementation returns zeros for <tt>*num_values</tt>, <tt>*num_indexes</tt> and <tt>*num_chars</tt>
	 * (i.e. there is no reduction operation performed).
	 */
	virtual void get_reduct_type_num_entries(
		int*   num_values
		,int*  num_indexes
		,int*  num_chars
		) const;
	///
	/** The default implementional just calls <tt>reduct_obj_create_raw()</tt>.
	 *
	 * This should be the only implementation needed by all subclasses.  This method
	 * simply sets up the <tt>reduct_obj</tt> encapulsation object so that it appears that
	 * it is a fully functional object.  In reality, <tt>reduct_obj</tt> contains a pointer
	 * to <tt>this</tt> and is not an independent object in that sense.  However, this method
	 * is no less useful but it could as well have been a nonmember function.  On the
	 * otherhand, making in a virtual member function allows an <tt>RTOp</tt> subclass to
	 * override it and monitor when the function is called for whatever reason.
	 */
	virtual void reduct_obj_create( ReductTargetT<Scalar>* reduct_obj ) const;
	///
	/** Apply the reduction/transformation operator to a set of sub-vectors.
	 *
	 * <tt>op(sub_vecs[],targ_sub_vecs[]),reduct_obj) -> targ_sub_vecs[],reduct_obj</tt>.
	 *
	 * This is the bread and butter of the whole design.  Through this method, a
	 * vector implementation applies a reduction/transformation operator to a
	 * set of sub-vectors.
	 *
	 *	@param	num_vecs
	 *				[in] Number of non-mutable sub-vectors in <tt>sub_vec[]</tt>.
	 *	@param	sub_vecs
	 *				[in] Array (length <tt>num_vecs</tt>) of non-mutable vectors to apply the
	 *				operator over.  The ordering of these sub-vectors
	 *				<tt>sub_vecs[k], for k = 0...num_vecs-1</tt>, is significant to the <tt>op</tt> object.
	 *             If <tt>num_vecs == 0</tt> then <tt>sub_vecs</tt> can be <tt>NULL</tt>.
	 *	@param	num_targ_vecs
	 *				[in] Number of mutable sub-vectors in <tt>targ_sub_vec[]</tt>.
	 *	@param	targ_sub_vecs
	 *				[in] Array (length <tt>num_targ_vecs</tt>) of mutable vectors to apply the
	 *				operator over and be mutated.  The ordering of these sub-vectors
	 *				<tt>targ_sub_vecs[k], for k = 0...num_targ_vecs-1</tt>, is significant to
	 *             the <tt>op</tt> object.  If <tt>num_targ_vecs == 0</tt> then <tt>targ_sub_vecs</tt> can be <tt>NULL</tt>.
	 *	@param	reduct_obj
	 *				[in/out] This reduction object must have been created by the
	 *              <tt>this->\ref reduct_obj_create_raw "reduct_obj_create_raw(&reduct_obj)"</tt>
	 *              method and it may have already passed through one or more other
	 *				reduction operations (accumulating the reductions along the way).
	 *              If <tt>this->\ref get_reduct_type_num_entries "get_reduct_type_num_entries(...)"</tt>
	 *              returns <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt>,
	 *              then <tt>reduct_obj</tt> should be set to <tt>RTOp_REDUCT_OBJ_NULL</tt>
	 *              and no reduction will be performed.
	 *
	 * Preconditions:<ul>
	 *	<li> <tt>num_vecs > 0 || num_targ_vecs > 0</tt>
	 *	<li> <tt>num_vecs > 0 || sub_vecs == NULL</tt>
	 *	<li> <tt>num_targ_vecs > 0 || targ_sub_vecs == NULL</tt>
	 *	<li> [<tt>num_vecs > 0</tt>] <tt>global_offset == sub_vecs[k].global_offset</tt>
	 *        , for <tt>k = 0,...,num_vecs-1</tt>
	 *	<li> [<tt>num_targ_vecs > 0</tt>] <tt>global_offset == targ_sub_vecs[k].global_offset</tt>
	 *        , for <tt>k = 0,...,num_targ_vecs-1</tt>
	 *	<li> [<tt>num_vecs > 0</tt>] <tt>sub_dim == sub_vecs[k].sub_dim</tt>
	 *       , for <tt>k = 0,...,num_vecs-1</tt>
	 *	<li> [<tt>num_targ_vecs > 0</tt>] <tt>sub_dim == targ_sub_vecs[k].sub_dim</tt>
	 *       , for <tt>k = 0,...,num_targ_vecs-1</tt>
	 *	</ul>
	 *
	 * If <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt> then the reduction operation will be accumlated as:
	 \verbatim

	 op(op(sub_vecs[],targ_sub_vecs[]),reduct_obj) -> reduct_obj
	 \endverbatim
	 * By allowing an in/out <tt>reduct_obj</tt> and an accumulation
	 * of the reduction, the maximum reuse of memory is achieved.
	 * If <tt>this->\ref reduct_obj_create_raw "reduct_obj_create_raw(&reduct_obj)"</tt>
	 * or <tt>this->\ref reduct_obj_reinit "reduct_obj_reinit(reduct_obj)"</tt> was called
	 * immediately before this function, then <tt>reduct_obj</tt> will of course only contain
	 * the reduction from this operation.
	 * 
	 *
	 * If <tt>num_vecs</tt> is incompatible with the underlying operator object then
	 * InvalidNumVecs is thrown.  If the sub-vectors are not compatible
	 * (i.e. <tt>global_offset</tt> and/or <tt>sub_dim</tt> not the same) then
	 * IncompatibleVecs is thrown.
	 */
	virtual void apply_op(
		const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
		,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
		,RTOp_ReductTarget reduct_obj
		) const = 0;
	///
	/** Reduce intermediate reduction target objects.
	 *
	 * The default implementation does not do anything (i.e. there is no reduction operation performed).
	 */
	virtual void reduce_reduct_objs(
		RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj
		) const;
	///
	/** Set a pointer to an MPI compatible reduction function.
	 *
	 * @param  reduct_op_func_ptr  [out] Pointer to an MPI compatible user defined
	 *                             reduction function.
	 *
	 * This function may set <tt>*reduct_op_func_ptr == NULL</tt> which is
	 * perfectly okay.  In this case, the vector implementation must still
	 * perform the reduction operation but without the explicit help of MPI.
	 *
	 * The default implmentation sets <tt>reduct_op_func_ptr==NULL</tt>.
	 * (i.e. Perform the reductions yourself).
	 *
	 */
    virtual void get_reduct_op( RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr ) const;

	/** @name Not for the casual C++ client.
	 *
	 * Use a <tt>ReductTargetT</tt> object initialized with <tt>this->reduct_obj_create()</tt>
	 * instead.
	 */
	//@{

	///
	/** Create an opaque reduction objet.
	 *
	 * The default implementation does nothing (i.e. there is no reduction operation performed).
	 */
	virtual void reduct_obj_create_raw( RTOp_ReductTarget* reduct_obj ) const;
	///
	/** Reinitialize an already created reduction object.
	 *
	 * The default implementation does nothing (i.e. there is no reduction operation performed).
	 *
	 * @param  reduct_obj  [in/out] Reduction object is reinitalized on output.
	 */
	virtual void reduct_obj_reinit( RTOp_ReductTarget reduct_obj ) const;
	///
	/** Free a previously created reduction object (i.e. there is no reduction operation performed).
	 *
	 * The default implementation does nothing (i.e. there is no reduction operation performed).
	 *
	 * @param  reduct_obj  [in/out] <tt>*reduct_obj</tt> is freed and then
	 *                     <tt>*reduct_obj == RTOp_REDUCT_OBJ_NULL</tt> is set on
	 *                     output.
	 */
	virtual void reduct_obj_free( RTOp_ReductTarget* reduct_obj ) const;
	///
	/** Extract the state of an already created reduction object.
	 *
	 * The default implementation does nothing (i.e. there is no reduction operation performed).
	 */
	virtual void extract_reduct_obj_state(
		const RTOp_ReductTarget   reduct_obj
		,int                      num_values
		,Scalar                   value_data[]
		,int                      num_indexes
		,RTOp_index_type          index_data[]
		,int                      num_chars
		,RTOp_char_type           char_data[]
		) const;
	///
	/** Load the state of an already created reduction object.
	 *
	 * The default implementation does nothing.
	 */
	virtual void load_reduct_obj_state(
		int                      num_values
		,const Scalar            value_data[]
		,int                     num_indexes
		,const RTOp_index_type   index_data[]
		,int                     num_chars
		,const RTOp_char_type    char_data[]
		,RTOp_ReductTarget       reduct_obj
		) const;

	//@}

	///
	/** Copy the state.
	 *
	 * This default version uses <tt>this->\ref extract_op_state "extract_op_state(...)"</tt>
	 * and <tt>this->\ref load_op_state "load_op_state(...)"</tt> to perform the copy.  No
	 * override should be needed by subclasses (unless a slightly more efficient implementation
	 * is desired).
	 */
	virtual RTOpT<Scalar>& operator=(const RTOpT<Scalar>& op);

}; // end class RTOp

///
/** Concrete encapulation class for an <tt>RTOp_ReductTarget</tt> object.
 *
 * This class simplifies the use of reduction objects with the
 * C++ interface to <tt>RTOpT</tt>.  Specifically, the underlying opaque reduction
 * objects are destroyed when objects of this class are destroyed.  It is
 * vitally important that the RTOp object that initialized objects of this
 * type are not destoryed first.  Every <tt>ReductTargetT</tt> object is
 * assicated with a <tt>RTOpT</tt> object and the <tt>%ReductTargetT</tt> object
 * derives all of its functionality from the <tt>%RTOpT</tt> object.
 *
 * Note that the purpose of this class is to abstract, not to encapsulate.
 * It is only designed to help in some memory management details and to provide
 * a more compact abstraction for reduction objects.
 */
template<class Scalar>
class ReductTargetT {
public:

	///
	/** Constructs the object given its operator object ready to use.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->op_ptr() == &this->op() == &op</tt>
	 * <li> <tt>this->obj()</tt> returns an already constructed <tt>RTOp_ReductTarget</tt> object
	 *      using <tt>op\ref RTOp::reduct_obj_create_raw ".reduct_obj_create_raw(&this->obj())"</tt>.
	 * </ul>
	 *
	 * @param  op  [in] The reduction/transformation operator object that owns this
	 *             reduction object.  Note that \c op must not be altered while <tt>this</tt>
	 *             is using it.
	 */
	ReductTargetT( const RTOpT<Scalar>& op );
	/// Calls <tt>this->\ref initialize "initialize(...)".
	ReductTargetT( const RTOpT<Scalar>* op = NULL, RTOp_ReductTarget obj = RTOp_REDUCT_OBJ_NULL, bool owns_obj = false );
	///
	/** Initialize <tt>this</tt> with a reduction/transformation operator object and a reduction object.
	 *
	 * @param  op       [in] Pointer to reduction/transformation object that <tt>this</tt> will be based on.
	 *                  Note that \c op can be \c NULL.
	 * @param  obj      [in] Pointer to the reduction object that <tt>this</tt> will encapsulate.
	 *                  Note that \c obj can be \c NULL.  Ignored if <tt>op == NULL</tt>.
	 * @param  owns_obj [in] If <tt>owns_obj == true</tt> then when <tt>this->free()</tt> is called
	 *                  \c obj will be freed by \c op.
	 *
	 * Note that <tt>this->free()</tt> is called to free <tt>this->obj()</tt> before <tt>this</tt> is 
	 * modified.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->op_ptr() == op</tt>
	 * <li> [<tt>op != NULL</tt>] <tt>this->obj() == obj</tt>
	 * <li> [<tt>op != NULL</tt>] <tt>this->owns_obj() == owns_obj</tt>
	 * <li> [<tt>op == NULL</tt>] <tt>this->obj() == RTOp_ReductTarget</tt>
	 * <li> [<tt>op == NULL</tt>] <tt>this->owns_obj() == true</tt>
	 * </ul>
	 */
	void initialize( const RTOpT<Scalar>* op = NULL, RTOp_ReductTarget obj = RTOp_REDUCT_OBJ_NULL, bool owns_obj = false );
	/// Calls <tt>this->free()</tt>
	~ReductTargetT();
	///
	/** Returns <tt>*this->op_ptr()</tt>.
	 *
	 * Warning!  Is an error to call if <tt>op_ptr() == NULL</tt>.
	 */
	const RTOpT<Scalar>& op() const;
	///
	/** Return a pointer to the reduction/transforamtion object that <tt>this</tt> is based on.
	 *
	 * Warning!  <tt>return</tt>can be <tt>NULL</tt>.
	 */
	const RTOpT<Scalar>* op_ptr() const;
	///
	/** Return the underlying <tt>RTOp_ReductTarget</tt> object that <tt>this</tt> abstracts!
	 */
	RTOp_ReductTarget& obj();
	///
	/** Return the underlying <tt>RTOp_ReductTarget</tt> object that <tt>this</tt> abstracts!
	 */
	const RTOp_ReductTarget& obj() const;
	/// Set whether <tt>this</tt> owns <tt>this->obj()</tt> in order to free it.
	void owns_obj(bool);
	/// Returns true if <tt>this</tt> owns <tt>this->obj()</tt> in order to free it.
	bool owns_obj() const;
	/** @name Methods derived from <tt>this->op</tt>.
	 *
	 * Note that all of these method calls will cause an exception if <tt>this->op_ptr() == NULL</tt>
	 * and may cause an exception if <tt>this->obj() == RTOp_REDUCT_OBJ_NULL</tt>.
	 */
	//@{
	/// Calls <tt>this->op()\ref RTOp::get_reduct_type_num_entries ".get_reduct_type_num_entries(,,,)"</tt>.
	void get_type_num_entries(
		int*   num_values
		,int*  num_indexes
		,int*  num_chars
		) const;
	/// Calls <tt>this->op()\ref RTOp::reduct_obj_reinit ".reduct_obj_reinit(&this->obj())"</tt>.
	void reinit();
	/// Calls <tt>this->op()\ref RTOp::reduct_obj_free ".reduct_obj_free(&this->obj())"</tt> if <tt>this->owns_obj() == true</tt>.
	void free();
	/// Calls <tt>this->op()\ref RTOp::extract_reduct_obj_state ".extract_reduct_obj_state(this->obj(),,,)"</tt>.
	void extract_state(
		int                       num_values
		,Scalar                   value_data[]
		,int                      num_indexes
		,RTOp_index_type          index_data[]
		,int                      num_chars
		,RTOp_char_type           char_data[]
		) const;
	/// Calls <tt>this->op()\ref RTOp::load_reduct_obj_state ".load_reduct_obj_state(,,,this->obj())"</tt>.
	void load_state(
		int                      num_values
		,const Scalar            value_data[]
		,int                     num_indexes
		,const RTOp_index_type   index_data[]
		,int                     num_chars
		,const RTOp_char_type    char_data[]
		);
	//@}

private:

	const RTOpT<Scalar>  *op_;
	RTOp_ReductTarget    obj_;
	bool                 owns_obj_;
	// not defined and not to be called
	ReductTargetT(const ReductTargetT<Scalar>&);
	ReductTargetT<Scalar>& operator=(const ReductTargetT<Scalar>&);

	friend class RTOpT<Scalar>;

}; // end class ReductTargetT

// /////////////////////////////////
// Inline member functions

// ReductTarget

template<class Scalar>
inline
ReductTargetT<Scalar>::ReductTargetT( const RTOpT<Scalar>& op )
	: op_(&op),obj_(RTOp_REDUCT_OBJ_NULL),owns_obj_(true)
{
	op_->reduct_obj_create(this);
}

template<class Scalar>
inline
ReductTargetT<Scalar>::ReductTargetT( const RTOpT<Scalar>* op, RTOp_ReductTarget obj, bool owns_obj )
	: op_(op), obj_(obj), owns_obj_(owns_obj)
{}

template<class Scalar>
inline
void ReductTargetT<Scalar>::initialize( const RTOpT<Scalar>* op, RTOp_ReductTarget obj, bool owns_obj )
{
	free();
	op_ = op;
	if(op_) {
		obj_       = obj;
		owns_obj_  = owns_obj;
	}
	else {
		obj_       = RTOp_REDUCT_OBJ_NULL;
		owns_obj_  = true;
	}
}

template<class Scalar>
inline
ReductTargetT<Scalar>::~ReductTargetT()
{
	free();
}

template<class Scalar>
inline
const RTOpT<Scalar>& ReductTargetT<Scalar>::op() const
{
	assert(op_);
	return *op_;
}

template<class Scalar>
inline
const RTOpT<Scalar>* ReductTargetT<Scalar>::op_ptr() const
{
	return op_;
}

template<class Scalar>
inline
RTOp_ReductTarget& ReductTargetT<Scalar>::obj()
{
	return obj_;
}

template<class Scalar>
inline
const RTOp_ReductTarget& ReductTargetT<Scalar>::obj() const
{
	return obj_;
}

template<class Scalar>
inline
void ReductTargetT<Scalar>::owns_obj(bool owns_obj)
{
	owns_obj_ = owns_obj;
}

template<class Scalar>
inline
bool ReductTargetT<Scalar>::owns_obj() const
{
	return owns_obj_;
}

template<class Scalar>
inline
void ReductTargetT<Scalar>::get_type_num_entries(
	int*   num_values
	,int*  num_indexes
	,int*  num_chars
	) const
{
	op().get_reduct_type_num_entries(num_values,num_indexes,num_chars);
}

template<class Scalar>
inline
void ReductTargetT<Scalar>::reinit()
{
	op().reduct_obj_reinit(obj_);
}

template<class Scalar>
inline
void ReductTargetT<Scalar>::free()
{
	if(obj() != RTOp_REDUCT_OBJ_NULL && owns_obj() == true)
		op().reduct_obj_free(&obj_);
}

template<class Scalar>
inline
void ReductTargetT<Scalar>::extract_state(
	int                       num_values
	,Scalar                   value_data[]
	,int                      num_indexes
	,RTOp_index_type          index_data[]
	,int                      num_chars
	,RTOp_char_type           char_data[]
	) const
{
	op().extract_reduct_obj_state(
		obj()
		,num_values,  value_data
		,num_indexes, index_data
		,num_chars,   char_data
		);
}

template<class Scalar>
inline
void ReductTargetT<Scalar>::load_state(
	int                      num_values
	,const Scalar            value_data[]
	,int                     num_indexes
	,const RTOp_index_type   index_data[]
	,int                     num_chars
	,const RTOp_char_type    char_data[]
	)
{
	op().load_reduct_obj_state(
		num_values ,  value_data
		,num_indexes, index_data
		,num_chars,   char_data
		,obj()
		);
}

} // end namespace RTOpPack

//@}

#endif // RTOPPACK_RTOP_CPP_DECL_HPP
