// //////////////////////////////////////////////////////////////////////////////////////
// RTOpUnaryFuncPtrDecl.hpp

#ifndef RTOPPACK_UNARY_FUNC_PTR_DECL_HPP
#define RTOPPACK_UNARY_FUNC_PTR_DECL_HPP

#include "RTOpCpp.hpp"

namespace RTOpPack {

///
/** <tt>RTOpT</tt> subclass for unary functions using a unary function pointer.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class RTOpUnaryFuncPtr : public RTOpT<Scalar> {
public:

	///
	typedef void (*unary_func_ptr_t) ( const Scalar x[], int x_dim, Scalar out[] );

	/// Construct to uninitialized
	RTOpUnaryFuncPtr();

	/// Calls <tt>initialize()</tt>
	RTOpUnaryFuncPtr(
		unary_func_ptr_t        unary_func_ptr
		,const std::string      &op_name
		);

	///
	/** Initialize.
	 *
	 * @param  unary_func_ptr
	 *               [in] Pointer to function that actually performs the unary operation.
	 * @param  op_name
	 *               [in] Name of the operation (for debugging mostly by clients)
	 *
	 * Preconditions:<ul>
	 * <li> <tt>unary_func_ptr != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 */
	void initialize(
		unary_func_ptr_t        unary_func_ptr
		,const std::string      &op_name
		);

	///
	/** Set uninitialized.
	 *
	 * @param  unary_func_ptr
	 *               [out] If <tt>unary_func_ptr!=NULL</tt> then <tt>*unary_func_ptr</tt>
	 *               is set to pointer to function that was passed in to <tt>initialize()</tt>.
	 * @param  op_name
	 *               [out] If <tt>op_name!=NULL</tt> then <tt>*op_name</tt>
	 *               is set to the operation name that was passed in to <tt>initialize()</tt>.
	 */
	void set_initialized(
		unary_func_ptr_t    *unary_func_ptr  = NULL
		,std::string        *op_name         = NULL
		);

	/** @name Overridden from RTOpT */
	//@{

	///
	const char* op_name() const;
	///
	void apply_op(
		const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
		,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
		,RTOp_ReductTarget reduct_obj
		) const;

	//@}

private:
	
	std::string        op_name_;
	unary_func_ptr_t   unary_func_ptr_;

	// Not defined and not to be called
	RTOpUnaryFuncPtr(const RTOpUnaryFuncPtr&);
	RTOpUnaryFuncPtr& operator=(const RTOpUnaryFuncPtr&);

};

} // end namespace RTOpPack

#endif // RTOPPACK_UNARY_FUNC_PTR_DECL_HPP
