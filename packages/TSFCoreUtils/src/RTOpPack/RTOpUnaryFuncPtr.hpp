// //////////////////////////////////////////////////////////////////////////////////////
// RTOpUnaryFuncPtr.hpp

#ifndef RTOPPACK_UNARY_FUNC_PTR_HPP
#define RTOPPACK_UNARY_FUNC_PTR_HPP

#include "RTOpUnaryFuncPtrDecl.hpp"

namespace RTOpPack {

template<class Scalar>
RTOpUnaryFuncPtr<Scalar>::RTOpUnaryFuncPtr()
{
	set_initialized();
}

template<class Scalar>
RTOpUnaryFuncPtr<Scalar>::RTOpUnaryFuncPtr(
	unary_func_ptr_t        unary_func_ptr
	,const std::string      &op_name
	)
{
	initialize(unary_func_ptr,op_name);
}

template<class Scalar>
void RTOpUnaryFuncPtr<Scalar>::initialize(
	unary_func_ptr_t        unary_func_ptr
	,const std::string      &op_name
	)
{
	TEST_FOR_EXCEPTION( unary_func_ptr==NULL, std::invalid_argument, "Error!" );
	unary_func_ptr_ = unary_func_ptr;
	op_name_ = op_name;
}

template<class Scalar>
void RTOpUnaryFuncPtr<Scalar>::set_initialized(
	unary_func_ptr_t    *unary_func_ptr
	,std::string        *op_name
	)
{
	if(unary_func_ptr) *unary_func_ptr = unary_func_ptr_;
	if(op_name) *op_name = op_name_;

	unary_func_ptr_ = NULL;
	op_name_ = "uninitialized()";
}

// Overridden from RTOpT

template<class Scalar>
const char* RTOpUnaryFuncPtr<Scalar>::op_name() const
{
	return op_name_.c_str();
}

template<class Scalar>
void RTOpUnaryFuncPtr<Scalar>::apply_op(
	const int   num_vecs,       const SubVectorT<Scalar>         sub_vecs[]
	,const int  num_targ_vecs,  const MutableSubVectorT<Scalar>  targ_sub_vecs[]
	,RTOp_ReductTarget reduct_obj
	) const
{
	TEST_FOR_EXCEPTION( num_vecs != 1 || sub_vecs == NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( num_targ_vecs != 1 || targ_sub_vecs == NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( reduct_obj != RTOp_REDUCT_OBJ_NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( sub_vecs[0].stride() != 1, std::invalid_argument, "Error, can't handle nonunit strides here!" );
	TEST_FOR_EXCEPTION( targ_sub_vecs[0].stride() != 1, std::invalid_argument, "Error, can't handle nonunit strides here!" );
	TEST_FOR_EXCEPTION( sub_vecs[0].subDim() != targ_sub_vecs[0].subDim(), std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( sub_vecs[0].globalOffset() != targ_sub_vecs[0].globalOffset(), std::invalid_argument, "Error!" );

	unary_func_ptr_( sub_vecs[0].values(), sub_vecs[0].subDim(), targ_sub_vecs[0].values() );

}

} // end namespace RTOpPack

#endif // RTOPPACK_UNARY_FUNC_PTR_HPP
