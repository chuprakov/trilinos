// /////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorBase.hpp

#ifndef TSFCORE_VECTOR_SERIAL_BASE_HPP
#define TSFCORE_VECTOR_SERIAL_BASE_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreSerialVectorBaseDecl.hpp"
#include "TSFCore_apply_op_helper.hpp"
#include "WorkspacePack.hpp"
#include "ThrowException.hpp"

namespace TSFCore {

template<class Scalar>
SerialVectorBase<Scalar>::SerialVectorBase()
	:in_applyOp_(false)
{}

template<class Scalar>
void SerialVectorBase<Scalar>::applyOp(
	const RTOpPack::RTOpT<Scalar>   &op
	,const size_t                   num_vecs
	,const Vector<Scalar>*          vecs[]
	,const size_t                   num_targ_vecs
	,Vector<Scalar>*                targ_vecs[]
	,RTOp_ReductTarget              reduct_obj
	,const Index                    first_ele
	,const Index                    sub_dim
	,const Index                    global_offset
	) const
{
#ifdef _DEBUG
	TSFCore::apply_op_validate_input(
		"SerialVectorBase::applyOp(...)"
		,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset
		);
	THROW_EXCEPTION(
		in_applyOp_, std::logic_error
		,"SerialVectorBase::applyOp(...): Error, something is not right here!" );
#endif
	in_applyOp_ = true;
	TSFCore::apply_op_serial(
		op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
		,first_ele,sub_dim,global_offset
		);
	in_applyOp_ = false;
}

} // end namespace TSFCore

#endif // TSFCORE_VECTOR_SERIAL_BASE_HPP
