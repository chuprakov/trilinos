// //////////////////////////////////////////////////////////////////////////////
// TSFCore_apply_op_helper_decl.hpp

#ifndef TSFCORE_APPLY_OP_HELPER_DECL_HPP
#define TSFCORE_APPLY_OP_HELPER_DECL_HPP

#include "TSFCoreTypes.hpp"
#include "RTOpCpp.hpp"

namespace TSFCore {

///
/** Validate the inputs to applyOp(...).
 *
 * Throws an exception if one of the preconditions is not met.
 *
 * ToDo: Finish documentation.
 */
template<class Scalar>
void apply_op_validate_input(
	const char                      func_name[]
	,const RTOpPack::RTOpT<Scalar>  &op
	,const size_t                   num_vecs
	,const Vector<Scalar>*          vecs[]
	,const size_t                   num_targ_vecs
	,Vector<Scalar>*                targ_vecs[]
	,RTOp_ReductTarget              reduct_obj
	,const Index                    first_ele
	,const Index                    sub_dim
	,const Index                    global_offset
	);

///
/** Implements reduction/transformation operators for any serial
 * vectors using just the public vector interface.
 *
 * Note that this function does not validate the input arguments so it is up to
 * the client to do that (i.e. by calling <tt>apply_op_validate_input()</tt>).
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
void apply_op_serial(
	const RTOpPack::RTOpT<Scalar>  &op
	,const size_t                  num_vecs
	,const Vector<Scalar>*         vecs[]
	,const size_t                  num_targ_vecs
	,Vector<Scalar>*               targ_vecs[]
	,RTOp_ReductTarget             reduct_obj
	,const Index                   first_ele
	,const Index                   sub_dim
	,const Index                   global_offset
	);

} // end namespace TSFCore

#endif // TSFCORE_APPLY_OP_HELPER_DECL_HPP
