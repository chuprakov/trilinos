// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// //////////////////////////////////////////////////////////////////////////////
// TSFCore_apply_op_helper_decl.hpp

#ifndef TSFCORE_APPLY_OP_HELPER_DECL_HPP
#define TSFCORE_APPLY_OP_HELPER_DECL_HPP

#include "TSFCoreTypes.hpp"
#include "RTOpPack_RTOpT.hpp"

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
	,const int                   num_vecs
	,const Vector<Scalar>*          vecs[]
	,const int                   num_targ_vecs
	,Vector<Scalar>*                targ_vecs[]
	,RTOpPack::ReductTarget         *reduct_obj
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
	,const int                  num_vecs
	,const Vector<Scalar>*         vecs[]
	,const int                  num_targ_vecs
	,Vector<Scalar>*               targ_vecs[]
	,RTOpPack::ReductTarget        *reduct_obj
	,const Index                   first_ele
	,const Index                   sub_dim
	,const Index                   global_offset
	);

} // end namespace TSFCore

#endif // TSFCORE_APPLY_OP_HELPER_DECL_HPP
