// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
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

#include "RTOpCppToC.hpp"
#include "RTOpCppC.hpp"
#include "RTOp_RTOp_C_Cpp.h"

namespace RTOpPack {

RTOpCppToC::RTOpCppToC( const RTOp &op_cpp )
	: op_cpp_(op_cpp), op_c_(NULL), op_vtbl_(NULL)
{
	if( const RTOpC *op_cpp_c = dynamic_cast<const RTOpC*>(&op_cpp_) ) {
		op_c_ = const_cast<RTOp_RTOp*>(&op_cpp_c->op());
	} else {
		op_vtbl_ = new RTOp_RTOp_vtbl_t;
	    RTOp_create_C_Cpp_vtbl(
			op_cpp_.get_op_create_func()
			,op_cpp_.get_op_free_func()
			,op_vtbl_
			);
		op_vtbl_->op_name    = op_cpp.op_name();
		op_c_                = new RTOp_RTOp;
		op_c_->vtbl          = op_vtbl_;
		op_c_->obj_data      = (void*)&op_cpp;
	}
}

RTOpCppToC::~RTOpCppToC()
{
	if(op_vtbl_) {
		RTOp_free_C_Cpp_vtbl( op_vtbl_ );
		delete op_vtbl_;
		delete op_c_;
	}
}

} // end namspace RTOpPack
