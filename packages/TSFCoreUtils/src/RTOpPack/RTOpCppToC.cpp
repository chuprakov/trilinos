// /////////////////////////////////////////////
// RTOpCppToC.cpp
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
		op_c_ = new RTOp_RTOp;
		op_c_->vtbl     = op_vtbl_;
		op_c_->obj_data = (void*)&op_cpp;
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
