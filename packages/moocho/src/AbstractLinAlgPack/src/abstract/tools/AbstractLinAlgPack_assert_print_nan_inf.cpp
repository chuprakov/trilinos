// ///////////////////////////////////////////////////////////
// assert_print_nan_inf.cpp
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

#include <ostream>
#include <iomanip>

#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/VectorWithOp.h"
#include "RTOpStdOpsLib/include/RTOp_ROp_find_nan_inf.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "RTOpPack/include/check_nan_inf.h"
#include "ThrowException.h"

namespace {

// Find a NaN or Inf element!
static RTOpPack::RTOpC          find_nan_inf_op;
static RTOpPack::ReductTarget   find_nan_inf_targ;

// Simple class for an object that will initialize the RTOp_Server.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// Operator and target obj for find_nan_inf
		if(0!=RTOp_ROp_find_nan_inf_construct(&find_nan_inf_op.op() ))
			assert(0);
		find_nan_inf_op.reduct_obj_create(&find_nan_inf_targ);
		if(0!=RTOp_Server_add_op_name_vtbl(
			   RTOp_ROp_find_nan_inf_name
			   ,&RTOp_ROp_find_nan_inf_vtbl
			   ))
			assert(0);
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

bool AbstractLinAlgPack::assert_print_nan_inf( const value_type& val, char name[]
	, bool throw_excpt, std::ostream* out )
{
	if( RTOp_is_nan_inf(val) ) {
		std::ostringstream omsg;
		omsg
			<< "The scalar \"" << name
			<< "\" = " << val << " is not a valid bounded number";
		if(out)
			*out << omsg.str() << std::endl;
		THROW_EXCEPTION(
			throw_excpt,NaNInfException
			,"assert_print_nan_inf(...) : Error, " << omsg.str() );
		return false;
	}
	return true;
}

bool AbstractLinAlgPack::assert_print_nan_inf( const VectorWithOp& v, char name[]
	, bool throw_excpt, std::ostream* out )
{
	find_nan_inf_targ.reinit();
	v.apply_reduction(find_nan_inf_op,0,NULL,0,NULL,find_nan_inf_targ.obj() );
	RTOp_ROp_find_nan_inf_reduct_obj_t
		ele =RTOp_ROp_find_nan_inf_val(find_nan_inf_targ.obj());
	if(out && ele.i) {
		*out
			<< "The vector \"" << name << "\" has the first following NaN or Inf element\n"
			<< name << "(" << ele.i << ") = " << ele.v0_i << std::endl;
	}
	THROW_EXCEPTION(
		ele.i && throw_excpt, NaNInfException
		,"assert_print_nan_inf(...) : Error, the vector named "
		<< name << " has at least one element which is NaN or Inf" );
	
	return ele.i == 0;
}
