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

#ifndef RTOP_FACTORY_TMPL_H
#define RTOP_FACTORY_TMPL_H

#include "RTOpCpp.hpp"

namespace RTOpPack {

///
/** Template class to create RTOp objects.
 *
 * For concrete C++ RTOp classes, this is the only
 * factory object needed as long as the class #Op# has a 
 * default constructor associated with it.  The class #Op#
 * must also allow #new# and #delete#.
 */
template<class Op>
class RTOpFactoryTmpl : public RTOpFactory
{
public:

	/** @name Overridden from RTOpFactory */
	//@{

	///
	op_create_func_t get_op_create_func() const;
	///
	op_free_func_t get_op_free_func() const;
	///
	op_ptr_t create_op() const;

	//@}

}; // class RTOpFactoryTmpl<...>

// //////////////////////////////////////////////////
// Template nonmember functions for C compatibility
//
// Pointers are taken from these member functions.
// It is imparative that the C++ compiler instantiate
// template functions when the address of the template
// function is taken.

template<class Op>
int op_create_func_tmpl(const RTOp_obj_type_vtbl_t* dummy1, const void* dummy2, void** op)
{
	*op = static_cast<RTOpPack::RTOp*>(new Op());
	return 0;
}

template<class Op>
int op_free_func_tmpl(const RTOp_obj_type_vtbl_t* dummy1, const void* dummy2, void** op)
{
    Op* _op = dynamic_cast<Op*>(static_cast<RTOpPack::RTOp*>(*op));
    assert(_op);
	delete _op;
    *op = NULL;
	return 0;
}

// //////////////////////////////////////
// Template member functions

template<class Op>
op_create_func_t RTOpFactoryTmpl<Op>::get_op_create_func() const
{
#if defined(_KAI_CXX) || defined(_CPQ_CXX)
//	return reinterpret_cast<op_create_func_t>(op_create_func_tmpl<Op>); // C++ compiler had better instantiate this function!
	return NULL; // Error!
#else
	return op_create_func_tmpl<Op>; // C++ compiler had better instantiate this function!
#endif
}

template<class Op>
op_free_func_t RTOpFactoryTmpl<Op>::get_op_free_func() const
{
#if defined(_KAI_CXX) || defined(_CPQ_CXX)
//	return reinterpret_cast<op_free_func_t>(op_free_func_tmpl<Op>); // C++ compiler had better instantiate this function!
	return NULL; // Error!
#else
	return op_free_func_tmpl<Op>; // C++ compiler had better instantiate this function!
#endif
}

template<class Op>
RTOpFactory::op_ptr_t RTOpFactoryTmpl<Op>::create_op() const
{
	namespace rcp = MemMngPack;
	return rcp::RefCountPtr<Op>(new Op());
}

} // end namespace RTOpPack

#endif // RTOP_FACTORY_TMPL_H
