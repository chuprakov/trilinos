// /////////////////////////////////////////////////////////////
// RTOpFactoryTmpl.hpp
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
