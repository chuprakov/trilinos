/* /////////////////////////////////////////////
// RTOp_RTOp_C_Cpp.h
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
*/

#ifndef RTOP_RTOP_C_CPP_H
#define RTOP_RTOP_C_CPP_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_RTOp_C_Cpp.h C interface for a C++ reduction/transformation operator.
 *
 * These functions encapsulate the population of a virtual function table for an C <tt>RTOp_RTOp</tt>
 * object that is based on an C++ <tt>RTOpPack::RTOp</tt> object.  In order to create a <tt>%RTOp_RTOp</tt>
 * object \c op_c given a <tt>%RTOpPack::RTOp</tt> object \c op_cpp a client could define a function like:
 \verbatim

 void create_C_op( const RTOpPack::RTOp& op_cpp, RTOp_RTOp* op_c )
 {
     RTOp_create_C_Cpp_vtbl(
         op_cpp.op_factory().get_op_create_func()
         ,op_cpp.op_factory().get_op_free_func()
		 ,op_c->vtbl
         );
     op_c->obj_data = (void*)&op_cpp;
 }
 \endverbatim
 * From the point of view of a C client, the \c op_c operator object will behave just as
 * if it had a native implementation in C.  The only draw back here is an extra level
 * of indirection but no big deal.
 *
 * To free the C operator object \c op_c, a function like the following could be defined:
 \verbatim

 void free_C_op( RTOp_RTOp* op_c )
 {
     RTOp_free_C_Cpp_vtbl( op_c->vtbl );
     op_c->obj_data = NULL;
 }
 \endverbatim
 *
 * Even though these functions have C linkage, they have C++ implementations.
 */
/*@{*/

typedef int (*RTOp_op_create_func_ptr_t)( const struct RTOp_obj_type_vtbl_t* vtbl, const void* dummy, void** op_data );
typedef int (*RTOp_op_free_func_ptr_t)( const struct RTOp_obj_type_vtbl_t* vtbl, const void* dummy, void** op_data );

/** Populate a virtual function table for <tt>RTOp_RTOp</tt> given the function <tt>op_create(...)</tt>
 * that will create an RTOpPack::RTOp object and <tt>op_free(...)</tt> that will free it.
 *
 * @param  op_create  [in] Pointer to function that will create <tt>RTOpPack::RTOp</tt> object.
 * @param  op_free    [in] Pointer to function that will free <tt>RTOpPack::RTOp</tt> object
 *                    created by <tt>op_create(...)</tt>.
 * @param  vtbl       [out] The virtual function table populated with virtual functions.
 *
 * Preconditions:<ul>
 * <li> <tt>op_create</tt> != NULL
 * <li> <tt>op_free</tt> != NULL
 * <li> <tt>vtbl->obj_data_vtbl == NULL</tt>
 * <li> <tt>vtbl->reduct_vtbl == NULL</tt>
 * <li> <tt>vtbl->apply == NULL</tt>
 * <li> <tt>vtbl->reduce_reduct_objs == NULL</tt>
 * <li> <tt>vtbl->get_reduct_op == NULL</tt>
 * </ul>
 */
int RTOp_create_C_Cpp_vtbl( 
	RTOp_op_create_func_ptr_t op_create, RTOp_op_free_func_ptr_t op_free
	,struct RTOp_RTOp_vtbl_t* vtbl );

/** Free any dynamically allocated data for a virtual function table populated by <tt>RTOp_create_C_Cpp_vtbl()</tt>.
 *
 * @param  vtbl       [in/out] On input, points to a virtual function table populated by <tt>RTOp_create_C_Cpp_tbl()</tt>.
 *                    On output, all members set to NULL.
 *
 * Preconditions:<ul>
 * <li> <tt>vtbl->obj_data_vtbl != NULL</tt>
 * <li> <tt>vtbl->reduct_vtbl != NULL</tt>
 * <li> <tt>vtbl->apply != NULL</tt>
 * </ul>
 *
 * Postconditions:<ul>
 * <li> <tt>vtbl->obj_data_vtbl == NULL</tt>
 * <li> <tt>vtbl->reduct_vtbl == NULL</tt>
 * <li> <tt>vtbl->apply == NULL</tt>
 * <li> <tt>vtbl->reduce_reduct_objs == NULL</tt>
 * <li> <tt>vtbl->get_reduct_op == NULL</tt>
 * </ul>
 */
int RTOp_free_C_Cpp_vtbl( struct RTOp_RTOp_vtbl_t* vtbl );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_RTOP_C_CPP_H */
