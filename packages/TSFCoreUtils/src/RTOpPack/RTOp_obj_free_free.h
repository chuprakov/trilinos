// /////////////////////////////////////////////////////////////////////////////
// RTOp_obj_free_free.h
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

#ifndef RTOP_OBJ_FREE_FREE_H
#define RTOP_OBJ_FREE_FREE_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

///
/** Definition of function to be used for </tt>obj_free<tt> in RTOp_RTOp_vtbl_t which
 * simply calls </tt>free(...)<tt>.
 *
 * @param  vtbl           [in] Totally ignored.
 * @param  instance_data  [in] Totally ignored.
 * @param  obj            [in/out]  If </tt>*obj != NULL<tt> on input then </tt>free(*obj)<tt> is called.
 *                        On output, </tt>*obj<tt> is set to NULL.
 *
 * @return Returns </tt>0<tt> for success.
 */
int RTOp_obj_free_free( const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data, void** obj );

#ifdef __cplusplus
}
#endif

#endif // RTOP_OBJ_FREE_FREE_H
