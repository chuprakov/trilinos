/* ///////////////////////////////////////// */
/* RTOp_obj_null_vtbl.h */
/* */
/* Copyright (C) 2001 Roscoe Ainsworth Bartlett */
/* */
/* This is free software; you can redistribute it and/or modify it */
/* under the terms of the "Artistic License" (see the web site */
/*   http://www.opensource.org/licenses/artistic-license.html). */
/* This license is spelled out in the file COPYING. */
/* */
/* This software is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* above mentioned "Artistic License" for more details. */

#ifndef RTOP_OBJ_NULL_VTBL_H
#define RTOP_OBJ_NULL_VTBL_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* */
/** Virtual function table for a null (none) object.
  *
  * The functions do the following:
  * <ul>
  *	<li> <tt>get_obj_type_num_entries(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>num_values</tt>    [out] Returns 0
  *		<li> <tt>num_indexes</tt>   [out] Returns 0
  *		<li> <tt>num_chars</tt>     [out] Returns 0
  *		</ul>
  *	<li> <tt>obj_create(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>obj</tt>           [out] set to NULL
  *		</ul>
  *	<li> <tt>obj_reinit(...)</tt>    Does nothing
  *	<li> <tt>obj_free(...)</tt>      Does nothing
  *	<li> <tt>extract_state(...)</tt> Does nothing
  *	<li> <tt>load_state(...)</tt>    Does nothing
  *	</ul>
  */
extern const struct RTOp_obj_type_vtbl_t   RTOp_obj_null_vtbl;

#ifdef __cplusplus
}
#endif

#endif /* RTOP_OBJ_NULL_VTBL_H */
