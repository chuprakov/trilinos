/* ///////////////////////////////////////////// */
/* RTOp_TOp_ele_wise_divide.h */
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

#ifndef RTOP_TOP_ELE_WISE_DIVIDE_H
#define RTOP_TOP_ELE_WISE_DIVIDE_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_ele_wise_divide.h Compute an element-wise division of two vectors.
 *
 * <tt>targ_vec[0](i) <- alpha * vec[0](i) / vec[1](i), for i = 1...n</tt>
 *
 * This operator is only admits dense vectors.
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_ele_wise_divide_vtbl;

/* Constructor */
int RTOp_TOp_ele_wise_divide_construct( RTOp_value_type alpha, struct RTOp_RTOp* op );

/* Set alpha */
int RTOp_TOp_ele_wise_divide_set_alpha( RTOp_value_type alpha, struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_ele_wise_divide_destroy( struct RTOp_RTOp* op );

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_TOP_ELE_WISE_DIVIDE_H */
