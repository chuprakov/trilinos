/* ///////////////////////////////////////////// */
/* RTOp_ROp_dot_prod.h */
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

#ifndef RTOP_ROP_DOT_PROD_H
#define RTOP_ROP_DOT_PROD_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_dot_prod.h Reduction operator for the dot product of two vectors.
  *
  * <tt>targ_obj <- sum( v[0](i) * v[1](i), i=1...n )</tt>
  *
  * This operator is defined to allow exactly two vector arguments
  * (<tt>num_vecs == 2</tt>) <tt>v[0]</tt>, <tt>v[1]</tt>, and it can handle
  * all combinations of sparse and dense vectors.
  */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_dot_prod_vtbl;

/* Constructor */
int RTOp_ROp_dot_prod_construct( struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_dot_prod_destroy( struct RTOp_RTOp* op );

/* Extract the value of the dot product */
RTOp_value_type RTOp_ROp_dot_prod_val(RTOp_ReductTarget targ_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_DOT_PROD_H */
