/* ///////////////////////////////////////////// */
/* RTOp_ROp_num_bounded.h */
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

#ifndef RTOP_ROP_NUM_BOUNDED_H
#define RTOP_ROP_NUM_BOUNDED_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_num_bounded.h Reduction operator counts the number of elements with finite bounds.
  *
  * <tt>targ_obj <- size( { i | v[0](i) > -inf_bnd || v[1](i) < +inf_bnd } )</tt>
  *
  * This operator is defined to allow exactly two vector arguments
  * (<tt>num_vecs == 2</tt>) <tt>v[0]</tt>, <tt>v[1]</tt> and can only handle dense vectors.
  */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_num_bounded_vtbl;

/* Constructor */
int RTOp_ROp_num_bounded_construct( RTOp_value_type inf_bnd, struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_num_bounded_destroy( struct RTOp_RTOp* op );

/* Reset inf_bnd */
int RTOp_ROp_num_bounded_set_inf_bnd( RTOp_value_type inf_bnd, struct RTOp_RTOp* op );

/* Extract the number of bounded variables */
RTOp_index_type RTOp_ROp_num_bounded_val(RTOp_ReductTarget targ_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_NUM_BOUNDED_H */
