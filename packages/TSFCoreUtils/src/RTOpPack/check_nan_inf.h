/* /////////////////////////////////////////////////////////////////////////// */
/* check_nan_inf.h */
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

#ifndef CHECK_NAN_INF_H
#define CHECK_NAN_INF_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file check_nan_inf.h Utility functions for determining if a number is not a regular floating point number.
 */
/*@{ */

extern RTOp_value_type
    RTOp_pos_inf,
	RTOp_neg_inf,
    RTOp_pos_nan,
    RTOp_neg_nan;

/* Return true if the number is +-NaN if this reprsentation exists. */
int RTOp_is_nan( RTOp_value_type val );

/* Return true if the number is +-infinity if this reprsentation exists. */
int RTOp_is_inf( RTOp_value_type val );

/* Return true if the number is +-NaN or +-infinity if this reprsentation exists. */
int RTOp_is_nan_inf( RTOp_value_type val );

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* CHECK_NAN_INF_H */
