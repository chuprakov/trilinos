// /////////////////////////////////////////////
// RTOp_TOp_sign.h

#ifndef RTOp_TOp_sign_H
#define RTOp_TOp_sign_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  // vtbl for reduction object data

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_sign.h
 *
 * This operator computes the sign of each element in v0 as:
 \verbatim
  z0(i) = sign(v0(i)), for i = 1...n
 \endverbatim
 * where
 \verbatim
                /  -1.0 : if alpha  < 0.0
 sign(alpha) =  |   0.0 : if alpha == 0.0
                \  +1.0 : if alpha  > 0.0
 \endverbatim
 */
//@{

/// Name of this transformation operator class
extern const char RTOp_TOp_sign_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_sign_vtbl;

/// Constructor
int RTOp_TOp_sign_construct(  struct RTOp_RTOp* op );

/// Destructor
int RTOp_TOp_sign_destroy( struct RTOp_RTOp* op );




//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOp_TOp_sign_H
