// /////////////////////////////////////////////
// ExampleNLPFirstOrderDirectRTOps.h

#ifndef EXAMPLE_NLP_FIRST_ORDER_DIRECT_RTOPS_H
#define EXAMPLE_NLP_FIRST_ORDER_DIRECT_RTOPS_H

#include "RTOpPack/include/RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \defgroup explnlp2_ops_grp Reduction/Transformation operators for example NLP subclass.
 *
 * These operator classes are used in the implementation of
 * <tt>\ref NLPInterfacePack::ExampleNLPFirstOrderDirect "ExampleNLPFirstOrderDirect"</tt>.
 * The binary code for these RTOp operator classes (as well as any others) must be loaded into
 * runtime environment in any process where a vector implementation must apply it.
 */
//@{

/** \defgroup explnlp2_eval_grp Evaluate the constraints for the example NLP.
 *
 * <tt>z[0](i) <- v[0](i) * (v[1](i) - 1) - 10 * v[1](i), for i = 1...n</tt>
 *
 * This operator is only admits dense vectors and is only defined for <tt>num_vecs == 2</tt>
 * and <tt>num_targ_vecs == 1</tt>.
 */
//@{

/// Name of this transformation operator class
extern const char RTOp_TOp_explnlp2_c_eval_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_explnlp2_c_eval_vtbl;

/// Constructor
int RTOp_TOp_explnlp2_c_eval_construct( struct RTOp_RTOp* op );

/// Destructor
int RTOp_TOp_explnlp2_c_eval_destroy( struct RTOp_RTOp* op );

//@}

/** \defgroup explnlp2_calc_py_D_grp Evaluate py = -inv(C)*c and/or D = inv(C)*N for the example %NLP.
 *
 * This operator performs the following:
 \verbatim

 task = 0 (py only, num_vecs = 2, num_targ_vecs = 1):
     py(i) <- c(i) / ( 1.0 - xD(i) ), i = 1...n
     where: xD = vec[0], c = vec[1], py = targ_vec[0]
 task = 1 (D only, num_vecs = 2, num_targ_vecs = 1):
     d(i) <-(xI(i) - 10.0) / (1.0 - xD(i)), i = 1...n
     where:  xD = vec[0], xI = vec[1], d = targ_vec[0]
 task = 2 (py and D, num_vecs = 3, num_targ_ves = 2)
     py(i) = c(i) / ( 1.0 - xD(i) ), i = 1...n
     d(i) = (xI(i) - 10.0) / (1.0 - xD(i)), i = 1...n
     where: xD = vec[0], xI = vec[1], c = vec[2], d = targ_vec[0], py = targ_vec[1]
 \endverbatim
 */
//@{

/// Name of this transformation operator class
extern const char RTOp_TOp_explnlp2_calc_py_D_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_explnlp2_calc_py_D_vtbl;

/// Constructor (task = 0, 1 or 2)
int RTOp_TOp_explnlp2_calc_py_D_construct( int task, struct RTOp_RTOp* op );

/// Set the task
int RTOp_TOp_explnlp2_calc_py_D_set_task( int task, struct RTOp_RTOp* op );

/// Destructor
int RTOp_TOp_explnlp2_calc_py_D_destroy( struct RTOp_RTOp* op );

//@}

//@}

#ifdef __cplusplus
}
#endif

#endif  // EXAMPLE_NLP_FIRST_ORDER_DIRECT_RTOPS_H
