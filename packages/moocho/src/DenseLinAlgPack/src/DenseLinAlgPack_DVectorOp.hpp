// //////////////////////////////////////////////////////////////////////////////////
// VectorOp.h
//

// See GenMatrixOp.h for description of naming convensions

#ifndef VECTOROP_H
#define VECTOROP_H

#include "VectorAssign.h"

/** @name {\bf Basic Vector Operation Functions (Level-1 BLAS)}.
  *
  * These are functions that perform basic operations with vectors such as element-wise
  * linear algebra operations (e.g. v1 = v2 + v3, v1 = sin(v2)) and other vector 
  * related functions.  The functions that have vectors as lhs arguments come in
  * two varieties: those with a Vector object as the lhs argument (v_lhs), those with
  * a VectorSlice object as a lhs argument (vs_lhs).  Having different functions
  * for Vector and VectorSlice objects as lhs arguments is important because
  * Vectors resize to the rhs expression while VectorSlice objects do not.
  *
  * Only VectorSlice objects are used as rhs arguments however.  When Vector objects
  * are used to call these fucntions as rhs arguments, the implicit type conversion
  * to a const temp VectorSlice will be performed to make the call work.
  * 
  * The implementations of these functions takes care of the following details:
  *
  * \begin{itemize}
  *	\item Resizing Vector LHS on assignment
  *	\item Test for aliasing of assign(...) but not other functions
  *	\item Check preconditions (sizes of arguments) if LINALGPACK_CHECK_RHS_SIZES is defined
  * \end{itemize}
  *
  * These functions share common behavior and precondtions which are listed below.
  *
  * Preconditions for functions with a VectorSlice object (vs_lhs) as a lhs argument
  * (e.g. vs_lhs = abs(vs_rhs), vs_lhs = vs_rhs1 + vs_rhs2).
  *	\begin{itemize}
  * \item #vs_lhs.size() ==# size of rhs expression  (throw #std::length_error#)
  * \end{itemize}
  *
  * Preconditions for functions with two VectorSlice objects (vs_rhs1, vs_rhs2) rhs arguments
  * (e.g. v_lhs = pow(vs_rhs1,vs_rhs2), result = trans(vs_rhs1,vs_rhs2)):
  *	\begin{itemize}
  * \item #vs_rhs1.size() == vs_rhs2.size()#  (throw #std::length_error#)
  * \end{itemize}
  *
  * Algebric functions are named according to the types of their arguments.  For example,
  * the function for the operation vs_lhs = vs_rhs1 - vs_rhs2 is named V_VmV(...).
  * For a description of this namming format see
  * \Ref{LinAlgOpPack}
  */

//@{
//		begin Basic Vector Operation Functions

namespace LinAlgPack {

/** @name {\bf Algebraic Functions}.
  *
  * The functions assign(...) are used by the implementation of the assignment operators for
  * Vector and VectorSlice and therefore the user can use the assignment operator to
  * perform the copies.
  */

//@{
//		begin Algebraic Functions


/// vs_lhs += alpha
void Vp_S(VectorSlice* vs_lhs, value_type alpha);
/// vs_lhs *= alpha (BLAS xSCAL)
void Vt_S(VectorSlice* vs_lhs, value_type alpha);
/// vs_lhs += alpha * vs_rhs (BLAS xAXPY)
void Vp_StV(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs);

/// v_lhs = alpha (elementwise)
//void assign(Vector* v_lhs, value_type alpha);
/// v_lhs = vs_rhs.
//void assign(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = vs_rhs1 + vs_rhs2
void V_VpV(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// v_lhs = vs_rhs1 - vs_rhs2
void V_VmV(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// v_lhs = - vs_rhs
void V_mV(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = alpha * vs_rhs
void V_StV(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs);

/// vs_lhs = alpha (elementwise)
//void assign(VectorSlice* vs_lhs, value_type alpha);
/// vs_lhs = vs_rhs
//void assign(VectorSlice* vs_lhs, const VectorSlice& vx_rhs);
/// vs_lhs = vs_rhs1 + vs_rhs2
void V_VpV(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// vs_lhs = vs_rhs1 - vs_rhs2
void V_VmV(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// vs_lhs = - vs_rhs
void V_mV(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = alpha * vs_rhs
void V_StV(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs);

///
/* Apply a plane (Givens) rotation.
 *
 * [  c  s ] * [ x' ] -> [ x' ]
 * [ -s  c ]   [ y' ]    [ y' ]
 *
 * See "Handbook for Matrix Computations" section 2.4
 */
void rot( const value_type c, const value_type s, VectorSlice* x, VectorSlice* y );

//		end Algebraic Functions
//@}

/** @name {\bf Elementwise Math Vector / VectorSlice Functions}. */

//@{
//		begin Elementsize Math Functions

/// vs_lhs = abs(vs_rhs)
void abs(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = asin(vs_rhs)
void asin(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = acos(vs_rhs)
void acos(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = atan(vs_rhs)
void atan(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = atan(vs_rhs1/vs_rhs2)
void atan2(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// vs_lhs = atan(vs_rhs/alpha)
void atan2(VectorSlice* vs_lhs, const VectorSlice& vs_rhs, value_type alpha);
/// vs_lhs = atan(alpha/vs_rhs)
void atan2(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs);
/// vs_lhs = cos(vs_rhs)
void cos(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = cosh(vs_rhs)
void cosh(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = exp(vs_rhs)
void exp(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = max(vs_rhs1,vs_rhs2)
void max(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// vs_lhs = max(alpha,vs_rhs)
void max(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs);
/// vs_lhs = min(vs_rhs1,vs_rhs2)
void min(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// vs_lhs = mim(alpha,vs_rhs)
void min(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs);
/// vs_lhs = pow(vs_rhs1,vs_rhs2)
void pow(VectorSlice* vs_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// vs_lhs = pow(vs_rhs,alpha)
void pow(VectorSlice* vs_lhs, const VectorSlice& vs_rhs, value_type alpha);
/// vs_lhs = pow(vs_rhs,n) 
void pow(VectorSlice* vs_lhs, const VectorSlice& vs_rhs, int n);
/// vs_lhs = pow(alpha,vs_rhs)
void pow(VectorSlice* vs_lhs, value_type alpha, const VectorSlice& vs_rhs);
/// vs_lhs = sqrt(vs_rhs)
void sqrt(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = sin(vs_rhs)
void sin(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = sinh(vs_rhs)
void sinh(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = tan(vs_rhs)
void tan(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);
/// vs_lhs = tanh(vs_rhs)
void tanh(VectorSlice* vs_lhs, const VectorSlice& vs_rhs);

/// v_lhs = abs(vs_rhs)
void abs(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = asin(vs_rhs)
void asin(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = acos(vs_rhs)
void acos(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = atan(vs_rhs)
void atan(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = atan(vs_rhs1/vs_rhs2)
void atan2(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// v_lhs = atan(vs_rhs/alpha)
void atan2(Vector* v_lhs, const VectorSlice& vs_rhs, value_type alpha);
/// v_lhs = atan(alpha/vs_rhs)
void atan2(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs);
/// v_lhs = cos(vs_rhs)
void cos(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = cosh(vs_rhs)
void cosh(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = exp(vs_rhs)
void exp(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = max(vs_rhs1,vs_rhs2)
void max(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// v_lhs = max(alpha,vs_rhs)
void max(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs);
/// v_lhs = min(vs_rhs1,vs_rhs2)
void min(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// v_lhs = mim(alpha,vs_rhs)
void min(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs);
/// v_lhs = pow(vs_rhs1,vs_rhs2)
void pow(Vector* v_lhs, const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// v_lhs = pow(vs_rhs,alpha)
void pow(Vector* v_lhs, const VectorSlice& vs_rhs, value_type alpha);
/// v_lhs = pow(vs_rhs,n) 
void pow(Vector* v_lhs, const VectorSlice& vs_rhs, int n);
/// v_lhs = pow(alpha,vs_rhs)
void pow(Vector* v_lhs, value_type alpha, const VectorSlice& vs_rhs2);
/// v_lhs = sqrt(vs_rhs)
void sqrt(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = sin(vs_rhs)
void sin(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = sinh(vs_rhs)
void sinh(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = tan(vs_rhs)
void tan(Vector* v_lhs, const VectorSlice& vs_rhs);
/// v_lhs = tanh(vs_rhs)
void tanh(Vector* v_lhs, const VectorSlice& vs_rhs);

//		end Elementsize Math Functions
//@}

/** @name {\bf Scalar Returning and Misc VectorSlice Functions}. */

//@{
//		begin Scalar Returning VectorSlice Functions}

/// result = vs_rhs1' * vs_rhs2 (BLAS xDOT)
value_type dot(const VectorSlice& vs_rhs1, const VectorSlice& vs_rhs2);
/// result = max(vs_rhs)
value_type max(const VectorSlice& vs_rhs);
/// result = min(vs_rhs)
value_type min(const VectorSlice& vs_rhs);
/// result = ||vs_rhs||1 (BLAS xASUM)
value_type norm_1(const VectorSlice& vs_rhs);
/// result = ||vs_rhs||2 (BLAS xNRM2)
value_type norm_2(const VectorSlice& vs_rhs);
/// result = ||vs_rhs||infinity (BLAS IxAMAX)
value_type norm_inf(const VectorSlice& vs_rhs);

// Misc. operations

/// swap(vs1, vs2). Swaps the contents of vs1 and vs2
void swap(VectorSlice* vs1, VectorSlice* vs2);

//		end Scalar Returning VectorSlice Functions
//@}

} // end namespace LinAlgPack

//		end Basic Vector Operation Functions
//@}

#endif // VECTOROP_H
