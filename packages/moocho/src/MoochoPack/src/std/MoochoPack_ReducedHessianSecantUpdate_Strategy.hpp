// ///////////////////////////////////////////////////////
// ReducedHessianSecantUpdate_Strategy.h

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_STRATEGY_H
#define REDUCED_HESSIAN_SECANT_UPDATE_STRATEGY_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

///
/** Strategy interface for performing secant updates {abstract}.
 *
 * This interface is used by the class \Ref{ReducedHessianSecantUpdateStd_Step}
 * to actually perform the secant updates.
 */
class ReducedHessianSecantUpdate_Strategy {
public:
	
    ///
	virtual ~ReducedHessianSecantUpdate_Strategy() {}

	///
	/** Perform the secant update.
	 *
	 * The function will update #rHL_k# so that #rHL_k * s_bfgs \approx y_bfgs#.
	 * Note that this post conditions for this function do not strictly require
	 * that the secant property #rHL_k * s_bfgs = y_bfgs# be satisfied.  This
	 * allows for more flexibility in how the update is perform.
	 *
	 * Preconditions\begin{itemize}
	 * \item #s_bfgs->size() == y_bfgs->size() == rHL_k->rows() == rHL_k->cols()# (throws ???)
	 * \end{itemize}
	 *
	 * @param s_bfgs        [in/w] Secant change vector on input.  May be modified as
	 *                      modified as workspace.
	 * @param y_bfgs        [in/w] Secant change vector on input.  May be modified as
 	 *                      modified as workspace.
	 * @param first_update  [in] If true then this is the first update after #rHL# was
	 *                      initialized to identity.  This is information that may be
	 *                      used in order to deliver a beter initial update.
	 * @param out           [out] Output stream journal data is written to.
	 * @param olevel        [in] Output level for printing to #out#
	 * @param algo          [in/out] The rSQPAlgo object.  This object can be queryed for
	 *                      information and also be called to redirect control (in which
	 *                      case this function should probably return false).
	 * @param s             [in/out] rSQPState object.  May be queried or modified if needed.
	 * @param rHL_k         [in/out] The matrix to be updated.  Note that #rHL_k# was already
	 *                      set to #rHL_km1# before this call was made.  Also, #rHL_k# will
	 *                      probably have to support the #MatrixSymSecantUpdateable# interface
	 *                      or an exception will be thrown.
	 * 
	 * @return Returns false if the algorithms path has been redirected through #algo#.
	 * Ohterwise, this function should return true.
	 */
	virtual bool perform_update(
		VectorSlice* s_bfgs, VectorSlice* y_bfgs, bool first_update
		,std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,MatrixWithOp *rHL_k
		) = 0;
	
	///
	/** This function will print a description of the computations and logic used
	 * in the update.
	 */
	virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class ReducedHessianSecantUpdate_Strategy

}  // end namespace ReducedSpaceSQPPack

#endif // REDUCED_HESSIAN_SECANT_UPDATE_STRATEGY_H
