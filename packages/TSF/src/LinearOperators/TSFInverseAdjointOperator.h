#ifndef TSFINVERSEADJOINTOPERATOR_H
#define TSFINVERSEADJOINTOPERATOR_H

#include "TSFDefs.h"
#include "TSFVectorSpace.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolver.h"

namespace TSF
{


  /** \ingroup LinearOperatorSubtypes
   * TSFInverseAdjointOperator
   */

  class TSFInverseAdjointOperator : public TSFLinearOperatorBase
    {
    public:
      /** */
      TSFInverseAdjointOperator(const TSFLinearOperator& op,
                                const TSFLinearSolver& solver = TSFLinearSolver());

      /** the usual virtual dtor */
      virtual ~TSFInverseAdjointOperator(){;}

      /**  */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;

      /**  */
      virtual void applyInverse(const TSFVector& in,
                                TSFVector& out) const ;

      /**  */
      virtual void applyInverse(const TSFLinearSolver& solver,
                                const TSFVector& in,
                                TSFVector& out) const ;

      /**  */
      virtual void applyAdjoint(const TSFVector& in,
                                TSFVector& out) const ;

      /**  */
      virtual void applyInverseAdjoint(const TSFVector& in,
                                       TSFVector& out) const ;
      /**  */
      virtual void applyInverseAdjoint(const TSFLinearSolver& solver,
                                       const TSFVector& in,
                                       TSFVector& out) const ;

      /**  */
      virtual void getInverse(const TSFLinearSolver& /* solver */,
                              const TSFLinearOperator& /* self */,
                              TSFLinearOperator& inv) const ;

      /** */
      virtual void getInverse(const TSFLinearOperator& /* self */,
                              TSFLinearOperator& inv) const ;

      /** */
      virtual void getAdjoint(const TSFLinearOperator& /* self */,
                              TSFLinearOperator& adj) const ;

      /** */
      virtual void getInverseAdjoint(const TSFLinearOperator& /* self */,
                                     TSFLinearOperator& invAdj) const ;

      /** */
      virtual void getInverseAdjoint(const TSFLinearSolver& /* solver */,
                                     const TSFLinearOperator& /* self */,
                                     TSFLinearOperator& invAdj) const ;
    protected:
      /** the forward operator */
      TSFLinearOperator op_;
      /** the solver  */
      TSFLinearSolver solver_;
    };
}

#endif
