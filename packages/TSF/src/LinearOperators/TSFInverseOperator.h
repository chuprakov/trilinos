#ifndef TSFINVERSEOPERATOR_H
#define TSFINVERSEOPERATOR_H

#include "TSFConfig.h"
#include "TSFVectorSpace.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolver.h"

namespace TSF
{	
	

	/** \ingroup LinearOperatorSubtypes
	 * TSFInverseOperator represents the inverse of some other operator. An 
	 * inverse operator object will contain an operator and a solver. The
	 * operator data member is the operator whose inverse this represents. The
	 * solver data member is the solver that will be used in applying the
	 * inverse. If the solver is null, the operator is assumed to have 
	 * self-contained ability to solve systems, as for a dense matrix that
	 * does solves by factoring and backsolves. 
	 */

	class TSFInverseOperator : public TSFLinearOperatorBase
		{
		public:
			/** */
			TSFInverseOperator(const TSFLinearOperator& op, 
												 const TSFLinearSolver& solver = TSFLinearSolver());

			/** the usual virtual dtor */
			virtual ~TSFInverseOperator(){;}

			/** apply this operator, which is applying the inverse of the 
			 * underlying operator. */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** apply the inverse of the inverse, which is the forward operator */
			virtual void applyInverse(const TSFVector& in,
																TSFVector& out) const ;

			/** apply the adjoint of the inverse, which is the inverse
			 * adjoint of the forward operator */
			virtual void applyAdjoint(const TSFVector& in,
																TSFVector& out) const ;

			/** apply the inverse adjoint of the inverse, which is the
			 * adjoint of the forward operator */
			virtual void applyInverseAdjoint(const TSFVector& in,
																TSFVector& out) const ;

			/** construct the inverse */
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
			/** the solver to be used when applying the inverse */
			TSFLinearSolver solver_;
		};
}

#endif
