#ifndef TSFADJOINTOPERATOR_H
#define TSFADJOINTOPERATOR_H

#include "TSFConfig.h"
#include "TSFVectorSpace.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolver.h"

namespace TSF
{	
	

	/** \ingroup LinearOperatorSubtypes
	 * TSFInverseOperator represents the adjoint of some other operator. 
	 */

	class TSFAdjointOperator : public TSFLinearOperatorBase
		{
		public:
			/** */
			TSFAdjointOperator(const TSFLinearOperator& op);

			/* the usual virtual dtor */
			virtual ~TSFAdjointOperator(){;}

			/** apply this operator, which is to apply the adjoint of the underlying
			 * operator */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** */
			virtual void applyInverse(const TSFVector& in,
																TSFVector& out) const ;
			/** */
			virtual void applyInverse(const TSFLinearSolver& solver,
																const TSFVector& in,
																TSFVector& out) const ;
			/** */
			virtual void applyAdjoint(const TSFVector& in,
																TSFVector& out) const ;

			/** */
			virtual void applyInverseAdjoint(const TSFVector& in,
																			 TSFVector& out) const ;	
			/** */
			virtual void applyInverseAdjoint(const TSFLinearSolver& solver,
																			 const TSFVector& in,
																			 TSFVector& out) const ;
			
			/** */
			virtual void getInverse(const TSFLinearSolver& solver,
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
			virtual void getInverseAdjoint(const TSFLinearSolver& solver,
																		 const TSFLinearOperator& /* self */,
																		 TSFLinearOperator& invAdj) const;
			/** */
			virtual bool isMatrixOperator() const;

			/** */
			TSFLinearOperator op() const { return op_; }

		protected:
			/** the operator */
			TSFLinearOperator op_;
		};
}

#endif
