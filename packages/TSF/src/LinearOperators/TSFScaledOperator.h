#ifndef TSFSCALEDOPERATOR_H
#define TSFSCALEDOPERATOR_H

#include "TSFConfig.h"
#include "TSFVectorSpace.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"


namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 * TSFScaledOperator is a TSFLinearOperator formed by multiplying 
	 * another TSFLinearOperator by a constant scaling factor. 
	 */

	class TSFScaledOperator : public TSFLinearOperatorBase
		{
		public:
			/** construct with a pair of operators and a boolean to indicate
			 * if addition or substraction is to be performed. */
			TSFScaledOperator(const TSFLinearOperator& op,
												const TSFReal& scale);
										 
			/** the usual virtual dtor */
			virtual ~TSFScaledOperator(){;}

			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** */
			virtual void applyAdjoint(const TSFVector& in, 
												 TSFVector& out) const ;

			/** */
			virtual void applyInverse(const TSFVector& in, 
												 TSFVector& out) const ;

		protected:
			TSFLinearOperator op_;
			TSFReal scale_;
		};
}

#endif
