#ifndef TSFSCALEDNONLINEAROPERATOR_H
#define TSFSCALEDNONLINEAROPERATOR_H

#include "TSFConfig.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFNonlinearOperator.h"

namespace TSF
{
	

	/** \ingroup NonlinearOperatorSubtypes
	 * TSFScaledNonlinearOperator is the product of a scalar and a nonlinear
	 * operator.
	 */

	class TSFScaledNonlinearOperator : public TSFNonlinearOperatorBase
		{
		public:
			/** construct with a scale factor and a nonlinear operator */
			TSFScaledNonlinearOperator(const TSFReal& scale,
																 const TSFNonlinearOperator& op); 
										 
			/** the usual virtual dtor */
			virtual ~TSFScaledNonlinearOperator(){;}

			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** get a linear operator representing the derivative */
			virtual TSFLinearOperator derivative(const TSFVector& evalPt) const ;

			
			/** write to a stream */
			virtual void print(ostream& os) const ;

		protected:
			TSFNonlinearOperator op_;
			TSFReal scale_;
		};
}

#endif
