#ifndef TSFADDVECTORNONLINEAROPERATOR_H
#define TSFADDVECTORNONLINEAROPERATOR_H

#include "TSFConfig.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFNonlinearOperator.h"

namespace TSF
{
	

	/** \ingroup NonlinearOperatorSubtypes
	 * TSFAddVectorNonlinearOperator is the sum (or difference) of a nonlinear 
	 * operator and a vector.
	 * Whether this object represents addition or subtraction is indicated
	 * by the subtraction argument to the ctor. 
	 */

	class TSFAddVectorNonlinearOperator : public TSFNonlinearOperatorBase
		{
		public:
			/** construct with an operator, a vector, and a bool to indicate
			 * if addition or substraction is to be performed. */
			TSFAddVectorNonlinearOperator(const TSFNonlinearOperator& op, 
																		const TSFVector& vec, 
																		bool subtraction = false);
										 
			/** the usual virtual dtor */
			virtual ~TSFAddVectorNonlinearOperator(){;}

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
			TSFVector vec_;
			bool subtraction_;
		};
}

#endif
