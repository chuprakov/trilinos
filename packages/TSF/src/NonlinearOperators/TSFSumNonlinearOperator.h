#ifndef TSFSUMNONLINEAROPERATOR_H
#define TSFSUMNONLINEAROPERATOR_H

#include "TSFConfig.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFNonlinearOperator.h"

namespace TSF
{
	

	/** \ingroup NonlinearOperatorSubtypes
	 * TSFSumNonlinearOperator is the sum (or difference) of two 
	 * nonlinear operators.
	 * Whether this object represents addition or subtraction is indicated
	 * by the subtraction argument to the ctor. 
	 */

	class TSFSumNonlinearOperator : public TSFNonlinearOperatorBase
		{
		public:
			/** construct with a pair of operators and a boolean to indicate
			 * if addition or substraction is to be performed. */
			TSFSumNonlinearOperator(const TSFNonlinearOperator& left, 
															const TSFNonlinearOperator& right, 
															bool subtraction = false);
										 
			/** the usual virtual dtor */
			virtual ~TSFSumNonlinearOperator(){;}

			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** get a linear operator representing the derivative */
			virtual TSFLinearOperator derivative(const TSFVector& evalPt) const ;

			
			/** write to a stream */
			virtual void print(ostream& os) const ;

		protected:
			TSFNonlinearOperator left_;
			TSFNonlinearOperator right_;
			bool subtraction_;
		};
}

#endif
