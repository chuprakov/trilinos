#ifndef TSFSUMOPERATOR_H
#define TSFSUMOPERATOR_H

#include "TSFConfig.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"

namespace TSF
{
	

	/** \ingroup LinearOperatorSubtypes
	 * TSFSumOperator is the sum (or difference) of two linear operators.
	 * Whether this object represents addition or subtraction is indicated
	 * by the subtraction argument to the ctor. 
	 */

	class TSFSumOperator : public TSFLinearOperatorBase
		{
		public:
			/** construct with a pair of operators and a boolean to indicate
			 * if addition or substraction is to be performed. */
			TSFSumOperator(const TSFLinearOperator& left, 
										 const TSFLinearOperator& right, 
										 bool subtraction = false);
										 
			/** the usual virtual dtor */
			virtual ~TSFSumOperator(){;}

			/** return domain space */
			const TSFVectorSpace& domain() const {return left_.domain();}
			/** return range space */
			const TSFVectorSpace& range() const {return left_.range();}

			/** apply operator to a vector in the domain space, returning a vector
			 * in the range space */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

		protected:
			TSFLinearOperator left_;
			TSFLinearOperator right_;
			bool subtraction_;
		};
}

#endif
