#ifndef TSFIDENTITYOPERATOR_H
#define TSFIDENTITYOPERATOR_H

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
	 * TSFIdentityOperator is the identity ("I") operator on a vector space. 
	 * The action of this operator, its adjoint, or its inverse, is to
	 * return the input vector. 
	 */

	class TSFIdentityOperator : public TSFLinearOperatorBase
		{
		public:
			/** The domain and range spaces for an identity operator
			 * are equivalent, so the ctor needs only a single space */
			TSFIdentityOperator(const TSFVectorSpace& space);
										 
			/* the usual virtual dtor */
			virtual ~TSFIdentityOperator(){;}

			/** apply returns the input vector */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;

			/** applyAdjoint returns the input vector */
			virtual void applyAdjoint(const TSFVector& in, 
																TSFVector& out) const ;

			/** applyInverse returns the input vector */
			virtual void applyInverse(const TSFVector& in, 
																TSFVector& out) const ;

		protected:
		};
}

#endif
