#ifndef TSFNONLINEAROPERATORBASE_H
#define TSFNONLINEAROPERATORBASE_H

#include "TSFConfig.h"
#include "TSFLinearOperator.h"
#include "TSFVectorSpace.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVector.h"
#include "TSFVectorBase.h"

namespace TSF
{
		
	

	
	/** \ingroup NonlinearOperator 
	 * Base class for nonlinear operator objects.
	 * 
	 */
	
	class TSFNonlinearOperatorBase
		{
		public:
			/** empty ctor constructs a null nonlinear operator. 
			 * This is primarily for
			 * use with templated container classes. */
			TSFNonlinearOperatorBase() {;}

			TSFNonlinearOperatorBase(const TSFVectorSpace& domain,
															 const TSFVectorSpace& range)
				: domain_(domain), range_(range) {;}

			/** the usual virtual dtor */
			virtual ~TSFNonlinearOperatorBase(){;}

			/** return domain space */
			const TSFVectorSpace& domain() const {return domain_;}
			/** return range space */
			const TSFVectorSpace& range() const {return range_;}

			/** apply the operator to a vector. */
			virtual void apply(const TSFVector& arg, TSFVector& out) const = 0 ;

			/** get a linear operator representing the derivative at an eval point */
			virtual TSFLinearOperator derivative(const TSFVector& evalPt) const ;

			/** write to a stream */
			virtual void print(ostream& os) const {} ;
		private:
			
			TSFVectorSpace domain_;
			TSFVectorSpace range_;
		};



}

#endif
