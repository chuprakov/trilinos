#ifndef SIMPLEQUADRATIC_H
#define SIMPLEQUADRATIC_H

#include "TSFConfig.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFNonlinearOperator.h"

namespace TSF
{
	

	class SimpleQuadratic : public TSFNonlinearOperatorBase
		{
		public:
			/** build a quadratic f(x) = a x^2 + b x + c  */
			SimpleQuadratic(double a, double b, double c);
			
			/** TUVD */
			virtual ~SimpleQuadratic(){;}

			/** */
			virtual void apply(const TSFVector& in, 
												 TSFVector& out) const ;


			/** */
			virtual TSFLinearOperator derivative(const TSFVector& evalPt) const ;

			/** write to a stream */
			virtual void print(ostream& os) const {os << "(SimpleQuadratic operator)";}
		private:
			TSFReal a_;
			TSFReal b_;
			TSFReal c_;
		};

}

#endif
