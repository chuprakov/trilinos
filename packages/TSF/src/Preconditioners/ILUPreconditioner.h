#ifndef ILUPRECONDITIONER_H
#define ILUPRECONDITIONER_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFLinearOperator.h"

namespace TSF
{
	using std::ostream;
	
	/** \ingroup LinearOperator
	 * Base class for preconditioners. A general preconditioner object
	 * is split into a left preconditioner M1^-1 and a right 
	 * preconditioner M2^-1. To solve A x = b, we define the auxiliary 
	 * system M2^-1 y = x, and solve M1^-1 A M2^-1 y = M1^-1 b to obtain y. 
	 * Having y, we can quickly recover x by applying M2^-1 to y. 
	 */
	
	class ILUPreconditioner : public TSFPreconditionerBase
		{
		public:
			/** empty ctor */
			ILUPreconditioner(){;}

			/** construct with approximate L and U */
			ILUPreconditioner(const TSFLinearOperator& approxL,
												const TSFLinearOperator& approxU);
			/** TUVC */
			virtual ~ILUPreconditioner(){;}

			/** Left preconditioner */
			virtual TSFLinearOperator left() const ;

			/** Right preconditioner. By default, this returns the identity
			 * operator */
			virtual TSFLinearOperator right() const ;

			/** print to a stream */
			virtual void print(ostream& os) const = 0 ;
		private:
		};


}

#endif
