#ifndef TSFDUMMYPRECONDITIONER_H
#define TSFDUMMYPRECONDITIONER_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFPreconditionerBase.h"

namespace TSF
{
	using std::string;
	
	/** \ingroup Preconditioner
	 * Dummy identity preconditioner
	 */
	
	class TSFDummyPreconditioner : public TSFPreconditionerBase
		{
		public:
			/** empty ctor */
			TSFDummyPreconditioner();
			/** TUVD */
			virtual ~TSFLeftPreconditioner(){;}

			/** Left preconditioner */
			virtual TSFLinearOperator left() const ;

			/** inform the world that the left preconditioner is trivial */
			virtual bool hasLeft() const {return false;}

			/** print to a string */
			virtual string toString() const ;
		private:
		};


}

#endif
