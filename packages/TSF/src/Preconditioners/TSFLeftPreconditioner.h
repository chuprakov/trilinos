#ifndef TSFLEFTPRECONDITIONER_H
#define TSFLEFTPRECONDITIONER_H

#include "TSFConfig.h"
#include "TSFVector.h"
#include "TSFPreconditionerBase.h"

namespace TSF
{
	using std::string;
	
	/** \ingroup Preconditioner
	 * Stores a linear operator M1^-1 that will be applied as 
	 * a left preconditioner.
	 */
	
	class TSFLeftPreconditioner : public TSFPreconditionerBase
		{
		public:
			/** empty ctor */
			TSFLeftPreconditioner(const TSFLinearOperator& leftOp);
			/** TUVD */
			virtual ~TSFLeftPreconditioner(){;}

			/** Left preconditioner */
			virtual TSFLinearOperator left() const {return leftOp_;}

			/** tell the world that we have a left preconditioner */
			virtual bool hasLeft() const {return true;}

			/** print to a string */
			virtual string toString() const ;
		private:
			TSFLinearOperator leftOp_;
		};


}

#endif
