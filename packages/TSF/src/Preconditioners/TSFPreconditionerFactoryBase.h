#ifndef TSFPRECONDITIONERFACTORYBASE_H
#define TSFPRECONDITIONERFACTORYBASE_H

#include "TSFConfig.h"
#include "TSFSmartPtr.h"
#include "TSFPreconditioner.h"
#include "TSFPreconditionerBase.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace TSF
{
	
	using std::string;
	using std::ostream;
	

	/** \ingroup Preconditioner 
	 * Base class for preconditioner factories
	 */

	class TSFPreconditionerFactoryBase
		{
		public: 
			/** empty ctor */
			TSFPreconditionerFactoryBase();
			/** TUVD */
			virtual ~TSFPreconditionerFactoryBase();

			/** create a concrete preconditioner */
			virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& A) const = 0 ;

			/** write to a string */
			virtual string toString() const = 0 ;

		private:
		};

}


#endif
