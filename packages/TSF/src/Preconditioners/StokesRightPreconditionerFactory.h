#ifndef STOKESRIGHTPRECONDITIONERFACTORY_H
#define STOKESRIGHTPRECONDITIONERFACTORY_H

#include "TSFConfig.h"
#include "TSFSmartPtr.h"
#include "TSFPreconditionerFactoryBase.h"
#include "TSFLinearSolver.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace TSF
{
	
	using std::string;
	using std::ostream;
	

	/** \ingroup Preconditioner 
	 * 
	 */

	class StokesRightPreconditionerFactory : public TSFPreconditionerFactoryBase
		{
		public: 
			/** */
			StokesRightPreconditionerFactory(const TSFLinearSolver& innerSolver);
			/** TUVD */
			virtual ~StokesRightPreconditionerFactory(){;}

			/** create a concrete preconditioner */
			virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& A) const;

			/** write to a string */
			virtual string toString() const;

		private:
			TSFLinearSolver innerSolver_;
		};

}


#endif
