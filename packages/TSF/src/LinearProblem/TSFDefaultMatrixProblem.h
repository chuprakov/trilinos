#ifndef TSFDEFAULTMATRIXPROBLEM_H
#define TSFDEFAULTMATRIXPROBLEM_H

#include "TSFConfig.h"
#include "TSFMatrixProblem.h"
#include "TSFVector.h"


namespace TSF
{
	using std::ostream;

	/** \ingroup CoreSubtypes
	 * TSFMatrixProblem
	 * 
	 */

	class TSFDefaultMatrixProblem : public TSFMatrixProblem
		{
		public:
			/** empty ctor only */
			TSFDefaultMatrixProblem(TSFMatrixOperator* matrix);

			/** virtual dtor */
			virtual ~TSFDefaultMatrixProblem(){;}

		protected:
			/** fill the matrix with values */
			virtual void fillMatrix() const ; 

			/** fill the RHS vector with values */
			virtual void fillRHS(TSFVector& rhs) const ; 

			/** fill the known solution vector with values. */
			virtual void fillKnownSolution(TSFVector& soln) const ; 

			/** create the update list */
			virtual TSFSmartPtr<TSFArray<int> > formUpdateList() const = 0 ;
			
			/** */
			virtual TSFSmartPtr<TSFArray<int> > formBandwidth() const ;
			
			/* */
			virtual void getRowValues(int row, TSFArray<int>& indices, 
																TSFArray<TSFReal>& values) const = 0;

			/* */
			virtual int getRowBandwidth(int row) const = 0 ;

			/* */
			virtual TSFReal getRHSValue(int row) const = 0 ;


			/* */
			virtual TSFReal getSolutionValue(int row) const = 0 ;
		};
}

#endif
