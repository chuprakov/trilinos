#ifndef TSFMATRIXFREEPROBLEM_H
#define TSFMATRIXFREEPROBLEM_H

#include "TSFConfig.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"


namespace TSF
{
	using std::ostream;

	/** \ingroup CoreSubtypes
	 * TSFMatrixProblem
	 * 
	 */

	class TSFMatrixFreeProblem : public TSFLinearProblemBase
		{
		public:
			/** empty ctor only */
			TSFMatrixFreeProblem();

			/** virtual dtor */
			virtual ~TSFMatrixFreeProblem(){;}

			/** returns a RHS vector of a type specified by the input
			 * vector space */
			virtual TSFVector getRHS(const TSFVectorSpace& space) const ;

			/** returns a solution vector of a type specified by the input
			 * vector space */
			virtual TSFVector getKnownSolution(const TSFVectorSpace& space) const ;

			/** */

		protected:
			virtual TSFVector formRHSVector(const TSFVectorSpace& space);

			virtual TSFVector formSolnVector(const TSFVectorSpace& space);

			virtual TSFReal getRHSRowValue(int row, 
																		 const TSFVectorSpace& space) const = 0 ;

			virtual TSFReal getSolutionRowValue(int row, 
																					const TSFVectorSpace& space) const  ;

			virtual TSFReal multiplyRow(
		};
}

#endif
