#ifndef TSFVECTORTYPEBASE_H
#define TSFVECTORTYPEBASE_H

#include "TSFConfig.h"
#include "TSFVectorSpace.h"
#include "TSFLinearSolver.h"

namespace TSF
{
	class TSFMatrixOperator;

	/** 
	 * Base class for vector types
	 */

	class TSFVectorTypeBase
		{
		public: 
			/** empty ctor */
			TSFVectorTypeBase(){;}
			/** TUVD */
			virtual ~TSFVectorTypeBase(){;}
	
			/** create a vector space */
			virtual TSFVectorSpace createSpace(int dimension) const ;

			/** create a vector space */
			virtual TSFVectorSpace createSpace(int dimension, int nLocal,
																				 int firstLocal) const ;

			/** create a vector space */
			virtual TSFVectorSpace createSpace(int dimension, int nLocal, 
																				 const int* localIndices) const ;	

			/** create a vector space */
			virtual TSFVectorSpace createSpace(int dimension, int nLocal,
																				 const int* localIndices,
																				 int nGhost,
																				 const int* ghostIndices) const ;

			/** Create a matrix useable with this vector type */
			virtual TSFMatrixOperator* createMatrix(const TSFVectorSpace& domain,
																							const TSFVectorSpace& range) const ;

			/** write to stream */
			virtual ostream& print(ostream& os) const = 0 ;

			/** return a default choice of linear solver that will work with this type */
			virtual TSFLinearSolver defaultSolver() const = 0 ;
			
		private:
		};
}

#endif
