#ifndef DENSESERIALVECTORTYPE_H
#define DENSESERIALVECTORTYPE_H

#include "TSFConfig.h"
#include "TSFVectorTypeBase.h"
#include "DenseSerialVectorSpace.h"
#include "TSFLinearSolver.h"

namespace TSF
{
	/** \ingroup VectorSpace
	 *
	 */

	class DenseSerialVectorType : public TSFVectorTypeBase
		{
		public: 
			/** empty ctor */
			DenseSerialVectorType(){;}
			/** TUVD */
			virtual ~DenseSerialVectorType(){;}
	
			/** create a vector space */
			virtual TSFVectorSpace createSpace(int dimension) const ;

			/** create a vector space */
			virtual TSFVectorSpace createSpace(int dimension, int nLocal,
																				 int firstLocal) const ;

			/** create a vector space */
			virtual TSFVectorSpace createSpace(int dimension, int nLocal, 
																				 const int* localIndices) const ;

			/** return a LAPACKGeneralMatrix */
			TSFMatrixOperator* createMatrix(const TSFVectorSpace& domain,
																			const TSFVectorSpace& range) const ;

			
			/** write to a stream */
			ostream& print(ostream& os) const {return os << "DenseSerialVectorType";}

			/** return a DirectSolver */
			virtual TSFLinearSolver defaultSolver() const ;
		private:
		};
}


#endif
