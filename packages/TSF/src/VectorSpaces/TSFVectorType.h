#ifndef TSFVECTORTYPE_H
#define TSFVECTORTYPE_H

#include "TSFConfig.h"
#include "TSFSmartPtr.h"
#include "TSFVectorTypeBase.h"
#include "TSFMatrixOperator.h"
#include <string>

namespace TSF
{
	using std::string;

	/** 
	 * Class TSFVectorType is a specification of the type of vector
	 * to be used in a given computation, without reference to the
	 * size or processor distribution of any vector instance. This gives us a way
	 * to specify that we will use, say, Petra or LAPACK before we have sufficient
	 * information to build an actual vector space. 
	 */

	class TSFVectorType
		{
		public: 
			/** empty ctor */
			TSFVectorType();
			/** construct with a pointer to a subtype */
			TSFVectorType(TSFVectorTypeBase* ptr);
	
			/** create a vector space of specified dimension */
			TSFVectorSpace createSpace(int dimension) const ;

			/** create a vector space of the specified dimension, with the specified 
			 * index range living on this processor. */
			TSFVectorSpace createSpace(int dimension, int nLocal, 
																 int firstLocal) const ;

			/** create a vector space of the specified dimension, with the specified 
			 * indices living on this processor. */
			TSFVectorSpace createSpace(int dimension, int nLocal, 
																 const int* localIndices) const ;

			/** create a vector space of the specified dimension, with the specified 
			 * indices owned by this processor and specified ghost indices viewed
			 * but not owned by this processor */
			TSFVectorSpace createSpace(int dimension, int nLocal,
																 const int* localIndices,
																 int nGhost,
																 const int* ghostIndices) const ;
			
			/** Create a matrix useable with this vector type */
			TSFMatrixOperator* createMatrix(const TSFVectorSpace& domain,
																			const TSFVectorSpace& range) const ;

			/** test equality with another vector type */
			bool operator==(const TSFVectorType& other) const ;

			/** write to a stream */
			ostream& print(ostream& os) const ;

			/** return a default choice of linear solver that will work with this type */
			TSFLinearSolver defaultSolver() const ;
			
		private:
			TSFSmartPtr<TSFVectorTypeBase> ptr_;
		};

	/** \relates TSFVectorType */
	inline ostream& operator<<(ostream& os, const TSFVectorType& f)
		{
			return f.print(os);
		}
}

#endif
