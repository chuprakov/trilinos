#ifndef DENSESERIALVECTORSPACE_H
#define DENSESERIALVECTORSPACE_H

#include "TSFVectorSpaceBase.h"

namespace TSF
{
	
	using std::string;

	/** \ingroup DenseSerial
	 * DenseSerialVectorSpace is a vector space that generates TSFSerialVector
	 * objects. 
	 */
	class DenseSerialVectorSpace : public TSFVectorSpaceBase
		{
		public:
			/** construct a DenseSerialVectorSpace of the given dimension */
			DenseSerialVectorSpace(int dim) : dim_(dim) {;}

			/** virtual dtor */
			virtual ~DenseSerialVectorSpace(){;}

			/** deep copy */
			virtual TSFVectorSpaceBase* deepCopy() const ;

			/** return dimension */
			virtual int dim() const {return dim_;}

			/** factory method */
			virtual TSFVectorBase* createMember(const TSFVectorSpace& handle) const ;

			/** write to stream */
			virtual ostream& print(ostream& os) const ;

		private:
			/* dimension of the space */
			int dim_;
		};
}
#endif
