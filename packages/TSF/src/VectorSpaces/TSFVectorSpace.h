#ifndef TSFVECTORSPACE_H
#define TSFVECTORSPACE_H

#include "TSFConfig.h"
#include "TSFVectorSpaceBase.h"
#include "TSFSmartPtr.h"

namespace TSF
{
	class TSFVector;
	class TSFVectorBase;
	using std::string;
	using std::ostream;

	/** \ingroup VectorSpace
	 * TSFVectorSpace represents a vector space, and acts as a factory
	 * class that can build vectors that are elements of that space. 
	 * TSFVectorSpace objects are also used to specify the domain and 
	 * range spaces of an operator. These domain and range spaces can produce
	 * appropriate vectors, or can be used to check the consistency of an 
	 * operation on a given vector. 
	 * 
	 */

	class TSFVectorSpace
		{
		public: 
			/** empty ctor */
			TSFVectorSpace();
			/** construct with a pointer to a subtype */
			TSFVectorSpace(TSFVectorSpaceBase* ptr);
	
			/** test equality between two spaces */
			bool operator==(const TSFVectorSpace& other) const ;
			/** test inequality of two spaces */
			bool operator!=(const TSFVectorSpace& other) const ;

			/** create a copy of this space */
			TSFVectorSpace deepCopy() const ;
	
			/** return the dimension of the space */
			int dim() const ;

			/** create a member vector */
			TSFVectorBase* createMember() const ;

			/** test whether the space contains a given vector */
			bool contains(const TSFVector& vec) const ;

            /** return true if ProductVectorSpace */
            bool isProductSpace() const;

			/** return the number of subblocks. */
			int numBlocks() const ;
			
			/** get the i-th subblock */
			TSFVectorSpace getBlock(int i) const ;

			/** write to stream */
			ostream& print(ostream& os) const ;

			/** access to raw pointer */
			const TSFVectorSpaceBase* ptr() const {return ptr_.get();}
			
		private:
			TSFSmartPtr<TSFVectorSpaceBase> ptr_;
		};

	inline ostream& operator<<(ostream& os, const TSFVectorSpace& space)
		{
			return space.print(os);
		}
}

#endif
