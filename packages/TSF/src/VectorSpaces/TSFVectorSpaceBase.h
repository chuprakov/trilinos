#ifndef TSFVECTORSPACEBASE_H
#define TSFVECTORSPACEBASE_H

#include "TSFConfig.h"
#include <string>

namespace TSF
{
	using std::string;
	using std::ostream;
	class TSFVectorSpace;
	class TSFVectorBase;
	


	/** \ingroup VectorSpaceSubtypes
	 * Base class for vector spaces 
	 */

	class TSFVectorSpaceBase
		{
		public: 
			/** the usual virtual dtor */
			virtual ~TSFVectorSpaceBase(){;}
	
			/** virtual copy ctor */
			virtual TSFVectorSpaceBase* deepCopy() const = 0 ;
	
			/** return dimension of space */
			virtual int dim() const = 0 ;

			/** create a vector that is a member of this space */
			virtual TSFVectorBase* createMember(const TSFVectorSpace& handle) const=0;
			
			/** test equality */
			virtual bool checkEquality(const TSFVectorSpaceBase* other) const ;

			/** */
			virtual int numBlocks() const {return 1;}

			/** */
			virtual void getBlock(int i, const TSFVectorSpace& self,
														TSFVectorSpace& sub) const ;

			/** write to stream */
			virtual ostream& print(ostream& os) const = 0 ;

		protected:
			static int getNewID();
		private:
			static int typeIDCounter_;
		};
}

#endif
