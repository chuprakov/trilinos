#ifndef TSFPRODUCTSPACE_H
#define TSFPRODUCTSPACE_H

#include "TSFDefs.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorSpace.h"
#include <string>
#include "TSFArray.h"

namespace TSF
{
  class TSFVectorBase;
  using std::string;

  using std::ostream;

  /** \ingroup VectorSpaceSubtypes
   * A product space represents a cartesian product of vector spaces
   */

  class TSFProductSpace : public TSFVectorSpaceBase
    {
    public:
      /** */
      TSFProductSpace(const TSFVectorSpace& space);

      /** */
      TSFProductSpace(const TSFVectorSpace& space0,
                      const TSFVectorSpace& space1);
      /** */
      TSFProductSpace(const TSFVectorSpace& space0,
                      const TSFVectorSpace& space1,
                      const TSFVectorSpace& space2);
      /** */
      TSFProductSpace(const TSFVectorSpace& space0,
                      const TSFVectorSpace& space1,
                      const TSFVectorSpace& space2,
                      const TSFVectorSpace& space3);

      /** */
      TSFProductSpace(const TSFArray<TSFVectorSpace>& spaces);

      /** the usual virtual dtor */
      virtual ~TSFProductSpace(){;}

      /** virtual copy ctor */
      virtual TSFVectorSpaceBase* deepCopy() const ;

      /** return dimension of space */
      virtual int dim() const ;

      /** create a vector that is a member of this space */
      virtual TSFVectorBase* createMember(const TSFVectorSpace& handle) const;

      /** test equality */
      virtual bool checkEquality(const TSFVectorSpaceBase* other) const ;

      /** returns true since this is a ProductSpace  */
      virtual bool isProductSpace() const
        {
          return true;
        }

      /** */
      virtual int numBlocks() const ;

      /** */
      virtual void getBlock(int i, const TSFVectorSpace& self,
                            TSFVectorSpace& sub) const ;
      /** write to stream */
      virtual void print(ostream& os) const ;
    protected:
      TSFArray<TSFVectorSpace> blocks_;
    };
}

#endif
