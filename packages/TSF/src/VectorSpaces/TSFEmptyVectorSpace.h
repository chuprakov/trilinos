#ifndef TSFEMPTYVECTORSPACE_H
#define TSFEMPTYVECTORSPACE_H

#include "TSFVectorSpaceBase.h"

namespace TSF
{
  using std::string;

  /** \ingroup VectorSpaceSubtypes
   * empty vector space, used only to initialize empty vector ctors.
   */

  class TSFEmptyVectorSpace : public TSFVectorSpaceBase
    {
    public:
      /** ctor */
      TSFEmptyVectorSpace(){;}
      /** the usual virtual dtor */
      virtual ~TSFEmptyVectorSpace(){;}

      /** virtual copy ctor */
      virtual TSFVectorSpaceBase* deepCopy() const ;

      /** return dimension = 0 */
      virtual int dim() const {return 0;}

      /** createMember() throws an exception */
      virtual TSFVectorBase* createMember(const TSFVectorSpace& handle) const ;

      /** write to stream */
      virtual void print(ostream& os) const {os << "TSFEmptyVectorSpace";}

    };
}

#endif
