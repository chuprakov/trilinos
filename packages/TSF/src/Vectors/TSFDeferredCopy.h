#ifndef TSFDEFERREDCOPY_H
#define TSFDEFERREDCOPY_H


#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include <typeinfo>


namespace TSF
{
  class TSFVector;
  class TSFVectorBase;

  /** TSFDeferredCopy holds a copy operation until it can be determined whether
   * space needs to be allocated for the copy */
  class TSFDeferredCopy
    {
    public:
      /** constructor */
      TSFDeferredCopy(const TSFSmartPtr<TSFVectorBase>& ptr);

      /** replace the argument with a clone of this vector */
      void createCopy(TSFVector& vec) const ;

      /** copy this vector into the argument */
      void copyInto(TSFVector& vec) const ;
    private:
      TSFSmartPtr<TSFVectorBase> ptr_;
    };
}
#endif
