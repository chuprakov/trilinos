#ifndef RIGHTBLOCKNSOPERATORSOURCE_H
#define RIGHTBLOCKNSOPERATORSOURCE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSource.h"
#include "TSFOperatorSourceBase.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace SPP
{
  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * TSFOperatorSource builds implementation-specific sources
   * (linear operators)
   * from an abstract specification.
   */

  class RightBlockNSOperatorSource : public TSFOperatorSourceBase
    {
    public:
      /** ctor */
      RightBlockNSOperatorSource();

      /** destructor */
      virtual ~RightBlockNSOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;


      /** write to a string */
      string toString() const ;

    private:
      TSFSmartPtr<TSFOperatorSourceBase> ptr_;
    };

  /** \relates TSFOperatorSource
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const TSFOperatorSource& x);

  /** \relates TSFOperatorSource
   * write to a string
   */
  string toString(const TSFOperatorSource& x);

}


#endif
