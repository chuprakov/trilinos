#ifndef SCHURFACTORYBASE_H
#define SCHURFACTORYBASE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSource.h"
#include "TSFOperatorSourceBase.h"
//#include "RightBlockNSOperatorSource.h"
//#include "KayLoghinRightOperatorSource.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace SPP
{

  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * Base class for SchurFactory factories
   */

  class SchurFactoryBase
    {
    public:
      /** empty ctor. Constructs a null vector */
      SchurFactoryBase();

      /** virtual destructor */
      virtual ~SchurFactoryBase();

      /** get a concrete linear operator */
      virtual TSFLinearOperator getSchurInvApprox(const TSFOperatorSource& opSrc) const = 0;

      /** get a concrete linear operator */
      // virtual TSFLinearOperator getSchurInvApprox(const KayLoghinRightOperatorSource& opSrc) const = 0;

      /** write to a string */
      virtual string toString() const = 0;

    private:
    };

}


#endif
