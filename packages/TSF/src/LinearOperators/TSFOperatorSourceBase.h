#ifndef TSFOPERATORSOURCEBASE_H
#define TSFOPERATORSOURCEBASE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace TSF
{

  using std::string;
  using std::ostream;


  /** \ingroup LinearOperator
   * Base class for TSFOperatorSource factories
   */
  // VEH

  class TSFOperatorSourceBase
    {
    public:
      /** empty ctor. Constructs a null vector */
      TSFOperatorSourceBase();

      /** virtual destructor */
      virtual ~TSFOperatorSourceBase();

      /** get a concrete linear operator */
      virtual TSFLinearOperator getOp() const = 0;

      /** write to a string */
      virtual string toString() const = 0;

    private:
    };

}


#endif
