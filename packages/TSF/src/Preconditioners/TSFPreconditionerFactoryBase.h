#ifndef TSFPRECONDITIONERFACTORYBASE_H
#define TSFPRECONDITIONERFACTORYBASE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFPreconditioner.h"
#include "TSFPreconditionerBase.h"
#include "TSFOperatorSource.h"
#include "TSFOperatorSourceBase.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace TSF
{

  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * Base class for preconditioner factories
   */

  class TSFPreconditionerFactoryBase
    {
    public:
      /** empty ctor */
      TSFPreconditionerFactoryBase();
      /** TUVD */
      virtual ~TSFPreconditionerFactoryBase();

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& A) const = 0 ;

      /* VEH Create a concrete preconditioner from an OperatorSource*/
      virtual TSFPreconditioner createPreconditioner(const TSFOperatorSource& S) const = 0 ;

      /** write to a string */
      virtual string toString() const = 0 ;

    private:
    };

}


#endif
