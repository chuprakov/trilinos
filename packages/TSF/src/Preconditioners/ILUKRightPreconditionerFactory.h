#ifndef ILUKRIGHTPRECONDITIONERFACTORY_H
#define ILUKRIGHTPRECONDITIONERFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFPreconditionerFactoryBase.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace TSF
{

  using std::string;
  using std::ostream;

  // VEH
  /** \ingroup Preconditioner
   * ILUKRightPreconditionerFactory builds a k-level incomplete factorization
   * (right) preconditioner for a linear operator.
   */

  class ILUKRightPreconditionerFactory : public TSFPreconditionerFactoryBase
    {
    public:
      /** construct with the number of levels of fill allowed in
       * the incomplete factorization */
      ILUKRightPreconditionerFactory(int fillLevels, int overlapFill=0);
      /** TUVD */
      virtual ~ILUKRightPreconditionerFactory(){;}

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& A) const;

      /** create a concrete preconditioner from an OperatorSource*/
      virtual TSFPreconditioner createPreconditioner(const TSFOperatorSource& S) const;

      /** write to a string */
      virtual string toString() const;

    private:
      int fillLevels_;
      int overlapFill_;
    };

}


#endif
