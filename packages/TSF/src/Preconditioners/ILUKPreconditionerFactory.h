#ifndef ILUKPRECONDITIONERFACTORY_H
#define ILUKPRECONDITIONERFACTORY_H

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


  /** \ingroup Preconditioner
   * ILUKPreconditionerFactory builds a k-level incomplete factorization
   * preconditioner for a linear operator.
   */

  class ILUKPreconditionerFactory : public TSFPreconditionerFactoryBase
    {
    public:
      /** construct with the number of levels of fill allowed in
       * the incomplete factorization */
      ILUKPreconditionerFactory(int fillLevels, int overlapFill=0);
      /** TUVD */
      virtual ~ILUKPreconditionerFactory(){;}

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& A) const;

      /** write to a string */
      virtual string toString() const;

    private:
      int fillLevels_;
      int overlapFill_;
    };

}


#endif
