#ifndef TSFPRECONDITIONERFACTORY_H
#define TSFPRECONDITIONERFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFPreconditioner.h"
#include "TSFPreconditionerBase.h"
#include "TSFPreconditionerFactoryBase.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace TSF
{

  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * TSFPreconditionerFactory builds an implementation-specific preconditioner
   * from an abstract specification.
   *
   */

  class TSFPreconditionerFactory
    {
    public:
      /** empty ctor. Constructs a null vector */
      TSFPreconditionerFactory();
      /** assume control of a pointer */
      TSFPreconditionerFactory(TSFPreconditionerFactoryBase* ptr);

      /** create a concrete preconditioner */
      TSFPreconditioner createPreconditioner(const TSFLinearOperator& A) const;

      /** write to a string */
      string toString() const ;

    private:
      TSFSmartPtr<TSFPreconditionerFactoryBase> ptr_;
    };

  /** \relates TSFPreconditionerFactory
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const TSFPreconditionerFactory& x);

  /** \relates TSFPreconditionerFactory
   * write to a string
   */
  string toString(const TSFPreconditionerFactory& x);

}


#endif
