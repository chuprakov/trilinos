#ifndef SCHURFACTORY_H
#define SCHURFACTORY_H

#include "SchurFactoryBase.h"
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

  using std::string;
  using std::ostream;
  using namespace TSF;

  /** \ingroup Preconditioner
   * SchurFactory builds an implementation-specific problem
   * from an abstract specification.
   *
   */

  class SchurFactory
    {
    public:
      /** empty ctor. Constructs a null vector */
      SchurFactory();

      /** assume control of a pointer */
      SchurFactory(SchurFactoryBase* ptr);

      /** get a concrete linear operator for Xinv */
      virtual TSFLinearOperator getSchurInvApprox(const TSFOperatorSource& opSrc) const;

      /** get a concrete linear operator for Xinv */
      //      TSFLinearOperator getSchurInvApprox(const KayLoghinRightOperatorSource& opSrc) const;

      /** write to a string */
      string toString() const ;

    private:
      TSFSmartPtr<SchurFactoryBase> ptr_;
    };

  /** \relates SchurFactory
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const SchurFactory& x);

  /** \relates SchurFactory
   * write to a string
   */
  string toString(const SchurFactory& x);

}


#endif
