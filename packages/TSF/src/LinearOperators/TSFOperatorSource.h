#ifndef TSFOPERATORSOURCE_H
#define TSFOPERATORSOURCE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSourceBase.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace TSF
{

  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * TSFOperatorSource builds implementation-specific sources
   * (linear operators)
   * from an abstract specification.
   */
  // VEH

  class TSFOperatorSource
    {
    public:
      /** empty ctor. Constructs a null vector */
      TSFOperatorSource();

      /** assume control of a pointer */
      TSFOperatorSource(TSFOperatorSourceBase* ptr);

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** write to a string */
      string toString() const ;

      /** read-only pointer access */
      const TSFOperatorSourceBase* ptr() const {return &(*ptr_);}
      /** read-write pointer access */
      TSFOperatorSourceBase* ptr() {return &(*ptr_);}


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
