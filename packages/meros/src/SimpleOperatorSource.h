  #ifndef SIMPLEOPERATORSOURCE_H
  #define SIMPLEOPERATORSOURCE_H

  #include "TSFDefs.h"
  #include "TSFSmartPtr.h"
  #include "TSFLinearOperator.h"
  #include "TSFOperatorSource.h"
  #include "RightBlockNSOperatorSource.h"
  #include "TSFArray.h"
  #include <iostream>
  #include <string>

  namespace SPP
  {
  using namespace TSF;
  using std::string;
  using std::ostream;

  // Need a better name for this class.

  class SimpleOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      /** ctor (no Mp) */
      SimpleOperatorSource(TSFLinearOperator& S);

      /** virtual destructor */
      virtual ~SimpleOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** get pressure Poisson operator */
      TSFLinearOperator getDinv() const;

      /** write to a string */
      string toString() const ;

    private:
      TSFLinearOperator S_;
      mutable TSFLinearOperator Dinv_;
      mutable bool hasDinv_;
    };
   }


 #endif
