  #ifndef SIMPLECOPERATORSOURCE_H
  #define SIMPLECOPERATORSOURCE_H

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

  class SimpleCOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      /** ctor (no Mp) */
      SimpleCOperatorSource(TSFLinearOperator& S);

      /** virtual destructor */
      virtual ~SimpleCOperatorSource(){;}

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
