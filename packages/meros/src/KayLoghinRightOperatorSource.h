#ifndef KAYLOGHINRIGHTOPERATORSOURCE_H
#define KAYLOGHINRIGHTOPERATORSOURCE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
//#include "KayLoghinSchurFactory.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace SPP
{
  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * KayLoghinRightOperatorSource
   *
   */

  // Need a better name for this class.

  class KayLoghinRightOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      // KayLoghinRightOperatorSource(TSFLinearOperator& S);

      /** ctor (no Mp) */
      KayLoghinRightOperatorSource(TSFLinearOperator& S,
                                   TSFLinearOperator& Fp,
                                   TSFLinearOperator& Ap);

      // /** ctor (with Mp) */
      // KayLoghinRightOperatorSource(TSFLinearOperator& S,
      //                             TSFLinearOperator& Fp,
      //                             TSFLinearOperator& Ap,
      //                             TSFLinearOperator& Mp);

      // /** ctor, we build Fp (no Mp) */
      // KayLoghinRightOperatorSource(TSFLinearOperator& S,
      //                             TSFLinearOperator& Ap);

      // I don't see how we can have an option where we build Fp but they supply Mp
      // No way to distinguish the constructors.

      /** virtual destructor */
      virtual ~KayLoghinRightOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** get pressure Poisson operator */
      TSFLinearOperator getAp() const;

      /** get pressure convection diffusion operator */
      TSFLinearOperator getFp() const;

      /** get mass matrix */
      TSFLinearOperator getMp() const;

      /** write to a string */
      string toString() const ;

    private:
      TSFLinearOperator S_;
      TSFLinearOperator Fp_;
      TSFLinearOperator Ap_;
      TSFLinearOperator Mp_;
    };


}


#endif
