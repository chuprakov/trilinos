#ifndef POISSONBOLTZMANN_H
#define POISSONBOLTZMANN_H

#include "TSFDefs.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFNonlinearOperator.h"

namespace TSF
{


  class PoissonBoltzmann : public TSFNonlinearOperatorBase
    {
    public:
      /** */
      PoissonBoltzmann(int n);

      /** TUVD */
      virtual ~PoissonBoltzmann(){;}

      /** */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;


      /** */
      virtual TSFLinearOperator derivative(const TSFVector& evalPt) const ;

      /** write to a stream */
      virtual void print(ostream& os) const {os << "(PoissonBoltzmann operator)";}
    private:
      int n_;
      TSFReal uRight_;
      TSFReal h_;
    };

}

#endif
