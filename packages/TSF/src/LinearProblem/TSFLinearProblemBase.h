#ifndef TSFLINEARPROBLEMBASE_H
#define TSFLINEARPROBLEMBASE_H

#include "TSFDefs.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"


namespace TSF
{
  using std::ostream;

  /** \ingroup CoreSubtypes
   * TSFLinearProblemBase is the base class for linear problems.
   *
   */

  class TSFLinearProblemBase
    {
    public:
      /** empty ctor only */
      TSFLinearProblemBase(){;}

      /** virtual dtor */
      virtual ~TSFLinearProblemBase(){;}

      /** returns the right-hand side vector b */
      virtual TSFVector getRHS() const ;

      /** returns a RHS vector of a type specified by the input
       * vector space */
      virtual TSFVector getRHS(const TSFVectorSpace& space) const = 0 ;

      /** By default, the solution is not known, so the base class
       * implementation throws an exception. */
      virtual TSFVector getKnownSolution(const TSFVectorSpace& space) const ;

      /** By default, the solution is not known, so the base class
       * implementation throws an exception. */
      virtual TSFVector getKnownSolution() const ;

      /** returns the linear operator A */
      virtual TSFLinearOperator getOperator() const = 0 ;

    private:
    };
}

#endif


