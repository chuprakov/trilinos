#ifndef MLSOLVERFACTORY_H
#define MLSOLVERFACTORY_H

#include "TSFDefs.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolverBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFParameterList.h"
#include "TSFParameterListImplem.h"

namespace SPP
{
  using namespace TSF;
  using std::ostream;

  /** \ingroup ConcreteSolverFactories
   *  Sets up default ML solver.
   *  createSolver makes an AZTECSolver object with ML preconditioner.
   *  If symmetric, CG preconditioned with ML
   *  else, GMRES preconditioned with ML.
   */

  class MLSolverFactory
    {
    public:
      /** */
      MLSolverFactory(bool isSymmetric);

      /** virtual destructor */
      virtual ~MLSolverFactory();


      virtual TSFLinearSolver createSolver(const TSFLinearOperator& op) const;


    protected:

    private:
      bool isSymmetric_;

    };

}


#endif
