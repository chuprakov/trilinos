#ifndef AZTECSOLVER_H
#define AZTECSOLVER_H

#include "TSFDefs.h"
#include "TSFLinearSolverBase.h"
#include "TSFTimeMonitor.h"
#include "TSFHashtable.h"
#include "TSFSmartPtr.h"



#include "AztecOO.h"

#define AZ_ml           101
#define AZ_ml_levels    102
#define AZ_ml_sym       103

namespace TSF
{
  using std::ostream;

  /** \ingroup ConcreteLinearSolvers
   * Aztec solver
   */

  class AZTECSolver : public TSFLinearSolverBase
    {
    public:
      /** construct a AZTEC solver with default options and parameters. */
      AZTECSolver();

      /** construct a AZTEC solver with the given options and parameters. */
      AZTECSolver(const TSFHashtable<int, int>& aztecOptions,
                  const TSFHashtable<int, double>& aztecParameters)
      /** TUVD */
      virtual ~AZTECSolver();

      /**
       * Solve the system with the given RHS, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      virtual bool solve(const TSFLinearOperator& op,
                         const TSFVector& rhs,
                         TSFVector& soln) const ;


    private:
      void setupML(const Epetra_RowMatrix* A);

      /** Aztec options */
      mutable TSFArray<int> options_;

      /** Aztec parameters */
      mutable TSFArray<double> parameters_;

      /** Flag indicating whether we are using ML preconditioning */
      bool useML_;

      /** Number of ML levels to use */
      int mlLevels_;

      /** whether ML should assume symmetric system */
      bool mlSymmetric_;

      /** ML preconditioner object */
      TSFSmartPtr<Epetra_Operator> prec_;

      /** Aztec status */
      TSFArray<double> aztec_status;

      /** Aztec proc_config */
      TSFArray<int> aztec_proc_config;
   };

}


#endif
