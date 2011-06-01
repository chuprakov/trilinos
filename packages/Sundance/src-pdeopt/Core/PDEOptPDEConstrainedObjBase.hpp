#ifndef PDEOPT_PDECONSTRAINEDOBJBASE_H
#define PDEOPT_PDECONSTRAINEDOBJBASE_H

#include "PlayaObjectiveBase.hpp"
#include "SundanceExpr.hpp"
#include "SundanceFunctional.hpp"

namespace Sundance
{
using namespace Playa;


/**
 * PDEConstrainedObj is a base class for objective functions of the
 * reduced-space variable where the constraint is a PDE in the
 * state variables. 
 * Objects of this type are suitable for use in adjoint gradient methods.
 *
 * One constructs such an objective function by giving it a Lagrangian in the
 * form of a %Sundance Functional, along with a specification of which functions
 * are to be regarded as the adjoint, state, and design variables. 
 */
class PDEConstrainedObjBase : public ObjectiveBase
{
public:
      
  /** Constructor */
  PDEConstrainedObjBase(
    const Functional& lagrangian,
    const Array<Expr>& stateVarVals,
    const Array<Expr>& adjointVarVals,
    const Expr& designVarVal,
    int verb=0);

  /** virtual dtor */
  virtual ~PDEConstrainedObjBase(){;}

  /** evaluate objective function and gradient */
  void  evalGrad(const Vector<double>& x, double& f, 
    Vector<double>& grad) const ;

  /** evaluate objective function without gradient. */
  void eval(const Vector<double>& x, double& f) const ;

  /** return an initial guess for the design vector  */
  Vector<double> getInit() const ;

  /** Hook for anything that needs to be done between the solution
   * of the state equations and the evaluation of the functional or
   * solution of adjoints. 
   *
   * Default is a no-op. */
  virtual void statePostprocCallback() const {;}

  /** */
  const Mesh& mesh() const {return Lagrangian_.mesh();}

  /** */
  int numFuncEvals() const {return numFuncEvals_;}

  /** Solve the state equations, followed by postprocessing.
   * At the end of this call, the system is ready for evaluation of
   * the objective function or solution of the adjoint equations. */
  virtual void solveState(const Vector<double>& x) const = 0 ;

  /** Solve the state equations, then do postprocessing,
   * then finally the adjoint equations in <b>reverse</b> order. 
   * At the end of this call, the system is ready for evaluation
   * of the objective function and its gradient. */
  virtual void solveStateAndAdjoint(const Vector<double>& x) const = 0 ;

  /** Set up the state and adjoint equations. This is left to the derived
   * class, because we can't know at this level whether the equations
   * are linear or nonlinear */
  virtual void initEquations(
    const Array<Expr>& stateVars,
    const Array<Expr>& adjointVars,
    const Array<Array<Expr> >& fixedVarsInStateEqns,
    const Array<Array<Expr> >& fixedVarsInStateEqnsVals,
    const Array<Array<Expr> >& fixedVarsInAdjointEqns,
    const Array<Array<Expr> >& fixedVarsInAdjointEqnsVals
    ) = 0 ;

  /** Set the scale of the Hessian */
  void setHScale(const double& H) {invHScale_ = 1.0/H;}

  /** return an initial approximation to the scale for the 
   * inverse of the Hessian */
  double getInvHScale() const {return invHScale_;}
  
protected:  
  /** */
  void init(
    const Array<Expr>& stateVars,
    const Array<Expr>& adjointVars,
    const Expr& designVar);

  /** Access to the design variable */
  Expr& designVarVal() const
    {return designVarVal_;}

  /** Access to the state variables */
  Expr& stateVarVals(int i) const
    {return stateVarVals_[i];}

  /** Access to the adjoint variables */
  Expr& adjointVarVals(int i) const
    {return adjointVarVals_[i];}

  const Functional& Lagrangian() const {return Lagrangian_;}
      
private:

  Functional Lagrangian_;

  mutable Expr designVarVal_;

  mutable Array<Expr> stateVarVals_;

  mutable Array<Expr> adjointVarVals_;

  FunctionalEvaluator fEval_;

  mutable int numFuncEvals_;

  mutable int numGradEvals_;

  double invHScale_;
};

}

#endif
