/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#include "TSFIfpackOperator.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFEpetraVector.hpp"

using namespace TSFExtended;
using namespace Teuchos;

IfpackOperator::IfpackOperator(const EpetraMatrix* A,
                               int fillLevels,
                               int overlapFill,
                               double relaxationValue,
                               double relativeThreshold,
                               double absoluteThreshold)
  : precondGraph_(),
    precond_(),
    domain_(A->domain()),
    range_(A->range())
{
  const Epetra_CrsMatrix* matrix = A->crsMatrix();

  const Epetra_CrsGraph& matrixGraph = matrix->Graph();
			
  Ifpack_IlukGraph* precondGraph 
    = new Ifpack_IlukGraph(matrixGraph, fillLevels, overlapFill);

  int ierr = precondGraph->ConstructFilledGraph();

  TEST_FOR_EXCEPTION(ierr != 0, runtime_error,
                     "IfpackOperator ctor: "
                     "precondGraph->ConstructFilledGraph() failed with ierr="
                     << ierr);

  Ifpack_CrsRiluk* precond = new Ifpack_CrsRiluk(*precondGraph);

  precond->SetRelaxValue(relaxationValue);
  precond->SetRelativeThreshold(relativeThreshold);
  precond->SetAbsoluteThreshold(absoluteThreshold);

  ierr = precond->InitValues(*matrix);

  TEST_FOR_EXCEPTION(ierr != 0, runtime_error,
                     "IfpackOperator ctor: "
                     "precond->InitValues() failed with ierr="
                     << ierr);

  ierr = precond->Factor();

  TEST_FOR_EXCEPTION(ierr != 0, runtime_error,
                     "IfpackOperator ctor: "
                     "precond->Factor() failed with ierr="
                     << ierr);
}


void IfpackOperator::apply(const TSFCore::ETransp            M_trans,
                           const TSFCore::Vector<double>    &x,
                           TSFCore::Vector<double>          *y,
                           const double            alpha,
                           const double            beta
                           ) const
{
  /* grab the epetra vector objects underlying the input and output vectors */
  const EpetraVector* epIn = dynamic_cast<const EpetraVector*>(&x);
  TEST_FOR_EXCEPTION(epIn == 0, runtime_error,
                     "IfpackOperator apply: input vector is "
                     "not an EpetraVector");

  const Epetra_Vector* in = epIn->epetra_vec().get();

  EpetraVector* epy = dynamic_cast<EpetraVector*>(y);
  TEST_FOR_EXCEPTION(epy == 0, runtime_error,
                     "IfpackOperator apply: output vector is "
                     "not an EpetraVector");

  Epetra_Vector* yy = epy->epetra_vec().get();

  /* if beta != 0, we have to create a temporary vector to hold the
   * intermediate result beta*y. */
  RefCountPtr<TSFCore::Vector<double> > tsfOut;
  Epetra_Vector* tmp;

  if (beta!=0.0) /* we need to have storage for the result of op*x */
    {
      tsfOut = y->space()->createMember();
      EpetraVector* ep = dynamic_cast<EpetraVector*>(tsfOut.get());
      tmp = ep->epetra_vec().get();
    }
  else /* we overwrite y with the application of the op */
    {
      tmp = yy;
    }

  /* ifpack's solve is logically const but declared non-const because
   *  internal data changes. So, do a const_cast. */
  Ifpack_CrsRiluk* p = const_cast<Ifpack_CrsRiluk*>(precond_.get());

  int ierr;

  /* do the solve (or transpose solve) */
  if (M_trans == TSFCore::NOTRANS)
    {
      ierr = p->Solve(false, *in, *tmp);
    }
  else
    {
      ierr = p->Solve(true, *in, *tmp);
    }

  /* if necessary, add beta*y */
  if (beta != 0.0)
    {
      /* Compute yy = alpha*tmp + beta*yy */
      yy->Update(alpha, *tmp, beta);
    }  
  else if (alpha != 1.0) /* compute yy = alpha*tmp */
    {
      /* Because beta != 0.0, an earlier conditional has set tmp=yy. 
       * We can therefore do the computation on tmp and it will 
       * modify the contents of yy as a side effect. */
      tmp->Scale(alpha);
    }

  /* At this point, the contents of y should be yy. We are done. */
}
  
