// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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
// ***********************************************************************
// @HEADER

#ifndef MEROS_PCD_PRECONDITIONERFACTORY_H
#define MEROS_PCD_PRECONDITIONERFACTORY_H

/*! \file Meros_PCDPreconditionerFactory.h
 *  \brief Factory for building pressure convection-diffusion block 
 *         preconditioners.
 */

#include "Teuchos_ParameterListAcceptor.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Meros_PCDOperatorSource.h"


namespace Meros
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using namespace Thyra;

  /*! \ingroup PreconditionerFactory
   *
   *  \brief Factory for building pressure convection-diffusion
   *  style block preconditioner. This class of preconditioners were
   *  originally proposed by Kay, Loghin, and Wathen (ref) and
   *  Silvester, Elman, Kay, and Wathen (ref).
   *  
   *  Meros 1.0 currently implements the PCD preconditioner,
   *  a.k.a. Fp preconditioner.
   *  
   *  The LDU factors of a saddle point system are given as follows:
   *  
   *  \f$
   *  \left[ \begin{array}{cc} A & B^T \\ B & C \end{array} \right]
   *  = \left[ \begin{array}{cc} I & \\ BF^{-1} & I \end{array} \right]
   *    \left[ \begin{array}{cc} F & \\  & -S \end{array} \right]
   *    \left[ \begin{array}{cc} I & F^{-1} B^T  \\  & I \end{array} \right],
   *  \f$
   *  
   *  where \f$S\f$ is the Schur complement \f$S = B F^{-1} B^T - C\f$.
   *  A pressure convection-diffusion style preconditioner is then given by
   *  
   *  \f$
   *  P^{-1} =
   *    \left[ \begin{array}{cc} F & B^T \\ & -\tilde S \end{array} \right]^{-1}
   *    = 
   *    \left[ \begin{array}{cc} F^{-1} &  \\  & I \end{array} \right]
   *    \left[ \begin{array}{cc} I & -B^T \\  & I \end{array} \right]
   *    \left[ \begin{array}{cc} I &  \\  & -\tilde S^{-1} \end{array} \right]
   *  \f$
   * 
   *  where for \f$\tilde S\f$ 
   *  is an approximation to the Schur complement S.
   * 
   *  To apply the above
   *  preconditioner, we need a linear solver on the (0,0) block
   *  and an approximation to the inverse of the Schur
   *  complement.
   *
   *  To build a concrete
   *  preconditioner object, we will also need a 2x2 block Thyra
   *  matrix or the 4 separate blocks as either Thyra or
   *  Epetra matrices.  If Thyra, assumes each block is a Thyra
   *  EpetraMatrix.
   */


  class PCDPreconditionerFactory 
    : public PreconditionerFactoryBase<double>
    {
    public:
      /** @name Constructors/initializers/accessors */
      //@{

      /** \brief Default constructor. */
      PCDPreconditionerFactory(); 
      

      /** \brief Constructor for Pressure Convection-Diffusion
       *  preconditioner factory. Takes an AztecOO parameter list for
       *  the F (convection-diffusion) solve and the Ap 
       *  (pressure Poisson) solve.  */
      PCDPreconditionerFactory(
         RCP<const LinearOpWithSolveFactoryBase<double> >   
	 const&  FSolveStrategy,
	 RCP<const LinearOpWithSolveFactoryBase<double> >   
	 const&  ApSolveStrategy
	 );
      
      /** \brief Constructor for Pressure Convection-Diffusion
       *  preconditioner factory. Takes an AztecOO parameter list for
       *  the F (convection-diffusion) solve the Ap (pressure Poisson)
       *  solve, and the Qp (pressure mass matrix) solve.  */
      PCDPreconditionerFactory(
	 RCP<const LinearOpWithSolveFactoryBase<double> > 
	 const& FSolveStrategy,
	 RCP<const LinearOpWithSolveFactoryBase<double> >
	 const& ApSolveStrategy,
	 RCP<const LinearOpWithSolveFactoryBase<double> >
	 const& QpSolveStrategy
	 );
      //@}

      /** @name Overridden from PreconditionerFactoryBase */
      //@{

      /** \brief Check that a <tt>LinearOperator</tt> object is compatible with
       * <tt>*this</tt> factory object.
       */
      bool isCompatible(const LinearOpSourceBase<double> &fwdOpSrc ) const;


      /** \brief Create an (uninitialized) <tt>LinearOperator</tt>
       * object to be initalized as the preconditioner later in
       * <tt>this->initializePrecOp()</tt>.
       *
       * Note that on output <tt>return->domain().get()==NULL</tt> may
       * be true which means that the operator is not fully
       * initialized.  In fact, the output operator object is not
       * guaranteed to be fully initialized until after it is passed
       * through <tt>this->initializePrecOp()</tt>.
       */
      RCP<PreconditionerBase<double> > createPrec() const;


      /** \brief Initialize the PCDPreconditioner object */
      void initializePrec(const RCP<const LinearOpSourceBase<double> >
			  &fwdOpSrc,
			  PreconditionerBase<double> *precOp,
			  const ESupportSolveUse 
			  supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED) const;



      /** \brief Uninitialize the PCDPreconditioner object */
      void uninitializePrec(PreconditionerBase<double> *prec,
			    RCP<const LinearOpSourceBase<double> > *fwdOpSrc,
			    ESupportSolveUse *supportSolveUse = NULL) const;


      //@}

      /** @name Overridden from ParameterListAcceptor */
      //@{
      
      /** \brief . */
      void setParameterList(RCP<ParameterList> const& paramList);
      /** \brief . */
      RCP<ParameterList> getNonconstParameterList();
      /** \brief . */
      RCP<ParameterList> unsetParameterList();
      /** \brief . */
      RCP<const ParameterList> getParameterList() const;
      /** \brief . */
      RCP<const Teuchos::ParameterList> getValidParameters() const;
      //@}

      /* Deal with all of the TSFHandleable features */
      /* GET_RCP(PreconditionerFactoryBase<double>); */

    private:
      mutable RCP<ParameterList>  validPL_;
      RCP<ParameterList>          paramList_;

      RCP<const LinearOpWithSolveFactoryBase<double> >   
	FSolveStrategy_;
      RCP<const LinearOpWithSolveFactoryBase<double> >   
	ApSolveStrategy_;
      RCP<const LinearOpWithSolveFactoryBase<double> >   
	QpSolveStrategy_;
    };

}  // namespace Meros

#endif  // MEROS_PCD_PRECONDITIONERFACTORY_H
