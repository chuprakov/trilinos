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

#ifndef TSFEPETRAVECTORSPACE_HPP
#define TSFEPETRAVECTORSPACE_HPP

#include "TSFConfigDefs.hpp"
#include "Epetra_Map.h"
#include "TSFHandleable.hpp"
#include "Thyra_ScalarProdVectorSpaceBase.hpp"
#include "Thyra_MPIVectorSpaceBase.hpp"


namespace TSFExtended
{
  using namespace Teuchos;
  using namespace Thyra;


  /**
   * Adaptor wrapping Epetra in the Thyra vector space system.
   * We derive from Thyra::ScalarProdVectorSpaceBase in order to
   * inherit defaults for scalar products and view creation.
   */
  class EpetraVectorSpace : virtual public MPIVectorSpaceBase<double>,
                            public Handleable<const VectorSpaceBase<double> >
  {
  public:
    GET_RCP(const Thyra::VectorSpaceBase<double>);

    /** */
    EpetraVectorSpace(const RefCountPtr<const Epetra_Map>& map);
    

    /** @name Overridden form Teuchos::Describable */
    //@{
    /** \brief . */
    std::string description() const;
    //@}

    /** @name Public overridden from VectorSpace */
    //@{
    /** \brief clone the space */
    Teuchos::RefCountPtr< const VectorSpaceBase<double> > clone() const;

    

    /** */
    const RefCountPtr<const Epetra_Map>& epetraMap() const 
    {return epetraMap_;}

    /** */
    MPI_Comm mpiComm() const {return mpiComm_;}

    /** */
    Index localSubDim() const {return localSubDim_;}

  protected:

    /** @name Protected overridden from VectorSpace */
    //@{
    /** \brief create a vector */
    RefCountPtr<VectorBase<double> > createMember() const;

    //@}
  private:
    /** */
    RefCountPtr<const Epetra_Map> epetraMap_;

    MPI_Comm mpiComm_;

    Index localSubDim_;
      
  };
  
}

#endif
