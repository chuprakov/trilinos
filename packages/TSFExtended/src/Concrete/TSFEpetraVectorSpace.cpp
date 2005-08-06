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

#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraVector.hpp"
#include "Teuchos_Utils.hpp"
#include "Epetra_MpiComm.h"

using namespace TSFExtended;
using namespace Teuchos;
using namespace Thyra;

EpetraVectorSpace::EpetraVectorSpace(const RefCountPtr<const Epetra_Map>& m)
  : MPIVectorSpaceBase<double>(),
    Handleable<const VectorSpaceBase<double> >(), 
    epetraMap_(m),
    mpiComm_(MPI_COMM_NULL),
    localSubDim_(epetraMap_->NumMyElements())
{
  const Epetra_Comm& comm = epetraMap_->Comm();
  const Epetra_MpiComm* epMPIComm = dynamic_cast<const Epetra_MpiComm*>(&comm);
  if (epMPIComm != 0) 
    {
#ifdef HAVE_MPI
      mpiComm_ = epMPIComm->GetMpiComm();
#else
      TEST_FOR_EXCEPTION(true, runtime_error,
                         "Epetra communicator is MPI_Comm, but MPI is not "
                         "enabled?!?!?");
#endif
    }
  else
    {
      mpiComm_ = MPI_COMM_NULL;
    }

  updateState(epetraMap_->NumGlobalElements());
}

// Overridden from VectorSpace

Teuchos::RefCountPtr<VectorBase<double> >
EpetraVectorSpace::createMember() const
{
  return rcp(new EpetraVector(rcp(this, false)));
}

Teuchos::RefCountPtr< const VectorSpaceBase<double> >
EpetraVectorSpace::clone() const
{
  return Teuchos::rcp(new EpetraVectorSpace(epetraMap_));
}



string EpetraVectorSpace::description() const
{
  return "EpetraVectorSpace[dim=" + Teuchos::toString(dim()) + ", local="
    + Teuchos::toString(localSubDim()) + "]";
}





