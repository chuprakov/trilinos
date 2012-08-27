/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK2_SPARSECONTAINER_DEF_HPP
#define IFPACK2_SPARSECONTAINER_DEF_HPP

#include "Ifpack2_SparseContainer_decl.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif

namespace Ifpack2 {

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::SparseContainer(const size_t NumRows, const size_t NumVectors) :
  NumRows_(NumRows),
  NumVectors_(NumVectors),
  IsInitialized_(false),
  IsComputed_(false)
{
#ifdef HAVE_MPI
  LocalComm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper((MPI_Comm)MPI_COMM_SELF)));
#else
  LocalComm_ = Teuchos::rcp( new Teuchos::SerialComm<int>() );
#endif
}

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::SparseContainer(const SparseContainer<MatrixType,InverseType>& rhs)
{
  // nobody should ever call this  
}

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::~SparseContainer()
{
  destroy();
}

//==============================================================================
template<class MatrixType, class InverseType>
size_t SparseContainer<MatrixType,InverseType>::getNumRows() const
{
  if (isInitialized()) return NumRows_;
  else return 0;
}

//==============================================================================
// Returns the number of vectors in X/Y.
template<class MatrixType, class InverseType>
size_t SparseContainer<MatrixType,InverseType>::getNumVectors() const
{
  return NumVectors_;
}
//==============================================================================
  // Sets the number of vectors for X/Y
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::setNumVectors(const size_t NumVectors_in)
{
  // NTS: We didn't implement this in Ifpack either.
  throw std::runtime_error("Ifpack2::SparseContainer: does not support setNumVectors.");
}

//==============================================================================
// Get the Y vector
template<class MatrixType, class InverseType>
const Teuchos::RCP<Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > & SparseContainer<MatrixType,InverseType>::getY()
{
  return Y_;
}

//==============================================================================
// Get the X vector
template<class MatrixType, class InverseType>
const Teuchos::RCP<Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > & SparseContainer<MatrixType,InverseType>::getX()
{
  return X_;
}

//==============================================================================
// Returns the ID associated to local row i. 
template<class MatrixType, class InverseType>
typename MatrixType::local_ordinal_type & SparseContainer<MatrixType,InverseType>::ID(const LocalOrdinal i)
{
  return GID_[i];
}

//==============================================================================
// Returns  true is the container has been successfully initialized.
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isInitialized() const
{
  return IsInitialized_;
}

//==============================================================================
// Returns  true is the container has been successfully computed.
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isComputed() const
{
  return IsComputed_;
}

//==============================================================================
// Sets all necessary parameters.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::setParameters(const Teuchos::ParameterList& List)
{
  List_ = List;
}
 
//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::initialize()
{
  if(IsInitialized_) destroy();
  IsInitialized_=false;

  // NTS: The use of LocalOrdinal/LocalOrdinal for the local matrices is somewhat dubious.
  // We should find a more general approach by letting the InverseType pick the ordinals for the
  // matrices I build.  Unfortunately Ifpack2::Preconditioner doesn't support that yet.

  Map_ = Teuchos::rcp( new Tpetra::Map<LocalOrdinal,LocalOrdinal,Node>(NumRows_,0,LocalComm_) );
  X_ = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,LocalOrdinal,Node>(Map_,NumVectors_) );
  Y_ = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,LocalOrdinal,Node>(Map_,NumVectors_) );
  GID_.resize(NumRows_);

  Matrix_ = Teuchos::rcp( new Tpetra::CrsMatrix<Scalar,LocalOrdinal,LocalOrdinal,Node>(Map_,0) );

  // create the inverse
  Inverse_ = Teuchos::rcp( new InverseType(Matrix_) );
  TEUCHOS_TEST_FOR_EXCEPTION( Inverse_ == Teuchos::null, std::runtime_error, "Ifpack2::SparseContainer::initialize error in inverse constructor.");  
  Inverse_->setParameters(List_);

  // Note from IFPACK: Call Inverse_->Initialize() in Compute(). This saves
  // some time, because the diagonal blocks can be extracted faster,
  // and only once.
  IsInitialized_ = true;
}

//==============================================================================
// Finalizes the linear system matrix and prepares for the application of the inverse.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::compute(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix)
{
  IsComputed_=false;
  TEUCHOS_TEST_FOR_EXCEPTION( !IsInitialized_, std::runtime_error, "Ifpack2::SparseContainer::compute please call initialize first.");  

  // extract the submatrices
  extract(Matrix);

  // initialize & compute the inverse operator
  Inverse_->initialize();
  Inverse_->compute();
  IsComputed_=true;
}

//==============================================================================
// Apply the inverse of the matrix to RHS, result is stored in LHS.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::apply()
{
 TEUCHOS_TEST_FOR_EXCEPTION( !IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::apply compute has not been called.");  
 Inverse_->apply(*X_, *Y_);
}

//==============================================================================
// Destroys all data.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//============================================================================== 
// Prints basic information on iostream. This function is used by operator<<
template<class MatrixType, class InverseType>
std::ostream& SparseContainer<MatrixType,InverseType>::print(std::ostream& os) const
{
  using std::endl;
  os << "================================================================================" << endl;
  os << "Ifpack2_SparseContainer" << endl;
  os << "Number of rows          = " << NumRows_ << endl;
  os << "Number of vectors       = " << NumVectors_ << endl;
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
  return os;
}


//==============================================================================
template<class MatrixType, class InverseType>
std::string SparseContainer<MatrixType,InverseType>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }

  oss << "}";
  return oss.str();
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl;
  os << "Ifpack2_SparseContainer" << endl;
  os << "Number of rows          = " << NumRows_ << endl;
  os << "Number of vectors       = " << NumVectors_ << endl;
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

//============================================================================== 
// Extract the submatrices identified by the ID set int ID().
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::extract(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix_in) 
{
  size_t MatrixInNumRows= Matrix_in->getNodeNumRows();

  // Sanity checking
  for (size_t j = 0 ; j < NumRows_ ; ++j) {
    TEUCHOS_TEST_FOR_EXCEPTION( GID_[j] < 0 || (size_t) GID_[j] >= MatrixInNumRows, std::runtime_error, "Ifpack2::SparseContainer::applyInverse compute has not been called.");  
  }  

  int Length = Matrix_in->getNodeMaxNumRowEntries();
  Teuchos::Array<Scalar>       Values;
  Teuchos::Array<LocalOrdinal> Indices;
  Teuchos::Array<Scalar>       Values_insert;
  Teuchos::Array<LocalOrdinal> Indices_insert;

  Values.resize(Length);
  Indices.resize(Length);
  Values_insert.resize(Length);
  Indices_insert.resize(Length);

  for (size_t j = 0 ; j < NumRows_ ; ++j) {
    LocalOrdinal LRID = ID(j);
    size_t NumEntries;

    Matrix_in->getLocalRowCopy(LRID,Indices(),Values(),NumEntries);

    size_t num_entries_found=0;
    for (size_t k = 0 ; k < NumEntries ; ++k) {
      LocalOrdinal LCID = Indices[k];

      // skip off-processor elements
      if ((size_t)LCID >= MatrixInNumRows) 
	continue;

      // for local column IDs, look for each ID in the list
      // of columns hosted by this object
      // FIXME: use STL
      int jj = -1;
      for (size_t kk = 0 ; kk < NumRows_ ; ++kk)
	if (ID(kk) == LCID)
	  jj = kk;

      if (jj != -1) {
	Indices_insert[num_entries_found] = jj;
	Values_insert[num_entries_found]  = Values[k];
	num_entries_found++;
      }

    }
    Matrix_->insertGlobalValues(j,Indices_insert(0,num_entries_found),Values_insert(0,num_entries_found));
  }

  Matrix_->fillComplete();
}


} // namespace Ifpack2
#endif // IFPACK2_SPARSECONTAINER_HPP
