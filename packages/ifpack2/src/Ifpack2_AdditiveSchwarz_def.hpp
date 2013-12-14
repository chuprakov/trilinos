/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_ADDITIVESCHWARZ_DEF_HPP
#define IFPACK2_ADDITIVESCHWARZ_DEF_HPP

#include "Ifpack2_AdditiveSchwarz_decl.hpp"

// AdditiveSchwarz uses OneLevelFactory to create a default inner
// preconditioner.
//
// FIXME (mfh 13 Dec 2013) For some inexplicable reason, I have to
// include the _decl and _def headers separately here; including just
// Ifpack2_Details_OneLevelFactory.hpp doesn't work.  It probably has
// something to do with ETI, but I don't fully understand what.
#include "Ifpack2_Details_OneLevelFactory_decl.hpp"
#include "Ifpack2_Details_OneLevelFactory_def.hpp"

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
#include "Xpetra_RowMatrix.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"
#include "Zoltan2_XpetraRowMatrixInput.hpp"
#include "Zoltan2_OrderingProblem.hpp"
#include "Zoltan2_OrderingSolution.hpp"
#endif

#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Ifpack2_LocalFilter_def.hpp"
#include "Ifpack2_OverlappingRowMatrix_def.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_ReorderFilter_def.hpp"
#include "Ifpack2_SingletonFilter_def.hpp"

#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

#include <locale> // std::toupper

namespace Ifpack2 {

namespace Details {
//! Map from an Ifpack2::Preconditioner subclass to its string name.
template<class PrecType>
class OneLevelPreconditionerNamer {
public:
  //! Name corresponding to Preconditioner subclass PrecType.
  static std::string name ();
};

//
// Partial specializations for each single-level preconditioner.
//

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Chebyshev<MatrixType> > {
public:
  static std::string name () {
    return "CHEBYSHEV";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Details::DenseSolver<MatrixType> > {
public:
  static std::string name () {
    return "DENSE";
  }
};

#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_AMESOS2)
template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Details::Amesos2Wrapper<MatrixType> > {
public:
  static std::string name () {
    return "AMESOS2";
  }
};
#endif

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Diagonal<MatrixType> > {
public:
  static std::string name () {
    return "DIAGONAL";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::ILUT<MatrixType> > {
public:
  static std::string name () {
    return "ILUT";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::Relaxation<MatrixType> > {
public:
  static std::string name () {
    return "RELAXATION";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::RILUK<MatrixType> > {
public:
  static std::string name () {
    return "RILUK";
  }
};

template<class MatrixType>
class OneLevelPreconditionerNamer< ::Ifpack2::IdentitySolver<MatrixType> > {
public:
  static std::string name () {
    return "IDENTITY";
  }
};

} // namespace Details


template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A) :
  Matrix_ (A),
  IsInitialized_(false),
  IsComputed_(false),
  IsOverlapping_(false),
  OverlapLevel_ (0),
  CombineMode_(Tpetra::ADD),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  ComputeCondest_(true),
  UseReordering_(false),
  ReorderingAlgorithm_("none"),
  UseSubdomain_(false),
  FilterSingletons_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0)
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::SerialComm;

  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::invalid_argument, "Ifpack2::AdditiveSchwarz "
    "constructor: The input matrix A must be nonnull.");

  RCP<const Teuchos::Comm<int> > comm = Matrix_->getComm ();
  RCP<const map_type> rowMap = Matrix_->getRowMap ();
  RCP<node_type> node = Matrix_->getNode ();
  const global_size_t INVALID =
    Teuchos::OrdinalTraits<global_size_t>::invalid ();

  // If there's only one process in the matrix's communicator,
  // then there's no need to compute overlap.
  if (comm->getSize () == 1) {
    OverlapLevel_ = 0;
    IsOverlapping_ = false;
  } else if (OverlapLevel_ != 0) {
    IsOverlapping_ = true;
  }

  if (OverlapLevel_ == 0) {
    const global_ordinal_type indexBase = rowMap->getIndexBase ();

    // FIXME (mfh 28 Sep 2013) I don't understand why this is called a
    // "serial Map."  It's the same Map as the input matrix's row Map!
    // It's also the same Map as "distributed Map"!  I would change it
    // myself, but I don't want to break anything, so I just
    // reformatted the code to comply better with Ifpack2 standards
    // and left the names alone.
    SerialMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));
    DistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));

    RCP<const SerialComm<int> > localComm (new SerialComm<int> ());

    LocalDistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeNumElements (),
                         indexBase, localComm, node));
  }


  // Set parameters to default values
  Teuchos::ParameterList plist;
  setParameters (plist);
}

template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A,
                 const int overlapLevel) :
  Matrix_ (A),
  IsInitialized_(false),
  IsComputed_(false),
  IsOverlapping_(false),
  OverlapLevel_ (overlapLevel),
  CombineMode_(Tpetra::ADD),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  ComputeCondest_(true),
  UseReordering_(false),
  ReorderingAlgorithm_("none"),
  UseSubdomain_(false),
  FilterSingletons_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0)
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::SerialComm;

  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::invalid_argument, "Ifpack2::AdditiveSchwarz "
    "constructor: The input matrix A must be nonnull.");

  RCP<const Teuchos::Comm<int> > comm = Matrix_->getComm ();
  RCP<const map_type> rowMap = Matrix_->getRowMap ();
  RCP<node_type> node = Matrix_->getNode ();
  const global_size_t INVALID =
    Teuchos::OrdinalTraits<global_size_t>::invalid ();

  // If there's only one process in the matrix's communicator,
  // then there's no need to compute overlap.
  if (comm->getSize () == 1) {
    OverlapLevel_ = 0;
    IsOverlapping_ = false;
  } else if (OverlapLevel_ != 0) {
    IsOverlapping_ = true;
  }

  if (OverlapLevel_ == 0) {
    const global_ordinal_type indexBase = rowMap->getIndexBase ();

    // FIXME (mfh 28 Sep 2013) I don't understand why this is called a
    // "serial Map."  It's the same Map as the input matrix's row Map!
    // It's also the same Map as "distributed Map"!  I would change it
    // myself, but I don't want to break anything, so I just
    // reformatted the code to comply better with Ifpack2 standards
    // and left the names alone.
    SerialMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));
    DistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));

    RCP<const SerialComm<int> > localComm (new SerialComm<int> ());

    LocalDistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeNumElements (),
                         indexBase, localComm, node));
  }

  // Set parameters to default values
  Teuchos::ParameterList plist;
  setParameters (plist);
}


template<class MatrixType,class LocalInverseType>
AdditiveSchwarz<MatrixType,LocalInverseType>::~AdditiveSchwarz () {}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type > >
AdditiveSchwarz<MatrixType,LocalInverseType>::getDomainMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::getDomainMap: "
    "The matrix A to precondition is null.");
  return Matrix_->getDomainMap ();
}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
AdditiveSchwarz<MatrixType,LocalInverseType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::getRangeMap: "
    "The matrix A to precondition is null.");
  return Matrix_->getRangeMap ();
}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > AdditiveSchwarz<MatrixType,LocalInverseType>::getMatrix() const
{
  return Matrix_;
}


template<class MatrixType,class LocalInverseType>
void
AdditiveSchwarz<MatrixType,LocalInverseType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const std::string timerName ("Ifpack2::AdditiveSchwarz::apply");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    const scalar_type ZERO = Teuchos::ScalarTraits<scalar_type>::zero ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! IsComputed_, std::runtime_error,
      "Ifpack2::AdditiveSchwarz::apply: "
      "isComputed() must be true before you may call apply().");

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
      "Ifpack2::AdditiveSchwarz::apply: "
      "X and Y must have the same number of columns.  X has "
      << X.getNumVectors() << " columns, but Y has " << Y.getNumVectors() << ".");

    TEUCHOS_TEST_FOR_EXCEPTION(
      Matrix_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::apply: "
      "The input matrix A is null, but the preconditioner says that it has "
      "been computed (isComputed() is true).  This should never happen.  "
      "Please report this bug to the Ifpack2 developers.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      Inverse_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::apply: "
      "Inverse_ is null, but the preconditioner says that it has been computed "
      "(isComputed() is true).  This should never happen.  "
      "Please report this bug to the Ifpack2 developers.");

    const size_t numVectors = X.getNumVectors ();

    RCP<MV> OverlappingX,OverlappingY,Xtmp;

    if (IsOverlapping_) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        OverlappingMatrix_.is_null (), std::logic_error,
        "Ifpack2::AdditiveSchwarz::apply: The overlapping matrix is null.  "
        "This should never happen if IsOverlapping_ is true.  "
        "Please report this bug to the Ifpack2 developers.");

      // Setup if we're overlapping
      OverlappingX = rcp (new MV (OverlappingMatrix_->getRowMap (), numVectors));
      OverlappingY = rcp (new MV (OverlappingMatrix_->getRowMap (), numVectors));
      // FIXME (mfh 28 Sep 2013) MV's constructor fills with zeros by default,
      // so there is no need to call putScalar().
      OverlappingY->putScalar (ZERO);
      OverlappingX->putScalar (ZERO);
      OverlappingMatrix_->importMultiVector (X, *OverlappingX, Tpetra::INSERT);
      // FIXME from Ifpack1: Will not work with non-zero starting solutions.
    }
    else {
      Xtmp = rcp (new MV (X));

      TEUCHOS_TEST_FOR_EXCEPTION(
        LocalDistributedMap_.is_null (), std::logic_error,
        "Ifpack2::AdditiveSchwarz::apply: "
        "LocalDistributedMap_ is null.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        DistributedMap_.is_null (), std::logic_error,
        "Ifpack2::AdditiveSchwarz::apply: "
        "DistributedMap_ is null.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        SerialMap_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::apply: "
        "SerialMap_ is null.");

      MV Serial (SerialMap_, numVectors);
      // Create Import object on demand, if necessary.
      if (SerialImporter_.is_null ()) {
        SerialImporter_ =
          rcp (new import_type (SerialMap_, Matrix_->getDomainMap ()));
      }
      Serial.doImport (*Xtmp, *SerialImporter_, Tpetra::INSERT);

      OverlappingX = rcp (new MV (LocalDistributedMap_, numVectors));
      OverlappingY = rcp (new MV (LocalDistributedMap_, numVectors));

      //OverlappingX->putScalar(0.0);
      //OverlappingY->putScalar(0.0);

      MV Distributed (DistributedMap_, numVectors);
      // Create Import object on demand, if necessary.
      if (DistributedImporter_.is_null ()) {
        DistributedImporter_ =
          rcp (new import_type (DistributedMap_, Matrix_->getDomainMap ()));
      }
      Distributed.doImport (*Xtmp, *DistributedImporter_, Tpetra::INSERT);

      // FIXME (mfh 28 Sep 2013) Please don't call replaceLocalValue()
      // for every entry.  It's important to understand how MultiVector
      // views work.
      Teuchos::ArrayRCP<const scalar_type> values = Distributed.get1dView();
      size_t index = 0;

      for (size_t v = 0; v < numVectors; v++) {
        for (size_t i = 0; i < Matrix_->getRowMap()->getNodeNumElements(); i++) {
          OverlappingX->replaceLocalValue(i, v, values[index]);
          index++;
        }
      }
    }

    if (FilterSingletons_) {
      // process singleton filter
      MV ReducedX (SingletonMatrix_->getRowMap (), numVectors);
      MV ReducedY (SingletonMatrix_->getRowMap (), numVectors);
      SingletonMatrix_->SolveSingletons (*OverlappingX, *OverlappingY);
      SingletonMatrix_->CreateReducedRHS (*OverlappingY, *OverlappingX, ReducedX);

      // process reordering
      if (! UseReordering_) {
        Inverse_->apply (ReducedX, ReducedY);
      }
      else {
        MV ReorderedX (ReducedX);
        MV ReorderedY (ReducedY);
        ReorderedLocalizedMatrix_->permuteOriginalToReordered (ReducedX, ReorderedX);
        Inverse_->apply (ReorderedX, ReorderedY);
        ReorderedLocalizedMatrix_->permuteReorderedToOriginal (ReorderedY, ReducedY);
      }

      // finish up with singletons
      SingletonMatrix_->UpdateLHS (ReducedY, *OverlappingY);
    }
    else {

      // process reordering
      if (! UseReordering_) {
        Inverse_->apply (*OverlappingX, *OverlappingY);
      }
      else {
        MV ReorderedX (*OverlappingX);
        MV ReorderedY (*OverlappingY);
        ReorderedLocalizedMatrix_->permuteOriginalToReordered (*OverlappingX, ReorderedX);
        Inverse_->apply (ReorderedX, ReorderedY);
        ReorderedLocalizedMatrix_->permuteReorderedToOriginal (ReorderedY, *OverlappingY);
      }
    }

    if (IsOverlapping_) {
      OverlappingMatrix_->exportMultiVector (*OverlappingY, Y, CombineMode_);
    }
    else {
      Teuchos::ArrayRCP<const scalar_type> values = OverlappingY->get1dView();
      size_t index = 0;

      // FIXME (mfh 28 Sep 2013) Please don't call replaceLocalValue()
      // for every entry.  It's important to understand how MultiVector
      // views work.
      for (size_t v = 0; v < numVectors; v++) {
        for (size_t i = 0; i < Matrix_->getRowMap()->getNodeNumElements(); i++) {
          Y.replaceLocalValue(i, v, values[index]);
          index++;
        }
      }
    }
  } // Stop timing here.

  ++NumApply_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ApplyTime_ = timer->totalElapsedTime ();
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::
setParameters (const Teuchos::ParameterList& plist)
{
  // mfh 18 Nov 2013: Ifpack2's setParameters() method passes in the
  // input list as const.  This means that we have to copy it before
  // validation or passing into setParameterList().
  List_ = plist;
  this->setParameterList (Teuchos::rcpFromRef (List_));
}



template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::
setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  using Tpetra::CombineMode;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntry;
  using Teuchos::ParameterEntryValidator;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::StringToIntegralParameterEntryValidator;

  if (plist.is_null ()) {
    // Assume that the user meant to set default parameters by passing
    // in an empty list.
    this->setParameterList (Teuchos::parameterList ());
  }
  // At this point, plist should be nonnull.
  TEUCHOS_TEST_FOR_EXCEPTION(
    plist.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setParameterList: plist is null.  This should never happen, since the "
    "method should have replaced a null input list with a nonnull empty list "
    "by this point.  Please report this bug to the Ifpack2 developers.");

  // try {
  //   List_.validateParameters (* getValidParameters ());
  // }
  // catch (std::exception& e) {
  //   std::cerr << "Ifpack2::AdditiveSchwarz::setParameterList: Validation failed with the following error message: " << e.what () << std::endl;
  //   throw e;
  // }

  // mfh 18 Nov 2013: Supplying the current value as the default value
  // when calling ParameterList::get() ensures "delta" behavior when
  // users pass in new parameters: any unspecified parameters in the
  // new list retain their values in the old list.  This preserves
  // backwards compatiblity with this class' previous behavior.  Note
  // that validateParametersAndSetDefaults() would have different
  // behavior: any parameters not in the new list would get default
  // values, which could be different than their values in the
  // original list.

  ComputeCondest_ = plist->get ("schwarz: compute condest", ComputeCondest_);

  bool gotCombineMode = false;
  try {
    CombineMode_ = getIntegralValue<Tpetra::CombineMode> (List_, "schwarz: combine mode");
    gotCombineMode = true;
  }
  catch (Teuchos::Exceptions::InvalidParameterName&) {
    // The caller didn't provide that parameter.  Just keep the
    // existing value of CombineMode_.
    gotCombineMode = true;
  }
  catch (Teuchos::Exceptions::InvalidParameterType&) {
    // The user perhaps supplied it as an Tpetra::CombineMode enum
    // value.  Let's try again (below).  If it doesn't succeed, we
    // know that the type is wrong, so we can let it throw whatever
    // exception it would throw.
  }
  // Try to get the combine mode as an integer.
  if (! gotCombineMode) {
    try {
      CombineMode_ = plist->get ("schwarz: combine mode", CombineMode_);
      gotCombineMode = true;
    }
    catch (Teuchos::Exceptions::InvalidParameterType&) {}
  }
  // Try to get the combine mode as a string.  If this works, use the
  // validator to convert to int.  This is painful, but necessary in
  // order to do validation, since the input list doesn't come with a
  // validator.
  if (! gotCombineMode) {
    const ParameterEntry& validEntry =
      getValidParameters ()->getEntry ("schwarz: combine mode");
    RCP<const ParameterEntryValidator> v = validEntry.validator ();
    typedef StringToIntegralParameterEntryValidator<CombineMode> vs2e_type;
    RCP<const vs2e_type> vs2e = rcp_dynamic_cast<const vs2e_type> (v, true);

    const ParameterEntry& inputEntry = plist->getEntry ("schwarz: combine mode");
    CombineMode_ = vs2e->getIntegralValue (inputEntry, "schwarz: combine mode");
    gotCombineMode = true;
  }
  (void) gotCombineMode; // forestall "set but not used" compiler warning

  OverlapLevel_ = plist->get ("schwarz: overlap level", OverlapLevel_);
  if (OverlapLevel_ != 0 && Matrix_->getComm ()->getSize () > 1) {
    IsOverlapping_ = true;
  }

  // Will we be doing reordering?  Unlike Ifpack, we'll use a
  // "schwarz: reordering list" to give to Zoltan2.
  UseReordering_ = plist->get ("schwarz: use reordering", UseReordering_);

#if !defined(HAVE_IFPACK2_XPETRA) || !defined(HAVE_IFPACK2_ZOLTAN2)
  TEUCHOS_TEST_FOR_EXCEPTION(
    UseReordering_, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
    "setParameters: You specified \"schwarz: use reordering\" = true.  "
    "This is only valid when Trilinos was built with Ifpack2, Xpetra, and "
    "Zoltan2 enabled.  Either Xpetra or Zoltan2 was not enabled in your build "
    "of Trilinos.");
#endif

  // FIXME (mfh 18 Nov 2013) Now would be a good time to validate the
  // "schwarz: reordering list" parameter list.  Currently, that list
  // gets extracted in setup().

  // Subdomain check
  if (plist->isParameter ("schwarz: subdomain id") && plist->get ("schwarz: subdomain id", -1) > 0) {
    UseSubdomain_ = true;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    UseSubdomain_, std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setParameters: You specified the \"schwarz: subdomain id\" parameter, "
    "with a value other than -1.  This parameter is not yet supported.");

  // if true, filter singletons. NOTE: the filtered matrix can still have
  // singletons! A simple example: upper triangular matrix, if I remove
  // the lower node, I still get a matrix with a singleton! However, filter
  // singletons should help for PDE problems with Dirichlet BCs.
  FilterSingletons_ = plist->get ("schwarz: filter singletons", FilterSingletons_);

  // FIXME (mfh 18 Nov 2013) If the inner solver exists, now might be
  // a good time to validate its parameters.
}



template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Teuchos::ParameterList>
AdditiveSchwarz<MatrixType,LocalInverseType>::
getValidParameters () const
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;

  if (validParams_.is_null ()) {
    const int overlapLevel = 0;
    const bool useReordering = false;
    const bool computeCondest = false;
    const bool filterSingletons = false;
    ParameterList reorderingSublist;
    reorderingSublist.set ("order_method", std::string ("rcm"));

    RCP<ParameterList> plist = parameterList ("Ifpack2::AdditiveSchwarz");

    Tpetra::setCombineModeParameter (*plist, "schwarz: combine mode");
    plist->set ("schwarz: overlap level", overlapLevel);
    plist->set ("schwarz: use reordering", useReordering);
    plist->set ("schwarz: reordering list", reorderingSublist);
    plist->set ("schwarz: compute condest", computeCondest);
    plist->set ("schwarz: filter singletons", filterSingletons);

    // FIXME (mfh 18 Nov 2013) Get valid parameters from inner solver.
    //
    // FIXME (mfh 18 Nov 2013) Get valid parameters from Zoltan2, if
    // Zoltan2 was enabled in the build.

    validParams_ = rcp_const_cast<const ParameterList> (plist);
  }
  return validParams_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  const std::string timerName ("Ifpack2::AdditiveSchwarz::initialize");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
    TEUCHOS_TEST_FOR_EXCEPTION(timer.is_null (), std::logic_error,
      "Timer is null; this shouldn't happen.");
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    IsInitialized_ = false;
    IsComputed_ = false;
    Condest_ = -Teuchos::ScalarTraits<magnitude_type>::one ();

    // compute the overlapping matrix if necessary
    if (IsOverlapping_) {
      if (UseSubdomain_) {
        const int sid = List_.get ("subdomain id", -1);
        OverlappingMatrix_ = rcp (new OverlappingRowMatrix<row_matrix_type> (Matrix_, OverlapLevel_, sid));
      } else {
        OverlappingMatrix_ = rcp (new OverlappingRowMatrix<row_matrix_type> (Matrix_, OverlapLevel_));
      }
    }

    // Setup
    setup ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      Inverse_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::initialize: "
      "Inverse_ is null after calling setup().  "
      "This should never happen.  "
      "Please report this bug to the Ifpack2 developers.");

    // Initialize subdomain solver.
    //
    // FIXME (mfh 28 Sep 2013) The "inverse" should have its own sublist
    // in the input ParameterList.  We shouldn't pass AdditiveSchwarz's
    // parameters directly to the "inverse."
    //
    // FIXME (mfh 28 Sep 2013) Why don't we call these methods in setup()?
    Inverse_->setParameters (List_);
    Inverse_->initialize ();
  } // Stop timing here.

  IsInitialized_ = true;
  ++NumInitialize_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  InitializeTime_ = timer->totalElapsedTime ();
}


template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isInitialized() const
{
  return IsInitialized_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::compute ()
{
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  if (! IsInitialized_) {
    initialize ();
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isInitialized (), std::logic_error, "Ifpack2::AdditiveSchwarz::compute: "
    "The preconditioner is not yet initialized, "
    "even though initialize() supposedly has been called.  "
    "This should never happen.  "
    "Please report this bug to the Ifpack2 developers.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    Inverse_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::compute: "
    "Inverse_ is null, but the preconditioner says that it has been initialized "
    "(isInitialized() is true).  This should never happen.  "
    "Please report this bug to the Ifpack2 developers.");

  const std::string timerName ("Ifpack2::AdditiveSchwarz::compute");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    IsComputed_ = false;
    Condest_ = -Teuchos::ScalarTraits<magnitude_type>::one ();
    Inverse_->compute ();
  } // Stop timing here.

  IsComputed_ = true;
  ++NumCompute_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ComputeTime_ = timer->totalElapsedTime ();
}

//==============================================================================
// Returns true if the  preconditioner has been successfully computed, false otherwise.
template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isComputed() const
{
  return IsComputed_;
}


template<class MatrixType,class LocalInverseType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
AdditiveSchwarz<MatrixType,LocalInverseType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &Matrix_in)
{
  // The preconditioner must have been computed in order to estimate
  // its condition number.
  if (! isComputed ()) {
    return -Teuchos::ScalarTraits<magnitude_type>::one ();
  }

  Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, Matrix_in);
  return Condest_;
}


template<class MatrixType,class LocalInverseType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
AdditiveSchwarz<MatrixType,LocalInverseType>::getCondEst() const
{
  return Condest_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumInitialize() const
{
  return NumInitialize_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumCompute() const
{
  return NumCompute_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumApply() const
{
  return NumApply_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getInitializeTime() const
{
  return InitializeTime_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getComputeTime() const
{
  return ComputeTime_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getApplyTime() const
{
  return ApplyTime_;
}


template<class MatrixType,class LocalInverseType>
std::string AdditiveSchwarz<MatrixType,LocalInverseType>::description() const
{
  using Teuchos::TypeNameTraits;

  std::ostringstream out;
  out << "Ifpack2::AdditiveSchwarz: {";
  out << "MatrixType: " << TypeNameTraits<MatrixType>::name ()
      << ", LocalInverseType: " << TypeNameTraits<LocalInverseType>::name ();
  if (this->getObjectLabel () != "") {
    out << ", Label: \"" << this->getObjectLabel () << "\"";
  }
  out << ", Initialized: " << (isInitialized () ? "true" : "false")
      << ", Computed: " << (isComputed () ? "true" : "false")
      << ", Overlap level: " << OverlapLevel_
      << ", Subdomain reordering: \"" << ReorderingAlgorithm_ << "\""
      << "}";
  return out.str ();
}


template<class MatrixType,class LocalInverseType>
void
AdditiveSchwarz<MatrixType,LocalInverseType>::
describe (Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::OSTab;
  using Teuchos::TypeNameTraits;
  using std::endl;

  const int myRank = Matrix_->getComm ()->getRank ();
  const int numProcs = Matrix_->getComm ()->getSize ();
  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl > Teuchos::VERB_NONE) {
    // describe() starts with a tab, by convention.
    OSTab tab0 (out);
    if (myRank == 0) {
      out << "Ifpack2::AdditiveSchwarz:";
    }
    OSTab tab1 (out);
    if (myRank == 0) {
      out << "MatrixType: " << TypeNameTraits<MatrixType>::name () << endl;
      out << "LocalInverseType: " << TypeNameTraits<LocalInverseType>::name () << endl;
      if (this->getObjectLabel () != "") {
        out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
      }

      out << "Overlap level: " << OverlapLevel_ << endl
          << "Combine mode: \"";
      if (CombineMode_ == Tpetra::INSERT) {
        out << "INSERT";
      } else if (CombineMode_ == Tpetra::ADD) {
        out << "ADD";
      } else if (CombineMode_ == Tpetra::REPLACE) {
        out << "REPLACE";
      } else if (CombineMode_ == Tpetra::ABSMAX) {
        out << "ABSMAX";
      } else if (CombineMode_ == Tpetra::ZERO) {
        out << "ZERO";
      }
      out << "\"" << endl;
    }

    if (Matrix_.is_null ()) {
      if (myRank == 0) {
        out << "Input matrix A: null" << endl;
      }
    }
    else {
      if (myRank == 0) {
        out << "Input matrix A:" << endl;
        std::flush (out);
      }
      Matrix_->getComm ()->barrier (); // wait for output to finish
      Matrix_->describe (out, Teuchos::VERB_LOW);
    }

    if (myRank == 0) {
      out << "Number of initialize calls: " << getNumInitialize () << endl
          << "Number of compute calls: " << getNumCompute () << endl
          << "Number of apply calls: " << getNumApply () << endl
          << "Total time in seconds for initialize: " << getInitializeTime () << endl
          << "Total time in seconds for compute: " << getComputeTime () << endl
          << "Total time in seconds for apply: " << getApplyTime () << endl;
    }

    if (Inverse_.is_null ()) {
      if (myRank == 0) {
        out << "Subdomain solver: null" << endl;
      }
    }
    else {
      if (vl < Teuchos::VERB_EXTREME) {
        if (myRank == 0) {
          out << "Subdomain solver: not null" << endl;
        }
      }
      else { // vl >= Teuchos::VERB_EXTREME
        for (int p = 0; p < numProcs; ++p) {
          if (p == myRank) {
            out << "Subdomain solver on Process " << myRank << ":";
            if (Inverse_.is_null ()) {
              out << "null" << endl;
            } else {
              out << endl;
              Inverse_->describe (out, vl);
            }
          }
          Matrix_->getComm ()->barrier ();
          Matrix_->getComm ()->barrier ();
          Matrix_->getComm ()->barrier (); // wait for output to finish
        }
      }
    }

    Matrix_->getComm ()->barrier (); // wait for output to finish
  }
}


template<class MatrixType,class LocalInverseType>
std::ostream& AdditiveSchwarz<MatrixType,LocalInverseType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getOverlapLevel() const
{
  return OverlapLevel_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::setup ()
{
#ifdef HAVE_MPI
  using Teuchos::MpiComm;
#endif // HAVE_MPI
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  typedef Xpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> XpetraMatrixType;
  typedef Xpetra::TpetraRowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> XpetraTpetraMatrixType;
#endif

  TEUCHOS_TEST_FOR_EXCEPTION(
    Matrix_.is_null (), std::runtime_error, "Ifpack2::AdditiveSchwarz::setup: "
    "The matrix A to precondition is null.");

  // Localized version of Matrix_ or OverlappingMatrix_.
  RCP<row_matrix_type> LocalizedMatrix;

  // The "most current local matrix."  At the end of this method, this
  // will be handed off to the inner solver.
  RCP<row_matrix_type> ActiveMatrix;

  // Create localized matrix.
  if (! OverlappingMatrix_.is_null ()) {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Ifpack2::AdditiveSchwarz::setup: subdomain code not yet supported.");
      //
      // FIXME (mfh 18 Nov 2013) btw what's the difference between
      // Ifpack_NodeFilter and Ifpack_LocalFilter?  The former's
      // documentation sounds a lot like what Ifpack2::LocalFilter
      // does.
      //
      //Ifpack2_NodeFilter *tt = new Ifpack2_NodeFilter(OverlappingMatrix_,nodeID); //FIXME
      //LocalizedMatrix = Teuchos::rcp(tt);
    }
    else
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (OverlappingMatrix_));
  }
  else {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Ifpack2::AdditiveSchwarz::setup: subdomain code not yet supported.");
    }
    else {
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (Matrix_));
    }
  }

  // Sanity check; I don't trust the logic above to have created LocalizedMatrix.
  TEUCHOS_TEST_FOR_EXCEPTION(
    LocalizedMatrix.is_null (), std::logic_error,
    "Ifpack2::AdditiveSchwarz::setup: LocalizedMatrix is null, after the code "
    "that claimed to have created it.  This should never be the case.  Please "
    "report this bug to the Ifpack2 developers.");

  // Mark localized matrix as active
  ActiveMatrix = LocalizedMatrix;

  // Singleton Filtering
  if (FilterSingletons_) {
    SingletonMatrix_ = rcp (new SingletonFilter<row_matrix_type> (LocalizedMatrix));
    ActiveMatrix = SingletonMatrix_;
  }

  // Do reordering
  if (UseReordering_) {
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
    // Unlike Ifpack, Zoltan2 does all the dirty work here.
    Teuchos::ParameterList zlist = List_.sublist ("schwarz: reordering list");

    // FIXME (mfh 18 Nov 2013) Shouldn't this come from the Zoltan2 sublist?
    ReorderingAlgorithm_ = List_.get<std::string> ("order_method", "rcm");
    XpetraTpetraMatrixType XpetraMatrix (ActiveMatrix);
    Zoltan2::XpetraRowMatrixAdapter<XpetraMatrixType> Zoltan2Matrix (rcpFromRef (XpetraMatrix));

    typedef Zoltan2::OrderingProblem<Zoltan2::XpetraRowMatrixAdapter<XpetraMatrixType> > ordering_problem_type;
#ifdef HAVE_MPI
    // Grab the MPI Communicator and build the ordering problem with that
    MPI_Comm myRawComm;

    RCP<const MpiComm<int> > mpicomm =
      rcp_dynamic_cast<const MpiComm<int> > (ActiveMatrix->getComm ());
    if (mpicomm == Teuchos::null) {
      myRawComm = MPI_COMM_SELF;
    } else {
      myRawComm = * (mpicomm->getRawMpiComm ());
    }
    ordering_problem_type MyOrderingProblem (&Zoltan2Matrix, &zlist, myRawComm);
#else
    ordering_problem_type MyOrderingProblem (&Zoltan2Matrix, &zlist);
#endif
    MyOrderingProblem.solve ();

    // Now create the reordered matrix & mark it as active

    typedef ReorderFilter<row_matrix_type> reorder_filter_type;
    typedef Zoltan2::OrderingSolution<global_ordinal_type, local_ordinal_type> ordering_solution_type;
    ReorderedLocalizedMatrix_ = rcp (new reorder_filter_type (ActiveMatrix, rcp (new ordering_solution_type (*MyOrderingProblem.getSolution ()))));

    ActiveMatrix = ReorderedLocalizedMatrix_;
#else
    // This is a logic_error, not a runtime_error, because
    // setParameters() should have excluded this case already.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Ifpack2::AdditiveSchwarz::setup: "
      "The Zoltan2 and Xpetra packages must be enabled in order "
      "to support reordering.");
#endif
  }

  innerMatrix_ = ActiveMatrix;

  TEUCHOS_TEST_FOR_EXCEPTION(
    innerMatrix_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setup: Inner matrix is null right before constructing inner solver.  "
    "Please report this bug to the Ifpack2 developers.");

  // Construct the inner solver if necessary.  We go through a bit
  // more trouble than usual to do so, because we want to exercise the
  // new setInnerPreconditioner feature.
  if (Inverse_.is_null ()) {
    // FIXME (mfh 13 Dec 2013) We want to get rid of the
    // LocalInverseType template parameter.  Soon, we will add an
    // "inner preconditioner" string parameter to the input
    // ParameterList.  For now, we map statically from
    // LocalInverseType to its string name, and use the string name to
    // create the inner preconditioner.
    const std::string innerPrecName =
      Details::OneLevelPreconditionerNamer<LocalInverseType>::name ();

    // The second null argument refers to the inner preconditioner's
    // matrix.  We don't know it yet so we leave it null.
    Details::OneLevelFactory<MatrixType> factory;
    RCP<prec_type> innerPrec = factory.create (innerPrecName, Teuchos::null);
    TEUCHOS_TEST_FOR_EXCEPTION(
      innerPrec.is_null (), std::logic_error,
      "Ifpack2::AdditiveSchwarz::setup: Failed to create inner preconditioner "
      "with name \"" << innerPrecName << "\".");

    // Make sure that the new inner solver knows how to have its matrix changed.
    typedef Details::CanChangeMatrix<row_matrix_type> can_change_type;
    can_change_type* innerSolver = dynamic_cast<can_change_type*> (&*innerPrec);
    TEUCHOS_TEST_FOR_EXCEPTION(
      innerSolver == NULL, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
      "setInnerPreconditioner: The input preconditioner does not implement the "
      "setMatrix() feature.  Only input preconditioners that inherit from "
      "Ifpack2::Details::CanChangeMatrix implement this feature.");

    Inverse_ = innerPrec;
  }

  setInnerPreconditioner (Inverse_);
}


template<class MatrixType, class LocalInverseType>
void AdditiveSchwarz<MatrixType, LocalInverseType>::
setInnerPreconditioner (const Teuchos::RCP<Preconditioner<scalar_type,
                                                          local_ordinal_type,
                                                          global_ordinal_type,
                                                          node_type> >& innerPrec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    innerPrec.is_null (), std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
    "setInnerPreconditioner: Inner preconditioner must be nonnull.");

  // Make sure that the new inner solver knows how to have its matrix changed.
  typedef Details::CanChangeMatrix<row_matrix_type> can_change_type;
  can_change_type* innerSolver = dynamic_cast<can_change_type*> (&*innerPrec);
  TEUCHOS_TEST_FOR_EXCEPTION(
    innerSolver == NULL, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
    "setInnerPreconditioner: The input preconditioner does not implement the "
    "setMatrix() feature.  Only input preconditioners that inherit from "
    "Ifpack2::Details::CanChangeMatrix implement this feature.");

  // Give the local matrix to the new inner solver.
  innerSolver->setMatrix (innerMatrix_);

  // Set the new inner solver.
  Inverse_ = innerPrec;
}



} // namespace Ifpack2

#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP
