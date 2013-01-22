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

//-----------------------------------------------------
// Ifpack2::KRYLOV is based on the Krylov iterative
// solvers in Belos.
// written by Paul Tsuji.
//-----------------------------------------------------

#ifndef IFPACK2_KRYLOV_DEF_HPP
#define IFPACK2_KRYLOV_DEF_HPP

namespace Ifpack2 {

//Definitions for the Krylov methods:
//==============================================================================
template <class MatrixType>
Krylov<MatrixType>::Krylov(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& A) :
  A_(A),
  Comm_(A->getRowMap()->getComm()),
  // Default values
  IterationType_(1),
  Iterations_(5),
  ResidualTolerance_(0.001),
  BlockSize_(1),
  ZeroStartingSolution_(true),
  PreconditionerType_(1),
  // General
  Condest_(-1.0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  Time_("Ifpack2::Krylov"),
  NumMyRows_(-1),
  NumGlobalNonzeros_(0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error,
      Teuchos::typeName(*this) << "::Krylov(): input matrix reference was null.");
}

//==========================================================================
template <class MatrixType>
Krylov<MatrixType>::~Krylov() {
}

//==========================================================================
template <class MatrixType>
void Krylov<MatrixType>::setParameters(const Teuchos::ParameterList& params) {
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  // Read in parameters
  Ifpack2::getParameter(params, "krylov: number of iterations",Iterations_);
  Ifpack2::getParameter(params, "krylov: residual tolerance",ResidualTolerance_);
  Ifpack2::getParameter(params, "krylov: block size",BlockSize_);
  Ifpack2::getParameter(params, "krylov: zero starting solution",ZeroStartingSolution_);
  Ifpack2::getParameter(params, "krylov: preconditioner type",PreconditionerType_);
  params_=params;
}

//==========================================================================
template <class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > &
Krylov<MatrixType>::getComm() const{
  return(Comm_);
}

//==========================================================================
template <class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Krylov<MatrixType>::getMatrix() const {
  return(A_);
}

//==========================================================================
template <class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >&
Krylov<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
template <class MatrixType>
const Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >&
Krylov<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==============================================================================
template <class MatrixType>
bool Krylov<MatrixType>::hasTransposeApply() const {
  return true;
}

//==========================================================================
template <class MatrixType>
int Krylov<MatrixType>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template <class MatrixType>
int Krylov<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template <class MatrixType>
int Krylov<MatrixType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template <class MatrixType>
double Krylov<MatrixType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class MatrixType>
double Krylov<MatrixType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class MatrixType>
double Krylov<MatrixType>::getApplyTime() const {
  return(ApplyTime_);
}

//=============================================================================
template<class MatrixType>
typename Krylov<MatrixType>::magnitude_type
Krylov<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix) {
  if (!isComputed()) { // cannot compute right now
    return(-1.0);
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -1.0) {
    Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);
  }
  return(Condest_);
}

//==========================================================================
template<class MatrixType>
void Krylov<MatrixType>::initialize() {
  // clear any previous allocation
  IsInitialized_ = false;
  IsComputed_ = false;
  Time_.start(true);
  // check only in serial
  TEUCHOS_TEST_FOR_EXCEPTION(Comm_->getSize() == 1 && A_->getNodeNumRows() != A_->getNodeNumCols(), std::runtime_error, "Ifpack2::Krylov::Initialize ERROR, matrix must be square");
  NumMyRows_ = A_->getNodeNumRows();
  // Belos parameter list
  belosList_ = Teuchos::rcp( new Teuchos::ParameterList("GMRES") );
  belosList_->set("Maximum Iterations",Iterations_);
  belosList_->set("Convergence Tolerance",ResidualTolerance_);
  //belosList_->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
  // nothing else to do here
  IsInitialized_ = true;
  ++NumInitialize_;
  Time_.stop();
  InitializeTime_ += Time_.totalElapsedTime();
}

//==========================================================================
template<class MatrixType>
void Krylov<MatrixType>::compute() {
  using Teuchos::as;
  //--------------------------------------------------------------------------
  // Ifpack2::Krylov
  // Computes necessary preconditioner info
  //--------------------------------------------------------------------------
  if (!isInitialized()) {
    initialize();
  }
  Time_.start(true);
  // Preconditioner choice
  if(PreconditionerType_==2) {
    ifpack2_prec_ = Teuchos::rcp( new Ifpack2::ILUT<MatrixType> (A_) );
  }
  else if(PreconditionerType_==3) {
    ifpack2_prec_ = Teuchos::rcp( new Ifpack2::Chebyshev<MatrixType> (A_) );
  }
  else {
    ifpack2_prec_ = Teuchos::rcp( new Ifpack2::Relaxation<MatrixType> (A_) );
  }
  ifpack2_prec_->initialize();
  ifpack2_prec_->setParameters(params_);
  ifpack2_prec_->compute();
  IsComputed_ = true;
  ++NumCompute_;
  Time_.stop();
  ComputeTime_ += Time_.totalElapsedTime();
}

//==========================================================================
template <class MatrixType>
void Krylov<MatrixType>::apply(const Tpetra::MultiVector<typename MatrixType::scalar_type,
			      typename MatrixType::local_ordinal_type,
			      typename MatrixType::global_ordinal_type,
			      typename MatrixType::node_type>& X,
			      Tpetra::MultiVector<typename MatrixType::scalar_type,
			      typename MatrixType::local_ordinal_type,
			      typename MatrixType::global_ordinal_type,
			      typename MatrixType::node_type>& Y,
			      Teuchos::ETransp mode,
			      typename MatrixType::scalar_type alpha,
			      typename MatrixType::scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
    "Ifpack2::Krylov::apply() ERROR, compute() hasn't been called yet.");

  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::Krylov::apply() ERROR, X.getNumVectors() != Y.getNumVectors().");

  Time_.start(true);

  // If X and Y are pointing to the same memory location,
  // we need to create an auxiliary vector, Xcopy
  Teuchos::RCP< const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Xcopy;
  if (X.getLocalMV().getValues() == Y.getLocalMV().getValues()) {
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(X) );
  }
  else {
    Xcopy = Teuchos::rcp( &X, false );
  }

  const Teuchos::RCP< Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Zeros = 
    Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> (A_->getRowMap(),X.getNumVectors()) );
  Zeros->putScalar((scalar_type) 0.0);

  // Belos Linear Problem
  Teuchos::RCP< Belos::LinearProblem<scalar_type,
    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
    Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > > belosProblem = 
    Teuchos::rcp( new Belos::LinearProblem<scalar_type,
		  Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
		  Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > (A_, Zeros, Xcopy) );
  belosProblem->setLeftPrec(ifpack2_prec_);
  belosProblem->setProblem(Zeros,Xcopy);
  // Belos Solver Manager
  Teuchos::RCP< Belos::SolverManager<scalar_type,
    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
    Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > > belosSolver = 
    Teuchos::rcp( new Belos::BlockGmresSolMgr<scalar_type,
		  Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
		  Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > (belosProblem, belosList_) );
  // iterative solve
  belosSolver->solve();
  Y=*Zeros;
  ++NumApply_;
  Time_.stop();
  ApplyTime_ += Time_.totalElapsedTime();
}

//=============================================================================
template <class MatrixType>
std::string Krylov<MatrixType>::description() const {
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
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}

//=============================================================================
template <class MatrixType>
void Krylov<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = Comm_->getRank();
  Teuchos::OSTab tab(out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if (vl != VERB_NONE && myImageID == 0) {
  }
}

}//namespace Ifpack2

#endif /* IFPACK2_KRYLOV_DEF_HPP */
