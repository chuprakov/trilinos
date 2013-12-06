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

//-----------------------------------------------------
// Ifpack2::SupportGraph is an implementation
// of Vaidya's maximum weight spanning tree preconditioner
//------------------------------------------------------

#ifndef IFPACK2_SUPPORTGRAPH_DEF_HPP
#define IFPACK2_SUPPORTGRAPH_DEF_HPP

// Ifpack2's CMake system should (and does) prevent Trilinos from
// attempting to build or install this class, if Boost is not enabled.
// We check for this case regardless, in order to catch any bugs that
// future development might introduce in the CMake scripts.

#ifdef HAVE_IFPACK2_BOOST
#  include <boost/graph/adjacency_list.hpp>
#  include <boost/graph/kruskal_min_spanning_tree.hpp>
#  include <boost/graph/prim_minimum_spanning_tree.hpp>
#  include <boost/config.hpp>
#else
#  error "Ifpack2::SupportGraph requires that Trilinos be built with Boost support."
#endif // HAVE_IFPACK2_BOOST

#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Heap.hpp"
#include "Ifpack2_LocalFilter.hpp"
#include "Ifpack2_Parameters.hpp"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>

namespace Ifpack2 {

template <class MatrixType>
SupportGraph<MatrixType>::
SupportGraph (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero()),
  Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one()),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one()),
  Randomize_ (1),
  NumForests_ (1),
  KeepDiag_ (1.0), // FIXME (mfh 14 Nov 2013) should be either scalar_type or magnitude_type
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false)
{}


template <class MatrixType>
SupportGraph<MatrixType>::~SupportGraph () {}


template <class MatrixType>
void SupportGraph<MatrixType>::
setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Default values of the various parameters.
  magnitude_type absThresh = STM::zero ();
  magnitude_type relThresh = STM::one ();

  try {
    absThresh = params.get<magnitude_type> ("fact: absolute threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    absThresh = as<magnitude_type> (params.get<double> ("fact: absolute threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relThresh = params.get<magnitude_type> ("fact: relative threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relThresh = as<magnitude_type> (params.get<double> ("fact: relative threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try{
    Randomize_ = params.get<int> ("MST: randomize");
  }
  catch (InvalidParameterName&) {
  }

  if (absThresh != -STM::one ()) {
    Athresh_ = absThresh;
  }
  if (relThresh != -STM::one ()) {
    Rthresh_ = relThresh;
  }

  try{
    NumForests_ = params.get<int> ("MST: forest number");
  }
  catch (InvalidParameterName&) {
  }

  // FIXME (mfh 14 Nov 2013) KeepDiag_ should be either scalar_type or
  // magnitude_type, considering how it is used in the code below.
  try{
    KeepDiag_ = params.get<double> ("MST: keep diagonal");
  }
  catch (InvalidParameterName&) {
  }
}



template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
SupportGraph<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::SupportGraph::getComm: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
SupportGraph<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
SupportGraph<MatrixType>::getDomainMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ILUT::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getDomainMap ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,
                               typename MatrixType::global_ordinal_type,
                               typename MatrixType::node_type> >
SupportGraph<MatrixType>::getRangeMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::ILUT::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getRangeMap ();
}


template <class MatrixType>
bool SupportGraph<MatrixType>::hasTransposeApply() const {
  return true;
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumInitialize() const
{
  return(NumInitialize_);
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumApply() const {
  return(NumApply_);
}


template <class MatrixType>
double SupportGraph<MatrixType>::getInitializeTime() const {
  return InitializeTime_;
}


template<class MatrixType>
double SupportGraph<MatrixType>::getComputeTime() const {
  return ComputeTime_;
}


template<class MatrixType>
double SupportGraph<MatrixType>::getApplyTime() const {
  return ApplyTime_;
}


template<class MatrixType>
typename SupportGraph<MatrixType>::magnitude_type
SupportGraph<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  if (! isComputed ()) {
    return -STM::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -STM::one ()) {
    Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
  }
  return Condest_;
}


template<class MatrixType>
void SupportGraph<MatrixType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.is_null () && A->getComm ()->getSize () == 1 &&
    A->getNodeNumRows () != A->getNodeNumCols (),
    std::runtime_error, "Ifpack2::ILUT::setMatrix: If A's communicator only "
    "contains one process, then A must be square.  Instead, you provided a "
    "matrix A with " << A->getNodeNumRows () << " rows and "
    << A->getNodeNumCols () << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous support
  // graph computation.
  IsInitialized_ = false;
  IsComputed_ = false;
  A_local_ = Teuchos::null;
  Support_ = Teuchos::null;
  solver_ = Teuchos::null;
  A_ = A;
}



template<class MatrixType>
void
SupportGraph<MatrixType>::findSupport ()
{
  // FIXME (mfh 14 Nov 2013) Please don't bring in all of Boost.  Only
  // add "using" declarations for the things you need.  That avoids
  // possible collisions and also improves compilation time.
  using namespace boost;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> crs_matrix_type;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> vec_type;
  typedef std::pair<int, int> E;

  // FIXME (mfh 14 Nov 2013) The convention for typedefs is lower case
  // with words separated by underscores, followed by "type".  Hence,
  // for example, "map_type", "row_matrix_type".
  typedef adjacency_list<vecS, vecS, undirectedS, no_property,
    property<edge_weight_t, magnitude_type> > Graph;
  typedef typename graph_traits<Graph>::edge_descriptor Edge;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;

  const scalar_type zero = STS::zero ();
  const scalar_type one = STS::one ();

  //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));

  size_t num_verts = A_local_->getNodeNumRows();
  // FIXME (mfh 14 Nov 2013) Don't use int.  getNodeNumEntries()
  // returns size_t, which could be larger than int.  This could cause
  // overflow.
  int num_edges  = (A_local_->getNodeNumEntries() - A_local_->getNodeNumDiags())/2;

  // Create data structures for the BGL code
  // and temp data structures for extraction
  E *edge_array = new E[num_edges];
  magnitude_type *weights = new magnitude_type[num_edges];

  size_t num_entries;
  size_t max_num_entries = A_local_->getNodeMaxNumRowEntries();

  // FIXME (mfh 14 Nov 2013) Use std::vector or Teuchos::Array, not raw pointers.
  scalar_type *valuestemp = new scalar_type[max_num_entries];

  local_ordinal_type *indicestemp = new local_ordinal_type[max_num_entries];
  magnitude_type * diagonal = new magnitude_type[num_verts];

  for (size_t i = 0; i < max_num_entries; ++i) {
    valuestemp[i] = zero;
    indicestemp[i] = 0;
  }

  Tpetra::ArrayView<scalar_type> values (valuestemp,
                                         sizeof (scalar_type) * max_num_entries);
  Tpetra::ArrayView<local_ordinal_type> indices (indicestemp,
                                                 sizeof (local_ordinal_type) * max_num_entries);

  // Extract from the tpetra matrix keeping only one edge per pair (assume symmetric)
  //
  // FIXME (mfh 14 Nov 2013) Don't mix size_t and int.  That's mixing
  // both different sizes and different signedness.
  //
  // FIXME (mfh 14 Nov 2013) Please don't use variable names like "k"
  // for variables that are clearly significant.  Give it a name which
  // indicates its purpose, e.g., "numNegEntries" (?).
  int k = 0;
  for (size_t i = 0; i < num_verts; ++i) {
    A_local_->getLocalRowCopy (i, indices, values, num_entries);
    for (size_t j = 0; j < num_entries; ++j) {

      if(i == Teuchos::as<size_t>(indices[j])) {
        diagonal[i] = values[j];
      }

      if((i < Teuchos::as<size_t>(indices[j])) && (values[j] < 0)) {
        edge_array[k] = E(i,indices[j]);
        weights[k] = values[j];

        if (Randomize_) {
          // Add small random pertubation.
          //
          // FIXME (mfh 14 Nov 2013) "Small" should be a function of scalar_type.
          // You could use Teuchos::ScalarTraits<scalar_type> (STS) for that.
          // For example: STS::squareroot (STS::eps ()).
          weights[k] *= (1.0 + 1e-8 * drand48());
        }

        k++;
      }
    }
  }

  // Create BGL graph
  Graph g(edge_array, edge_array + num_edges, weights, num_verts);
  typedef typename property_map <Graph, edge_weight_t>::type type;
  type weight = get (edge_weight, g);
  std::vector<Edge> spanning_tree;

  // Run Kruskal, actually maximal weight ST since edges are negative
  kruskal_minimum_spanning_tree (g, std::back_inserter (spanning_tree));

  // Create array to store the exact number of non-zeros per row
  Teuchos::ArrayRCP<size_t> NumNz (num_verts, 1);

  typedef typename std::vector<Edge>::iterator edge_iterator;

  // Find the degree of all the vertices
  for (edge_iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
    local_ordinal_type localsource = source (*ei,g);
    local_ordinal_type localtarget = target (*ei,g);

    // We only want upper triangular entries, might need to swap
    if (localsource > localtarget) {
      localsource = target (*ei, g);
      localtarget = source (*ei, g);
    }

    NumNz[localsource] = NumNz[localsource] + 1;
  }


  // Create an stl vector of stl vectors to hold indices and values
  //
  // FIXME (mfh 14 Nov 2013) Don't do this.  It's perfectly valid to
  // make an std::vector<std::vector<local_ordinal_type> >, for
  // example.  If you insist on using raw pointers, just use them
  // everywhere.
  std::vector<local_ordinal_type> **Indices = new std::vector<local_ordinal_type>*[num_verts];
  std::vector<magnitude_type> **Values = new std::vector<magnitude_type>*[num_verts];

  for (size_t i = 0; i < num_verts; ++i) {
    // FIXME (mfh 14 Nov 2013) Don't do this.  It's perfectly valid to
    // make an std::vector<std::vector<local_ordinal_type> >.
    std::vector<local_ordinal_type> *temp = new std::vector<local_ordinal_type>(NumNz[i],0);
    std::vector<magnitude_type> *temp2 = new std::vector<scalar_type>(NumNz[i],0);

    Indices[i] = temp;
    Values[i] = temp2;
  }

  // The local ordering might be different from the global
  // ordering and we need the local number of non-zeros per row
  // to correctly allocate the preconditioner matrix memory
  Teuchos::ArrayRCP<size_t> localnumnz (num_verts, 1);

  for (size_t i = 0; i < num_verts; ++i) {
    // FIXME (mfh 14 Nov 2013) What is this??? Why are you calling the
    // at() method?  Why can't you just use operator[]???
    Indices[i]->at(0) = i;
  }


  // Add each spanning forest (tree) to the support graph and
  // remove it from original graph
  for (int i = 0; i < NumForests_; ++i) {
    // If a tree has already been added then we need to rerun Kruskall and
    // update the arrays containing size information
    if (i > 0) {
      spanning_tree.clear ();
      kruskal_minimum_spanning_tree (g, std::back_inserter (spanning_tree));

      for (edge_iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
        // FIXME (mfh 14 Nov 2013) Use += here.
        NumNz[source(*ei,g)] = NumNz[source(*ei,g)] + 1;
      }

      // FIXME (mfh 14 Nov 2013) Are you sure that all this resizing
      // is a good idea?
      for (size_t i = 0; i < num_verts; ++i) {
        Indices[i]->resize (NumNz[i]);
        Values[i]->resize (NumNz[i]);
      }
    }

    for (edge_iterator ei = spanning_tree.begin ();
         ei != spanning_tree.end (); ++ei) {
      local_ordinal_type localsource = source (*ei, g);
      local_ordinal_type localtarget = target (*ei, g);

      if (localsource > localtarget) {
        localsource = target(*ei,g);
        localtarget = source(*ei,g);
      }

      // FIXME (mfh 14 Nov 2013) Why are you calling the at() method???
      // Why can't you just use operator[]???

      // Assume standard Laplacian with constant row-sum.
      // Edge weights are negative, so subtract to make diagonal positive
      Values[localtarget]->at(0) = Values[localtarget]->at(0) - weight[*ei];
      Values[localsource]->at(0) = Values[localsource]->at(0) - weight[*ei];

      Indices[localsource]->at(localnumnz[localsource]) = localtarget;
      Values[localsource]->at(localnumnz[localsource]) = weight[*ei];
      localnumnz[localsource] = localnumnz[localsource] + 1;

      remove_edge (*ei,g);
    }
  }

  // Set diagonal to weighted average of Laplacian preconditioner
  // and the original matrix

  // First compute the "diagonal surplus" (in the original input matrix)
  // If input is a (pure, Dirichlet) graph Laplacian , this will be 0
  vec_type ones (A_local_->getDomainMap ());
  vec_type surplus (A_local_->getRangeMap ());

  ones.putScalar (one);
  A_local_->apply (ones, surplus);

  Teuchos::ArrayRCP<const scalar_type> surplusaccess = surplus.getData(0);

  for (size_t i = 0; i < num_verts; ++i) {
    if (surplusaccess[i] > 0) {
      // FIXME (mfh 14 Nov 2013) Why are you calling the at() method???
      // Why can't you just use operator[]???
      Values[i]->at(0) += surplusaccess[i];
    }

    // If the original diagonal is less than the row sum then we aren't going to use it
    // regardless of the diagonal option, shouldn't happen for proper Laplacian
    //
    // FIXME (mfh 14 Nov 2013) Why are you calling the at() method???
    // Why can't you just use operator[]???
    if (diagonal[i] < Values[i]->at(0)) {
      diagonal[i] = Values[i]->at(0);
    }

    // FIXME (mfh 14 Nov 2013) Why are you calling the at() method???
    // Why can't you just use operator[]???
    //
    // FIXME (mfh 14 Nov 2013) You went through the trouble of
    // defining scalar_type one at the top of this method; why are you
    // using "1." here then?
    Values[i]->at(0) = KeepDiag_*diagonal[i] + (1.-KeepDiag_) * Values[i]->at(0);

    // Modify the diagonal with user specified scaling
    //
    // FIXME (mfh 14 Nov 2013) Why are you calling the at() method???
    // Why can't you just use operator[]???
    if (Rthresh_) {
      Values[i]->at(0) *= Rthresh_;
    }
    if (Athresh_) {
      Values[i]->at(0) += Athresh_;
    }
  }

  // Create the CrsMatrix for the support graph
  Support_ = rcp (new crs_matrix_type (A_local_->getRowMap (),
                                       A_local_->getColMap (),
                                       localnumnz, Tpetra::StaticProfile));

  // Fill in the matrix with the stl vectors for each row
  for (size_t i = 0; i < num_verts; ++i) {
    Teuchos::ArrayView<local_ordinal_type>
      IndicesInsert (*Indices[Teuchos::as<local_ordinal_type>(i)]);
    Teuchos::ArrayView<scalar_type>
      ValuesInsert (*Values[Teuchos::as<local_ordinal_type>(i)]);
    Support_->insertLocalValues (i, IndicesInsert, ValuesInsert);
  }

  Support_->fillComplete();

  // Clean up all the memory allocated
  //
  // FIXME (mfh 14 Nov 2013) Please use std::vector or Teuchos::Array
  // instead of raw arrays.
  delete edge_array;
  delete weights;
  delete diagonal;
  delete Values;
  delete Indices;
}

template<class MatrixType>
Teuchos::RCP<const typename SupportGraph<MatrixType>::row_matrix_type>
SupportGraph<MatrixType>::
makeLocalFilter (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A->getComm ()->getSize () > 1) {
    return Teuchos::rcp (new LocalFilter<MatrixType> (A));
  } else {
    return A;
  }
}

template<class MatrixType>
void SupportGraph<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  // Create a timer for this method, if it doesn't exist already.
  // TimeMonitor::getNewCounter registers the timer, so that
  // TimeMonitor's class methods like summarize() will report the
  // total time spent in successful calls to this method.
  const std::string timerName ("Ifpack2::SupportGraph::initialize");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::SupportGraph::initialize: "
      "The matrix to precondition is null.  Please call setMatrix() with a "
      "nonnull input before calling this method.");

    // Clear any previous computations.
    IsInitialized_ = false;
    IsComputed_ = false;
    A_local_ = Teuchos::null;
    Support_ = Teuchos::null;
    solver_ = Teuchos::null;

    A_local_ = makeLocalFilter (A_); // Compute the local filter.
    findSupport (); // Compute the support.

    // Set up the solver and compute the symbolic factorization.
    solver_ = Amesos2::create<MatrixType, MV> ("amesos2_cholmod", Support_);
    solver_->symbolicFactorization ();

    IsInitialized_ = true;
    ++NumInitialize_;
  } // Stop timing here.

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  InitializeTime_ = timer->totalElapsedTime ();
}



template<class MatrixType>
void SupportGraph<MatrixType>::compute () {
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  // Don't count initialization in the compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  // Create a timer for this method, if it doesn't exist already.
  // TimeMonitor::getNewCounter registers the timer, so that
  // TimeMonitor's class methods like summarize() will report the
  // total time spent in successful calls to this method.
  const std::string timerName ("Ifpack2::SupportGraph::compute");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    Teuchos::TimeMonitor timeMon (*timer);
    solver_->numericFactorization ();
    IsComputed_ = true;
    ++NumCompute_;
  } // Stop timing here.

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ComputeTime_ = timer->totalElapsedTime ();
}


template <class MatrixType>
void
SupportGraph<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,
                                 local_ordinal_type,
                                 global_ordinal_type,
                                 node_type>& X,
       Tpetra::MultiVector<scalar_type,
                           local_ordinal_type,
                           global_ordinal_type,
                           node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  typedef scalar_type DomainScalar;
  typedef scalar_type RangeScalar;
  typedef Tpetra::MultiVector<DomainScalar, local_ordinal_type,
    global_ordinal_type, node_type> MV;

  RCP<FancyOStream> out = getFancyOStream (rcpFromRef (std::cout));

  // Create a timer for this method, if it doesn't exist already.
  // TimeMonitor::getNewCounter registers the timer, so that
  // TimeMonitor's class methods like summarize() will report the
  // total time spent in successful calls to this method.
  const std::string timerName ("Ifpack2::SupportGraph::apply");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    Teuchos::TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed (), std::runtime_error,
      "Ifpack2::SupportGraph::apply: You must call compute() to compute the incomplete "
      "factorization, before calling apply().");

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::SupportGraph::apply: X and Y must have the same number of columns.  "
      "X has " << X.getNumVectors () << " columns, but Y has "
      << Y.getNumVectors () << " columns.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      beta != STS::zero (), std::logic_error,
      "Ifpack2::SupportGraph::apply: This method does not currently work when beta != 0.");

    // If X and Y are pointing to the same memory location,
    // we need to create an auxiliary vector, Xcopy
    RCP<const MV> Xcopy;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      Xcopy = rcp (new MV (X));
    }
    else {
      Xcopy = rcpFromRef (X);
    }

    if (alpha != STS::one ()) {
      Y.scale (alpha);
    }

    RCP<MV> Ycopy = rcpFromRef(Y);

    solver_->setB (Xcopy);
    solver_->setX (Ycopy);

    solver_->solve ();
  } // Stop timing here.

  ++NumApply_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ApplyTime_ = timer->totalElapsedTime ();
}


template <class MatrixType>
std::string SupportGraph<MatrixType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status: [initialized, computed]";
    }
    else {
      oss << "{status: [initialized, not computed]";
    }
  }
  else {
    oss << "{status: [not initialized, not computed]";
  }

  if (A_.is_null ()) {
    oss << ", A_: null";
  }
  else {
    oss << ", A_: nonnull, "
        << ", global number of rows: " << A_->getGlobalNumRows ()
        << ", global number of columns: " << A_->getGlobalNumCols ()
        << "}";
  }
  return oss.str();
}


template <class MatrixType>
void SupportGraph<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
  Teuchos::OSTab tab (out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if(vl != VERB_NONE && getComm()->getRank() == 0) {
    out << this->description() << endl;
    out << endl;
    out << "==============================================================================="
        << endl;
    out << "Absolute threshold: " << getAbsoluteThreshold () << endl;
    out << "Relative threshold: " << getRelativeThreshold () << endl;

    out << "Condition number estimate: " << Condest_ << endl;

    if (isComputed ()) {
      out << "Number of nonzeros in A: " << A_->getGlobalNumEntries() << endl;
      out << "Number of nonzeros in A_local: " << A_local_->getGlobalNumEntries() << endl;
      out << "Number of edges in support graph: "
          << Support_->getGlobalNumEntries () - Support_->getGlobalNumDiags () << endl;

      const double popFrac =
        static_cast<double> (Support_->getGlobalNumEntries () - Support_->getGlobalNumDiags ()) /
        ((A_->getGlobalNumEntries () - A_->getGlobalNumDiags ()) / 2.0);

      out << "Fraction of off diagonals of supportgraph/off diagonals of original: "
          << popFrac << endl;
    }
    out << endl;
    out << "Phase           # calls    Total Time (s) " << endl;
    out << "------------    -------    ---------------" << endl;
    out << "initialize()    " << setw(7) << getNumInitialize() << "    "
        << setw(15) << getInitializeTime() << endl;
    out << "compute()       " << setw(7) << getNumCompute()    << "    "
        << setw(15) << getComputeTime()    << endl;
    out << "apply()         " << setw(7) << getNumApply()      << "    "
        << setw(15) << getApplyTime()      << endl;
    out << "==============================================================================="
        << endl;
    out << endl;

    solver_->printTiming (out, verbLevel);
  }
}


}//namespace Ifpack2

#endif /* IFPACK2_SUPPORTGRAPH_DEF_HPP */

