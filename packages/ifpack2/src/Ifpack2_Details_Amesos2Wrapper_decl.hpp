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

/// \file Ifpack2_Details_Amesos2Wrapper_decl.hpp
/// \brief Declaration of wrapper for Amesos2 solvers
///
/// \warning This header file and its contents is an implementation
///   detail of Ifpack2.  Users of Ifpack2 must NOT depend on this
///   header file continuing to exist, on its name, or on any of its
///   contents.

#ifndef IFPACK2_DETAILS_AMESOS2WRAPPER_DECL_HPP
#define IFPACK2_DETAILS_AMESOS2WRAPPER_DECL_HPP

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Details_CanChangeMatrix.hpp>

#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_AMESOS2)
#include <Amesos2_config.h>
#include <Amesos2.hpp>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Ifpack2 {
namespace Details {

/// \class Amesos2Wrapper
/// \tparam MatrixType A specialization of Tpetra::CrsMatrix.
///
/// This class computes a sparse direct factorization of a local
/// matrix using Amesos2.
///
/// The "local matrix" is the square diagonal block of the matrix
/// owned by the calling process.  Thus, if the input matrix is
/// distributed over multiple MPI processes, this preconditioner is
/// equivalent to nonoverlapping additive Schwarz domain decomposition
/// over the MPI processes, with Amesos2 as the subdomain Wrapper on
/// each process.
///
/// \warning This class is an implementation detail of Ifpack2.  Users
///   must not rely on this class.  It may go away or its interface
///   may change at any time.
///
/// \warning (mfh 10 Dec 2013) I strongly object to the need for this
///   class.  I will leave it in place because users need a way to
///   specify a sparse direct solver as a subdomain solver for
///   AdditiveSchwarz and SupportGraph.  It should be removed as soon
///   as we have a Stratimikos-like central factory solution for
///   solvers in the Tpetra-based solver stack.
///
/// \warning \c MatrixType <i>must</i> be a specialization of
///   Tpetra::CrsMatrix.  It may <i>not</i> just be a specialization
///   of Tpetra::RowMatrix.  This is a requirement of the interface of
///   Amesos2.
///
/// \note This class does <i>not</i> apply a LocalFilter to the input
///   matrix A.  This is unnecessary for subdomain solvers in
///   AdditiveSchwarz, because AdditiveSchwarz already applies a
///   LocalFilter to the matrix it passes to its subdomain solver.
///   Furthermore, some Amesos2 sparse factorizations do support MPI
///   parallelism, so it is reasonable to leave this option open to
///   users, rather than forcing a LocalFilter.
template<class MatrixType>
class Amesos2Wrapper :
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type>,
    virtual public Ifpack2::Details::CanChangeMatrix<Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                                       typename MatrixType::local_ordinal_type,
                                                                       typename MatrixType::global_ordinal_type,
                                                                       typename MatrixType::node_type> >
{
public:
  //! \name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;


  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;


  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;


  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

  //! Type of the Tpetra::RowMatrix specialization that this class uses.
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;

  //! Type of the Tpetra::Map specialization that this class uses.
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  //@}
  //! \name Constructors and destructor
  //@{

  /// \brief Constructor
  ///
  /// \param A [in] The sparse matrix to factor, as a
  ///   Tpetra::RowMatrix.  (Tpetra::CrsMatrix inherits from this, so
  ///   you may use a Tpetra::CrsMatrix here instead.)
  explicit Amesos2Wrapper (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor
  virtual ~Amesos2Wrapper();

  //@}
  //! \name Methods for setting up and computing the factorization
  //@{

  /// \brief Set parameters.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Compute the preordering and symbolic factorization of the matrix.
  ///
  /// You must call this method before you may call compute() or
  /// apply(), under the following conditions:
  /// <ol>
  /// <li> If you have not yet called initialize() before on this instance </li>
  /// <li> If you have just called setMatrix() with a nonnull matrix </li>
  /// <li> If the structure of the sparse matrix has changed </li>
  /// </ol>
  /// Please also see the documentation of compute() for the
  /// conditions under which you must call compute() before calling
  /// apply().
  void initialize ();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return IsInitialized_;
  }

  /// \brief Compute the numeric factorization of the matrix.
  ///
  /// You must call this method before you may call apply(), under the
  /// following conditions:
  /// <ol>
  /// <li> If you have not yet called compute() before on this instance </li>
  /// <li> If the values in the sparse matrix has changed </li>
  /// </ol>
  /// Please also see the documentation of initialize() for the
  /// conditions under which you must call initialize() before calling
  /// compute() or apply().
  void compute();

  //! True if compute() completed successfully, else false.
  inline bool isComputed() const {
    return IsComputed_;
  }

  //@}
  //! \name Implementation of Ifpack2::Details::CanChangeMatrix
  //@{

  /// \brief Change the matrix to be preconditioned.
  ///
  /// \param A [in] The new matrix.
  ///
  /// \post <tt>! isInitialized ()</tt>
  /// \post <tt>! isComputed ()</tt>
  ///
  /// Calling this method resets the preconditioner's state.  After
  /// calling this method with a nonnull input, you must first call
  /// initialize() and compute() (in that order) before you may call
  /// apply().
  ///
  /// You may call this method with a null input.  If A is null, then
  /// you may not call initialize() or compute() until you first call
  /// this method again with a nonnull input.  This method invalidates
  /// any previous factorization whether or not A is null, so calling
  /// setMatrix() with a null input is one way to clear the
  /// preconditioner's state (and free any memory that it may be
  /// using).
  ///
  /// The new matrix A need not necessarily have the same Maps or even
  /// the same communicator as the original matrix.
  virtual void
  setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

  //@}
  //! \name Implementation of Tpetra::Operator
  //@{

  /// \brief Apply the preconditioner to X, resulting in Y.
  ///
  /// If \f$M^{-1}\f$ represents the action of the preconditioner,
  /// then this routine computes \f$Y := beta \cdot Y + alpha \cdot
  /// Op(M^{-1}) X\f$, where \f$Op(M^{-1})\f$ can be either
  /// \f$M^{-1}\f$, \f$M^{-T}\f$, or \f$M^{-H}\f$, depending on \c
  /// mode.
  ///
  /// \param X [in] Input multivector; "right-hand side" of the solve.
  /// \param Y [out] Output multivector; result of the solve.
  /// \param mode [in] Whether to solve with the original matrix, its
  ///   transpose, or its conjugate transpose.
  /// \param alpha [in] Scaling factor for the result of applying the
  ///   preconditioner.
  /// \param alpha [in] Scaling factor for Y.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
         Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
         scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const map_type> getDomainMap () const;

  //! Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const map_type> getRangeMap () const;

  //! Whether this object's apply() method can apply the transpose (or conjugate transpose, if applicable).
  bool hasTransposeApply () const;

  //@}
  //! \name Mathematical functions
  //@{

  /// \brief Compute the condition number estimate and return its value.
  ///
  /// \warning This method is DEPRECATED.  It was inherited from
  ///   Ifpack, and Ifpack never clearly stated what this method
  ///   computes.  Furthermore, Ifpack's method just estimates the
  ///   condition number of the matrix A, and ignores the
  ///   preconditioner -- which is probably not what users thought it
  ///   did.  If there is sufficient interest, we might reintroduce
  ///   this method with a different meaning and a better algorithm.
  virtual magnitude_type TEUCHOS_DEPRECATED
  computeCondEst (CondestType CT = Cheap,
                  local_ordinal_type MaxIters = 1550,
                  magnitude_type Tol = 1e-9,
                  const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Matrix_in = Teuchos::null);

  /// \brief Return the computed condition number estimate, or -1 if not computed.
  ///
  /// \warning This method is DEPRECATED.  See warning for computeCondEst().
  virtual magnitude_type TEUCHOS_DEPRECATED getCondEst () const {
    return Condest_;
  }

  //! The input matrix's communicator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

  //! The input matrix; the matrix to be preconditioned.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  //! The total number of successful calls to initialize().
  int getNumInitialize () const;

  //! The total number of successful calls to compute().
  int getNumCompute () const;

  //! The total number of successful calls to apply().
  int getNumApply () const;

  //! The total time in seconds spent in successful calls to initialize().
  double getInitializeTime () const;

  //! The total time in seconds spent in successful calls to compute().
  double getComputeTime () const;

  //! The total time in seconds spent in successful calls to apply().
  double getApplyTime () const;

  //@}
  //! \name Implementation of Teuchos::Describable
  //@{

  //! A one-line description of this object.
  std::string description () const;

  //! Print the object with some verbosity level to the given output stream.
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef typename Teuchos::Array<local_ordinal_type>::size_type size_type;

  //! Type of the Tpetra::MultiVector specialization that this class uses.
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;

  //! Copy constructor (declared private and undefined; may not be used)
  Amesos2Wrapper (const Amesos2Wrapper<MatrixType>& RHS);

  //! operator= (declared private and undefined; may not be used)
  Amesos2Wrapper<MatrixType>& operator= (const Amesos2Wrapper<MatrixType>& RHS);

  //! Amesos2 solver; it contains the factorization of the matrix A_.
  Teuchos::RCP<Amesos2::Solver<MatrixType, MV> > amesos2solver_;

  //! The matrix to be preconditioned.
  Teuchos::RCP<const MatrixType> A_;

  //@}
  // \name Parameters (set by setParameters())
  //@{

  //@}
  // \name Other internal data
  //@{

  //! Condition number estimate (DEPRECATED; DO NOT USE)
  magnitude_type Condest_;
  //! Total time in seconds for all successful calls to initialize().
  double InitializeTime_;
  //! Total time in seconds for all successful calls to compute().
  double ComputeTime_;
  //! Total time in seconds for all successful calls to apply().
  mutable double ApplyTime_;
  //! The number of successful calls to initialize().
  int NumInitialize_;
  //! The number of successful call to compute().
  int NumCompute_;
  //! The number of successful call to apply().
  mutable int NumApply_;
  //! \c true if \c this object has been initialized
  bool IsInitialized_;
  //! \c true if \c this object has been computed
  bool IsComputed_;
  //@}
}; // class Amesos2Wrapper

} // namespace Details
} // namespace Ifpack2

#endif // HAVE_IFPACK2_AMESOS2

#endif // IFPACK2_DETAILS_AMESOS2WRAPPER_DECL_HPP
