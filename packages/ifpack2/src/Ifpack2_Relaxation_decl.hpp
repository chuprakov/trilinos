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

#ifndef IFPACK2_RELAXATION_DECL_HPP
#define IFPACK2_RELAXATION_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Parameters.hpp"

#include <Tpetra_Vector.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <string>
#include <iostream>
#include <sstream>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Ifpack2 {
enum RelaxationType {
  JACOBI,
  GS,
  SGS
};

/** \class Relaxation
\brief Relaxation preconditioners for Tpetra::RowMatrix and Tpetra::CrsMatrix sparse matrices.
\author Michael A. Heroux (Sandia)
\tparam MatrixType A specialization of Tpetra::CrsMatrix.

\section Ifpack_Relaxation_Summary Summary

This class implements several different relaxation preconditioners for
Tpetra::RowMatrix and Tpetra::CrsMatrix.  Relaxation is derived from
Preconditioner, which is itself derived from Tpetra::Operator.
Therefore this object can be used as a preconditioner everywhere an
apply() method is required in the preconditioning step.

This class implements the following relaxation methods:
- Jacobi
- Gauss-Seidel
- Symmetric Gauss-Seidel

All methods allow you to set an optional damping parameter.  The
"Gauss-Seidel" methods technically only perform Gauss-Seidel within an
MPI process, but Jacobi between processes.  To compensate, these
methods include an "L1" option, which can improve convergence by
weighting contributions near process boundaries differently.  For more
details, please refer to the following publication:

A. H. Baker, R. D. Falgout, T. V. Kolev, and U. M. Yang.
"Multigrid Smoothers for Ultraparallel Computing."
<i>SIAM J. Sci. Comput.</i>, Vol. 33, No. 5. (2011),
pp. 2864-2887.

\section Ifpack_Relaxation_Performance Performance

Jacobi will always use your matrix's native sparse matrix-vector
multiply kernel.  This should give good performance, since we have
spent a lot of effort tuning Tpetra's kernels.  Depending on the
Kokkos Node type of your Tpetra matrix, it may also exploit threads
for additional parallelism within each MPI process.  In contrast,
Gauss-Seidel and symmetric Gauss-Seidel are intrinsically sequential
methods within an MPI process.  This prevents us from exposing more
parallelism via threads.  The difference should become more apparent
as your code moves away from a "one MPI process per core" model, to a
"one MPI process per socket or node" model, assuming that you are
using a threaded Kokkos Node type.

Relaxation works with any Tpetra::RowMatrix.  If your
Tpetra::RowMatrix happens to be a Tpetra::CrsMatrix, the Gauss-Seidel
and symmetric Gauss-Seidel relaxations may be able to exploit this for
better performance.  You normally don't have to do anything to figure
this out (we test via \c dynamic_cast), but it may help to use a
Tpetra::CrsMatrix specialization as the \c MatrixType template
parameter, rather than a Tpetra::RowMatrix specialization.  (This
matters if you are using a nondefault value of the fifth template
parameter of Tpetra::CrsMatrix.)

\section Ifpack_Relaxation_Algorithms Algorithms

We now briefly describe the relaxation algorithms this class
implements.  Consider a linear system of type
\f[
A x = b,
\f]
where \f$A\f$ is a square matrix, and \f$x, b\f$ are two compatible
vectors. We begin with the decomposition
\f[
A = D - E - F
\f]
where \f$D\f$ is the diagonal of \f$A\f$, \f$-E\f$ is the strict lower
part of \f$A\f$, and \f$-F\f$ is the strict upper part of \f$A\f$.  We
assume that the diagonal entries of \f$A\f$ are all nonzero.

Given an starting solution \f$x_0\f$, an iteration of the (damped) Jacobi
method can be written in matrix form as follows:
\f[
x_{k+1} = \omega D^{-1}(E + F) x_k + D_{-1}b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter.

Users may specify the number of sweeps (\f$k_{max}\f$) and the damping
parameter \f$\omega \f$. If only one sweep is used, then this class
simply applies the inverse of the diagonal of \f$A\f to the input
vector.

Given a starting solution \f$x_0\f$, an iteration of the (damped)
Gauss-Seidel method can be written in matrix form as follows:
\f[
(D - E) x_{k+1} = \omega F x_k + b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter. Equivalently,
the Gauss-Seidel preconditioner can be defined as
\f[
P_{GS}^{-1} = (D - E)^{-1}.
\f]
Clearly, the roles of \f$E\f$ and \f$F\f$ can be interchanged.  Users
may interchange \f$E\f$ and \f$F\f$ by setting the "relaxation: backward mode"
option.

For a list of supported parameters, please refer to the documentation
of the setParameters() method.
*/
template<class MatrixType>
class Relaxation : virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {

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

  //@}
  //! \name Constructors and destructors
  //@{

  /// \brief Constructor.
  ///
  /// \param Matrix [in] The matrix for which to make the constructor.
  ///   Tpetra::RowMatrix is the base class of Tpetra::CrsMatrix, so
  ///   you may give either a Tpetra::RowMatrix or a Tpetra::CrsMatrix
  ///   here.
  ///
  /// The results of apply() are undefined if you change the diagonal
  /// entries of the sparse matrix after invoking this constructor.
  /// In particular, the compute() method may extract the diagonal
  /// entries and precompute their inverses, in order to speed up
  /// Gauss-Seidel or to implement the L1 version of various
  /// relaxation methods.  If you plan to change the diagonal entries
  /// of the matrix after making a Relaxation instance with that
  /// matrix, you must destroy the old Relaxation instance and create
  /// a new one after changing the diagonal entries.
  ///
  /// The "explicit" keyword just means that you must invoke the
  /// Relaxation constructor explicitly; you aren't allowed to use it
  /// as an implicit conversion ("cast").  For example, you may do
  /// this (namespaces and Tpetra template parameters omitted for
  /// brevity):
  /// \code
  /// RCP<const CrsMatrix<...> > A = ...;
  /// Relaxation<CrsMatrix<...> > R (A);
  /// \endcode
  /// but you may not do this:
  /// \code
  /// void foo (const Relaxation<CrsMatrix<...> >& R);
  ///
  /// RCP<const CrsMatrix<...> > A = ...;
  /// foo (A);
  /// \endcode
  explicit Relaxation(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Matrix);

  //! Destructor.
  virtual ~Relaxation();

  //@}
  // \name Preconditioner computation methods
  //@{

  /// \brief Set the relaxation / preconditioner parameters.
  ///
  /// \warning All parameters are case sensitive.  We make no attempt
  ///   to check the spelling of parameter names in your input list.
  ///
  /// The "relaxation: type" (string) parameter sets the relaxation /
  /// preconditioner method you want to use.  It currently accepts the
  /// following values (the default is "Jacobi"):
  /// - "Jacobi"
  /// - "Gauss-Seidel"
  /// - "Symmetric Gauss-Seidel"
  ///
  /// The "relaxation: sweeps" (int) parameter sets the number of
  /// sweeps, that is, the number of times to apply the relaxation on
  /// each invocation of apply().  The default number of sweeps is 1.
  ///
  /// The "relaxation: damping factor" (scalar_type -- the type of the
  /// entries of the matrix) parameter is the value of the damping
  /// factor \f$\omega \f$.  The main documentation of this class
  /// explains how we use this value.  The default value is 1.0.
  ///
  /// The "relaxation: min diagonal value" (scalar_type) parameter
  /// limits how close to zero the diagonal elements of the matrix are
  /// allowed to be.  If the magnitude of a diagonal element of the
  /// matrix is less than the magnitude of this value, then we set
  /// that diagonal element to this value.  (We don't actually modify
  /// the matrix; we just remember the diagonal values.)  The use of
  /// magnitude rather than the value itself makes this well defined
  /// if scalar_type is complex.  The default value of this parameter
  /// is zero, meaning that we do not impose a minimum diagonal value
  /// by default.
  ///
  /// The "relaxation: zero starting solution" (bool) parameter
  /// governs whether or not we use the existing values in the output
  /// multivector Y when applying the relaxation.  Its default value
  /// is true, meaning that we fill Y with zeros before applying
  /// relaxation sweeps.  If false, we use the existing values in Y.
  ///
  /// If the "relaxation: backward mode" (bool) parameter is true, we
  /// perform Gauss-Seidel in reverse mode.  The default value is
  /// false, meaning that we do forward-mode Gauss-Seidel.  This only
  /// affects standard Gauss-Seidel, not symmetric Gauss-Seidel.
  ///
  /// The last two parameters govern the L1 variant of Gauss-Seidel.
  /// The "relaxation: use l1" (bool) parameter, if true, turns on the
  /// L1 variant.  (In "l1", the first character is a lower-case L,
  /// and the second character is the numeral 1 (one).)  This
  /// parameter's value is false by default.  The "relaxation: l1 eta"
  /// (magnitude_type) parameter is the \f$\eta \f$ parameter
  /// associated with that method; its default value is 1.5.  Recall
  /// that "magnitude_type" is the type of the absolute value of a
  /// scalar_type value.  This is the same as scalar_type for
  /// real-valued floating-point types (like \c float and \c double).
  /// If scalar_type is <tt>std::complex<T></tt> for some type \c T,
  /// then magnitude_type is \c T.
  void setParameters(const Teuchos::ParameterList& params);

  //! Initialize the preconditioner.
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  //! Compute the preconditioner for the specified matrix, diagonal perturbation thresholds and relaxation parameters.
  void compute();

  //! Return true if compute() has been called.
  inline bool isComputed() const {
    return(IsComputed_);
  }

  //@}
  //! \name Implementation of the Tpetra::Operator interface
  //@{

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// This method computes Y = beta*Y + alpha*M*X, where M*X
  /// represents the action of the preconditioner on the input
  /// multivector X.
  ///
  /// \param X [in] The multivector input of the preconditioner.
  /// \param Y [in/out] The multivector output of the preconditioner.
  /// \param mode [in] Whether to apply the transpose or conjugate
  ///   transpose of the preconditioner.  Not all preconditioners
  ///   support options other than the default (no transpose); please
  ///   call hasTransposeApply() to determine whether nondefault
  ///   options are supported.
  /// \param alpha [in] Scaling factor for the preconditioned input.
  /// \param beta [in] Scaling factor for the output.
  void apply(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
             Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                 scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& getRangeMap() const;

  //! Whether apply() and applyMat() let you apply the transpose or conjugate transpose.
  bool hasTransposeApply() const;

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// This method computes Y = M*X, where M*X represents the action of
  /// the preconditioner on the input multivector X.
  ///
  /// \param X [in] The multivector input of the preconditioner.
  /// \param Y [in/out] The multivector output of the preconditioner.
  /// \param mode [in] Whether to apply the transpose or conjugate
  ///   transpose of the preconditioner.  Not all preconditioners
  ///   support options other than the default (no transpose); please
  ///   call hasTransposeApply() to determine whether nondefault
  ///   options are supported.
  void applyMat(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Mathematical functions
  //@{

  /// \brief Computes and returns the estimated condition number.
  ///
  /// We use an iterative process to estimate the condition number.
  /// You can control the number of iterations we use, the iteration
  /// tolerance, and how hard to work at estimating.
  ///
  /// \param CondestType [in] How hard to work at estimating the
  ///   condition number.  \c Cheap means not very hard.
  /// \param MaxIters [in] Maximum number of iterations for estimating
  ///   the condition number.
  /// \param Tol [in] Iteration tolerance.
  magnitude_type computeCondEst(CondestType CT = Cheap,
                               local_ordinal_type MaxIters = 1550,
                               magnitude_type Tol = 1e-9,
                               const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix = Teuchos::null);

  //@}
  //! \name Attribute accessor methods
  //@{

  /// \brief The computed estimated condition number, or -1 if not previously computed.
  ///
  /// If you haven't yet called computeCondEst(), then this method returns -1.
  magnitude_type getCondEst() const;

  //! The communicator over which the matrix and vectors are distributed.
  const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! The matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > getMatrix() const;

  //! Total number of floating-point operations over all calls to compute().
  double getComputeFlops() const;

  //! Total number of floating-point operations over all calls to apply().
  double getApplyFlops() const;

  //! Total number of calls to initialize().
  int getNumInitialize() const;

  //! Total number of calls to compute().
  int getNumCompute() const;

  //! Total number of calls to apply().
  int getNumApply() const;

  //! Total time in seconds spent in all calls to initialize().
  double getInitializeTime() const;

  //! Total time in seconds spent in all calls to compute().
  double getComputeTime() const;

  //! Total time in seconds spent in all calls to apply().
  double getApplyTime() const;

  //@}
  //! @name Implementation of Teuchos::Describable interface
  //@{

  //! A simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to a Teuchos::FancyOStream.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:

  //! @name Unimplemented methods that you are syntactically forbidden to call.
  //@{

  //! Copy constructor (not implemented; you are not allowed to call this).
  Relaxation (const Relaxation<MatrixType>& RHS);

  //! Assignment operator (not implemented; you are not allowed to call this).
  Relaxation<MatrixType>& operator= (const Relaxation<MatrixType>& RHS);

  //@}
  //! @name Internal methods
  //@{

  //! Apply Jacobi to X, returning the result in Y.
  void ApplyInverseJacobi(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Gauss-Seidel to X, returning the result in Y.
  void ApplyInverseGS(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Gauss-Seidel for a Tpetra::RowMatrix specialization.
  void ApplyInverseGS_RowMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply Gauss-Seidel for a Tpetra::CrsMatrix specialization.
  void ApplyInverseGS_CrsMatrix(
        const MatrixType& A,
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply symmetric Gauss-Seidel to X, returning the result in Y.
  void ApplyInverseSGS(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply symmetric Gauss-Seidel for a Tpetra::RowMatrix specialization.
  void ApplyInverseSGS_RowMatrix(
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //! Apply symmetric Gauss-Seidel for a Tpetra::CrsMatrix specialization.
  void ApplyInverseSGS_CrsMatrix(
        const MatrixType& A,
        const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y) const;

  //@}
  /// \name Internal data and parameters
  //@{

  //! The matrix for which to construct the preconditioner or smoother.
  const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > A_;
  //! Communicator over whose processes the matrix and vectors are distributed.
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;
  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! Importer for parallel Gauss-Seidel and symmetric Gauss-Seidel.
  Teuchos::RCP<const Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type> > Importer_;
  //! Contains the diagonal elements of \c Matrix.
  mutable Teuchos::RCP<Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Diagonal_;
  //! How many times to apply the relaxation per apply() call.
  int NumSweeps_;
  //! Which relaxation method to use.
  int PrecType_;
  //! Minimum diagonal value
  scalar_type MinDiagonalValue_;
  //! Damping factor
  scalar_type DampingFactor_;
  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  //! If true, do backward-mode Gauss-Seidel.
  bool DoBackwardGS_;
  //! If true, do the L1 version of Jacobi, Gauss-Seidel, or symmetric Gauss-Seidel.
  bool DoL1Method_;
  //! Eta parameter for modified L1 method
  magnitude_type L1Eta_;
  //! Condition number estimate
  magnitude_type Condest_;
  //! If \c true, the preconditioner has been initialized successfully.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;
  //! The number of successful calls to initialize().
  int NumInitialize_;
  //! the number of successful calls to compute().
  int NumCompute_;
  //! The number of successful calls to apply().
  mutable int NumApply_;
  //! Total time in seconds for all successful calls to initialize().
  double InitializeTime_;
  //! Total time in seconds for all successful calls to compute().
  double ComputeTime_;
  //! Total time in seconds for all successful calls to apply().
  mutable double ApplyTime_;
  //! The total number of floating-point operations for all successful calls to compute().
  double ComputeFlops_;
  //! The total number of floating-point operations for all successful calls to apply().
  mutable double ApplyFlops_;
  //! Number of local rows in the sparse matrix.
  size_t NumMyRows_;
  //! Number of global rows in the sparse matrix.
  global_size_t NumGlobalRows_;
  //! Number of global nonzeros in the sparse matrix.
  global_size_t NumGlobalNonzeros_;

  //@}

}; //class Relaxation

}//namespace Ifpack2

#endif // IFPACK2_RELAXATION_DECL_HPP

