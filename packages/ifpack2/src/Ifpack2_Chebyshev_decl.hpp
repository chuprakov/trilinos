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

#ifndef IFPACK2_CHEBYSHEV_DECL_HPP
#define IFPACK2_CHEBYSHEV_DECL_HPP


#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Parameters.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <iostream>
#include <string>
#include <sstream>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Ifpack2 {

/// \class Chebyshev 
/// \brief Diagonally scaled Chebyshev iteration for Tpetra sparse matrices.
/// \tparam MatrixType A specialization of Tpetra::CrsMatrix.
///
/// \section Ifpack_Chebyshev_Summary Summary
///
/// This class implements a Chebyshev polynomial preconditioner or
/// smoother for a Tpetra sparse matrix.  Given a matrix A, it applies
/// Chebyshev iteration to the left-scaled matrix \f$D^{-1} A\f$,
/// where D = diag(A) is the matrix of the diagonal entries of A.
/// This class accepts either a Tpetra::RowMatrix or a
/// Tpetra::CrsMatrix; its template parameter must be a specialization
/// of Tpetra::CrsMatrix.
///
/// Chebyshev is derived from Preconditioner, which itself is derived
/// from Tpetra::Operator.  Therefore, a Chebyshev instance may be
/// used as an operator in any code that invokes the operator as
/// apply().
///
/// \warning Our implementation currently <i>only</i> works with a
///   real symmetric positive definite (SPD) matrix.  Results for
///   matrices that are not SPD, or for complex-valued Scalar types,
///   are not defined.
///
/// \section Ifpack_Chebyshev_Algorithm Algorithm
///
/// Given a matrix A, a right-hand side X, and an initial guess Y,
/// this class performs Chebyshev iteration using the left-scaled
/// matrix \f$D^{-1} A\f$, where D is the matrix of the diagonal
/// elements of A.  (You may control left scaling yourself if you
/// wish, by providing an optional vector of the entries of
/// \f$D^{-1}\f$.)  While Chebyshev iteration works for any matrix, we
/// have chosen only to allow real-valued, symmetric positive definite
/// matrices.
///
/// Chebyshev iteration was originally intended as an iterative solver
/// for linear systems.  See the following publication:
///
/// T. Manteuffel, "The Tchebychev iteration for nonsymmetric linear
/// systems," Numer. Math., 28 (1977), pp. 307-327.
///
/// It also works as a smoother for algebraic multigrid, which is the
/// target use case of this implementation.
///
/// \section Ifpack_Chebyshev_Eig Eigenvalue bounds
///
/// We require that the input matrix A be real-valued and symmetric
/// positive definite.  Thus, all of its eigenvalues must lie in a
/// positive interval on the real line.  Furthermore, if D is the
/// matrix of the diagonal elements of A, then the same is true of
/// \f$D^{-1} A\f$.  
///
/// Suppose \f$[\lambda_{min}, \lambda_{max}]\f$ is the interval of
/// the eigenvalues of \f$D^{-1} A\f$.  We require users to give us at
/// least (an estimate of) the maximum eigenvalue \f$\lambda_{max}\f$.
/// They may optionally also give us the (estimated) ratio \f$\eta =
/// \lambda_{max} / \lambda_{min}\f$, or (an estimate of) the minimum
/// eigenvalue \f$\lambda_{min}\f$.  The \f$\eta\f$ parameter
/// corresponds to the "smoother: Chebyshev alpha" parameter of ML.
/// (We use "eta" instead of "alpha" to avoid confusion with the
/// "alpha" argument of the apply() method of Tpetra::Operator.)
///
/// When using Chebyshev iteration by itself to solve linear systems,
/// it is important to have good estimates of both the minimum and
/// maximum eigenvalues.  However, when using a small number of
/// Chebyshev iterations as a smoother in multigrid, the maximum
/// eigenvalue estimate is more important.  (The point of a smoother
/// is to smooth out the high-frequency components of the error, that
/// is, those that correspond to the largest eigenvalues.  The coarser
/// grids below the current grid will take care of the lower-frequency
/// components of the error.)  This is why we use a ratio \f$\eta =
/// \lambda_{max} / \lambda_{min}$, rather than requiring a guess for
/// $\lambda_{min}$.  In fact, we only use \f$\lambda_{min}\f$ for
/// error checking, not when determining the Chebyshev coefficients.
/// Often, if users give us \f$\lambda_{max}\f$, our default value of
/// \f$\eta\f$ suffices.
///
/// Underestimating \f$\lambda_{min}\f$ may make Chebyshev fail to
/// converge, or fail to reduce the highest-frequency components of
/// the error, if used as a smoother.  Thus, we always multiply the
/// given \f$\lambda_{min}\f$ by a small factor (1.1).  This heuristic
/// accounts for the fact that typical methods for estimating extremal
/// eigenvalues (like Lanczos or CG) underestimate them.
///
/// Some Chebyshev implementations attempt to estimate the eigenvalue
/// interval automatically.  We do not attempt to do this.  Users who
/// want estimates of both the smallest and largest eigenvalues should
/// run a few iterations of Lanczos or CG.  For just the largest
/// eigenvalue, a few iterations of the power method may be enough.
/// For an example of a Chebyshev implementation that estimates
/// eigenvalue bounds, see Steve Ashby's CHEBYCODE:
///
/// S. ASHBY, "CHEBYCODE: A Fortran implementation of Manteuffel's
/// adaptive Chebyshev algorithm," Tech. Rep. UIUCDCS-R-85-1203,
/// University of Illinois, 1985.
///
/// \section Ifpack_Chebyshev_Params Setting parameters
///
/// Call the setParameters() method to give this instance your
/// estimates of \f$\lambda_{max}\f$ and \f$\eta = \lambda_{max} /
/// \lambda_{min}\f$, as well as to set other options controlling the
/// behavior of Chebyshev iteration.  The documentation of
/// setParameters() lists all the parameters that this class accepts.
/// Where possible, we list comparable parameters in the Ifpack
/// package and the ML multigrid package.
///
/// \section Ifpack_Chebyshev_Performance Performance
///
/// Chebyshev should spend most of its time in Tpetra's native sparse
/// matrix-vector multiply kernel.  This should give good performance,
/// since we have spent a lot of effort tuning that kernel.  Depending
/// on the Kokkos Node type of your Tpetra matrix, the kernel may also
/// exploit threads for additional parallelism within each MPI process
/// ("hybrid parallelism" a.k.a. "MPI + X").  If your application
/// depends on hybrid parallelism for performance, you should favor
/// Chebyshev smoothers whenever possible over "serial within a
/// process" smoothers like Gauss-Seidel or SOR (Symmetric
/// Over-Relaxation).
///
/// \section Ifpack_Chebyshev_History History
///
/// The original implementation of this class was an adaptation of
/// ML's ML_Cheby routine.  The original author was Ulrich Hetmaniuk,
/// a Sandia employee in what was then (2006) Org 1416.  Ifpack2 has
/// seen significant development since then.
template<class MatrixType>
class Chebyshev : 
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
					   typename MatrixType::local_ordinal_type,
					   typename MatrixType::global_ordinal_type,
					   typename MatrixType::node_type> 
{
public:
  //! \name Typedefs
  //@{ 
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  //@}
  // \name Constructors and destructors
  //@{

  /// \brief Constructor.
  ///
  /// \param A [in] The sparse matrix to which to apply Chebyshev
  ///   iteration.  The matrix A must be square, and its domain Map
  ///   and range Map must be the same.  The latter means that the
  ///   vectors x and y in the sparse matrix-vector product y = A*x
  ///   must both have the same distribution over process(es).
  ///
  /// We currently also require that the row Map and the range Map of
  /// A be the same.  This is an arbitrary requirement that we plan to
  /// relax eventually.  However, respecting this requirement will
  /// reduce set-up costs.
  ///
  /// The constructor will only check the requirements on the various
  /// Maps of A if the CMake configuration option
  /// <tt>Teuchos_ENABLE_DEBUG</tt> was set to <tt>ON</tt> before
  /// building Trilinos.  The checks require \f$O(1)\f$ global
  /// reductions over all processes in A's communicator, so we prefer
  /// to avoid them if we can.
  explicit Chebyshev (const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A);

  //! Destructor.
  virtual ~Chebyshev ();

  //@}
  //! \name Preconditioner computation methods
  //@{ 

  /// \brief Set parameters for the preconditioner.
  ///
  /// The following parameters control the Chebyshev coefficients, the
  /// number of iterations, various other properties of the algorithm,
  /// and error checking.  Where applicable, for each parameter, we
  /// list the name of its corresponding ML parameter.
  /// - "chebyshev: max eigenvalue" (\c Scalar): (An estimate of) the
  ///   largest eigenvalue \f$\lambda_{max}\f$ of the matrix.  You
  ///   should always provide this value, since otherwise Chebyshev
  ///   will not be an effective smoother.
  /// - "chebyshev: ratio eigenvalue" (\c Scalar): (Estimated) ratio
  ///   \f$\eta\f$ between the largest and smallest eigenvalue.  The
  ///   default value is 30.0.  We use this to define the interval on
  ///   the real line that determines the Chebyshev coefficients.
  /// - "chebyshev: min eigenvalue" (\c Scalar): (An estimate of) the
  ///   smallest eigenvalue \f$\lambda_{min}\f$ of the matrix.  This
  ///   parameter is optional and is only used to check whether the
  ///   input matrix is the identity matrix.  We do <i>not</i> use
  ///   this to define the interval on the real line that determines
  ///   the Chebyshev coefficients.
  /// - "chebyshev: degree" (\c int): The polynomial degree; the
  ///   number of times that apply() invokes sparse matrix-vector
  ///   multiply.
  /// - "chebyshev: min diagonal value" (\c Scalar): Threshold for
  ///   diagonal values.  Values smaller than this are not inverted.
  /// - "chebyshev: zero starting solution" (\c bool): If true, the
  ///   input vector(s) is/are filled with zeros on entry to the
  ///   apply() method.
  /// - "chebyshev: operator inv diagonal" (<tt>Tpetra::Vector</tt>).
  ///   A (raw) pointer to the inverse of the diagonal entries of the
  ///   matrix, stored as a <tt>Tpetra::Vector<Scalar, LocalOrdinal,
  ///   GlobalOrdinal, Node></tt>.  If not provided, we compute this
  ///   ourselves from the matrix provided to the constructor.  If
  ///   provided here, we make a deep copy.
  ///
  /// The following list maps from an ML parameter to its
  /// corresponding Ifpack2 parameter.
  /// - "smoother: Chebyshev alpha": "chebyshev: ratio eigenvalue"
  /// - "smoother: sweeps": "chebyshev: degree"
  ///
  /// ML does not have a parameter corresponding to "chebyshev: max
  /// eigenvalue", because ML estimates the spectral radius
  /// automatically.  Ifpack2 does not do this automatically, but a
  /// multigrid package such as MueLu that uses Ifpack2 could do this.
  /// Similarly, ML does not have a parameter corresponding to
  /// "chebyshev: min eigenvalue".
  ///
  /// The following list maps from an Ifpack parameter to its
  /// corresponding Ifpack2 parameter.  Many of the parameters have
  /// the same names, in which case we simply write <i>same</i>.
  /// - "chebyshev: max eigenvalue": same
  /// - "chebyshev: ratio eigenvalue": same
  /// - "chebyshev: min eigenvalue": same
  /// - "chebyshev: degree": same
  /// - "relaxation: sweeps": "chebyshev: degree"
  /// - "chebyshev: min diagonal value": same
  /// - "relaxation: min diagonal value": "chebyshev: min diagonal value"
  /// - "chebyshev: zero starting solution": same
  /// - "relaxation: zero starting solution": "chebyshev: zero starting solution"
  /// - "chebyshev: operator inv diagonal": same
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Initialize the preconditioner.
  ///
  /// The compute() method will call initialize() automatically if it
  /// has not yet been called, so you do not normally need to call
  /// this.  However, it is correct to call initialize() yourself, and
  /// compute() will not call it again if it already has been called.
  void initialize();

  /// Whether the preconditioner has been successfully initialized
  /// (by calling initialize()).
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  /// \brief Compute the preconditioner.
  ///
  /// You must call this method before calling apply().  If any of the
  /// values of the diagonal entries of the matrix have changed, or if
  /// the structure of the matrix has changed, you must call this
  /// method again before calling apply().
  ///
  /// This method will call initialize() if it has not already been
  /// called.  However, you may call initialize() before calling this
  /// method if you wish.
  void compute();

  /// Whether the preconditioner has been successfully computed
  /// (by calling compute()).
  inline bool isComputed() const {
    return IsComputed_;
  }

  //@}
  //! @name Implementation of Tpetra::Operator
  //@{ 

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// This method actually computes Y = beta*Y + alpha*(M*X), where
  /// M*X represents the result of Chebyshev iteration on X, using the
  /// matrix Op(A).  Op(A) is A if <tt>mode</tt> is
  /// <tt>Teuchos::NO_TRANS</tt>, \f$A^T\f$ if <tt>mode</tt> is
  /// <tt>Teuchos::TRANS</tt>, and \f$A^H\f$ (the Hermitian transpose)
  /// if <tt>mode</tt> is <tt>Teuchos::CONJ_TRANS</tt>.
  ///
  /// Since this class currently requires A to be real and symmetric
  /// positive definite, setting <tt>mode</tt> should not affect the
  /// result.
  ///
  /// \param X [in] A (multi)vector to which to apply the preconditioner.
  /// \param Y [in/out] A (multi)vector containing the result of
  ///   applying the preconditioner to X.
  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                 Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! The Tpetra::Map representing the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const;

  //! The Tpetra::Map representing the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const;

  //! Whether it's possible to apply the transpose of this operator.
  bool hasTransposeApply() const;

  /// \brief Compute Y = Op(A)*X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  ///
  /// \param X [in] Input (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, X must be in the
  ///   domain Map of the matrix A.  Otherwise, X must be in the range
  ///   Map of A.
  /// \param Y [out] Output (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, Y must be in the
  ///   range Map of the matrix A.  Otherwise, Y must be in the domain
  ///   Map of A.
  /// \param mode [in] Whether to apply the matrix A, its transpose
  ///   \f$A^T\f$, or its conjugate transpose \f$A^H\f$.  Op(A) is A
  ///   if <tt>mode</tt> is <tt>Teuchos::NO_TRANS</tt>, \f$A^T\f$ if
  ///   <tt>mode</tt> is <tt>Teuchos::TRANS</tt>, and \f$A^H\f$ (the
  ///   Hermitian transpose) if <tt>mode</tt> is
  ///   <tt>Teuchos::CONJ_TRANS</tt>.
  ///
  /// Since this class currently requires A to be real and symmetric
  /// positive definite, setting <tt>mode</tt> should not affect the
  /// result.
  void applyMat(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                      Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Mathematical functions
  //@{

  //! Compute and return the estimated condition number.
  magnitudeType 
  computeCondEst (CondestType CT = Cheap, 
		  LocalOrdinal MaxIters = 1550,
		  magnitudeType Tol = 1e-9,
		  const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix = Teuchos::null);

  //@}
  //! \name Attribute accessor methods
  //@{ 

  //! The estimated condition number, or -1.0 if it has not yet been computed.
  magnitudeType getCondEst() const;

  //! The communicator over which the matrix is distributed.
  const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! The matrix for which this is a preconditioner.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getMatrix() const;

  //! The total number of floating-point operations taken by all calls to compute().
  double getComputeFlops() const;

  //! The total number of floating-point operations taken by all calls to apply().
  double getApplyFlops() const;

  //! The total number of successful calls to initialize().
  int getNumInitialize() const;

  //! The total number of successful calls to compute().
  int getNumCompute() const;

  //! The total number of successful calls to apply().
  int getNumApply() const;

  //! The total time spent in all calls to initialize().
  double getInitializeTime() const;

  //! The total time spent in all calls to compute().
  double getComputeTime() const;

  //! The total time spent in all calls to apply().
  double getApplyTime() const;

  //@}
  //! @name Implementation of Teuchos::Describable 
  //@{

  //! A simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to a Teuchos::FancyOStream.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}
  //! \name Utility methods
  //@{

  //! Simple power method to compute lambda_max.
  static void PowerMethod(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator,
                         const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal,
                         const int MaximumIterations, 
                         Scalar& LambdaMax);

  //! Not currently implemented: Use CG to estimate lambda_min and lambda_max.
  static void CG(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator, 
                const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal, 
                const int MaximumIterations, 
                Scalar& lambda_min, Scalar& lambda_max);

  //@}

private:
  
  //! Copy constructor (use is syntactically forbidden)
  Chebyshev (const Chebyshev<MatrixType>&);

  //! Assignment operator (use is syntactically forbidded)
  Chebyshev<MatrixType>& operator= (const Chebyshev<MatrixType>&);

  //@}
  //! \name The sparse matrix and related data
  //@{

  //! The matrix to be preconditioned
  const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! The communicator over which the matrix is distributed
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;
  //! The inverse of the diagonal elements of the matrix
  mutable Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > InvDiagonal_;
  //! Number of local rows in the matrix.
  size_t NumMyRows_;
  //! Number of global rows in the matrix.
  global_size_t NumGlobalRows_;
  //! Number of global nonzeros in the matrix.
  global_size_t NumGlobalNonzeros_;

  //@}
  //! \name Algorithmic parameters (set via setParameters())
  //@{

  //! The number of iterations to apply; the degree of the Chebyshev polynomial.
  int PolyDegree_;
  //! Estimate of the ratio LambdaMax_ / LambdaMin_.
  magnitudeType EigRatio_;
  //! Approximation of the smallest eigenvalue.
  Scalar LambdaMin_;
  //! Approximation of the largest eigenvalue.
  Scalar LambdaMax_;
  //! Minimum allowed value on the diagonal of the matrix.
  Scalar MinDiagonalValue_;
  /// If \c true, then the starting solution is always the zero vector.
  bool ZeroStartingSolution_;

  //@}
  //! \name Other internal state
  //@{

  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! The estimated condition number.
  magnitudeType Condest_;
  //! If \c true, initialize() has completed successfully.
  bool IsInitialized_;
  //! If \c true, compute() has completed successfully.
  bool IsComputed_;
  //! The total number of successful calls to initialize().
  int NumInitialize_;
  //! The total number of successful calls to compute().
  int NumCompute_;
  /// \brief The total number of successful calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable int NumApply_;
  //! The total time in seconds over all calls to initialize().
  double InitializeTime_;
  //! The total time in seconds over all calls to compute().
  double ComputeTime_;
  /// \brief The total time in seconds over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyTime_;
  //! The total number of floating-point operations over all calls to compute().
  double ComputeFlops_;
  /// \brief The total number of floating-point operations over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyFlops_;

  //@}
}; // class Chebyshev

}//namespace Ifpack2

#endif // IFPACK2_CHEBYSHEV_DECL_HPP

