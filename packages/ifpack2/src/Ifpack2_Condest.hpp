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

#ifndef IFPACK2_CONDEST_HPP
#define IFPACK2_CONDEST_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_CondestType.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include <Teuchos_Ptr.hpp>

namespace Ifpack2 {

/// \fn Condest
/// \brief Estimate the condition number of the matrix.
///
/// The template parameters of this function are the same and occur in
/// the same order as the template parameters for the Preconditioner
/// class.
///
/// \param TIFP [in] The Ifpack2 preconditioner.  We need this if
///   <tt>matrix_in</tt> is null.
///
/// \param CT [in] The method to use for computing the condition
///   number estimate.  Currently, the only supported option is
///   <tt>Cheap</tt>.  Unsupported options will throw
///   <tt>std::logic_error</tt>.
///
/// \param MaxIters [in] The number of iterations used to compute the
///   condition number estimate.  Currently, this parameter is ignored.
///
/// \param Tol [in] The convergence tolerance used to compute the
///   condition number estimate.  Currently, this parameter is
///   ignored.
///
/// \param matrix_in [in] Pointer to a Tpetra::RowMatrix.  If nonnull,
///   estimate the condition number of this matrix.  If null, estimate
///   the condition number of TIFP's matrix (as returned by its
///   getMatrix() method).
///
/// The "Cheap" condition number estimate is just \f$\max_i |y_i|\f$,
/// where \f$y = A*[1, \dots, 1]^T\f$.  That is, if the input matrix
/// is \f$A\f$, we multiply it on the right by a vector of ones, and
/// return the infinity norm (maximum absolute value) of the result.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Condest (const Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>& TIFP,
	 const Ifpack2::CondestType CT,
	 const int MaxIters = 1550,
	 const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& Tol = Teuchos::as<Scalar> (1e-9),
	 const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& matrix_in = Teuchos::null)
{
  using Teuchos::Ptr;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> row_matrix_type;
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vec_type;

  MT condNumEst = -STS::one ();

  // Users may either provide a matrix for which to estimate the
  // condition number, or use the Preconditioner's built-in matrix.
  Ptr<const row_matrix_type> matrix = matrix_in;
  if (matrix_in == Teuchos::null) {
    matrix = TIFP.getMatrix ().ptr ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      matrix == Teuchos::null,
      std::logic_error,
      "Ifpack2::Condest: Both the input matrix (matrix_in) and the Ifpack2 preconditioner's matrix are null, so we have no matrix with which to compute a condition number estimate.  This probably indicates a bug in Ifpack2, since no Ifpack2::Preconditioner subclass should accept a null matrix.");
  }

  if (CT == Ifpack2::Cheap) {
#ifdef HAVE_TEUCHOS_DEBUG
    {
      MT infNorm = STS::zero ();

      // mfh 30 Jan 2013: Make sure that A is not bogus.
      std::cerr << "*** Ifpack2::Condest ***" << std::endl;

      const MT frobNormA = matrix->getFrobeniusNorm ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        Teuchos::ScalarTraits<MT>::isnaninf (frobNormA),
        std::runtime_error,
        "Ifpack2::Condest: Frobenius norm of the (original) matrix is " 
        << frobNormA << ", which is Inf or NaN.");
      vec_type X (matrix->getDomainMap ());
      vec_type Y (matrix->getRangeMap ());
      X.putScalar (STS::one ());
      Y.putScalar (STS::zero ());
      matrix->apply (X, Y);
      infNorm = Y.normInf ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        Teuchos::ScalarTraits<MT>::isnaninf (infNorm),
        std::runtime_error,
        "Ifpack2::Condest: inf-norm of A*[1,...,1]^T is " 
        << infNorm << ", which is Inf or NaN.");

      X.putScalar (STS::zero ());
      Y.putScalar (STS::zero ());
      matrix->apply (X, Y);
      infNorm = Y.normInf ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        Teuchos::ScalarTraits<MT>::isnaninf (infNorm),
        std::runtime_error,
        "Ifpack2::Condest: inf-norm of A*[0,...,0]^T is " 
        << infNorm << ", which is Inf or NaN.");

      X.putScalar (STS::zero ());
      Y.putScalar (STS::zero ());
      TIFP.apply (X, Y);
      infNorm = Y.normInf ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        Teuchos::ScalarTraits<MT>::isnaninf (infNorm),
        std::runtime_error,
        "Ifpack2::Condest: The result of applying the preconditioner "
	"to the zero vector has infinity norm " << infNorm << ", which "
	"is Inf or NaN.");
    }
#endif // HAVE_TEUCHOS_DEBUG
    vec_type ones (TIFP.getDomainMap ()); // Vector of ones
    ones.putScalar (STS::one ());
    vec_type onesResult (TIFP.getRangeMap ()); // A*ones
    onesResult.putScalar (STS::zero ());
    TIFP.apply (ones, onesResult);
    condNumEst = onesResult.normInf (); // max (abs (A*ones))
    TEUCHOS_TEST_FOR_EXCEPTION(
      STM::isnaninf (condNumEst),
      std::runtime_error,
      "Ifpack2::Condest: $\\|A*[1, ..., 1]^T\\|_{\\infty}$ = " << condNumEst << " is NaN or Inf.");
  } else if (CT == Ifpack2::CG) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Ifpack2::Condest: Condition number estimation using CG is not currently supported.");
  } else if (CT == Ifpack2::GMRES) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Ifpack2::Condest: Condition number estimation using GMRES is not currently supported.");
  }
  return condNumEst;
}

}//namespace Ifpack2

#endif // IFPACK2_CONDEST_HPP

