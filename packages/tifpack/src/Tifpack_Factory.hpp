/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef TIFPACK_FACTORY_HPP
#define TIFPACK_FACTORY_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_Relaxation.hpp"
#include "Tifpack_Diagonal.hpp"
#include "Tifpack_Chebyshev.hpp"
#include "Tifpack_RILUK.hpp"
#include "Tifpack_ILUT.hpp"

/** Classes and functions for templated preconditioning.  */
namespace Tifpack {

/** \brief Return true if the specified precondtioner type supports
 * unsymmetric matrices. */
bool supportsUnsymmetric(const std::string& prec_type);

//! A factory class to create Tifpack preconditioners.
/*!
Tifpack::Factory contains just one method: create().
Using Tifpack::Factory::create(), users can easily create a variety of
Tifpack preconditioners. 

create requires 3 arguments:
- a string, indicating the preconditioner to be built;
- a pointer to a Tpetra::RowMatrix, representing the matrix
  to be used to define the preconditioner;
- an optional integer (defaulted to 0), that specifies the amount of
  overlap among the processes.

The first argument can assume the following values:
- \c "DIAGONAL"  : returns an instance of Tifpack::Diagonal.
- \c "RELAXATION"  : returns an instance of Tifpack::Relaxation.
- \c "CHEBYSHEV"   : returns an instance of Tifpack::Chebyshev (overlap is ignored).
- \c "ILUT"        : returns an instance of Tifpack::ILUT.
- \c "RILUK"       : returns an instance of Tifpack::RILUK.
- otherwise, create() returns Teuchos::null.


<P> The following fragment of code shows the
basic usage of this class.
\code
#include "Tifpack_Factory.hpp"

...
typedef double Scalar;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TCrsMatrix;
typedef Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> TPrecond;
...
Tifpack::Factory<Scalar,LocalOrdinal,GlobalOrdinal,Node> Factory;

Teuchos::RCP<TCrsMatrix> A; // A is fillComplete()'d.
std::string PrecType = "ILUT"; // use incomplete LU on each process
Teuchos::RCP<TPrecond> Prec = Factory.create(PrecType, A);
assert (Prec != Teuchos::null);

Teuchos::ParameterList List;
List.set("fact: ilut level-of-fill", 5.0); // use ILUT(fill=5, drop=0)

Prec->SetParameters(List);
Prec->Initialize();
Prec->Compute();

// now Prec can be used as a preconditioner

\endcode

*/
class Factory {
public:

  /** \brief Creates an instance of Tifpack_Preconditioner given the string
   * name of the preconditioner type (throws exception if given unrecognized name).
   *
   * \param PrecType (In) - String name of preconditioner type to be created. 
   *
   * \param Matrix (In) - Matrix used to define the preconditioner
   *
   * \param overlap (In) - specified overlap, defaulted to 0.
   *
   * Returns <tt>0</tt> if the preconditioner with that input name does not
   * exist.  Otherwise, return a newly created preconditioner object.  Note
   * that the client is responsible for calling <tt>delete</tt> on the
   * returned object once it is finished using it!
   */
  template<class MatrixType>
  static
  Teuchos::RCP<Tifpack::Preconditioner<typename MatrixType::scalar_type,
                                       typename MatrixType::local_ordinal_type,
                                       typename MatrixType::global_ordinal_type,
                                       typename MatrixType::node_type> >
  create(const std::string& prec_type,
         const Teuchos::RCP<const MatrixType>& matrix,
         const int overlap = 0);
};

/////////////////////////////////////////////
/////////////////////////////////////////////

template<class MatrixType>
Teuchos::RCP<Tifpack::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Factory::create(const std::string& prec_type,
                const Teuchos::RCP<const MatrixType>& matrix,
                const int overlap)
{
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename MatrixType::mat_vec_type LocalMatVec;
  typedef typename MatrixType::mat_solve_type LocalMatSolve;
  (void)overlap;
  Teuchos::RCP<Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec;

  if (prec_type == "ILUT") {
    prec = Teuchos::rcp(new Tifpack::ILUT<MatrixType>(matrix));
  }
  else if (prec_type == "RILUK") {
    prec = Teuchos::rcp(new Tifpack::RILUK<MatrixType>(matrix));
  }
  else if (prec_type == "RELAXATION") {
    prec = Teuchos::rcp(new Tifpack::Relaxation<MatrixType>(matrix));
  }
  else if (prec_type == "CHEBYSHEV") {
    prec = Teuchos::rcp(new Tifpack::Chebyshev<MatrixType>(matrix));
  }
  else if (prec_type == "DIAGONAL") {
    prec = Teuchos::rcp(new Tifpack::Diagonal<MatrixType>(matrix));
  }
  else {
    std::ostringstream os;
    os << "Tifpack::Factory::Create ERROR, invalid preconditioner type ("
       << prec_type << ")";
    std::string str = os.str();
    throw std::runtime_error(str);
  }
  return prec;
}

} //namespace Tifpack

#endif
