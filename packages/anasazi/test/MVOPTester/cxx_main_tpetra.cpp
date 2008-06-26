//@HEADER
// ************************************************************************
// 
//
//                 Anasazi: Block Eigensolvers Package
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
// ************************************************************************
//@HEADER
//
//  This test instantiates the Anasazi classes using a complex scalar type
//  and checks functionality.
//

#include "AnasaziConfigDefs.hpp"
#include "Tpetra_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_MultiVector.hpp"

typedef int OrdinalType;
#if 1
#ifdef HAVE_COMPLEX
typedef std::complex<double> ScalarType;
#elif HAVE_COMPLEX_H
typedef ::complex<double> ScalarType;
#endif
#else
typedef double ScalarType;
#endif

#include "AnasaziOperator.hpp"
#include "AnasaziMultiVec.hpp"

namespace Anasazi {

template<class OrdinalType, class ScalarType>
class TpetraOperator : public Anasazi::Operator<ScalarType>
{
public:
  TpetraOperator(Teuchos::RCP<Tpetra::CisMatrix<OrdinalType, ScalarType> > Op) :
    Op_(Op)
  {}

  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
			    Anasazi::MultiVec<ScalarType>& Y) const
  {
    Tpetra::MultiVector<OrdinalType, ScalarType>* MyX;
    MyX = dynamic_cast<Tpetra::MultiVector<OrdinalType, ScalarType>*>(const_cast<Anasazi::MultiVec<ScalarType>*>(&X)); 
    assert (MyX != 0);
    
    Tpetra::MultiVector<OrdinalType, ScalarType>* MyY;
    MyY = dynamic_cast<Tpetra::MultiVector<OrdinalType, ScalarType>*>(&Y); 
    assert (MyY != 0);

    for (int i = 0 ; i < X.GetNumberVecs() ; ++i)
      Op_->apply(*(MyX->GetVector(i)), *(MyY->GetVector(i)), false);

    return(Anasazi::Ok);
  }

private:
  Teuchos::RCP<Tpetra::CisMatrix<OrdinalType, ScalarType> > Op_;

}; // class TpetraOperator

class TpetraMultiVec : public Anasazi::MultiVec<ScalarType>, 
                       public Tpetra::MultiVector<OrdinalType, ScalarType>
{
public:
  //@{ \name Constructors and Destructors
  //! Basic constructor
  TpetraMultiVec(const Tpetra::VectorSpace<OrdinalType, ScalarType>& vectorSpace, const int numvecs) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(vectorSpace, numvecs)
  {}

  // This is a deep copy
  TpetraMultiVec(const Tpetra::VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                 std::vector<Tpetra::Vector<OrdinalType, ScalarType> const *> list) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(vectorSpace, list)
  {}

  // This is a shallow copy
  TpetraMultiVec(const Tpetra::VectorSpace<OrdinalType, ScalarType>& vectorSpace, 
                 std::vector<Teuchos::RCP<Tpetra::Vector<OrdinalType, ScalarType> > > list) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(vectorSpace, list)
  {}

  //! Copy constructor.
  TpetraMultiVec(const Tpetra::MultiVector<OrdinalType, ScalarType>& rhs) :
    Tpetra::MultiVector<OrdinalType, ScalarType>(rhs)
  {}

  //! Destructor
  virtual ~TpetraMultiVec() {};

  //@}
  //@{ \name Creation methods

  /*! \brief Creates a new empty TpetraMultiVec containing \c numvecs columns.

    \returns Pointer to an TpetraMultiVec
  */
  MultiVec<ScalarType> * Clone (const int numvecs) const
  {
    return(new TpetraMultiVec(vectorSpace(), numvecs));
  }

  /*! \brief Creates a new TpetraMultiVec and copies contents of \c *this into
      the new vector (deep copy).

    \returns Pointer to an TpetraMultiVec
  */	
  MultiVec<ScalarType> * CloneCopy () const
  {
    return(new TpetraMultiVec(*this));
  }

  /*! \brief Creates a new TpetraMultiVec and copies the selected contents of \c *this 
      into the new vector (deep copy).  
      
    The copied vectors from \c *this are indicated by the \c index.size() indices in \c index.

    \returns Pointer to an TpetraMultiVec
  */
  MultiVec<ScalarType> * CloneCopy ( const std::vector<int>& index ) const
  {
    std::vector<Tpetra::Vector<OrdinalType, ScalarType> const*> list(index.size());
    for (int i = 0 ; i < index.size() ; ++i)
      list[i] = this->GetVector(index[i]);

    return(new TpetraMultiVec(vectorSpace(), list));
  }
    
  /*! \brief Creates a new TpetraMultiVec that shares the selected contents of \c *this.
      
    The index of the \c numvecs vectors shallow copied from \c *this are indicated by the
    indices given in \c index.
    
    \returns Pointer to an TpetraMultiVec
  */
  MultiVec<ScalarType> * CloneView ( const std::vector<int>& index )
  {
    std::vector<Teuchos::RCP<Tpetra::Vector<OrdinalType, ScalarType> > > list(index.size());
    for (int i = 0 ; i < index.size() ; ++i)
      list[i] = this->GetRCP(index[i]);

    return(new TpetraMultiVec(vectorSpace(), list));
  }

  //@}
  //@{ \name Attribute methods	

  //! Obtain the vector length of *this.
  int GetNumberVecs () const { return getNumVectors(); }

  //! Obtain the number of vectors in *this.
  int GetVecLength () const { return getNumGlobalEntries(); }

  //@}

  //@{ \name Update methods
  /*! \brief Update \c *this with \f$\alpha AB + \beta (*this)\f$.
   */
  void MvTimesMatAddMv(const ScalarType alpha, const MultiVec<ScalarType>& A, 
                       const Teuchos::SerialDenseMatrix<int,ScalarType>& B, const ScalarType beta )
  {
    // FIXME: only works in serial
    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert(MyA!=NULL);

    if (this == MyA)
    {
      // If this == A, then need additional storage ...
      // This situation is a bit peculiar but it may be required by
      // certain algorithms.
      
      std::vector<ScalarType> tmp(getNumVectors());

      for (int i = 0 ; i < getNumMyEntries() ; ++i)
      {
	for (int v = 0; v < getNumVectors() ; ++v) tmp[v] = (*MyA)(i, v);

        for (int v = 0 ; v < B.numCols() ; ++v)
        {
          (*this)(i, v) *= beta; 
          ScalarType res = 0.0;
          for (int j = 0 ; j < A.GetNumberVecs() ; ++j)
          {
            res +=  tmp[j] * B(j, v);
          }

          (*this)(i, v) += alpha * res;
        }
      }
    }
    else
    {
      for (int i = 0 ; i < getNumMyEntries() ; ++i)
      {
        for (int v = 0 ; v < B.numCols() ; ++v)
        {
          (*this)(i, v) *= beta; 
          ScalarType res = 0.0;
          for (int j = 0 ; j < A.GetNumberVecs() ; ++j)
          {
            res +=  (*MyA)(i, j) * B(j, v);
          }

          (*this)(i, v) += alpha * res;
        }
      }
    }
  }

  /*! \brief Replace \c *this with \f$\alpha A + \beta B\f$.
   */
  void MvAddMv (const ScalarType alpha, const MultiVec<ScalarType>& A, const ScalarType beta,
                const MultiVec<ScalarType>& B)
  {
    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    TpetraMultiVec* MyB;
    MyB = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(B)); 
    assert (MyB != 0);
    
    for (int v = 0 ; v < getNumVectors() ; ++v)
      for (int i = 0 ; i < getNumMyEntries() ; ++i)
	(*this)(i, v) = alpha * (*MyA)(i, v) + beta * (*MyB)(i, v);
  }

  /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$\alpha A^T(*this)\f$.
  */
  void MvTransMv (const ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B, Anasazi::ConjType = Anasazi::CONJ ) const
  {
    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetVecLength() == getNumGlobalEntries());
    assert (getNumVectors() == B.numCols());
    assert (A.GetNumberVecs() == B.numRows());
    
    // FIXME: DO ONLY HALF OF THIS...
    for (int v = 0 ; v < A.GetNumberVecs() ; ++v)
    {
      for (int w = 0 ; w < getNumVectors() ; ++w)
      {
        B(v, w) = alpha * MyA->GetVector(v)->dotProduct(*(GetVector(w)));
      }
    }
  }

  /*! \brief Compute a vector \c b where the components are the individual dot-products, i.e. \f$ b[i] = A[i]^T(this[i])\f$ where \c A[i] is the i-th column of \c A.
  */
  void MvDot (const MultiVec<ScalarType>& A, std::vector<ScalarType> &b,
              Anasazi::ConjType = Anasazi::CONJ) const
  {
    const TpetraMultiVec* MyA;
    MyA = dynamic_cast<const TpetraMultiVec*>(&A); 
    assert (MyA != 0);
    
    assert (getNumVectors() == (int)bsize());
    assert (getNumVectors() == A.GetNumberVecs());
    assert (getNumGlobalEntries() == A.GetVecLength());
    
    // hack here, it is not so good
    //vector<ScalarType> b2(getNumVectors());

    this->dotProduct(*MyA, &(b[0]));

    //for (int i = 0 ; i < getNumVectors() ; ++i)
      //b[i] = Teuchos::ScalarTraits<ScalarType>::magnitude(b2[i]);
  }

  //@}
  //@{ \name Norm method

  /*! \brief Compute the 2-norm of each individual vector of \c *this.  
    Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
    */
  void MvNorm (std::vector<Teuchos::ScalarTraits<ScalarType>::magnitudeType > &normvec) const 
  {
    assert (getNumVectors() == (int)normvec.size());
    
    norm2(&(normvec[0]));
  }
  //@}

  //@{ \name Initialization methods
  /*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  

    The \c numvecs vectors in \c A are copied to a subset of vectors in \c *this
    indicated by the indices given in \c index.
  */
  void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index )
  {
    TpetraMultiVec* MyA;
    MyA = dynamic_cast<TpetraMultiVec*>(&const_cast<Anasazi::MultiVec<ScalarType> &>(A)); 
    assert (MyA != 0);
    
    assert (A.GetNumberVecs() >= (int)index.size());
    assert (A.GetVecLength() == getNumGlobalEntries());
    
    for (unsigned int v = 0 ; v < index.size() ; ++v)
    {
      for (int i = 0 ; i < getNumMyEntries() ; ++i)
	(*this)(i, index[v])  = (*MyA)(i, v);
    }
  }

  /*! \brief Fill the vectors in \c *this with random numbers.
   */
  void MvRandom()
  {
    setAllToRandom();
  }

  /*! \brief Replace each element of the vectors in \c *this with \c alpha.
   */
  void MvInit (const ScalarType alpha)
  {
    setAllToScalar(alpha);
  }

  //@}
  //@{ \name Print method.
  /*! \brief Print \c *this TpetraMultiVec.
   */
  void MvPrint( ostream& os ) const 
  {
    Print(); // FIXME
  }
  //@}


}; // class TpetraMultiVec
}; // namespace Anasazi


// Get zero and one for the OrdinalType

OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

// Get zero and one for the ScalarType

ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

Teuchos::RCP<Anasazi::TpetraOperator<OrdinalType, ScalarType> >
  CreateMatrix(Tpetra::VectorSpace<OrdinalType, ScalarType>& vectorSpace,
               const ScalarType& a, const ScalarType& b, const ScalarType& c)
{
  Teuchos::RCP<Tpetra::CisMatrix<OrdinalType,ScalarType> > matrix =
    Teuchos::rcp(new Tpetra::CisMatrix<OrdinalType, ScalarType>(vectorSpace));

  OrdinalType NumMyElements = vectorSpace.elementSpace().getNumMyElements();
  OrdinalType NumGlobalElements = vectorSpace.elementSpace().getNumGlobalElements();
  std::vector<OrdinalType> const& MyGlobalElements = vectorSpace.elementSpace().getMyGlobalElements();

  for (OrdinalType LID = OrdinalZero ; LID < NumMyElements ; ++LID)
  {
    OrdinalType GID = MyGlobalElements[LID];
    if (GID != OrdinalZero)
      matrix->submitEntry(Tpetra::Add, GID, b, GID - 1);
    if (GID != NumGlobalElements - 1)
      matrix->submitEntry(Tpetra::Add, GID, c, GID + 1);
    matrix->submitEntry(Tpetra::Add, GID, a, GID);
  }

  matrix->fillComplete();

  return(Teuchos::rcp(new Anasazi::TpetraOperator<OrdinalType, ScalarType>(matrix)));
}




// =========== //
// main driver //
// =========== //

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockJacobiDavidson.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_RCP.hpp"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  // Creates a vector of size `length', then set the elements values.
  
  OrdinalType length    = OrdinalOne * 100;
  OrdinalType indexBase = OrdinalZero;

  // 1) Creation of a platform
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // 2) We can now create a space:

  Tpetra::ElementSpace<OrdinalType> elementSpace(length, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> 
    vectorSpace(elementSpace, platformV);


  typedef Anasazi::MultiVec<ScalarType> MV;        
  typedef Anasazi::Operator<ScalarType> OP;        
  typedef Anasazi::TpetraOperator<OrdinalType, ScalarType> TOP;        
  typedef Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

  Teuchos::RCP<TOP> A = CreateMatrix(vectorSpace, 4.0 * ScalarOne, -ScalarOne, -ScalarOne);
  Teuchos::RCP<TOP> B = CreateMatrix(vectorSpace, 2.0 * ScalarOne, -ScalarOne, -ScalarOne);
  Teuchos::RCP<TOP> K = CreateMatrix(vectorSpace, ScalarOne, ScalarZero, ScalarZero);

  Anasazi::ReturnType returnCode = Anasazi::Ok;	

  Teuchos::ScalarTraits<ScalarType>::magnitudeType TOL = 1e-4;

  int MAXITER = 500;
  int NEV = 3; 
  int SMIN = 20;
  int SMAX = 30;
  // FIXME: if I set BLOCKSIZE = 2 then nothing good happens...
  int BLOCKSIZE = 1;

  ScalarType TARGET = 1.5;

  // ================== //
  // Sets up the solver //
  // ================== //


  Teuchos::RCP<MV> ivec = Teuchos::rcp(new Anasazi::TpetraMultiVec(vectorSpace, BLOCKSIZE));        
  ivec->MvRandom();

  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<ScalarType, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<ScalarType, MV, OP>(A, B, ivec) );

  MyProblem->SetPrec(K); 

  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true); 

  // Set the number of eigenvalues requested
  MyProblem->SetNEV(NEV);

  // Inform the eigenproblem that you are finishing passing it information
  MyProblem->SetProblem();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Anasazi::OutputManager<ScalarType> > MyOM =
    Teuchos::rcp(new Anasazi::OutputManager<ScalarType>(Comm.getMyImageID()));
  MyOM->SetVerbosity(Anasazi::FinalSummary);	

  // Create a sort manager
  Teuchos::RCP<Anasazi::BasicSort<ScalarType, MV, OP> > MySM =
    Teuchos::rcp(new Anasazi::BasicSort<ScalarType, MV, OP>("SM"));


  // Create parameter list to pass into solver
  // FIXME: ADD PARAMTERS
  Teuchos::ParameterList MyPL;
  MyPL.set("BlockSize", BLOCKSIZE);
  MyPL.set("SMIN", SMIN);
  MyPL.set("SMAX", SMAX);
  MyPL.set("Max Iters", MAXITER);
  MyPL.set("Tol", TOL);
  MyPL.set("Target", TARGET);

  // Initialize the Block Jacobi-Davidson solver
  Anasazi::BlockJacobiDavidson<ScalarType, MV, OP> MySolver(MyProblem, MySM, MyOM, MyPL);
                           
  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();

  MySolver.currentStatus();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

}
