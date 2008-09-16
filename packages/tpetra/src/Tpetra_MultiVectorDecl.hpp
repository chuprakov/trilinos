// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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
// ***********************************************************************
// @HEADER

// FINISH: some of these arrayview objects should be something else, like Ptr

#ifndef TPETRA_MULTIVECTOR_DECL_HPP
#define TPETRA_MULTIVECTOR_DECL_HPP

#include <Teuchos_CompObject.hpp>
#include <Teuchos_Object.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of MultiVectorData, needed to prevent circular inclusions
  template<typename Ordinal, typename Scalar> class MultiVectorData;
#endif

  /*! multivector */
  template<class Ordinal, class Scalar>
  class MultiVector : public Teuchos::CompObject, public DistObject<Ordinal,Scalar>{

    public:

    /* FINISH: what with these?
       void replaceMap (const Map<Ordinal> &map);
       void reduce ();
     */

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Basic MultiVector constuctor.
    MultiVector(const Map<Ordinal> &map, Ordinal numVectors, bool zeroOut=true);

    //! MultiVector copy constructor.
    MultiVector(const MultiVector<Ordinal,Scalar> &source);

    //! Set multi-vector values from two-dimensional array. (copy)
    MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Scalar> &A, Ordinal LDA, Ordinal numVectors);

    //! Set multi-vector values from array of pointers. (copy)
    MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &arrayOfArrays, Ordinal numVectors);

    //! MultiVector destructor.
    virtual ~MultiVector();

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified (GlobalRow, VectorIndex) location with ScalarValue.
    void replaceGlobalValue(Ordinal globalRow, Ordinal vectorIndex, const Scalar &value);

    //! Adds ScalarValue to existing value at the specified (GlobalRow, VectorIndex) location.
    void sumIntoGlobalValue(Ordinal globalRow, Ordinal vectorIndex, const Scalar &value);

    //! Replace current value at the specified (MyRow, VectorIndex) location with ScalarValue.
    void replaceMyValue(Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue);

    //! Adds ScalarValue to existing value at the specified (MyRow, VectorIndex) location.
    void sumIntoMyValue(Ordinal MyRow, Ordinal VectorIndex, const Scalar &ScalarValue);

    //! Initialize all values in a multi-vector with constant value.
    void putScalar(const Scalar &ScalarConstant);

    //! Set multi-vector values to random numbers.
    void random();

    //@} 

    //! @name Extraction methods
    //@{ 

    /*
    //! Returns a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Ordinal,Scalar> > subCopy(const Teuchos::Range1D &colRng) const;
    */

    // Returns a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Ordinal,Scalar> > subCopy(const Teuchos::ArrayView<const Teuchos_Index> &cols) const;

    /*
    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Ordinal,Scalar> > subView(const Teuchos::Range1D &colRng);
    */

    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Ordinal,Scalar> > subView(const Teuchos::ArrayView<const Teuchos_Index> &cols);

    /*
    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Ordinal,Scalar> > subViewConst(const Teuchos::Range1D &colRng) const;
    */

    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Ordinal,Scalar> > subViewConst(const Teuchos::ArrayView<const Teuchos_Index> &cols) const;

    //! Return multi-vector values in user-provided two-dimensional array.
    void extractCopy(const Teuchos::ArrayView<Scalar> &A, Ordinal &MyLDA) const;

    //! Return multi-vector values in user-provided array of pointers.
    void extractCopy(Teuchos::ArrayView<Teuchos::ArrayView<Scalar> > arrayOfArrays) const;

    //! Return non-const pointers to multi-vector values in user-provided two-dimensional array.
    void extractView(Teuchos::ArrayRCP<Scalar> &A, Ordinal &MyLDA);

    //! Return const pointers to multi-vector values in user-provided two-dimensional array.
    void extractConstView(Teuchos::ArrayRCP<const Scalar> &A, Ordinal &MyLDA) const;

    //! Return non-const pointers to multi-vector values in user-provided array of pointers.
    void extractView(Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > &arrayOfArrays);

    //! Return const pointers to multi-vector values in user-provided array of pointers.
    void extractConstView(Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > &arrayOfArrays) const;

    //@} 

    //! @name Mathematical methods
    //@{ 

    // FINISH: expand documentation of these functions

    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    void dot(const MultiVector<Ordinal,Scalar> &A, const Teuchos::ArrayView<Scalar> &dots) const;

    //! Puts element-wise absolute values of input Multi-vector in target, this = abs(A).
    void abs(const MultiVector<Ordinal,Scalar> &A);

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<Ordinal,Scalar> &A);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const Scalar &alpha);

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    void scale(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A);

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta);

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta, const MultiVector<Ordinal,Scalar> &B, const Scalar &gamma);

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute 2-norm of each vector in multi-vector.
    void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    void normWeighted(const MultiVector<Ordinal,Scalar> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute minimum value of each vector in multi-vector.
    void minValue(Teuchos::Array<Scalar> &mins) const;

    //! Compute maximum value of each vector in multi-vector.
    void maxValue(Teuchos::Array<Scalar> &maxs) const;

    //! Compute mean (average) value of each vector in multi-vector.
    void meanValue(Teuchos::Array<Scalar> &means) const;

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const MultiVector<Ordinal,Scalar> &B, const Scalar &beta);

    //! Multiply a MultiVector with another, element-by-element: this(i,j) = beta*this(i,j) + alpha*A(i,j)*B(i,j)
    void multiply(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const MultiVector<Ordinal,Scalar> &B, const Scalar &beta);

    //! Multiply a MultiVector by the reciprocal of another, element-by-element. this(i,j) = beta*this(i,j) + alpha*B(i,j)/A(i,j)
    void reciprocalMultiply(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const MultiVector<Ordinal,Scalar> &B, const Scalar &beta);

    //@} 

    //! @name Overloaded operators
    //@{ 

    //! = Operator.
    /*! \param In 
      A - Multivector to copy
     */
    MultiVector<Ordinal,Scalar>& operator=(const MultiVector<Ordinal,Scalar> &source);

    //! Vector access function.
    /*! ArrayRCP to the local values in the ith vector of this multi-vector.
     */
    const Teuchos::ArrayRCP<Scalar> & operator[](Ordinal i);

    //! Vector access function.
    /** ArrayRCP to the local values in the ith vector of this multi-vector.
     */
    Teuchos::ArrayRCP<const Scalar> operator[](Ordinal i) const;

    //! Vector access function.
    Vector<Ordinal,Scalar> & operator()(Ordinal i);

    //! Vector access function.
    const Vector<Ordinal,Scalar> & operator() (Ordinal i) const;

    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    Ordinal numVectors() const;

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    Ordinal myLength() const;

    //! Returns the global vector length of vectors in the multi-vector.
    Ordinal globalLength() const;

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    Ordinal stride() const;

    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    bool constantStride() const;

    //@} 

    //! @name I/O methods
    //@{ 

    //! Print method.
    void print(std::ostream &os) const;
    void printValues(std::ostream &os) const;

    //@} 

    //! @name Expert-only unsupported methods
    //@{ 

/*
    //! Reset the view of an existing multivector to point to new user data.
    void resetView(const Teuchos::ArrayRCP<const Teuchos::ArrayRCP<Scalar> > &arrayOfArrays);

    //! Get pointer to MultiVector values.
    const Teuchos::ArrayRCP<Scalar> & values();

    //! Get pointer to MultiVector values.
    const Teuchos::ArrayRCP<const Scalar> & valuesConst() const;

    //! Get pointer to individual vector pointers.
    const Teuchos::ArrayRCP<const Teuchos::ArrayRCP<Scalar> > & pointers();

    //! Get pointer to individual vector pointers.
    const Teuchos::ArrayRCP<const Teuchos::ArrayRCP<const Scalar> > & pointersConst() const;
*/

    //@}

    protected:

    Teuchos::RCP<MultiVectorData<Ordinal,Scalar> > MVData_;

    // Advanced MultiVector constuctor for creating views.
    MultiVector(const Map<Ordinal> &map, Teuchos::RCP<MultiVectorData<Ordinal,Scalar> > &mvdata);

    // four functions needed for DistObject derivation
    bool checkSizes(const DistObject<Ordinal,Scalar> & sourceObj);

    void copyAndPermute(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numSameIDs,
               Ordinal numPermuteIDs,
               const Teuchos::ArrayView<const Ordinal> &permuteToLIDs,
               const Teuchos::ArrayView<const Ordinal> &permuteFromLIDs);

    void packAndPrepare(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numExportIDs,
               const Teuchos::ArrayView<const Ordinal> &exportLIDs,
               const Teuchos::ArrayView<Scalar> &exports,
               Ordinal &packetSize,
               Distributor<Ordinal> &distor);
  
    void unpackAndCombine(Ordinal numImportIDs,
               const Teuchos::ArrayView<const Ordinal> &importLIDs,
               const Teuchos::ArrayView<const Scalar> &imports,
               Distributor<Ordinal> &distor,
               CombineMode CM);


  }; // class MultiVector

} // namespace Tpetra


#endif // TPETRA_MULTIVECTOR_DECL_HPP
