#ifndef CTHULHU_MULTIVECTOR_HPP
#define CTHULHU_MULTIVECTOR_HPP

#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>

#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultArithmetic.hpp>

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_DistObject.hpp"
//#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Map.hpp"

#include "Cthulhu_Import.hpp"
#include "Cthulhu_Export.hpp"
#include "Cthulhu_CombineMode.hpp"

namespace Cthulhu {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Vector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class Vector;
#endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class MultiVector : public DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! MultiVector destructor.
    virtual ~MultiVector() {  }

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Initialize all values in a multi-vector with specified value.
    virtual void putScalar(const Scalar &value) =0;

    //! Set multi-vector values to random numbers.
    virtual void randomize() =0;

    //! Set seed for Random function.
    /** Note: not a method of the Tpetra interface. Added for MueLu. */
    virtual void setSeed(unsigned int seed) =0;

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    virtual Teuchos::ArrayRCP<const Scalar> getData(size_t j) const =0;

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    virtual Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j) =0;

    //@}

    //! @name Mathematical methods
    //@{ 

    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    virtual void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Teuchos::ArrayView<Scalar> &dots) const =0;

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    virtual void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) =0;

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    virtual void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) =0;

    //! Scale the current values of a multi-vector, this = alpha*this.
    virtual void scale(const Scalar &alpha) =0;

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    virtual void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta) =0;

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    virtual void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &gamma) =0;

    //! Compute 1-norm of each vector in multi-vector.
    virtual void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;

    //! Compute 2-norm of each vector in multi-vector.
    virtual void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;

    //! Compute Inf-norm of each vector in multi-vector.
    virtual void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    virtual void normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;

    //! Compute mean (average) value of each vector in multi-vector.
    virtual void meanValue(const Teuchos::ArrayView<Scalar> &means) const =0;

    // Added, not present in Tpetra
    //! Compute max value of each vector in multi-vector.
    virtual void maxValue(const Teuchos::ArrayView<Scalar> &maxs) const =0;

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    virtual void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta) =0;

    //! Element-wise multiply of a Vector A with a MultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
    virtual void elementWiseMultiply(Scalar scalarAB,
                                     const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,
                                     const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                                     Scalar scalarThis) =0;
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    virtual size_t getNumVectors() const =0;

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    virtual size_t getLocalLength() const =0;

    //! Returns the global vector length of vectors in the multi-vector.
    virtual global_size_t getGlobalLength() const =0;

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{
    /** \brief Return a simple one-line description of this object. */
    virtual std::string description() const =0;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;

    //@}

  }; // class MultiVector

} // namespace Cthulhu

#define CTHULHU_MULTIVECTOR_SHORT
#endif // CTHULHU_MULTIVECTOR_DECL_HPP
