#ifndef CTHULHU_TPETRAMULTIVECTOR_HPP
#define CTHULHU_TPETRAMULTIVECTOR_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include "Cthulhu_Vector.hpp"
#include "Cthulhu_MultiVector.hpp"

#include "Cthulhu_TpetraMap.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

#include "Cthulhu_Exceptions.hpp"
#include "Cthulhu_Debug.hpp"

namespace Cthulhu {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of TpetraVector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class TpetraVector;
#endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class TpetraMultiVector : public virtual Cthulhu::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

    // The following typedef are used by the CTHULHU_DYNAMIC_CAST() macro.
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Basic TpetraMultiVector constuctor.
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t NumVectors, bool zeroOut=true) {
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMap, map, tMap, "Cthulhu::TpetraMultiVector constructors only accept Cthulhu::TpetraMap as input arguments.");
      vec_ = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tMap->getTpetra_Map(), NumVectors, zeroOut));
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! TpetraMultiVector copy constructor.
    TpetraMultiVector(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source){ CTHULHU_DEBUG_ME; }
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Scalar> &A, size_t LDA, size_t NumVectors) {
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMap, map, tMap, "Cthulhu::TpetraMultiVector constructors only accept Cthulhu::TpetraMap as input arguments.");
      vec_ = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tMap->getTpetra_Map(), A, LDA, NumVectors));
    } 

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, size_t NumVectors) { 
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const TpetraMap, map, tMap, "Cthulhu::TpetraMultiVector constructors only accept Cthulhu::TpetraMap as input arguments.");
      vec_ = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tMap->getTpetra_Map(), ArrayOfPtrs, NumVectors));
    } 
  
    TpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &vec) : vec_(vec) { CTHULHU_DEBUG_ME; } //TODO removed const

    //! TpetraMultiVector destructor.
    virtual ~TpetraMultiVector() { CTHULHU_DEBUG_ME; }

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->replaceGlobalValue(globalRow, vectorIndex, value); }

    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->sumIntoGlobalValue(globalRow, vectorIndex, value); }

    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->replaceLocalValue(myRow, vectorIndex, value); }

    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { CTHULHU_DEBUG_ME; vec_->sumIntoLocalValue(myRow, vectorIndex, value); }

    //! Initialize all values in a multi-vector with specified value.
    inline void putScalar(const Scalar &value) { CTHULHU_DEBUG_ME; vec_->putScalar(value); }

    //! Set multi-vector values to random numbers.
    inline void randomize() { CTHULHU_DEBUG_ME; vec_->randomize(); }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Replace the underlying Map with a compatible one.
    inline void replaceMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map) { CTHULHU_DEBUG_ME; vec_->replaceMap(map); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Instruct a local (non-distributed) TpetraMultiVector to sum values across all nodes.
    inline void reduce() { CTHULHU_DEBUG_ME; vec_->reduce(); }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! = Operator.
    /*! \param In A - Multivector to copy
     */
    inline TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& operator=(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) { CTHULHU_DEBUG_ME; return vec_->(source); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //@}

    //! @name Data Copy and View get methods
    /** These methods are used to get the data underlying the TpetraMultiVector. They return data in one of three forms: 
        - a TpetraMultiVector with a subset of the columns of the target TpetraMultiVector
        - a raw C pointer or array of raw C pointers
        - one of the Teuchos memory management classes
        Not all of these methods are valid for a particular TpetraMultiVector. For instance, calling a method that accesses a 
        view of the data in a 1-D format (i.e., get1dView) requires that the target TpetraMultiVector has constant stride.
    */
    //@{
#ifdef CTHULHU_NOT_IMPLEMENTED

    //! Returns a MultiVector with copies of selected columns.
    inline Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::Range1D &colRng) const { CTHULHU_DEBUG_ME; return vec_->subCopy(colRng); }

    //! Returns a TpetraMultiVector with copies of selected columns.
    inline Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::ArrayView<const size_t> &cols) const { CTHULHU_DEBUG_ME; return vec_->subCopy(cols); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::Range1D &colRng) const { CTHULHU_DEBUG_ME; return vec_->subView(colRng); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::ArrayView<const size_t> &cols) const { CTHULHU_DEBUG_ME; return vec_->subView(cols); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::Range1D &colRng) { CTHULHU_DEBUG_ME; return vec_->subViewNonConst(colRng); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::ArrayView<const size_t> &cols) { CTHULHU_DEBUG_ME; return vec_->subViewNonConst(cols); }

    //! \brief Returns a const MultiVector view of a subset of rows.
    /** 
        Returns a const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetView(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset) const { CTHULHU_DEBUG_ME; return vec_->offsetView(subMap, offset); }

    //! \brief Returns a non-const MultiVector view of a subset of rows.
    /** 
        Returns a non-const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetViewNonConst(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset) { CTHULHU_DEBUG_ME; return vec_->offsetViewNonConst(subMap, offset); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Const Vector access function.
    inline Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVector(size_t j) const { 
      CTHULHU_DEBUG_ME; 
      const RCP<Tpetra::Vector<Scalar,LocalOrdinal, GlobalOrdinal, Node> > v = vec_->getVectorNonConst(j);

      rcp(new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(v)); 
      return Teuchos::null;//rcp(new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(v)); 
    }

    //! Vector access function. //TODO see getVector
    inline Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVectorNonConst(size_t j) { CTHULHU_DEBUG_ME; return vec_->getVectorNonConst(j); }

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<const Scalar> getData(size_t j) const { CTHULHU_DEBUG_ME; return vec_->getData(j); }

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j) { CTHULHU_DEBUG_ME; return vec_->getDataNonConst(j); }

    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    inline void get1dCopy(Teuchos::ArrayView<Scalar> A, size_t LDA) const { CTHULHU_DEBUG_ME; vec_->get1dCopy(A, LDA); }

    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    inline void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const { CTHULHU_DEBUG_ME; vec_->get2dCopy(ArrayOfPtrs); }

    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
    inline Teuchos::ArrayRCP<const Scalar> get1dView() const { CTHULHU_DEBUG_ME; return vec_->get1dView(); }

    //! Return const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView() const { CTHULHU_DEBUG_ME; return vec_->get2dView(); }

    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<Scalar> get1dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->(); }
    inline Teuchos::ArrayRCP<Scalar> get1dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->get1dViewNonConst(); }

    //! Return non-const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->get2dViewNonConst(); }

    //! Return a const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline const Kokkos::MultiVector<Scalar,Node> & getLocalMV() const { CTHULHU_DEBUG_ME; return vec_->getLocalMV(); }

    //! Return a non-const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline Kokkos::MultiVector<Scalar,Node> & getLocalMVNonConst() { CTHULHU_DEBUG_ME; return vec_->getLocalMVNonConst(); }

    //@}

    //! @name Mathematical methods
    //@{ 

    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    inline void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Teuchos::ArrayView<Scalar> &dots) const { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, A, tA, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      vec_->dot(*tA.getTpetra_MultiVector(), dots); 
    }

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    inline void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, A, tA, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      vec_->abs(*tA.getTpetra_MultiVector()); 
    }

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    inline void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, A, tA, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      vec_->reciprocal(*tA.getTpetra_MultiVector()); 
    }

    //! Scale the current values of a multi-vector, this = alpha*this.
    inline void scale(const Scalar &alpha) { 
      CTHULHU_DEBUG_ME;
      vec_->scale(alpha); 
    }

    //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
    inline void scale(Teuchos::ArrayView<const Scalar> alpha) { 
      CTHULHU_DEBUG_ME; 
      vec_->scale(alpha); 
    }

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    inline void scale(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, A, tA, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      vec_->scale(alpha, *tA.getTpetra_MultiVector()); 
    }

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    inline void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta) { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, A, tA, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      vec_->update(alpha, *tA.getTpetra_MultiVector(), beta); 
    }

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    inline void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &gamma) { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, A, tA, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, B, tB, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      vec_->update(alpha, *tA.getTpetra_MultiVector(), beta, *tB.getTpetra_MultiVector(), gamma); 
    }

    //! Compute 1-norm of each vector in multi-vector.
    inline void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { 
      CTHULHU_DEBUG_ME; 
      vec_->norm1(norms); 
    }

    //! Compute 2-norm of each vector in multi-vector.
    inline void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { 
      CTHULHU_DEBUG_ME; 
      vec_->norm2(norms); 
    }

    //! Compute Inf-norm of each vector in multi-vector.
    inline void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { 
      CTHULHU_DEBUG_ME; 
      vec_->normInf(norms); 
    }

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    inline void normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { 
      CTHULHU_DEBUG_ME; 
      CTHULHU_DYNAMIC_CAST(const TpetraMultiVector, weights, tWeights, "This Cthulhu::TpetraMultiVector method only accept Cthulhu::TpetraMultiVector as input arguments.");
      vec_->normWeighted(*tWeights.getTpetra_MultiVector(), norms);
    }

    //! Compute mean (average) value of each vector in multi-vector.
    inline void meanValue(const Teuchos::ArrayView<Scalar> &means) const { CTHULHU_DEBUG_ME; vec_->meanValue(means); }

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    inline void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta) { CTHULHU_DEBUG_ME; vec_->multiply(transA, transB, alpha, A, B, beta); }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Element-wise multiply of a Vector A with a TpetraMultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
    inline void elementWiseMultiply(Scalar scalarAB, const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, Scalar scalarThis) { CTHULHU_DEBUG_ME; vec_->elementWiseMultiply(scalarAB, A , B , scalarThis); }
#endif // CTHULHU_NOT_IMPLEMENTED
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    inline size_t getNumVectors() const { CTHULHU_DEBUG_ME; return vec_->getNumVectors(); }

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    inline size_t getLocalLength() const { CTHULHU_DEBUG_ME; return vec_->getLocalLength(); }

    //! Returns the global vector length of vectors in the multi-vector.
    inline global_size_t getGlobalLength() const { CTHULHU_DEBUG_ME; return vec_->getGlobalLength(); }

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    inline size_t getStride() const { CTHULHU_DEBUG_ME; return vec_->getStride(); }

    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    inline bool isConstantStride() const { CTHULHU_DEBUG_ME; return vec_->isConstantStride(); }

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { CTHULHU_DEBUG_ME; return vec_->description(); }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { CTHULHU_DEBUG_ME; vec_->describe(out, verbLevel); }

    //@}

    RCP< Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_MultiVector() const { CTHULHU_DEBUG_ME; return vec_; }
    
  private:
    RCP< Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > vec_;


  }; // class TpetraMultiVector

#ifdef CTHULHU_NOT_IMPLEMENTED
  /** \brief Non-member function to create a TpetraMultiVector from a specified Map.
      \relates TpetraMultiVector
  */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createTpetraMultiVector(const Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t numVectors) { CTHULHU_DEBUG_ME; }
#endif // CTHULHU_NOT_IMPLEMENTED
  
} // namespace Cthulhu

#endif // CTHULHU_MULTIVECTOR_DECL_HPP
