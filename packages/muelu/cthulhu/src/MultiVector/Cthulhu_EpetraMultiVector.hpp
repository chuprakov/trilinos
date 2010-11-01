#ifndef CTHULHU_EPETRAMULTIVECTOR_DECL_HPP
#define CTHULHU_EPETRAMULTIVECTOR_DECL_HPP

#include "Cthulhu_Classes.hpp" //TMP

#include "Cthulhu_MultiVector.hpp"

#include "Cthulhu_EpetraMap.hpp" 
#include "Epetra_MultiVector.h"

#include "Cthulhu_Debug.hpp"

namespace Cthulhu {

  // #ifndef DOXYGEN_SHOULD_SKIP_THIS
  //   // forward declaration of Vector, needed to prevent circular inclusions
  //   template<class S, class LO, class GO, class N> class Vector;
  // #endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
    This class is templated on \c double, \c int and \c GlobalOrdinal. 
    The \c int type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c int type.
  */
  class EpetraMultiVector : public Cthulhu::MultiVector<double,int,int> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 
    //! Basic EpetraMultiVector constuctor.
    EpetraMultiVector(const Teuchos::RCP<const Map<int,int> > &map, size_t NumVectors, bool zeroOut=true) {
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, map, eMap, "Cthulhu::TpetraMultiVector constructors only accept Cthulhu::TpetraMap as input arguments.");
      vec_ = rcp(new Epetra_MultiVector(eMap->getEpetra_Map(), NumVectors, zeroOut));
    }

#ifdef CTHULHU_NOT_IMPLEMENTED    
    //! EpetraMultiVector copy constructor.
    EpetraMultiVector(const EpetraMultiVector<double,int,int,Node> &source){ CTHULHU_DEBUG_ME; } 
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    EpetraMultiVector(const Teuchos::RCP<const Map<int,int,Node> > &map, const Teuchos::ArrayView<const double> &A, size_t LDA, size_t NumVectors) { CTHULHU_DEBUG_ME; } 
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    EpetraMultiVector(const Teuchos::RCP<const Map<int,int,Node> > &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const double> > &ArrayOfPtrs, size_t NumVectors) { CTHULHU_DEBUG_ME; } 
#endif // CTHULHU_NOT_IMPLEMENTED
  
    EpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector> &vec) : vec_(vec) { CTHULHU_DEBUG_ME; }

    //! EpetraMultiVector destructor.
    virtual ~EpetraMultiVector() { CTHULHU_DEBUG_ME; }

    //@}

    //! @name Post-construction modification routines
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(int globalRow, size_t vectorIndex, const double &value) { CTHULHU_DEBUG_ME; vec_->replaceGlobalValue(globalRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(int globalRow, size_t vectorIndex, const double &value) { CTHULHU_DEBUG_ME; vec_->sumIntoGlobalValue(globalRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(int myRow, size_t vectorIndex, const double &value) { CTHULHU_DEBUG_ME; vec_->replaceLocalValue(myRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(int myRow, size_t vectorIndex, const double &value) { CTHULHU_DEBUG_ME; vec_->sumIntoLocalValue(myRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Initialize all values in a multi-vector with specified value.
    inline void putdouble(const double &value) { CTHULHU_DEBUG_ME; vec_->putdouble(value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Set multi-vector values to random numbers.
    inline void randomize() { CTHULHU_DEBUG_ME; vec_->Random(); }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Replace the underlying Map with a compatible one.
    inline void replaceMap(const Teuchos::RCP<const Map<int,int,Node> > &map) { CTHULHU_DEBUG_ME; vec_->replaceMap(map); }
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Instruct a local (non-distributed) EpetraMultiVector to sum values across all nodes.
    inline void reduce() { CTHULHU_DEBUG_ME; vec_->reduce(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! = Operator.
    /*! \param In A - Multivector to copy
     */
    inline EpetraMultiVector<double,int,int,Node>& operator=(const EpetraMultiVector<double,int,int,Node> &source) { CTHULHU_DEBUG_ME; return vec_->(source); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //@}

    //! @name Data Copy and View get methods
    /** These methods are used to get the data underlying the EpetraMultiVector. They return data in one of three forms: 
        - a EpetraMultiVector with a subset of the columns of the target EpetraMultiVector
        - a raw C pointer or array of raw C pointers
        - one of the Teuchos memory management classes
        Not all of these methods are valid for a particular EpetraMultiVector. For instance, calling a method that accesses a 
        view of the data in a 1-D format (i.e., get1dView) requires that the target EpetraMultiVector has constant stride.
    */
    //@{
#ifdef CTHULHU_NOT_IMPLEMENTED

    //! Returns a MultiVector with copies of selected columns.
    inline Teuchos::RCP<MultiVector<double,int,int,Node> > subCopy(const Teuchos::Range1D &colRng) const { CTHULHU_DEBUG_ME; return vec_->subCopy(colRng); }

    //! Returns a EpetraMultiVector with copies of selected columns.
    inline Teuchos::RCP<MultiVector<double,int,int,Node> > subCopy(const Teuchos::ArrayView<const size_t> &cols) const { CTHULHU_DEBUG_ME; return vec_->subCopy(cols); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<double,int,int,Node> > subView(const Teuchos::Range1D &colRng) const { CTHULHU_DEBUG_ME; return vec_->subView(colRng); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<double,int,int,Node> > subView(const Teuchos::ArrayView<const size_t> &cols) const { CTHULHU_DEBUG_ME; return vec_->subView(cols); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<double,int,int,Node> > subViewNonConst(const Teuchos::Range1D &colRng) { CTHULHU_DEBUG_ME; return vec_->subViewNonConst(colRng); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<double,int,int,Node> > subViewNonConst(const Teuchos::ArrayView<const size_t> &cols) { CTHULHU_DEBUG_ME; return vec_->subViewNonConst(cols); }

    //! \brief Returns a const MultiVector view of a subset of rows.
    /** 
        Returns a const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<const MultiVector<double,int,int,Node> > offsetView(const Teuchos::RCP<const Map<int,int,Node> > &subMap, size_t offset) const { CTHULHU_DEBUG_ME; return vec_->offsetView(subMap, offset); }

    //! \brief Returns a non-const MultiVector view of a subset of rows.
    /** 
        Returns a non-const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<MultiVector<double,int,int,Node> > offsetViewNonConst(const Teuchos::RCP<const Map<int,int,Node> > &subMap, size_t offset) { CTHULHU_DEBUG_ME; return vec_->offsetViewNonConst(subMap, offset); }

    //! Const Vector access function.
    inline Teuchos::RCP<const Vector<double,int,int,Node> > getVector(size_t j) const { CTHULHU_DEBUG_ME; return vec_->getVector(j); }

    //! Vector access function.
    inline Teuchos::RCP<Vector<double,int,int,Node> > getVectorNonConst(size_t j) { CTHULHU_DEBUG_ME; return vec_->getVectorNonConst(j); }
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<const double> getData(size_t j) const { CTHULHU_DEBUG_ME; return vec_->getData(j); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<double> getDataNonConst(size_t j) { 
      CTHULHU_DEBUG_ME; 

      double ** arrayOfPointers;
      
      vec_->ExtractView(&arrayOfPointers);
     
      double * data = arrayOfPointers[j];
      int localLength = vec_->MyLength();
      
      return ArrayRCP<double>(data, 0, localLength, false); // not ownership
    }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    inline void get1dCopy(Teuchos::ArrayView<double> A, size_t LDA) const { CTHULHU_DEBUG_ME; vec_->get1dCopy(A, LDA); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    inline void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<double> > ArrayOfPtrs) const { CTHULHU_DEBUG_ME; vec_->get2dCopy(ArrayOfPtrs); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
    inline Teuchos::ArrayRCP<const double> get1dView() const { CTHULHU_DEBUG_ME; return vec_->get1dView(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > get2dView() const { CTHULHU_DEBUG_ME; return vec_->get2dView(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<double> get1dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->(); }
    inline Teuchos::ArrayRCP<double> get1dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->get1dViewNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return non-const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > get2dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->get2dViewNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return a const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline const Kokkos::MultiVector<double,Node> & getLocalMV() const { CTHULHU_DEBUG_ME; return vec_->getLocalMV(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return a non-const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline Kokkos::MultiVector<double,Node> & getLocalMVNonConst() { CTHULHU_DEBUG_ME; return vec_->getLocalMVNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@}

    //! @name Mathematical methods
    //@{ 
#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    inline void dot(const EpetraMultiVector<double,int,int,Node> &A, const Teuchos::ArrayView<double> &dots) const { CTHULHU_DEBUG_ME; vec_->dot(A, dots); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    inline void abs(const EpetraMultiVector<double,int,int,Node> &A) { CTHULHU_DEBUG_ME; vec_->abs(A); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    inline void reciprocal(const EpetraMultiVector<double,int,int,Node> &A) { CTHULHU_DEBUG_ME; vec_->reciprocal(A); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Scale the current values of a multi-vector, this = alpha*this.
    inline void scale(const double &alpha) { CTHULHU_DEBUG_ME; vec_->scale(alpha); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
    inline void scale(Teuchos::ArrayView<const double> alpha) { CTHULHU_DEBUG_ME; vec_->scale(alpha); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    inline void scale(const double &alpha, const EpetraMultiVector<double,int,int,Node> &A) { CTHULHU_DEBUG_ME; vec_->scale(A); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    inline void update(const double &alpha, const EpetraMultiVector<double,int,int,Node> &A, const double &beta) { CTHULHU_DEBUG_ME; vec_->update(A, beta); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    inline void update(const double &alpha, const EpetraMultiVector<double,int,int,Node> &A, const double &beta, const EpetraMultiVector<double,int,int,Node> &B, const double &gamma) { CTHULHU_DEBUG_ME; vec_->update(A, beta); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Compute 1-norm of each vector in multi-vector.
    inline void norm1(const Teuchos::ArrayView<typename Teuchos::doubleTraits<double>::magnitudeType> &norms) const { CTHULHU_DEBUG_ME; vec_->norm1(norms); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Compute 2-norm of each vector in multi-vector.
    inline void norm2(const Teuchos::ArrayView<typename Teuchos::doubleTraits<double>::magnitudeType> &norms) const { CTHULHU_DEBUG_ME; vec_->norm2(norms); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Compute Inf-norm of each vector in multi-vector.
    inline void normInf(const Teuchos::ArrayView<typename Teuchos::doubleTraits<double>::magnitudeType> &norms) const { CTHULHU_DEBUG_ME; vec_->normInf(norms); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    inline void normWeighted(const EpetraMultiVector<double,int,int,Node> &weights, const Teuchos::ArrayView<typename Teuchos::doubleTraits<double>::magnitudeType> &norms) const { CTHULHU_DEBUG_ME; vec_->normWeighted(weights, norms); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Compute mean (average) value of each vector in multi-vector.
    inline void meanValue(const Teuchos::ArrayView<double> &means) const { CTHULHU_DEBUG_ME; vec_->meanValue(means); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    inline void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const double &alpha, const EpetraMultiVector<double,int,int,Node> &A, const EpetraMultiVector<double,int,int,Node> &B, const double &beta) { CTHULHU_DEBUG_ME; vec_->multiply(transA, transB, alpha, A, B, beta); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
    inline void elementWiseMultiply(double scalarAB, const Vector<double,int,int,Node> &A, const EpetraMultiVector<double,int,int,Node> &B, double scalarThis) { CTHULHU_DEBUG_ME; vec_->elementWiseMultiply(scalarAB, A , B , scalarThis); }
#endif // CTHULHU_NOT_IMPLEMENTED
    //@} 

    //! @name Attribute access functions
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the number of vectors in the multi-vector.
    inline size_t getNumVectors() const { CTHULHU_DEBUG_ME; return vec_->getNumVectors(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    inline size_t getLocalLength() const { CTHULHU_DEBUG_ME; return vec_->getLocalLength(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the global vector length of vectors in the multi-vector.
    inline global_size_t getGlobalLength() const { CTHULHU_DEBUG_ME; return vec_->getGlobalLength(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    inline size_t getStride() const { CTHULHU_DEBUG_ME; return vec_->getStride(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    inline bool isConstantStride() const { CTHULHU_DEBUG_ME; return vec_->isConstantStride(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { CTHULHU_DEBUG_ME; return vec_->description(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { CTHULHU_DEBUG_ME; vec_->describe(out, verbLevel); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@}

    RCP< Epetra_MultiVector > getEpetra_MultiVector() const { CTHULHU_DEBUG_ME; return vec_; }
    
  private:
    RCP< Epetra_MultiVector > vec_;


  }; // class EpetraMultiVector

#ifdef CTHULHU_NOT_IMPLEMENTED
  /** \brief Non-member function to create a EpetraMultiVector from a specified Map.
      \relates EpetraMultiVector
  */
  template <class double, class int, class int, class Node>
  Teuchos::RCP< EpetraMultiVector<double,int,int,Node> >
  createEpetraMultiVector(const Teuchos::RCP< const Map<int,int,Node> > &map, size_t numVectors) { CTHULHU_DEBUG_ME; }
#endif // CTHULHU_NOT_IMPLEMENTED
  
} // namespace Cthulhu

#endif // CTHULHU_MULTIVECTOR_DECL_HPP
