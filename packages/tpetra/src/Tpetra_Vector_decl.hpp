// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_VECTOR_DECL_HPP
#define TPETRA_VECTOR_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_MultiVector_decl.hpp"

namespace Tpetra {

/// \class Vector
/// \brief A distributed dense vector.
///
/// This class inherits from MultiVector, and has the same template
/// parameters.  A Vector is a special case of a MultiVector that has
/// only one vector (column).  It may be used wherever a MultiVector
/// may be used.
template<class Scalar,
         class LocalOrdinal=int,
         class GlobalOrdinal=LocalOrdinal,
         class Node=KokkosClassic::DefaultNode::DefaultNodeType>
class Vector : public MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
  // need this so that MultiVector::operator() can call Vector's private view constructor
  friend class MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;


public:
  typedef typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot_type dot_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  explicit Vector(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, bool zeroOut=true);

  //! Vector copy constructor.
  Vector(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);

  //! \brief Set vector values from an existing array (copy)
  Vector(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const ArrayView<const Scalar> &A);

  //! Destructor.
  virtual ~Vector();

  //@}
  //! @name Post-construction modification routines
  //@{

  //!Create a cloned Vector for a different node type
  template <class Node2>
  RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
  clone(const RCP<Node2> &node2);


  //! Replace current value at the specified location with specified value.
  /** \pre \c globalRow must be a valid global element on this node, according to the row map.
   */
  void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

  //! Adds specified value to existing value at the specified location.
  /** \pre \c globalRow must be a valid global element on this node, according to the row map.
   */
  void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value);

  //! Replace current value at the specified location with specified values.
  /** \pre \c localRow must be a valid local element on this node, according to the row map.
   */
  void replaceLocalValue(LocalOrdinal myRow, const Scalar &value);

  //! Adds specified value to existing value at the specified location.
  /** \pre \c localRow must be a valid local element on this node, according to the row map.
   */
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value);

  //@}

  //! @name Extraction methods
  //@{

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dCopy; // overloading, not hiding
  //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
  void get1dCopy(ArrayView<Scalar> A) const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDataNonConst; // overloading, not hiding
  //! View of the local values of this vector.
  Teuchos::ArrayRCP<Scalar> getDataNonConst()     { return getDataNonConst(0); }

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getData; // overloading, not hiding
  //! Const view of the local values of this vector.
  Teuchos::ArrayRCP<const Scalar> getData() const { return getData(0); }

  //@}

  //! @name Mathematical methods
  //@{

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot; // overloading, not hiding
  //! Computes dot product of this Vector against input Vector x.
  Scalar dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1; // overloading, not hiding
  //! Return 1-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2; // overloading, not hiding
  //! Compute 2-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf; // overloading, not hiding
  //! Compute Inf-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normWeighted; // overloading, not hiding
  //! Compute Weighted 2-norm (RMS Norm) of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType normWeighted(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights) const;

  using MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue; // overloading, not hiding
  //! Compute mean (average) value of this Vector.
  Scalar meanValue() const;

  //@}
  /// \name Methods to get a row subset view of a Vector (as a Vector).
  ///
  /// These methods hide those of the same names in MultiVector.  This
  /// is necessary because users who use Vectors would very much like
  /// a view of a Vector to be a Vector, rather than a MultiVector.
  /// Not only do these methods in MultiVector return MultiVector
  /// pointers, those pointers can't even be dynamic_cast to Vector,
  /// because the returned objects are not Vectors.
  ///
  /// There are two options to make this work:
  /// <ol>
  /// <li> Make these methods virtual, and override them in Vector. </li>
  /// <li> Hide these methods in MultiVector (that return MultiVector)
  ///      with methods in Vector that return Vector. </li>
  /// </ol>
  /// The first option only solves the dynamic_cast problem mentioned
  /// above.  The methods still would have to return MultiVector
  /// pointers.  Thus, we prefer the second option.  Even though the
  /// methods in Vector must therefore hide those in MultiVector,
  /// Vector is a MultiVector, so the returned objects from the Vector
  /// methods are still compatible with MultiVector operations.
  //@{

  /// \brief Return a const view of a subset of rows.
  ///
  /// Return a const (nonmodifiable) view of this Vector consisting of
  /// a subset of the rows, as specified by an offset and a subset Map
  /// of this Vector's current row Map.  If you want X1 or X2 to be
  /// nonconst (modifiable) views, use offsetViewNonConst() with the
  /// same arguments.  "View" means "alias": if the original (this)
  /// MultiVector's data change, the view will see the changed data.
  ///
  /// \param subMap [in] The row Map for the new Vector.
  ///   This must be a subset Map of this MultiVector's row Map.
  /// \param offset [in] The local row offset at which to start the view.
  ///
  /// Suppose that you have a Vector X, and you want to view X, on all
  /// processes in X's (MPI) communicator, as split into two row
  /// blocks X1 and X2.  One could express this in Matlab notation as
  /// X = [X1; X2], except that here, X1 and X2 are views into X,
  /// rather than copies of X's data.  This method assumes that the
  /// <i>local</i> indices of X1 and X2 are each contiguous, and that
  /// the local indices of X2 follow those of X1.  If that is not the
  /// case, you cannot use views to divide X into blocks like this;
  /// you must instead use the Import or Export functionality, which
  /// copies the relevant rows of X.
  ///
  /// Here is how you would construct the views X1 and X2.
  /// \code
  /// using Teuchos::RCP;
  /// typedef Tpetra::Map<LO, GO, Node> map_type;
  /// typedef Tpetra::Vector<Scalar, LO, GO, Node> V;
  ///
  /// V X (...); // the input Vector
  /// // ... fill X with data ...
  ///
  /// // Map that on each process in X's communicator,
  /// // contains the global indices of the rows of X1.
  /// RCP<const map_type> map1 (new map_type (...));
  /// // Map that on each process in X's communicator,
  /// // contains the global indices of the rows of X2.
  /// RCP<const map_type> map2 (new map_type (...));
  ///
  /// // Create the first view X1.  The second argument, the offset,
  /// // is the index of the local row at which to start the view.
  /// // X1 is the topmost block of X, so the offset is zero.
  /// RCP<const V> X1 = X.offsetView (map1, 0);
  ///
  /// // Create the second view X2.  X2 is directly below X1 in X,
  /// // so the offset is the local number of rows in X1.  This is
  /// // the same as the local number of entries in map1.
  /// RCP<const V> X2 = X.offsetView (map2, X1->getLocalLength ());
  /// \endcode
  ///
  /// It is legal, in the above example, for X1 or X2 to have zero
  /// local rows on any or all process(es).  In that case, the
  /// corresponding Map must have zero local entries on that / those
  /// process(es).  In particular, if X2 has zero local rows on a
  /// process, then the corresponding offset on that process would be
  /// the number of local rows in X (and therefore in X1) on that
  /// process.  This is the only case in which the sum of the local
  /// number of entries in \c subMap (in this case, zero) and the \c
  /// offset may equal the number of local entries in <tt>*this</tt>.
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
              size_t offset) const;

  /// \brief Return a nonconst view of a subset of rows.
  ///
  /// Return a nonconst (modifiable) view of this Vector consisting of
  /// a subset of the rows, as specified by an offset and a subset Map
  /// of this Vector's current row Map.  If you want X1 or X2 to be
  /// const (nonmodifiable) views, use offsetView() with the same
  /// arguments.  "View" means "alias": if the original (this)
  /// Vector's data change, the view will see the changed data, and if
  /// the view's data change, the original Vector will see the changed
  /// data.
  ///
  /// \param subMap [in] The row Map for the new Vector.  This must be
  ///   a subset Map of this Vector's row Map.
  /// \param offset [in] The local row offset at which to start the view.
  ///
  /// See the documentation of offsetView() for a code example and
  /// an explanation of edge cases.
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  offsetViewNonConst (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap,
                      size_t offset);

  //@}
  //! @name Implementation of the Teuchos::Describable interface
  //@{

  //! A simple one-line description of this object.
  virtual std::string description() const;

  //! Print the object with some verbosity level to a FancyOStream.
  virtual void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}

protected:

  template <class S,class LO,class GO,class N>
  friend RCP< Vector<S,LO,GO,N> >
  createVectorFromView(const RCP<const Map<LO,GO,N> > &,const ArrayRCP<S> &);

  // view constructor, sitting on user allocated data, only for CPU nodes
  // and his non-member constructor friend
  Vector (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
          const ArrayRCP<Scalar> &view,
          EPrivateHostViewConstructor /* dummy tag */);

  //! Advanced constructor accepting parallel buffer view, used by MultiVector to break off Vector objects
  Vector (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
          const ArrayRCP<Scalar> & data);

  //! Advanced constructor accepting parallel buffer view, used by MultiVector to break off Vector objects
  Vector (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
          const ArrayRCP<Scalar> & data,
          EPrivateComputeViewConstructor /* dummy tag */);

  typedef KokkosClassic::MultiVector<Scalar,Node> KMV;
  typedef KokkosClassic::DefaultArithmetic<KMV>   MVT;
}; // class Vector

/** \brief Non-member function to create a Vector from a specified Map.
    \relatesalso Vector
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
createVector (const RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> > &map)
{
  const bool DO_INIT_TO_ZERO = true;
  return rcp (new Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> (map, DO_INIT_TO_ZERO));
}

//! \brief Non-member function to create a Vector with view semantics using user-allocated data.
/*! This use case is not supported for all nodes. Specifically, it is not typically supported for accelerator-based nodes like KokkosClassic::ThrustGPUNode.
    \relatesalso Vector
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
createVectorFromView (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
                      const ArrayRCP<Scalar> &view)
{
  return rcp(
    // this is a protected constructor, but we are friends
    new Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(
      map,
      // this will fail to compile for unsupported node types
      Tpetra::details::ViewAccepter<Node>::template acceptView<Scalar>(view),
      HOST_VIEW_CONSTRUCTOR)
  );
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class Node2>
RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::clone(const RCP<Node2> &node2){
        typedef Map<LocalOrdinal,GlobalOrdinal,Node2> Map2;
        typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node2> V2;
        Teuchos::ArrayRCP<const Scalar> V_view = this->getData();
        Teuchos::RCP<const Map2> clonedMap = this->getMap()->template clone<Node2> (node2);
        Teuchos::RCP<V2> clonedV = Teuchos::rcp(new V2(clonedMap));
        Teuchos::ArrayRCP<Scalar> clonedV_view = clonedV->getDataNonConst();
        clonedV_view.deepCopy(V_view());
        clonedV_view = Teuchos::null;
        return clonedV;
  }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node >
  createCopy( const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& src);
} // namespace Tpetra

// Include KokkosRefactor partial specialisation if enabled
#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Tpetra_KokkosRefactor_Vector_decl.hpp"
#endif

#endif // TPETRA_VECTOR_DECL_HPP
