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

#ifndef TPETRA_EXPERIMENTAL_BLOCKMULTIVECTOR_DECL_HPP
#define TPETRA_EXPERIMENTAL_BLOCKMULTIVECTOR_DECL_HPP

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Experimental_BlockView.hpp>

namespace Tpetra {
namespace Experimental {

/// \class BlockMultiVector
/// \brief MultiVector for multiple degrees of freedom per mesh point
/// \author Mark Hoemmen
///
/// \tparam Scalar The type of each entry of the vector
/// \tparam LO The local ordinal type
/// \tparam GO The global ordinal type
/// \tparam Node The Kokkos Node type
///
/// BlockMultiVector is like ::Tpetra::MultiVector, but its interface
/// supports multiple degrees of freedom per mesh point.  You can
/// specify a mesh point by its local or global index, and read or
/// write the values at that point.  Every mesh point must have the
/// same number of degrees of freedom.  We call the number of degrees
/// of freedom per mesh point the <i>block size</i>.
///
/// BlockMultiVector has <i>view semantics</i>.  This means that the
/// copy constructor and assignment operator (operator=) all do
/// shallow copies.  If you want to create a deep copy, call
/// createCopy(); if you want to copy deeply between two
/// BlockMultiVector objects, call deep_copy().
///
/// The replace*, sumInto*, and get* methods all are meant to be
/// thread-safe.  They don't throw exceptions on error; instead, they
/// return an error code.  They are marked "const" because Kokkos
/// requires this in order to use them in parallel kernels.  "Const"
/// is legitimate here, even though some of these methods may modify
/// entries of the vector.  This is because they do not change any
/// constants or pointers (e.g., they do not reallocate memory).
///
/// BlockMultiVector stores entries in a column corresponding to
/// degrees of freedom for the same mesh point contiguously.  This is
/// not strictly necessary; we could generalize so that they are
/// stored in a strided fashion, or even take away the layout
/// assumption and only allow access to entries by copy (e.g.,
/// getLocalRowCopy() instead of getLocalRowView()).
///
/// Here is how you are supposed to use this class:
///
///   1. Fill the BlockMultiVector.
///   2. Do global assembly with BlockMultiVector, if necessary.  (It
///      implements Tpetra::DistObject, so you can use it with Import
///      or Export.)
///   3. Call getMultiVectorView() to get a Tpetra::MultiVector which
///      views the BlockMultiVector's data.
///   4. Give the Tpetra::MultiVector to Trilinos' solvers.
///
/// Note that BlockMultiVector is not a subclass of
/// ::Tpetra::MultiVector.  I did this on purpose, for simplicity.
/// You should really think of BlockMultiVector as an in-place
/// reinterpretation of a ::Tpetra::MultiVector.
///
/// This is the first Tpetra class to assume view semantics of all
/// Tpetra objects.  In particular, it assumes that Map has have view
/// semantics, so that its copy constructor and assignment operator do
/// shallow copies, and its default constructor (that takes no
/// arguments) compiles and succeeds, creating a meaningful empty
/// object.
///
/// Design philosophy of the new block objects:
///
///   1. Each mesh point has the same number of degrees of freedom
///   2. A BlockMultiVector might <i>have</i> a MultiVector, but <i>is
///      not</i> a MultiVector (i.e., is not a subclass of MultiVector)
///
/// Point 1 lets us omit the whole untested and probably broken
/// infrastructure of Tpetra::BlockMap.  Users fill by mesh points,
/// not degrees of freedom.  Thus, they mainly care about the
/// distribution of mesh points, not so much about what GID we assign
/// to which degree of freedom.  The latter is just another Map from
/// their perspective, if they care at all.  Preconditioners that
/// respect the block structure also just want the mesh Map.
/// Iterative linear solvers treat the matrix and preconditioner as
/// black-box operators.  This means that they only care about the
/// domain and range point Maps.
///
/// The latter motivates Point 2.  BlockMultiVector views a
/// MultiVector, and iterative solvers (and users) can access the
/// MultiVector and work with it directly, along with its (point) Map.
/// It doesn't make sense for BlockMultiVector to implement
/// MultiVector, because the desired fill interfaces of the two
/// classes are different.
template<class Scalar,
         class LO = int,
         class GO = LO,
         class Node = KokkosClassic::DefaultNode::DefaultNodeType>
class BlockMultiVector :
#if 0 && defined(TPETRA_USE_KOKKOS_DISTOBJECT)
    public Tpetra::DistObjectKA<Scalar, LO, GO, Node>
#else
    public Tpetra::DistObject<Scalar, LO, GO, Node>
#endif // TPETRA_USE_KOKKOS_DISTOBJECT
{
private:
#if 0 && defined(TPETRA_USE_KOKKOS_DISTOBJECT)
  typedef Tpetra::DistObjectKA<Scalar, LO, GO, Node> dist_object_type;
#else
  typedef Tpetra::DistObject<Scalar, LO, GO, Node> dist_object_type;
#endif // TPETRA_USE_KOKKOS_DISTOBJECT

  typedef Teuchos::ScalarTraits<Scalar> STS;

protected:
  typedef Scalar packet_type;

public:
  //! \name Typedefs to facilitate template metaprogramming.
  //@{

  //! The type of entries in the vector.
  typedef Scalar scalar_type;
  //! The type of local indices.
  typedef LO local_ordinal_type;
  //! The type of global indices.
  typedef GO global_ordinal_type;
  //! The Kokkos Node type.
  typedef Node node_type;

  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type,
                      node_type> map_type;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> mv_type;

  typedef LittleVector<Scalar, LO> little_vec_type;
  typedef LittleVector<const Scalar, LO> const_little_vec_type;

  //@}
  //! \name Constructors
  //@{

  /// \brief Constructor that takes a mesh Map, a block size, and a
  ///   number of vectors (columns).
  ///
  /// \param meshMap [in] Map that describes the distribution of mesh
  ///   points (rather than the distribution of unknowns for those
  ///   mesh points).
  ///
  /// \param blockSize [in] The number of degrees of freedom per mesh
  ///   point.  We assume that this is the same for all mesh points in
  ///   the above Map.
  ///
  /// \param numVecs [in] Number of vectors (columns) in the MultiVector.
  ///
  /// The <i>mesh Map</i> describes the distribution of mesh points.
  /// Its corresponding <i>point Map</i> describes the distribution of
  /// degrees of freedom corresponding to those mesh points.  If you
  /// have already computed the point Map corresponding to the above
  /// mesh Map, then it is more efficient to call the four-argument
  /// constructor below, that takes both the mesh Map and the point
  /// Map.
  ///
  /// There are two ways to get the point Map corresponding to a given
  /// mesh Map and block size.  You may either call the class method
  /// makePointMap(), or you may call this two-argument constructor,
  /// and then call getPointMap().
  ///
  /// The point Map enables reinterpretation of a BlockMultiVector as
  /// a standard Tpetra::MultiVector.  This lets users solve linear
  /// systems with Trilinos' solvers and preconditioners, that expect
  /// vectors as Tpetra::MultiVector instances.
  BlockMultiVector (const map_type& meshMap,
                    const LO blockSize,
                    const LO numVecs);

  /// \brief Constructor that takes a mesh Map, a point Map, a block
  ///   size, and a number of vectors (columns).
  ///
  /// See the documentation of the three-argument constructor above.
  BlockMultiVector (const map_type& meshMap,
                    const map_type& pointMap,
                    const LO blockSize,
                    const LO numVecs);

  /// \brief View an existing MultiVector.
  ///
  /// \param X_mv [in/out] The MultiVector to view.  It MUST have view
  ///   semantics; otherwise this constructor throws.  Its Map must be
  ///   the same (in the sense of isSameAs) as the point Map
  ///   corresponding to the given mesh Map and block size.
  ///
  /// \param meshMap [in] The mesh Map to use for interpreting the
  ///   given MultiVector (in place) as a BlockMultiVector.
  ///
  /// \param blockSize [in] The number of degrees of freedom per mesh
  ///   point.  We assume that this is the same for all mesh points.
  BlockMultiVector (const mv_type& X_mv,
                    const map_type& meshMap,
                    const LO blockSize);

  /// \brief Default constructor.
  ///
  /// Creates an empty BlockMultiVector (with zero rows and columns,
  /// and block size zero).
  BlockMultiVector ();

  //@}
  //! \name Access to Maps, the block size, and a MultiVector view.
  //@{

  /// \brief Create and return the point Map corresponding to the
  ///   given mesh Map and block size.
  ///
  /// This is a class ("static") method so that you can make and reuse
  /// a point Map for creating different BlockMultiVector instances,
  /// using the more efficient four-argument constructor.
  static map_type
  makePointMap (const map_type& meshMap, const LO blockSize);

  /// \brief Get this BlockMultiVector's (previously computed) point Map.
  ///
  /// It is always valid to call this method.  A BlockMultiVector
  /// always has a point Map.  We do not compute the point Map lazily.
  map_type getPointMap () const {
    return pointMap_;
  }

  //! Get the number of degrees of freedom per mesh point.
  LO getBlockSize () const {
    return blockSize_;
  }

  //! Get the number of columns (vectors) in the BlockMultiVector.
  LO getNumVectors () const {
    return static_cast<LO> (mv_.getNumVectors ());
  }

  /// \brief Get a ::Tpetra::MultiVector that views this BlockMultiVector's data.
  ///
  /// This is how you can give a BlockMultiVector to Trilinos' solvers
  /// and preconditioners.
  mv_type getMultiVectorView ();

  //@}
  //! \name Coarse-grained operations
  //@{

  //! Fill all entries with the given value \c val.
  void putScalar (const Scalar& val);

  //! Multiply all entries in place by the given value \c val.
  void scale (const Scalar& val);

  //@}
  //! \name Fine-grained data access
  //@{

  /// \brief Replace all values at the given mesh point, using local
  ///   row and column indices.
  ///
  /// \param localRowIndex [in] Local index of the mesh point.
  /// \param colIndex [in] Column (vector) to modify.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given local index of the mesh point
  ///   is invalid on the calling process.
  ///
  /// \note This method, the other "replace" and "sumInto" methods,
  ///   and the view methods, are marked const.  This is because they
  ///   do not change pointers.  They do, of course, change the values
  ///   in the BlockMultiVector, but that does not require marking the
  ///   methods as nonconst.
  bool replaceLocalValues (const LO localRowIndex, const LO colIndex, const Scalar vals[]) const;

  /// \brief Replace all values at the given mesh point, using a global index.
  ///
  /// \param globalRowIndex [in] Global index of the mesh point.
  /// \param colIndex [in] Column (vector) to modify.
  /// \param vals [in] Input values with which to sum into whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given global index of the mesh point
  ///   is invalid on the calling process.
  bool replaceGlobalValues (const GO globalRowIndex, const LO colIndex, const Scalar vals[]) const;

  /// \brief Sum into all values at the given mesh point, using a local index.
  ///
  /// \param localRowIndex [in] Local index of the mesh point.
  /// \param colIndex [in] Column (vector) to modify.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given local index of the mesh point
  ///   is invalid on the calling process.
  bool sumIntoLocalValues (const LO localRowIndex, const LO colIndex, const Scalar vals[]) const;

  /// \brief Sum into all values at the given mesh point, using a global index.
  ///
  /// \param globalRowIndex [in] Global index of the mesh point.
  /// \param colIndex [in] Column (vector) to modify.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given global index of the mesh point
  ///   is invalid on the calling process.
  bool sumIntoGlobalValues (const GO globalRowIndex, const LO colIndex, const Scalar vals[]) const;

  /// \brief Get a writeable view of the entries at the given mesh
  ///   point, using a local index.
  ///
  /// \param localRowIndex [in] Local index of the mesh point.
  /// \param colIndex [in] Column (vector) to modify.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given local index of the mesh point
  ///   is invalid on the calling process.
  bool getLocalRowView (const LO localRowIndex, const LO colIndex, Scalar*& vals) const;

  /// \brief Get a writeable view of the entries at the given mesh
  ///   point, using a global index.
  ///
  /// \param globalRowIndex [in] Global index of the mesh point.
  /// \param colIndex [in] Column (vector) to modify.
  /// \param vals [in] Input values with which to replace whatever
  ///   existing values are at the mesh point.
  ///
  /// \return true if successful, else false.  This method will
  ///   <i>not</i> succeed if the given global index of the mesh point
  ///   is invalid on the calling process.
  bool getGlobalRowView (const GO globalRowIndex, const LO colIndex, Scalar*& vals) const;

  /// \brief Get a view of the degrees of freedom at the given mesh point.
  ///
  /// \warning This method's interface may change or disappear at any
  ///   time.  Please do not rely on it in your code yet.
  ///
  /// The preferred way to refer to little_vec_type is to get it from
  /// BlockMultiVector's typedef.  This is because different
  /// specializations of BlockMultiVector reserve the right to use
  /// different types to implement little_vec_type.  This gives us a
  /// porting strategy to move from "classic" Tpetra to the Kokkos
  /// refactor version.
  little_vec_type
  getLocalBlock (const LO localRowIndex, const LO colIndex) const;
  //@}

protected:
  /// \brief \name Implementation of DistObject (or DistObjectKA).
  ///
  /// The methods here implement Tpetra::DistObject or
  /// Tpetra::DistObjectKA, depending on a configure-time option.
  /// They let BlockMultiVector participate in Import and Export
  /// operations.  Users don't have to worry about these methods.
  //@{

  virtual bool checkSizes (const Tpetra::SrcDistObject& source);

#if 0 && defined(TPETRA_USE_KOKKOS_DISTOBJECT)
  virtual void
  copyAndPermute (const Tpetra::SrcDistObject& source,
                  size_t numSameIDs,
                  const Kokkos::View<const LO*, device_type> &permuteToLIDs,
                  const Kokkos::View<const LO*, device_type> &permuteFromLIDs);

  virtual void
  packAndPrepare (const Tpetra::SrcDistObject& source,
                  const Kokkos::View<const LO*, device_type>& exportLIDs,
                  Kokkos::View<packet_type*, device_type>& exports,
                  const Kokkos::View<size_t*, device_type>& numPacketsPerLID,
                  size_t& constantNumPackets,
                  Tpetra::Distributor& distor);

  virtual void
  unpackAndCombine (const Kokkos::View<const LO*, device_type>& importLIDs,
                    const Kokkos::View<const packet_type*, device_type>& imports,
                    const Kokkos::View<size_t*, device_type>& numPacketsPerLID,
                    size_t constantNumPackets,
                    Tpetra::Distributor& distor,
                    Tpetra::CombineMode CM);
#else
  virtual void
  copyAndPermute (const Tpetra::SrcDistObject& source,
                  size_t numSameIDs,
                  const Teuchos::ArrayView<const LO>& permuteToLIDs,
                  const Teuchos::ArrayView<const LO>& permuteFromLIDs);

  virtual void
  packAndPrepare (const Tpetra::SrcDistObject& source,
                  const Teuchos::ArrayView<const LO>& exportLIDs,
                  Teuchos::Array<packet_type>& exports,
                  const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                  size_t& constantNumPackets,
                  Tpetra::Distributor& distor);

  virtual void
  unpackAndCombine (const Teuchos::ArrayView<const LO> &importLIDs,
                    const Teuchos::ArrayView<const packet_type> &imports,
                    const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    Tpetra::Distributor& distor,
                    Tpetra::CombineMode CM);
#endif // TPETRA_USE_KOKKOS_DISTOBJECT
  //@}

private:
  /// \brief Mesh Map given to constructor.
  ///
  /// This is stored by value, not as a Teuchos::RCP, because the
  /// latter is not thread-safe.  I would like GID->LID lookups to be
  /// thread-safe.
  map_type meshMap_;

  //! The point Map (describing the distribution of degrees of freedom).
  map_type pointMap_;

  //! The Tpetra::MultiVector used to represent the data.
  mv_type mv_;

  /// \brief Raw pointer to the above Tpetra::MultiVector's data.
  ///
  /// Keeping this is a temporary measure that ensures that the
  /// replace, sumInto, and view methods are thread-safe.  (Updates to
  /// ArrayRCP's reference count are not thread safe.)  Eventually,
  /// this will become a Kokkos::View.  It's safe to keep this,
  /// because the data belong to the Vector, which we keep for the
  /// lifetime of the BlockMultiVector.
  Scalar* mvData_;

  //! The number of degrees of freedom per mesh point.
  LO blockSize_;

  //! Implementation of replaceLocalValues; does not check localRowIndex.
  void
  replaceLocalValuesImpl (const LO localRowIndex,
                          const LO colIndex,
                          const Scalar vals[]) const;
  //! Implementation of sumIntoLocalValues; does not check localRowIndex.
  void
  sumIntoLocalValuesImpl (const LO localRowIndex,
                          const LO colIndex,
                          const Scalar vals[]) const;

  static Teuchos::RCP<const mv_type>
  getMultiVectorFromSrcDistObject (const Tpetra::SrcDistObject&);

  static Teuchos::RCP<const BlockMultiVector<Scalar, LO, GO, Node> >
  getBlockMultiVectorFromSrcDistObject (const Tpetra::SrcDistObject& src);
};

} // namespace Experimental
} // namespace Tpetra

#endif // TPETRA_EXPERIMENTAL_BLOCKMULTIVECTOR_DECL_HPP
