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

#ifndef TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP
#define TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP

#include <Kokkos_DefaultArithmetic.hpp>
#include <Teuchos_as.hpp>
#include <Tpetra_Util.hpp>
#include <Tpetra_Vector.hpp>

#include <Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_KokkosRefactor_MultiVector_decl.hpp"
#endif

#include <KokkosCompat_View.hpp>
#include <Kokkos_MV.hpp>

namespace Tpetra {

namespace { // (anonymous)
  template<class DstType, class SrcType, class DeviceType,bool DstConstStride,bool SrcConstStride>
  struct DeepCopySelectedVectors {
    typedef DeviceType device_type;
    DstType dst;
    SrcType src;
    Kokkos::View<int*,DeviceType> whichVectorDst;
    Kokkos::View<int*,DeviceType> whichVectorSrc;
    int n;

    DeepCopySelectedVectors (DstType dst_,
                             SrcType src_,
                             Kokkos::View<int*,DeviceType> whichVectorDst_,
                             Kokkos::View<int*,DeviceType> whichVectorSrc_):
      dst (dst_),
      src (src_),
      whichVectorDst (whichVectorDst_),
      whichVectorSrc (whichVectorSrc_),
      n (whichVectorSrc_.dimension_0 ())
    {}

    void KOKKOS_INLINE_FUNCTION operator() (int i) const {
      if (DstConstStride) {
        if (SrcConstStride) {
          for (int j = 0; j < n; ++j) {
            dst(i,j) = src(i,j);
          }
        } else {
          for (int j = 0; j < n; ++j) {
            dst(i,j) = src(i,whichVectorSrc(j));
          }
        }
      } else {
        if (SrcConstStride) {
          for (int j = 0; j < n; ++j) {
            dst(i,whichVectorDst(j)) = src(i,j);
          }
        } else {
          for (int j = 0; j < n; ++j) {
            dst(i,whichVectorDst(j)) = src(i,whichVectorSrc(j));
          }
        }
      }
    }
  };
} // namespace (anonymous)


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  vectorIndexOutOfRange(size_t VectorIndex) const {
    return (VectorIndex < 1 && VectorIndex != 0) || VectorIndex >= getNumVectors();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector () :
    base_type (Teuchos::rcp (new Map<LocalOrdinal, GlobalOrdinal, node_type> ())),
    lclMV_ (this->getMap ()->getNode ())
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               size_t NumVectors,
               bool zeroOut) : /* default is true */
    base_type (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;

    (void) zeroOut; // View allocation does first touch automatically.

    TEUCHOS_TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument,
      "Tpetra::MultiVector::MultiVector(): NumVectors must be strictly positive.");
    const size_t myLen = getLocalLength();
    if (myLen > 0) {
      RCP<Node> node = map->getNode();
      // On host-type Kokkos Nodes, allocBuffer() just calls the
      // one-argument version of arcp to allocate memory.  This should
      // not fill the memory by default, otherwise we would lose the
      // first-touch allocation optimization.
      //ArrayRCP<Scalar> data = node->template allocBuffer<Scalar>(myLen*NumVectors);

      // Allocate a DualView from new Kokkos, wrap its device data into an ArrayRCP
      view_ = dual_view_type("MV::dual_view",myLen,NumVectors);
      ArrayRCP<Scalar> data = Kokkos::Compat::persistingView(view_.d_view);
      // getting stride of view: if second dimension is 0 stride might
      // be 0, so take view_dimension instead
      size_t stride[8];
      view_.stride (stride);
      const size_t LDA = view_.dimension_1() > 1?stride[1]:view_.dimension_0();
      MVT::initializeValues(lclMV_,myLen,NumVectors,data,LDA);
    }
    else {
      view_ = dual_view_type("MV::dual_view",0,NumVectors);
      MVT::initializeValues(lclMV_,0,NumVectors,Teuchos::null,0);
    }

    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& source) :
    base_type (source),
    lclMV_ (MVT::getNode (source.lclMV_)),
    view_ (source.view_),
    origView_ (source.origView_),
    whichVectors_ (source.whichVectors_)
  {
    using Teuchos::ArrayRCP;

    const LocalOrdinal myLen = source.getLocalLength ();
    const size_t numVecs = source.getNumVectors ();
    ArrayRCP<Scalar> data = (myLen > 0) ?
      Kokkos::Compat::persistingView (view_.d_view) :
      Teuchos::null;
    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();

    // This just sets the dimensions, pointer, and stride of lclMV_.
    // This is only a shallow copy.
    MVT::initializeValues (lclMV_, myLen, numVecs, data, LDA);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& source,
               const Teuchos::DataAccess copyOrView) :
    base_type (source),
    lclMV_ (MVT::getNode (source.lclMV_)),
    view_ (source.view_),
    origView_ (source.origView_),
    whichVectors_ (source.whichVectors_)
  {
    if (copyOrView == Teuchos::Copy) {
      // Reuse the conveniently already existing function that creates
      // a deep copy.
      MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,node_type> cpy =
        createCopy (source);
      this->lclMV_ = cpy.lclMV_;
      this->view_ = cpy.view_;
      this->origView_ = cpy.origView_;
      this->whichVectors_ = cpy.whichVectors_;
    }
    else if (copyOrView == Teuchos::View) {
      const LocalOrdinal myLen = source.getLocalLength ();
      const size_t numVecs = source.getNumVectors ();
      Teuchos::ArrayRCP<Scalar> data = (myLen > 0) ?
        Kokkos::Compat::persistingView (view_.d_view) :
        Teuchos::null;
      // Get stride of view: if second dimension is 0, the
      // stride might be 0, so take view_dimension instead.
      size_t stride[8];
      view_.stride (stride);
      const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();

      // This just sets the dimensions, pointer, and stride of lclMV_.
      // This is only a shallow copy.
      MVT::initializeValues (lclMV_, myLen, numVecs, data, LDA);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::MultiVector copy constructor: "
        "The second argument 'copyOrView' has an invalid value " << copyOrView
        << ".  Valid values include Teuchos::Copy = " << Teuchos::Copy <<
        " and Teuchos::View = " << Teuchos::View << ".");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const dual_view_type& view) :
    base_type (map),
    lclMV_ (map->getNode ()),
    view_ (view),
    origView_ (view)
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "Tpetra::MultiVector(map,view)";

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();
    const size_t numVecs = view_.dimension_1 ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(numVecs < 1, std::invalid_argument,
      ": numVecs must be strictly positive, but you specified numVecs = "
      << numVecs << ".");
    const size_t myLen = getLocalLength ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");

    // This is a shallow copy into the KokkosClassic::MultiVector.
    MVT::initializeValues (lclMV_, myLen, numVecs,
                           Kokkos::Compat::persistingView (view_.d_view), LDA);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const dual_view_type& view,
               const dual_view_type& origView) :
    base_type (map),
    lclMV_ (map->getNode ()),
    view_ (view),
    origView_ (origView)
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "Tpetra::MultiVector(map,view,origView)";

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();
    const size_t numVecs = view_.dimension_1 ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(numVecs < 1, std::invalid_argument,
      ": numVecs must be strictly positive, but you specified numVecs = "
      << numVecs << ".");
    const size_t myLen = getLocalLength ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");

    const size_t origNumRows = origView.dimension_0 ();
    const size_t origNumCols = origView.dimension_1 ();
    // This is a shallow copy into the KokkosClassic::MultiVector.
    MVT::initializeValues (lclMV_, myLen, numVecs,
                           Kokkos::Compat::persistingView (view_.d_view), LDA,
                           origNumRows, origNumCols);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const dual_view_type& view,
               const Teuchos::ArrayView<const size_t>& whichVectors) :
    base_type (map),
    lclMV_ (map->getNode ()),
    view_ (view),
    origView_ (view),
    whichVectors_ (whichVectors.begin (), whichVectors.end ())
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "MultiVector(map,view,whichVectors)";

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();
    size_t numVecs = view_.dimension_1 ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(numVecs < 1, std::invalid_argument,
      ": numVecs must be strictly positive, but you specified numVecs = "
      << numVecs << ".");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");

    if (whichVectors.size () == 1) {
      // If whichVectors has only one entry, we don't need to bother
      // with nonconstant stride.  Just shift the view over so it
      // points to the desired column.
      //
      // NOTE (mfh 10 May 2014) This is a special case where we set
      // origView_ just to view that one column, not all of the
      // original columns.  This ensures that the use of origView_ in
      // offsetView works correctly.
      view_ = subview<dual_view_type> (view_, ALL (), whichVectors[0]);
      origView_ = subview<dual_view_type> (origView_, ALL (), whichVectors[0]);
      numVecs = 1;
      // whichVectors_.size() == 0 means "constant stride."
      whichVectors_.clear ();
    }
    // This is a shallow copy into the KokkosClassic::MultiVector.
    // KokkosClassic::MultiVector doesn't know about whichVectors_.
    MVT::initializeValues (lclMV_, myLen, numVecs,
                           Kokkos::Compat::persistingView (view_.d_view), LDA);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const dual_view_type& view,
               const dual_view_type& origView,
               const Teuchos::ArrayView<const size_t>& whichVectors) :
    base_type (map),
    lclMV_ (map->getNode ()),
    view_ (view),
    origView_ (origView),
    whichVectors_ (whichVectors.begin (), whichVectors.end ())
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "MultiVector(map,view,whichVectors)";

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();
    size_t numVecs = view_.dimension_1 ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(numVecs < 1, std::invalid_argument,
      ": numVecs must be strictly positive, but you specified numVecs = "
      << numVecs << ".");
    const size_t myLen = getLocalLength();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < myLen, std::invalid_argument,
      ": LDA must be large enough to accomodate the local entries.");

    if (whichVectors.size () == 1) {
      // If whichVectors has only one entry, we don't need to bother
      // with nonconstant stride.  Just shift the view over so it
      // points to the desired column.
      //
      // NOTE (mfh 10 May 2014) This is a special case where we set
      // origView_ just to view that one column, not all of the
      // original columns.  This ensures that the use of origView_ in
      // offsetView works correctly.
      view_ = subview<dual_view_type> (view_, ALL (), whichVectors[0]);
      origView_ = subview<dual_view_type> (origView_, ALL (), whichVectors[0]);
      numVecs = 1;
      // whichVectors_.size() == 0 means "constant stride."
      whichVectors_.clear ();
    }

    const size_t origNumRows = origView.dimension_0 ();
    const size_t origNumCols = origView.dimension_1 ();
    // This is a shallow copy into the KokkosClassic::MultiVector.
    // KokkosClassic::MultiVector doesn't know about whichVectors_.
    MVT::initializeValues (lclMV_, myLen, numVecs,
                           Kokkos::Compat::persistingView (view_.d_view), LDA,
                           origNumRows, origNumCols);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const Scalar>& data,
               const size_t LDA,
               const size_t numVecs) :
    base_type (map),
    lclMV_ (map->getNode ())
  {
    // Deep copy constructor, constant stride (NO whichVectors_).
    // There is no need for a deep copy constructor with nonconstant stride.

    const char tfecfFuncName[] = "MultiVector(map,data,LDA,numVecs)";
    const size_t numRows = this->getLocalLength ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(LDA < numRows, std::runtime_error,
      ": LDA = " << LDA << " < numRows = " << numRows << ".");

    view_ = dual_view_type ("MV::view_", numRows, numVecs);
    for (size_t i = 0; i < numRows; ++i) {
      for (size_t j = 0; j < numVecs; ++j) {
        view_.h_view(i,j) = data[j*LDA+i];
      }
    }
    view_.template modify<typename dual_view_type::host_mirror_device_type> ();

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    view_.stride (stride);
    const size_t stride_1 = (view_.dimension_1 () > 1) ?
      stride[1] : view_.dimension_0 ();
    MVT::initializeValues (lclMV_, numRows, numVecs,
                           Kokkos::Compat::persistingView (view_.d_view),
                           stride_1);
    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
               const Teuchos::ArrayView<const ArrayView<const Scalar> >& ArrayOfPtrs,
               const size_t NumVectors) :
    base_type (map),
    lclMV_ (map->getNode ())
  {
    using Teuchos::as;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    const char tfecfFuncName[] = "MultiVector(map,ArrayOfPtrs,NumVectors)";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      NumVectors < 1 || NumVectors != as<size_t>(ArrayOfPtrs.size()),
      std::runtime_error,
      ": ArrayOfPtrs.size() must be strictly positive and as large as ArrayOfPtrs.");
    const size_t myLen = getLocalLength ();
    view_ = dual_view_type ("MV::view_", myLen, NumVectors);

    // TODO: write a functor and use parallel_for.

    // Get stride of view: if second dimension is 0, the
    // stride might be 0, so take view_dimension instead.
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();

    for (size_t i = 0; i < myLen; ++i) {
      for (size_t j = 0; j < NumVectors; ++j) {
        view_.h_view(i,j) = ArrayOfPtrs[j][i];
      }
    }
    view_.template modify<typename dual_view_type::t_host::device_type> ();
    MVT::initializeValues (lclMV_, myLen, NumVectors, Kokkos::Compat::persistingView (view_.d_view), LDA);

    origView_ = view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  ~MultiVector () {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  isConstantStride () const {
    return whichVectors_.empty ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getLocalLength () const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getNodeNumElements ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getGlobalLength () const
  {
    if (this->getMap ().is_null ()) { // possible, due to replaceMap().
      return static_cast<size_t> (0);
    } else {
      return this->getMap ()->getGlobalNumElements ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getStride () const
  {
    if (isConstantStride ()) {
      // Get stride of view: if second dimension is 0, the
      // stride might be 0, so take view_dimension instead.
      size_t stride[8];
      view_.stride (stride);
      const size_t LDA = (view_.dimension_1 () > 1) ? stride[1] : view_.dimension_0 ();
      return LDA;
    }
    else {
      return static_cast<size_t> (0);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  checkSizes (const SrcDistObject& sourceObj)
  {
    // Check whether the source object is a MultiVector.  If not, then
    // we can't even compare sizes, so it's definitely not OK to
    // Import or Export from it.
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> this_type;
    const this_type* src = dynamic_cast<const this_type*> (&sourceObj);
    if (src == NULL) {
      return false;
    } else {
      // The target of the Import or Export calls checkSizes() in
      // DistObject::doTransfer().  By that point, we've already
      // constructed an Import or Export object using the two
      // multivectors' Maps, which means that (hopefully) we've
      // already checked other attributes of the multivectors.  Thus,
      // all we need to do here is check the number of columns.
      return src->getNumVectors () == this->getNumVectors ();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  constantNumberOfPackets () const {
    return this->getNumVectors ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  copyAndPermuteNew (
    const SrcDistObject& sourceObj,
    size_t numSameIDs,
    const Kokkos::View<const LocalOrdinal*, device_type> &permuteToLIDs,
    const Kokkos::View<const LocalOrdinal*, device_type> &permuteFromLIDs)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    //typedef typename ArrayView<const LocalOrdinal>::size_type size_type; // unused
    const char tfecfFuncName[] = "copyAndPermute";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
      ": permuteToLIDs and permuteFromLIDs must have the same size."
      << std::endl << "permuteToLIDs.size() = " << permuteToLIDs.size ()
      << " != permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    const size_t numCols = this->getNumVectors ();
    //const size_t stride = MVT::getStride (lclMV_);

    // TODO (mfh 15 Sep 2013) When we replace
    // KokkosClassic::MultiVector with a Kokkos::View, there are two
    // ways to copy the data:
    //
    // 1. Get a (sub)view of each column and call deep_copy on that.
    // 2. Write a custom kernel to copy the data.
    //
    // The first is easier, but the second might be more performant in
    // case we decide to use layouts other than LayoutLeft.  It might
    // even make sense to hide whichVectors_ in an entirely new layout
    // for Kokkos Views.

    // Copy rows [0, numSameIDs-1] of the local multivectors.
    //
    // For GPU Nodes: All of this happens using device pointers; this
    // does not require host views of either source or destination.
    //
    // Note (ETP 2 Jul 2014)  We need to always copy one column at a
    // time, even when both multivectors are constant-stride, since
    // deep_copy between strided subviews with more than one column
    // doesn't currently work.
    if (numSameIDs > 0) {
      const std::pair<size_t, size_t> rows( 0, numSameIDs );
      for (size_t j = 0; j < numCols; ++j) {
        const size_t dstCol =
          isConstantStride() ? j : whichVectors_[j];
        const size_t srcCol =
          sourceMV.isConstantStride() ? j : sourceMV.whichVectors_[j];
        dual_view_type dst_j =
          Kokkos::subview<dual_view_type>( view_, rows, dstCol );
        dual_view_type src_j =
          Kokkos::subview<dual_view_type>( sourceMV.view_, rows, srcCol );
        Kokkos::deep_copy( dst_j, src_j ); // Copy src_j into dest_j
      }
    }

    // For the remaining GIDs, execute the permutations.  This may
    // involve noncontiguous access of both source and destination
    // vectors, depending on the LID lists.
    //
    // FIXME (mfh 20 June 2012) For an Export with duplicate GIDs on
    // the same process, this merges their values by replacement of
    // the last encountered GID, not by the specified merge rule
    // (such as ADD).

    // If there are no permutations, we are done
    if (permuteFromLIDs.size() == 0 || permuteToLIDs.size() == 0)
      return;

    if (this->isConstantStride ()) {
      KokkosRefactor::Details::permute_array_multi_column(
        getKokkosView(),
        sourceMV.getKokkosView(),
        permuteToLIDs,
        permuteFromLIDs,
        numCols);
    }
    else {
      KokkosRefactor::Details::permute_array_multi_column_variable_stride(
        getKokkosView(),
        sourceMV.getKokkosView(),
        permuteToLIDs,
        permuteFromLIDs,
        getKokkosViewDeepCopy<device_type> (whichVectors_ ()),
        getKokkosViewDeepCopy<device_type> (sourceMV.whichVectors_ ()),
        numCols);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  packAndPrepareNew (
    const SrcDistObject& sourceObj,
    const Kokkos::View<const LocalOrdinal*, device_type> &exportLIDs,
    Kokkos::View<Scalar*, device_type> &exports,
    const Kokkos::View<size_t*, device_type> &numExportPacketsPerLID,
    size_t& constantNumPackets,
    Distributor & /* distor */ )
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    //typedef Array<size_t>::size_type size_type; // unused

    // If we have no exports, there is nothing to do
    if (exportLIDs.size () == 0) {
      return;
    }

    // We've already called checkSizes(), so this cast must succeed.
    const MV& sourceMV = dynamic_cast<const MV&> (sourceObj);

    // We don't need numExportPacketsPerLID; forestall "unused
    // variable" compile warnings.
    (void) numExportPacketsPerLID;

    /* The layout in the export for MultiVectors is as follows:
       exports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    // FIXME (mfh 15 Sep 2013) Would it make sense to rethink the
    // packing scheme in the above comment?  The data going to a
    // particular process must be contiguous, of course, but those
    // data could include entries from multiple LIDs.  DistObject just
    // needs to know how to index into that data.  Kokkos is good at
    // decoupling storage intent from data layout choice.

    const size_t numCols = sourceMV.getNumVectors ();

    // This spares us from needing to fill numExportPacketsPerLID.
    // Setting constantNumPackets to a nonzero value signals that
    // all packets have the same number of entries.
    constantNumPackets = numCols;

    const size_t numExportLIDs = exportLIDs.size ();
    const size_t newExportsSize = numCols * numExportLIDs;
    if (exports.size () != newExportsSize) {
      Kokkos::Compat::realloc (exports, newExportsSize);
    }

    if (numCols == 1) { // special case for one column only
      // MultiVector always represents a single column with constant
      // stride, but it doesn't hurt to implement both cases anyway.
      //
      // ETP:  I'm not sure I agree with the above statement.  Can't a single-
      // column multivector be a subview of another multi-vector, in which case
      // sourceMV.whichVectors_[0] != 0 ?  I think we have to handle that case
      // separately.
      if (sourceMV.isConstantStride ()) {
        KokkosRefactor::Details::pack_array_single_column(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          0);
      }
      else {
        KokkosRefactor::Details::pack_array_single_column(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          sourceMV.whichVectors_[0]);
      }
    }
    else { // the source MultiVector has multiple columns
      if (sourceMV.isConstantStride ()) {
        KokkosRefactor::Details::pack_array_multi_column(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          numCols);
      }
      else {
        KokkosRefactor::Details::pack_array_multi_column_variable_stride(
          exports,
          sourceMV.getKokkosView (),
          exportLIDs,
          getKokkosViewDeepCopy<device_type> (sourceMV.whichVectors_ ()),
          numCols);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  unpackAndCombineNew (
    const Kokkos::View<const LocalOrdinal*, device_type> &importLIDs,
    const Kokkos::View<const Scalar*, device_type> &imports,
    const Kokkos::View<size_t*, device_type> &numPacketsPerLID,
    size_t constantNumPackets,
    Distributor & /* distor */,
    CombineMode CM)
  {
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    const char tfecfFuncName[] = "unpackAndCombine";

    // If we have no imports, there is nothing to do
    if (importLIDs.size () == 0) {
      return;
    }

    // We don't need numPacketsPerLID; forestall "unused variable"
    // compile warnings.
    (void) numPacketsPerLID;

    /* The layout in the export for MultiVectors is as follows:
       imports = { all of the data from row exportLIDs.front() ;
                   ....
                   all of the data from row exportLIDs.back() }
      This doesn't have the best locality, but is necessary because
      the data for a Packet (all data associated with an LID) is
      required to be contiguous. */

    const size_t numVecs = getNumVectors ();

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (imports.size()) != getNumVectors()*importLIDs.size(),
      std::runtime_error,
      ": 'imports' buffer size must be consistent with the amount of data to "
      "be sent.  " << std::endl << "imports.size() = " << imports.size()
      << " != getNumVectors()*importLIDs.size() = " << getNumVectors() << "*"
      << importLIDs.size() << " = " << getNumVectors() * importLIDs.size()
      << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (constantNumPackets) == as<size_t> (0), std::runtime_error,
      ": constantNumPackets input argument must be nonzero.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (numVecs) != as<size_t> (constantNumPackets),
      std::runtime_error, ": constantNumPackets must equal numVecs.");
#endif // HAVE_TPETRA_DEBUG

    if (numVecs > 0 && importLIDs.size () > 0) {
      //const size_t myStride = MVT::getStride (lclMV_);

      // NOTE (mfh 10 Mar 2012, 24 Mar 2014) If you want to implement
      // custom combine modes, start editing here.  Also, if you trust
      // inlining, it would be nice to condense this code by using a
      // binary function object f in the pack functors.
      if (CM == INSERT || CM == REPLACE) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::InsertOp(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<device_type>(whichVectors_ ()),
            KokkosRefactor::Details::InsertOp(),
            numVecs);
        }
      }
      else if (CM == ADD) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::AddOp(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<device_type>(whichVectors_ ()),
            KokkosRefactor::Details::AddOp(),
            numVecs);
        }
      }
      else if (CM == ABSMAX) {
        if (isConstantStride()) {
          KokkosRefactor::Details::unpack_array_multi_column(
            getKokkosView(),
            imports,
            importLIDs,
            KokkosRefactor::Details::AbsMaxOp(),
            numVecs);
        }
        else {
          KokkosRefactor::Details::unpack_array_multi_column_variable_stride(
            getKokkosView(),
            imports,
            importLIDs,
            getKokkosViewDeepCopy<device_type>(whichVectors_ ()),
            KokkosRefactor::Details::AbsMaxOp(),
            numVecs);
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          CM != ADD && CM != REPLACE && CM != INSERT && CM != ABSMAX,
          std::invalid_argument, ": Invalid CombineMode: " << CM << ".  Valid "
          "CombineMode values are ADD, REPLACE, INSERT, and ABSMAX.");
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  inline size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getNumVectors () const
  {
    if (isConstantStride ()) {
      return static_cast<size_t> (view_.dimension_1 ());
    }
    else {
      return static_cast<size_t> (whichVectors_.size ());
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A,
       const Teuchos::ArrayView<dot_type>& dots) const
  {
    typedef typename dual_view_type::host_mirror_device_type host_mirror_device_type;
    typedef Kokkos::View<dot_type*, host_mirror_device_type, Kokkos::MemoryUnmanaged> host_dots_view_type;
    typedef Kokkos::View<dot_type*, device_type> dev_dots_view_type;

    const size_t numDots = static_cast<size_t> (dots.size ());
    host_dots_view_type dotsHostView (dots.getRawPtr (), numDots);
    dev_dots_view_type dotsDevView ("MV::dot tmp", numDots);
    this->dot (A, dotsDevView); // Do the computation on the device.
    Kokkos::deep_copy (dotsHostView, dotsDevView); // Bring back result to host
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& A,
       const Kokkos::View<dot_type*, device_type>& dots) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    // View of a MultiVector's local data (all columns).
    typedef typename dual_view_type::t_dev mv_view_type;
    // View of a single column of a MultiVector's local data.
    //
    // FIXME (mfh 14 Jul 2014) It would be better to get this typedef
    // from mv_view_type itself, in case the layout changes.
    typedef Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type> vec_view_type;
    // View of all the dot product results.
    typedef Kokkos::View<dot_type*, device_type> dots_view_type;
    // Scalar view; view of a single dot product result.
    typedef Kokkos::View<dot_type, device_type> dot_view_type;
    const char tfecfFuncName[] = "Tpetra::MultiVector::dot";

    // FIXME (mfh 11 Jul 2014) These exception tests may not
    // necessarily be thrown on all processes consistently.  We should
    // instead pass along error state with the inner product.  We
    // could do this by setting an extra slot to
    // Kokkos::Details::ArithTraits<dot_type>::one() on error.  The
    // final sum should be
    // Kokkos::Details::ArithTraits<dot_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != A.getLocalLength (), std::runtime_error,
      ": MultiVectors do not have the same local length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getNumVectors () != A.getNumVectors (), std::runtime_error,
      ": MultiVectors must have the same number of columns (vectors).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      dots.dimension_0 () < getNumVectors (), std::runtime_error, ": dots."
      "dimension_0() must be at least as large as the number of columns "
      "(vectors) in *this and A.");

    // We're computing using the device's data, so we need to make
    // sure first that the device is in sync with the host.
    A.view_.template sync<DeviceType> ();
    view_.template sync<DeviceType> ();

    // All the "min"s here ensure that incorrect input won't segfault.
    const size_t numVecs = std::min (getNumVectors (), A.getNumVectors ());
    const size_t lclNumRows = std::min (getLocalLength (), A.getLocalLength ());
    const size_t numDots = std::min (dots.dimension_0 (), numVecs);

    // In case the input dimensions don't match, make sure that we
    // don't overwrite memory that doesn't belong to us, by using
    // subset views with the minimum dimensions over all input.
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numDots);
    dots_view_type theDots = subview<dots_view_type> (dots, colRng);
    mv_view_type X = subview<mv_view_type> (view_.d_view, rowRng, colRng);
    mv_view_type Y = subview<mv_view_type> (A.view_.d_view, rowRng, colRng);

    // FIXME (mfh 14 Jul 2014) How come ALL() works as the first
    // argument, but not a row range?  The first line below doesn't
    // compile, but the second line does.  See
    // kokkos/core/unit_test/TestViewAPI.hpp, in particular
    // run_test_vector(), for an example of allowed subview arguments.
    //
    //vec_view_type X_0 = subview<vec_view_type> (X, rowRng, static_cast<size_t> (0));
    //vec_view_type X_0 = subview<vec_view_type> (X, ALL (), static_cast<size_t> (0));

    if (numDots == 1) {
      // Special case 1: Both MultiVectors only have a single column.
      // The single-vector dot product kernel may be more efficient.
      const size_t ZERO = static_cast<size_t> (0);
      vec_view_type X_k = subview<vec_view_type> (X, ALL (), ZERO);
      vec_view_type Y_k = subview<vec_view_type> (Y, ALL (), ZERO);
      dot_view_type dot_k = subview<dot_view_type> (theDots, ZERO);
      Kokkos::VecDotFunctor<vec_view_type> f (X_k, Y_k, dot_k);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else if (isConstantStride () && A.isConstantStride ()) {
      // Special case 2: Both MultiVectors have constant stride.
      Kokkos::MultiVecDotFunctor<mv_view_type> f (X, Y, theDots);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else {
      // FIXME (mfh 14 Jul 2014) This does a kernel launch for every
      // column.  It might be better to have a kernel that does the
      // work all at once.  On the other hand, we don't prioritize
      // performance of MultiVector views of noncontiguous columns.
      for (size_t k = 0; k < numDots; ++k) {
        const size_t X_col = isConstantStride () ? k : whichVectors_[k];
        const size_t Y_col = A.isConstantStride () ? k : A.whichVectors_[k];
        vec_view_type X_k = subview<vec_view_type> (X, ALL (), X_col);
        vec_view_type Y_k = subview<vec_view_type> (Y, ALL (), Y_col);
        dot_view_type dot_k = subview<dot_view_type> (theDots, k);
        Kokkos::VecDotFunctor<vec_view_type> f (X_k, Y_k, dot_k);
        Kokkos::parallel_reduce (lclNumRows, f);
      }
    }

    // If the MultiVectors are distributed over multiple processes,
    // sum the results across processes.  We assume that the MPI
    // implementation can read from and write to device memory.
    //
    // replaceMap() may have removed some processes.  Those processes
    // have a null Map.  They must not participate in any collective
    // operations.  We ask first whether the Map is null, because
    // isDistributed() defers that question to the Map.  We still
    // compute and return local dot products for processes not
    // participating in collective operations; those probably don't
    // make any sense, but it doesn't hurt to do them, since it's
    // illegal to call dot() on those processes anyway.
    if (! this->getMap ().is_null () && this->isDistributed ()) {
      // MPI doesn't allow aliasing of arguments, so we have to make a
      // copy of the local sum.
      dots_view_type lclDots ("MV::dot lcl", numDots);
      Kokkos::deep_copy (lclDots, theDots);
      const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
      const dot_type* const lclSum = lclDots.ptr_on_device ();
      dot_type* const gblSum = theDots.ptr_on_device ();
      reduceAll<int, dot_type> (comm, REDUCE_SUM, static_cast<int> (numDots),
                                lclSum, gblSum);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  norm2 (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef typename dual_view_type::host_mirror_device_type host_mirror_device_type;
    typedef Kokkos::View<mag_type*, host_mirror_device_type, Kokkos::MemoryUnmanaged> host_norms_view_type;
    typedef Kokkos::View<mag_type*, device_type> dev_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    dev_norms_view_type normsDevView ("MV::norm2 tmp", numNorms);
    this->norm2 (normsDevView); // Do the computation on the device.
    Kokkos::deep_copy (normsHostView, normsDevView); // Bring back result to host
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  norm2 (const Kokkos::View<mag_type*, device_type>& norms) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    // View of a MultiVector's local data (all columns).
    typedef typename dual_view_type::t_dev mv_view_type;
    // View of a single column of a MultiVector's local data.
    //
    // FIXME (mfh 14 Jul 2014) It would be better to get this typedef
    // from mv_view_type itself, in case the layout changes.
    typedef Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type> vec_view_type;
    // View of all the norm results.
    typedef Kokkos::View<mag_type*, device_type> norms_view_type;
    // Scalar view; view of a single norm result.
    typedef Kokkos::View<mag_type, device_type> norm_view_type;
    const char tfecfFuncName[] = "Tpetra::MultiVector::norm2";

    // FIXME (mfh 11 Jul 2014) These exception tests may not
    // necessarily be thrown on all processes consistently.  We should
    // instead pass along error state with the inner product.  We
    // could do this by setting an extra slot to
    // Kokkos::Details::ArithTraits<mag_type>::one() on error.  The
    // final sum should be
    // Kokkos::Details::ArithTraits<mag_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      norms.dimension_0 () < getNumVectors (), std::runtime_error, ": "
      "norms.dimension_0() must be at least as large as the number of "
      "columns (vectors) in *this.");

    // We're computing using the device's data, so we need to make
    // sure first that the device is in sync with the host.
    view_.template sync<DeviceType> ();

    // All the "min"s here ensure that incorrect input won't segfault.
    const size_t numVecs = getNumVectors ();
    const size_t lclNumRows = getLocalLength ();
    const size_t numNorms =
      std::min (static_cast<size_t> (norms.dimension_0 ()), numVecs);

    // In case the input dimensions don't match, make sure that we
    // don't overwrite memory that doesn't belong to us, by using
    // subset views with the minimum dimensions over all input.
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numNorms);
    norms_view_type theNorms = subview<norms_view_type> (norms, colRng);
    mv_view_type X = subview<mv_view_type> (view_.d_view, rowRng, colRng);

    if (numNorms == 1) {
      // Special case 1: The MultiVector only has a single column.
      // The single-vector norm kernel may be more efficient.
      const size_t ZERO = static_cast<size_t> (0);
      vec_view_type X_k = subview<vec_view_type> (X, ALL (), ZERO);
      norm_view_type norm_k = subview<norm_view_type> (theNorms, ZERO);
      Kokkos::VecNorm2SquaredFunctor<vec_view_type> f (X_k, norm_k);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else if (isConstantStride ()) {
      // Special case 2: The MultiVector has constant stride.
      Kokkos::MultiVecNorm2SquaredFunctor<mv_view_type> f (X, theNorms);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else {
      // FIXME (mfh 14 Jul 2014) This does a kernel launch for every
      // column.  It might be better to have a kernel that does the
      // work all at once.  On the other hand, we don't prioritize
      // performance of MultiVector views of noncontiguous columns.
      for (size_t k = 0; k < numNorms; ++k) {
        const size_t X_col = isConstantStride () ? k : whichVectors_[k];
        vec_view_type X_k = subview<vec_view_type> (X, ALL (), X_col);
        norm_view_type norm_k = subview<norm_view_type> (theNorms, k);
        Kokkos::VecNorm2SquaredFunctor<vec_view_type> f (X_k, norm_k);
        Kokkos::parallel_reduce (lclNumRows, f);
      }
    }

    // If the MultiVectors are distributed over multiple processes,
    // sum the results across processes.  We assume that the MPI
    // implementation can read from and write to device memory.
    //
    // replaceMap() may have removed some processes.  Those processes
    // have a null Map.  They must not participate in any collective
    // operations.  We ask first whether the Map is null, because
    // isDistributed() defers that question to the Map.  We still
    // compute and return local norms for processes not participating
    // in collective operations; those probably don't make any sense,
    // but it doesn't hurt to do them, since it's illegal to call
    // norm2() on those processes anyway.
    if (! this->getMap ().is_null () && this->isDistributed ()) {
      // MPI doesn't allow aliasing of arguments, so we have to make a
      // copy of the local sum.
      norms_view_type lclNorms ("MV::norm2 lcl", numNorms);
      Kokkos::deep_copy (lclNorms, theNorms);
      const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
      const mag_type* const lclSum = lclNorms.ptr_on_device ();
      mag_type* const gblSum = theNorms.ptr_on_device ();
      reduceAll<int, mag_type> (comm, REDUCE_SUM, static_cast<int> (numNorms),
                                lclSum, gblSum);

      // Replace the norm-squared results with their square roots in
      // place, to get the final output.  If the device memory and the
      // host memory are the same, it probably doesn't pay to launch a
      // parallel kernel for that, since there isn't enough
      // parallelism for the typical MultiVector case.
      typedef typename device_type::host_mirror_device_type host_mirror_device_type;
      const bool inHostMemory = Kokkos::Impl::is_same<typename device_type::memory_space,
        typename host_mirror_device_type::memory_space>::value;
      if (inHostMemory) {
        for (size_t j = 0; j < numNorms; ++j) {
          theNorms(j) = Kokkos::Details::ArithTraits<mag_type>::sqrt (theNorms(j));
        }
      }
      else {
        // There's not as much parallelism now, but that's OK.  The
        // point of doing parallel dispatch here is to keep the norm
        // results on the device, thus avoiding a copy to the host and
        // back again.
        Kokkos::SquareRootFunctor<norms_view_type> f (theNorms);
        Kokkos::parallel_for (numNorms, f);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  normWeighted (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& weights,
                const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const
  {
    // KR FIXME MVT::WeightedNorm (might already exist).
    //
    // KR FIXME Overload this method to take a View.

    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::arcp_const_cast;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    using Teuchos::ScalarTraits;
    typedef Teuchos::ScalarTraits<Scalar> SCT;
    typedef typename SCT::magnitudeType Mag;
    const char tfecfFuncName[] = "normWeighted";

    const Mag OneOverN = Teuchos::ScalarTraits<Mag>::one () / static_cast<Mag> (getGlobalLength ());
    bool OneW = false;
    const size_t numVecs = this->getNumVectors();
    if (weights.getNumVectors() == 1) {
      OneW = true;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(weights.getNumVectors() != numVecs, std::runtime_error,
        ": MultiVector of weights must contain either one vector or the same number of vectors as this.");
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !this->getMap()->isCompatible(*weights.getMap()), std::runtime_error,
      ": MultiVectors do not have compatible Maps:" << std::endl
      << "this->getMap(): " << std::endl << *this->getMap()
      << "weights.getMap(): " << std::endl << *weights.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != weights.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(as<size_t>(norms.size()) != numVecs, std::runtime_error,
      ": norms.size() must be as large as the number of vectors in *this.");
    typedef Kokkos::View<Scalar*,DeviceType> view_type;
    view_.template sync<DeviceType>();
    weights.view_.template sync<DeviceType>();
    if (isConstantStride ()) {
      if(OneW) {
        view_type weights_0 = subview<view_type> (weights.view_.d_view, ALL (), 0);
        Kokkos::MV_DotWeighted (&norms[0], weights_0, view_.d_view, getLocalLength ());
      } else
        Kokkos::MV_DotWeighted (&norms[0], weights.view_.d_view , view_.d_view, getLocalLength ());
    }
    else {
      // FIXME (mfh 11 Mar 2014) Once we have strided Views, we won't
      // have to write the explicit for loop over columns any more.
      if(OneW) {
        view_type weights_0 = subview<view_type> (weights.view_.d_view, ALL (), 0);
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t curCol = whichVectors_[k];
          view_type vector_k = subview<view_type> (view_.d_view, ALL (), curCol);
          norms[k] = Kokkos::V_DotWeighted (weights_0, vector_k,getLocalLength());
        }
      } else {
        for (size_t k = 0; k < numVecs; ++k) {
          const size_t curCol = whichVectors_[k];
          view_type weights_k = subview<view_type> (weights.view_.d_view, ALL (), curCol);
          view_type vector_k = subview<view_type> (view_.d_view, ALL (), curCol);
          norms[k] = Kokkos::V_DotWeighted (weights_k, vector_k,getLocalLength());
        }
      }
    }
    if (this->isDistributed ()) {
      Array<Mag> lnorms(norms);
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int>(numVecs),
                 lnorms.getRawPtr (), norms.getRawPtr ());
    }
    for (typename ArrayView<Mag>::iterator n = norms.begin(); n != norms.begin()+numVecs; ++n) {
      *n = ScalarTraits<Mag>::squareroot(*n * OneOverN);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  norm1 (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef typename dual_view_type::host_mirror_device_type host_mirror_device_type;
    typedef Kokkos::View<mag_type*, host_mirror_device_type, Kokkos::MemoryUnmanaged> host_norms_view_type;
    typedef Kokkos::View<mag_type*, device_type> dev_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    dev_norms_view_type normsDevView ("MV::norm1 tmp", numNorms);
    this->norm1 (normsDevView); // Do the computation on the device.
    Kokkos::deep_copy (normsHostView, normsDevView); // Bring back result to host
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  norm1 (const Kokkos::View<mag_type*, device_type>& norms) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    // View of a MultiVector's local data (all columns).
    typedef typename dual_view_type::t_dev mv_view_type;
    // View of a single column of a MultiVector's local data.
    //
    // FIXME (mfh 14 Jul 2014) It would be better to get this typedef
    // from mv_view_type itself, in case the layout changes.
    typedef Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type> vec_view_type;
    // View of all the norm results.
    typedef Kokkos::View<mag_type*, device_type> norms_view_type;
    // Scalar view; view of a single norm result.
    typedef Kokkos::View<mag_type, device_type> norm_view_type;
    const char tfecfFuncName[] = "Tpetra::MultiVector::norm1";

    // FIXME (mfh 11 Jul 2014) These exception tests may not
    // necessarily be thrown on all processes consistently.  We should
    // instead pass along error state with the inner product.  We
    // could do this by setting an extra slot to
    // Kokkos::Details::ArithTraits<mag_type>::one() on error.  The
    // final sum should be
    // Kokkos::Details::ArithTraits<mag_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      norms.dimension_0 () < getNumVectors (), std::runtime_error, ": "
      "norms.dimension_0() must be at least as large as the number of "
      "columns (vectors) in *this.");

    // We're computing using the device's data, so we need to make
    // sure first that the device is in sync with the host.
    view_.template sync<DeviceType> ();

    // All the "min"s here ensure that incorrect input won't segfault.
    const size_t numVecs = getNumVectors ();
    const size_t lclNumRows = getLocalLength ();
    const size_t numNorms =
      std::min (static_cast<size_t> (norms.dimension_0 ()), numVecs);

    // In case the input dimensions don't match, make sure that we
    // don't overwrite memory that doesn't belong to us, by using
    // subset views with the minimum dimensions over all input.
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numNorms);
    norms_view_type theNorms = subview<norms_view_type> (norms, colRng);
    mv_view_type X = subview<mv_view_type> (view_.d_view, rowRng, colRng);

    if (numNorms == 1) {
      // Special case 1: The MultiVector only has a single column.
      // The single-vector norm kernel may be more efficient.
      const size_t ZERO = static_cast<size_t> (0);
      vec_view_type X_k = subview<vec_view_type> (X, ALL (), ZERO);
      norm_view_type norm_k = subview<norm_view_type> (theNorms, ZERO);
      Kokkos::VecNorm1Functor<vec_view_type> f (X_k, norm_k);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else if (isConstantStride ()) {
      // Special case 2: The MultiVector has constant stride.
      Kokkos::MultiVecNorm1Functor<mv_view_type> f (X, theNorms);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else {
      // FIXME (mfh 14 Jul 2014) This does a kernel launch for every
      // column.  It might be better to have a kernel that does the
      // work all at once.  On the other hand, we don't prioritize
      // performance of MultiVector views of noncontiguous columns.
      for (size_t k = 0; k < numNorms; ++k) {
        const size_t X_col = isConstantStride () ? k : whichVectors_[k];
        vec_view_type X_k = subview<vec_view_type> (X, ALL (), X_col);
        norm_view_type norm_k = subview<norm_view_type> (theNorms, k);
        Kokkos::VecNorm1Functor<vec_view_type> f (X_k, norm_k);
        Kokkos::parallel_reduce (lclNumRows, f);
      }
    }

    // If the MultiVectors are distributed over multiple processes,
    // sum the results across processes.  We assume that the MPI
    // implementation can read from and write to device memory.
    //
    // replaceMap() may have removed some processes.  Those processes
    // have a null Map.  They must not participate in any collective
    // operations.  We ask first whether the Map is null, because
    // isDistributed() defers that question to the Map.  We still
    // compute and return local norms for processes not participating
    // in collective operations; those probably don't make any sense,
    // but it doesn't hurt to do them, since it's illegal to call
    // norm1() on those processes anyway.
    if (! this->getMap ().is_null () && this->isDistributed ()) {
      // MPI doesn't allow aliasing of arguments, so we have to make a
      // copy of the local sum.
      norms_view_type lclNorms ("MV::norm1 lcl", numNorms);
      Kokkos::deep_copy (lclNorms, theNorms);
      const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
      const mag_type* const lclSum = lclNorms.ptr_on_device ();
      mag_type* const gblSum = theNorms.ptr_on_device ();
      reduceAll<int, mag_type> (comm, REDUCE_SUM, static_cast<int> (numNorms),
                                lclSum, gblSum);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  normInf (const Teuchos::ArrayView<mag_type>& norms) const
  {
    typedef typename dual_view_type::host_mirror_device_type host_mirror_device_type;
    typedef Kokkos::View<mag_type*, host_mirror_device_type, Kokkos::MemoryUnmanaged> host_norms_view_type;
    typedef Kokkos::View<mag_type*, device_type> dev_norms_view_type;

    const size_t numNorms = static_cast<size_t> (norms.size ());
    host_norms_view_type normsHostView (norms.getRawPtr (), numNorms);
    dev_norms_view_type normsDevView ("MV::normInf tmp", numNorms);
    this->normInf (normsDevView); // Do the computation on the device.
    Kokkos::deep_copy (normsHostView, normsDevView); // Bring back result to host
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
              Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  normInf (const Kokkos::View<mag_type*, device_type>& norms) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    // View of a MultiVector's local data (all columns).
    typedef typename dual_view_type::t_dev mv_view_type;
    // View of a single column of a MultiVector's local data.
    //
    // FIXME (mfh 14 Jul 2014) It would be better to get this typedef
    // from mv_view_type itself, in case the layout changes.
    typedef Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type> vec_view_type;
    // View of all the norm results.
    typedef Kokkos::View<mag_type*, device_type> norms_view_type;
    // Scalar view; view of a single norm result.
    typedef Kokkos::View<mag_type, device_type> norm_view_type;
    const char tfecfFuncName[] = "Tpetra::MultiVector::normInf";

    // FIXME (mfh 11 Jul 2014) These exception tests may not
    // necessarily be thrown on all processes consistently.  We should
    // instead pass along error state with the inner product.  We
    // could do this by setting an extra slot to
    // Kokkos::Details::ArithTraits<mag_type>::one() on error.  The
    // final sum should be
    // Kokkos::Details::ArithTraits<mag_type>::zero() if not error.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      norms.dimension_0 () < getNumVectors (), std::runtime_error, ": "
      "norms.dimension_0() must be at least as large as the number of "
      "columns (vectors) in *this.");

    // We're computing using the device's data, so we need to make
    // sure first that the device is in sync with the host.
    view_.template sync<DeviceType> ();

    // All the "min"s here ensure that incorrect input won't segfault.
    const size_t numVecs = getNumVectors ();
    const size_t lclNumRows = getLocalLength ();
    const size_t numNorms =
      std::min (static_cast<size_t> (norms.dimension_0 ()), numVecs);

    // In case the input dimensions don't match, make sure that we
    // don't overwrite memory that doesn't belong to us, by using
    // subset views with the minimum dimensions over all input.
    const std::pair<size_t, size_t> rowRng (0, lclNumRows);
    const std::pair<size_t, size_t> colRng (0, numNorms);
    norms_view_type theNorms = subview<norms_view_type> (norms, colRng);
    mv_view_type X = subview<mv_view_type> (view_.d_view, rowRng, colRng);

    // All of the functors set the local infinity-norm to zero if the
    // MultiVector has zero rows.  Thus, any MPI processes with zero
    // rows don't contribute to the global maximum.
    if (numNorms == 1) {
      // Special case 1: The MultiVector only has a single column.
      // The single-vector norm kernel may be more efficient.
      const size_t ZERO = static_cast<size_t> (0);
      vec_view_type X_k = subview<vec_view_type> (X, ALL (), ZERO);
      norm_view_type norm_k = subview<norm_view_type> (theNorms, ZERO);
      Kokkos::VecNormInfFunctor<vec_view_type> f (X_k, norm_k);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else if (isConstantStride ()) {
      // Special case 2: The MultiVector has constant stride.
      Kokkos::MultiVecNormInfFunctor<mv_view_type> f (X, theNorms);
      Kokkos::parallel_reduce (lclNumRows, f);
    }
    else {
      // FIXME (mfh 15 Jul 2014) This does a kernel launch for every
      // column.  It might be better to have a kernel that does the
      // work all at once.  On the other hand, we don't prioritize
      // performance of MultiVector views of noncontiguous columns.
      for (size_t k = 0; k < numNorms; ++k) {
        const size_t X_col = isConstantStride () ? k : whichVectors_[k];
        vec_view_type X_k = subview<vec_view_type> (X, ALL (), X_col);
        norm_view_type norm_k = subview<norm_view_type> (theNorms, k);
        Kokkos::VecNormInfFunctor<vec_view_type> f (X_k, norm_k);
        Kokkos::parallel_reduce (lclNumRows, f);
      }
    }

    // If the MultiVectors are distributed over multiple processes,
    // compute the max of the results across processes.  We assume
    // that the MPI implementation can read from and write to device
    // memory.
    //
    // replaceMap() may have removed some processes.  Those processes
    // have a null Map.  They must not participate in any collective
    // operations.  We ask first whether the Map is null, because
    // isDistributed() defers that question to the Map.  We still
    // compute and return local norms for processes not participating
    // in collective operations; those probably don't make any sense,
    // but it doesn't hurt to do them, since it's illegal to call
    // normInf() on those processes anyway.
    if (! this->getMap ().is_null () && this->isDistributed ()) {
      // MPI doesn't allow aliasing of arguments, so we have to make a
      // copy of the local sum.
      norms_view_type lclNorms ("MV::normInf lcl", numNorms);
      Kokkos::deep_copy (lclNorms, theNorms);
      const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
      const mag_type* const lclSum = lclNorms.ptr_on_device ();
      mag_type* const gblSum = theNorms.ptr_on_device ();
      reduceAll<int, mag_type> (comm, REDUCE_MAX, static_cast<int> (numNorms),
                                lclSum, gblSum);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  meanValue (const Teuchos::ArrayView<Scalar> &means) const
  {
    // KR FIXME MVT::Sum (might already exist).
    //
    // KR FIXME Overload this method to take a View.

    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::arcp_const_cast;
    using Teuchos::as;
    using Teuchos::arcp_const_cast;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;

    using Kokkos::ALL;
    using Kokkos::subview;

    typedef Teuchos::ScalarTraits<Scalar> SCT;

    const size_t numVecs = getNumVectors();

    TEUCHOS_TEST_FOR_EXCEPTION(as<size_t>(means.size()) != numVecs, std::runtime_error,
      "Tpetra::MultiVector::meanValue(): means.size() must be as large as the number of vectors in *this.");
    // compute local components of the means
    // sum these across all nodes
    view_.template sync<DeviceType>();
    if (isConstantStride ()) {
      Kokkos::MV_Sum (&means[0], view_.d_view, getLocalLength ());
    }
    else {
      // FIXME (mfh 11 Mar 2014) Once we have strided Views, we won't
      // have to write the explicit for loop over columns any more.
      for (size_t k = 0; k < numVecs; ++k) {
        typedef Kokkos::View<Scalar*,DeviceType> view_type;
        const size_t curCol = whichVectors_[k];
        view_type vector_k = subview<view_type> (view_.d_view, ALL (), curCol);
        means[k] = Kokkos::V_Sum (vector_k);
      }
    }
    if (this->isDistributed()) {
      Array<Scalar> lmeans(means);
      // only combine if we are a distributed MV
      reduceAll (*this->getMap ()->getComm (), REDUCE_SUM, as<int> (numVecs),
                 lmeans.getRawPtr (), means.getRawPtr ());
    }
    // mfh 12 Apr 2012: Don't take out the cast from the ordinal type
    // to the magnitude type, since operator/ (std::complex<T>, int)
    // isn't necessarily defined.
    const Scalar OneOverN =
      SCT::one() / as<typename SCT::magnitudeType> (getGlobalLength ());
    for (typename ArrayView<Scalar>::iterator i = means.begin();
         i != means.begin() + numVecs;
         ++i)
    {
      (*i) = (*i) * OneOverN;
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  randomize ()
  {
    // KR FIXME *sigh* MVT::Random.  Hack around this for now by
    // filling on the host.

    if (isConstantStride ()) {
      MVT::Random (lclMV_);
    }
    else {
      const size_t numVecs = this->getNumVectors ();
      KMV v (MVT::getNode (lclMV_));
      Teuchos::ArrayRCP<Scalar> vj;
      for (size_t j = 0; j < numVecs; ++j) {
        vj = MVT::getValuesNonConst (lclMV_, whichVectors_[j]);
        MVT::initializeValues (v, MVT::getNumRows (lclMV_), 1, vj, MVT::getStride (lclMV_));
        MVT::Random (v);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  putScalar (const Scalar& alpha)
  {
    const size_t numVecs = getNumVectors();
    if (isConstantStride ()) {
      Kokkos::Impl::ViewFill<typename dual_view_type::t_dev> (view_.d_view, alpha);
    }
    else {
      typedef Kokkos::View<Scalar*, DeviceType> view_type;
      for (size_t k = 0; k < numVecs; ++k) {
        view_type vector_k = Kokkos::subview<view_type> (view_.d_view, Kokkos::ALL (), whichVectors_[k]);
        Kokkos::Impl::ViewFill<view_type> (vector_k, alpha);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  replaceMap (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& newMap)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::Comm;
    using Teuchos::RCP;

    // mfh 28 Mar 2013: This method doesn't forget whichVectors_, so
    // it might work if the MV is a column view of another MV.
    // However, things might go wrong when restoring the original
    // Map, so we don't allow this case for now.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! this->isConstantStride (), std::logic_error,
      "Tpetra::MultiVector::replaceMap: This method does not currently work "
      "if the MultiVector is a column view of another MultiVector (that is, if "
      "isConstantStride() == false).");

    // Case 1: current Map and new Map are both nonnull on this process.
    // Case 2: current Map is nonnull, new Map is null.
    // Case 3: current Map is null, new Map is nonnull.
    // Case 4: both Maps are null: forbidden.
    //
    // Case 1 means that we don't have to do anything on this process,
    // other than assign the new Map.  (We always have to do that.)
    // It's an error for the user to supply a Map that requires
    // resizing in this case.
    //
    // Case 2 means that the calling process is in the current Map's
    // communicator, but will be excluded from the new Map's
    // communicator.  We don't have to do anything on the calling
    // process; just leave whatever data it may have alone.
    //
    // Case 3 means that the calling process is excluded from the
    // current Map's communicator, but will be included in the new
    // Map's communicator.  This means we need to (re)allocate the
    // local (KokkosClassic::)MultiVector if it does not have the right
    // number of rows.  If the new number of rows is nonzero, we'll
    // fill the newly allocated local data with zeros, as befits a
    // projection operation.
    //
    // The typical use case for Case 3 is that the MultiVector was
    // first created with the Map with more processes, then that Map
    // was replaced with a Map with fewer processes, and finally the
    // original Map was restored on this call to replaceMap.

#ifdef HAVE_TEUCHOS_DEBUG
    // mfh 28 Mar 2013: We can't check for compatibility across the
    // whole communicator, unless we know that the current and new
    // Maps are nonnull on _all_ participating processes.
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   origNumProcs == newNumProcs && ! this->getMap ()->isCompatible (*map),
    //   std::invalid_argument, "Tpetra::MultiVector::project: "
    //   "If the input Map's communicator is compatible (has the same number of "
    //   "processes as) the current Map's communicator, then the two Maps must be "
    //   "compatible.  The replaceMap() method is not for data redistribution; "
    //   "use Import or Export for that purpose.");

    // TODO (mfh 28 Mar 2013) Add compatibility checks for projections
    // of the Map, in case the process counts don't match.
#endif // HAVE_TEUCHOS_DEBUG

    if (this->getMap ().is_null ()) { // current Map is null
      // If this->getMap() is null, that means that this MultiVector
      // has already had replaceMap happen to it.  In that case, just
      // reallocate the DualView with the right size.

      TEUCHOS_TEST_FOR_EXCEPTION(
        newMap.is_null (), std::invalid_argument,
        "Tpetra::MultiVector::replaceMap: both current and new Maps are null.  "
        "This probably means that the input Map is incorrect.");

      // Case 3: current Map is null, new Map is nonnull.
      // Reallocate the DualView with the right dimensions.
      const size_t newNumRows = newMap->getNodeNumElements ();
      const size_t origNumRows = view_.dimension_0 ();
      const size_t numCols = this->getNumVectors ();

      if (origNumRows != newNumRows || view_.dimension_1 () != numCols) {
        view_ = dual_view_type ("MV::dual_view", newNumRows, numCols);

        // KR FIXME (mfh 10 Jul 2014) Once MultiVector no longer needs
        // lclMV_, remove the code below.
        Teuchos::ArrayRCP<Scalar> data =
          Kokkos::Compat::persistingView (view_.d_view);
        const size_t stride = newNumRows;
        MVT::initializeValues (lclMV_, newNumRows, numCols, data, stride);
      }
    }
    else if (newMap.is_null ()) { // Case 2: current Map is nonnull, new Map is null
      // I am an excluded process.  Reinitialize my data so that I
      // have 0 rows.  Keep the number of columns as before.
      const size_t newNumRows = static_cast<size_t> (0);
      const size_t numCols = this->getNumVectors ();
      view_ = dual_view_type ("MV::dual_view", newNumRows, numCols);

      // KR FIXME (mfh 10 Jul 2014) Once MultiVector no longer needs
      // lclMV_, remove the code below.
      Teuchos::ArrayRCP<Scalar> data =
        Kokkos::Compat::persistingView (view_.d_view);
      const size_t stride = newNumRows;
      MVT::initializeValues (lclMV_, newNumRows, numCols, data, stride);
    }

    this->map_ = newMap;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  scale (const Scalar &alpha)
  {
    using Kokkos::ALL;
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    // NOTE: can't substitute putScalar(0.0) for scale(0.0), because
    //       the former will overwrite NaNs present in the MultiVector, while the
    //       semantics of this call require multiplying them by 0, which IEEE requires to be NaN
    const size_t numVecs = getNumVectors();
    if (alpha == Teuchos::ScalarTraits<Scalar>::one()) {
      // do nothing
    }
    else if (isConstantStride ()) {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      Kokkos::MV_MulScalar (view_.d_view, alpha, view_.d_view);
    }
    else {
      typedef Kokkos::View<Scalar*, DeviceType> view_type;

      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = Kokkos::subview<view_type> (view_.d_view, ALL (), this_col);
        Kokkos::V_MulScalar (vector_k, alpha, vector_k);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  scale (Teuchos::ArrayView<const Scalar> alphas)
  {
    using Kokkos::ALL;
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (alphas.size ()) != numVecs, std::invalid_argument,
      "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as "
      "the number of vectors in *this.");

    const size_t myLen = view_.dimension_0();
    if (myLen == 0) {
      return;
    }

    if (isConstantStride ()) {
      Kokkos::DualView<Scalar*,device_type> k_alphas("Alphas::tmp",alphas.size());
      for(int i=0; i<alphas.size(); i++)
         k_alphas.h_view(i) = alphas[i];
      k_alphas.template modify<host_mirror_device_type>();
      k_alphas.template sync<device_type>();
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      Kokkos::MV_MulScalar (view_.d_view, k_alphas.d_view, view_.d_view);
    }
    else {
      typedef Kokkos::View<Scalar*, DeviceType> view_type;

      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = Kokkos::subview<view_type> (view_.d_view, ALL (), this_col);
        Kokkos::V_MulScalar (vector_k, alphas[k], vector_k);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  scale (const Kokkos::View<const Scalar*, device_type> alphas)
  {
    using Kokkos::ALL;
    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;

    const size_t numVecs = this->getNumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (alphas.dimension_0()) != numVecs, std::invalid_argument,
      "Tpetra::MultiVector::scale(alphas): alphas.size() must be as large as "
      "the number of vectors in *this.");

    const size_t myLen = view_.dimension_0();
    if (myLen == 0) {
      return;
    }

    if (isConstantStride ()) {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      Kokkos::MV_MulScalar (view_.d_view, alphas, view_.d_view);
    }
    else {
      typedef Kokkos::View<Scalar*, DeviceType> view_type;

      typename Kokkos::View<const Scalar*, device_type>::HostMirror h_alphas = Kokkos::create_mirror_view(alphas);
      Kokkos::deep_copy(h_alphas,alphas);
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = Kokkos::subview<view_type> (view_.d_view, ALL (), this_col);
        Kokkos::V_MulScalar (vector_k, alphas(k), vector_k);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  scale (const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A)
  {
    // KR FIXME Need MVT::Scale: Kernels for this already exist.
    //
    // *this := alpha * A

    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;

    const char tfecfFuncName[] = "scale(alpha,A)";

    const size_t numVecs = getNumVectors();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != numVecs, std::runtime_error,
      ": MultiVectors must have the same number of vectors.");
    if (isConstantStride() && A.isConstantStride()) {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      Kokkos::MV_MulScalar (view_.d_view, alpha, A.view_.d_view);
    }
    else {
      typedef Kokkos::View<Scalar*, DeviceType> view_type;

      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      A.view_.template sync<DeviceType>();
      A.view_.template modify<DeviceType>();
      for (size_t k=0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = Kokkos::subview<view_type> (view_.d_view, Kokkos::ALL (), this_col);
        const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
        view_type vector_Ak = Kokkos::subview<view_type> (A.view_.d_view, Kokkos::ALL (), A_col);
        Kokkos::V_MulScalar (vector_k, alpha, vector_Ak);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  reciprocal (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > &A)
  {
    // KR FIXME Need MVT::Recip (assign the elementwise reciprocal of
    // A to *this).

    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    const char tfecfFuncName[] = "reciprocal";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
      ": MultiVectors must have the same number of vectors.");
    const size_t numVecs = getNumVectors();
    try {
      if (isConstantStride() && A.isConstantStride()) {
        view_.template sync<DeviceType>();
        view_.template modify<DeviceType>();
        Kokkos::MV_Reciprocal(view_.d_view,A.view_.d_view);
      }
      else {
        typedef Kokkos::View<Scalar*, DeviceType> view_type;

        view_.template sync<DeviceType>();
        view_.template modify<DeviceType>();
        A.view_.template sync<DeviceType>();
        A.view_.template modify<DeviceType>();
        for (size_t k=0; k < numVecs; ++k) {
          const size_t this_col = isConstantStride () ? k : whichVectors_[k];
          view_type vector_k = Kokkos::subview<view_type> (view_.d_view, Kokkos::ALL (), this_col);
          const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
          view_type vector_Ak = Kokkos::subview<view_type> (A.view_.d_view, Kokkos::ALL (), A_col);
          Kokkos::V_Reciprocal(vector_k, vector_Ak);
        }
      }
    }
    catch (std::runtime_error &e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true,std::runtime_error,
          ": caught exception from Kokkos:" << std::endl
          << e.what() << std::endl);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  abs (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A)
  {
    // KR FIXME Need MVT::Abs (assign the elementwise absolute value
    // of A to *this).

    using Teuchos::arcp_const_cast;
    using Teuchos::ArrayRCP;
    const char tfecfFuncName[] = "abs";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getLocalLength() != A.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(A.getNumVectors() != this->getNumVectors(), std::runtime_error,
      ": MultiVectors must have the same number of vectors.");

    const size_t numVecs = getNumVectors();
    if (isConstantStride() && A.isConstantStride()) {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      Kokkos::MV_Abs(view_.d_view,A.view_.d_view);
    }
    else {
      typedef Kokkos::View<Scalar*, DeviceType> view_type;

      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      A.view_.template sync<DeviceType>();
      A.view_.template modify<DeviceType>();
      for (size_t k=0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = Kokkos::subview<view_type> (view_.d_view, Kokkos::ALL (), this_col);
        const size_t A_col = isConstantStride () ? k : A.whichVectors_[k];
        view_type vector_Ak = Kokkos::subview<view_type> (A.view_.d_view, Kokkos::ALL (), A_col);
        Kokkos::V_Abs(vector_k, vector_Ak);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  update (const Scalar& alpha,
          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
          const Scalar& beta)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::View<Scalar*, DeviceType> view_type;
    const char tfecfFuncName[] = "update";

    // this = beta*this + alpha*A
    // must support case where &this == &A
    // can't short circuit on alpha==0.0 or beta==0.0, because 0.0*NaN != 0.0

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != A.getLocalLength (), std::invalid_argument,
      ": The input MultiVector A has " << A.getLocalLength () << " local "
      "row(s), but this MultiVector has " << getLocalLength () << " local "
      "row(s).");

    const size_t numVecs = getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != numVecs, std::invalid_argument,
      ": The input MultiVector A has " << A.getNumVectors () << " column(s), "
      "but this MultiVector has " << numVecs << " column(s).");

    if (isConstantStride () && A.isConstantStride ()) {
      Kokkos::MV_Add (view_.d_view, alpha, A.view_.d_view, beta, view_.d_view,
                      getLocalLength ());
    }
    else {
      // TODO: make sure it only uses LocalLength for add.
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        const size_t A_col = A.isConstantStride () ? k : A.whichVectors_[k];
        Kokkos::V_Add (subview<view_type> (view_.d_view, ALL (), this_col),
                       alpha,
                       subview<view_type> (A.view_.d_view, ALL (), A_col),
                       beta,
                       subview<view_type> (view_.d_view, ALL (), this_col));
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  update (const Scalar& alpha,
          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
          const Scalar& beta,
          const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& B,
          const Scalar& gamma)
  {
    using Kokkos::ALL;
    using Kokkos::V_Add;
    using Kokkos::subview;
    typedef Kokkos::View<Scalar*, DeviceType> view_type;
    const char tfecfFuncName[] = "update";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != A.getLocalLength (), std::invalid_argument,
      ": The input MultiVector A has " << A.getLocalLength () << " local "
      "row(s), but this MultiVector has " << getLocalLength () << " local "
      "row(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength () != B.getLocalLength (), std::invalid_argument,
      ": The input MultiVector B has " << B.getLocalLength () << " local "
      "row(s), but this MultiVector has " << getLocalLength () << " local "
      "row(s).");
    const size_t numVecs = getNumVectors ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      A.getNumVectors () != numVecs, std::invalid_argument,
      ": The input MultiVector A has " << A.getNumVectors () << " column(s), "
      "but this MultiVector has " << numVecs << " column(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      B.getNumVectors () != numVecs, std::invalid_argument,
      ": The input MultiVector B has " << B.getNumVectors () << " column(s), "
      "but this MultiVector has " << numVecs << " column(s).");

    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero ();
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one ();

    if (isConstantStride() && A.isConstantStride () && B.isConstantStride ()) {
      if (gamma == zero) {
        Kokkos::MV_Add (view_.d_view, alpha, A.view_.d_view, beta,
                        B.view_.d_view);
      } else {
        Kokkos::MV_Add (view_.d_view, alpha, A.view_.d_view, gamma, view_.d_view);
        Kokkos::MV_Add (view_.d_view, beta, B.view_.d_view, one, view_.d_view);
      }
    } else { // some input (or *this) is not constant stride
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        const size_t A_col = A.isConstantStride () ? k : A.whichVectors_[k];
        const size_t B_col = B.isConstantStride () ? k : B.whichVectors_[k];
        if (gamma == zero) {
          // TODO: make sure it only uses LocalLength for add.
          V_Add (subview<view_type> (view_.d_view, ALL (), this_col),
                 alpha,
                 subview<view_type> (A.view_.d_view, ALL (), A_col),
                 beta,
                 subview<view_type> (B.view_.d_view, ALL (), B_col));
        } else {
          V_Add (subview<view_type> (view_.d_view, ALL (), this_col),
                 alpha,
                 subview<view_type> (A.view_.d_view, ALL (), A_col),
                 gamma,
                 subview<view_type> (view_.d_view, ALL (), this_col));
          V_Add (subview<view_type> (view_.d_view, ALL (), this_col),
                 beta,
                 subview<view_type> (B.view_.d_view, ALL (), B_col),
                 one,
                 subview<view_type> (view_.d_view, ALL (), this_col));
        }
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getData (size_t j) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef typename dual_view_type::host_mirror_device_type host_type;
    typedef typename dual_view_type::t_host host_view_type;

    // NOTE (mfh 09 2014) Any MultiVector method that called the
    // (classic) Kokkos Node's viewBuffer or viewBufferNonConst
    // methods always implied a device->host synchronization.  Thus,
    // we synchronize here as well, though might want to change this
    // in the future.
    view_.template sync<host_type> ();

    // Get a host view of the entire MultiVector's data.
    host_view_type hostView = view_.template view<host_type> ();
    // Get a subview of column j.
    host_view_type hostView_j;
    if (isConstantStride ()) {
      hostView_j = subview<host_view_type> (hostView, ALL (), j);
    } else {
      hostView_j = subview<host_view_type> (hostView, ALL (), whichVectors_[j]);
    }

    // Wrap up the subview of column j in an ArrayRCP<const scalar_type>.
    Teuchos::ArrayRCP<scalar_type> dataAsArcp =
      Kokkos::Compat::persistingView (hostView_j, 0, getLocalLength ());
    return Teuchos::arcp_const_cast<const scalar_type> (dataAsArcp);

    // Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
    // return node->template viewBuffer<Scalar> (getLocalLength (),
    //                                           getSubArrayRCP (MVT::getValues (lclMV_), j));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getDataNonConst (size_t j)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef typename dual_view_type::host_mirror_device_type host_type;
    typedef typename dual_view_type::t_host host_view_type;

    // NOTE (mfh 09 2014) Any MultiVector method that called the
    // (classic) Kokkos Node's viewBuffer or viewBufferNonConst
    // methods always implied a device->host synchronization.  Thus,
    // we synchronize here as well, though might want to change this
    // in the future.
    view_.template sync<host_type> ();

    // Get a host view of the entire MultiVector's data.
    host_view_type hostView = view_.template view<host_type> ();
    // Get a subview of column j.
    host_view_type hostView_j;
    if (isConstantStride ()) {
      hostView_j = subview<host_view_type> (hostView, ALL (), j);
    } else {
      hostView_j = subview<host_view_type> (hostView, ALL (), whichVectors_[j]);
    }

    // FIXME (mfh 10 May 2014) Calling getDataNonConst() implies that
    // the user plans to modify the values in the MultiVector, so we
    // should call modify on the view here.

    // Wrap up the subview of column j in an ArrayRCP<const scalar_type>.
    Teuchos::ArrayRCP<scalar_type> dataAsArcp =
      Kokkos::Compat::persistingView (hostView_j, 0, getLocalLength ());
    return dataAsArcp;

    // Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
    // return node->template viewBufferNonConst<Scalar> (KokkosClassic::ReadWrite,
    //                                                   getLocalLength (),
    //                                                   getSubArrayRCP (MVT::getValuesNonConst (lclMV_), j));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >&
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  operator= (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& source)
  {
    // KR FIXME This is just the assignment operator.

    if (this != &source) {
      base_type::operator= (source);
      //
      // operator= implements view semantics (shallow copy).
      //

      // OK: Kokkos::View operator= also implements view semantics.
      view_ = source.view_;
      origView_ = source.origView_;

      // NOTE (mfh 24 Mar 2014) Christian wrote here that assigning
      // whichVectors_ is "probably not ok" (probably constitutes deep
      // copy).  I would say that it's OK, because whichVectors_ is
      // immutable (from the user's perspective); it's analogous to
      // the dimensions or stride.  Once we make whichVectors_ a
      // Kokkos::View instead of a Teuchos::Array, all debate will go
      // away and we will unquestionably have view semantics.
      whichVectors_ = source.whichVectors_;

      const LocalOrdinal myLen = source.getLocalLength ();
      const size_t numVecs = source.getNumVectors ();
      Teuchos::ArrayRCP<Scalar> data = (myLen > 0) ?
        Kokkos::Compat::persistingView (view_.d_view) :
        Teuchos::null;
      const size_t stride = (myLen > 0) ?
        static_cast<size_t> (myLen) :
        static_cast<size_t> (0);
      // This just sets the dimensions, pointer, and stride of lclMV_.
      MVT::initializeValues (lclMV_, myLen, numVecs, data, stride);
    }
    return *this;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subCopy (const Teuchos::ArrayView<const size_t> &cols) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef typename dual_view_type::host_mirror_device_type
      host_mirror_device_type;
    typedef typename dual_view_type::t_host host_view_type;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;

    // Sync the source MultiVector (*this) to host first.  Copy it to
    // the output View on host, then sync the output View (only) to
    // device.  Doing copies on host saves us the trouble of copying
    // whichVecsSrc and whichVecsDst over to the device.
    view_.template sync<host_mirror_device_type> ();

    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const std::pair<size_t, size_t> rowRange (0, numRows);
    const std::pair<size_t, size_t> colRange (0, numCols);
    const size_t numColsToCopy = static_cast<size_t> (cols.size ());

    // Create a DualView which will be a contiguously stored deep copy of this MV's view.
    dual_view_type dstView ("MV::dual_view", numRows, numColsToCopy);
    Kokkos::View<int*, host_mirror_device_type> whichVecsDst ("whichVecsDst", numColsToCopy);
    Kokkos::View<int*, host_mirror_device_type> whichVecsSrc ("whichVecsSrc", numColsToCopy);

    if (! this->isConstantStride ()) {
      for (size_t j = 0; j < numColsToCopy; ++j) {
        whichVecsSrc(j) = static_cast<int> (this->whichVectors_[cols[j]]);
      }
    }
    else {
      for (size_t j = 0; j < numColsToCopy; ++j) {
        whichVecsSrc(j) = static_cast<int> (cols[j]);
      }
    }
    for (size_t j = 0; j < numColsToCopy; ++j) {
      whichVecsDst(j) = static_cast<int> (j);
    }

    //
    // Do the copy on host first.
    //
    view_.template modify<host_mirror_device_type> ();
    host_view_type srcView =
      Kokkos::subview<host_view_type> (view_.h_view, rowRange, colRange);
    DeepCopySelectedVectors<host_view_type,
      host_view_type, host_mirror_device_type,
      false, false> f (dstView.h_view, srcView, whichVecsDst, whichVecsSrc);
    Kokkos::parallel_for (numRows, f);

    // Sync the output DualView (only) back to device.
    dstView.template sync<device_type> ();

    // Create and return a MultiVector using the new DualView.
    return rcp (new MV (this->getMap (), dstView));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subCopy (const Teuchos::Range1D &colRng) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef typename dual_view_type::host_mirror_device_type
      host_mirror_device_type;
    typedef typename dual_view_type::t_host host_view_type;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;

    // Sync the source MultiVector (*this) to host first.  Copy it to
    // the output View on host, then sync the output View (only) to
    // device.  Doing copies on host saves us the trouble of copying
    // whichVecsSrc and whichVecsDst over to the device.
    view_.template sync<host_mirror_device_type> ();

    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const std::pair<size_t, size_t> rowRange (0, numRows);
    const std::pair<size_t, size_t> colRange (0, numCols);

    // Range1D is an inclusive range.  If the upper bound is less then
    // the lower bound, that signifies an invalid range, which we
    // interpret as "copy zero columns."
    const size_t numColsToCopy = (colRng.ubound () >= colRng.lbound ()) ?
      static_cast<size_t> (colRng.size ()) :
      static_cast<size_t> (0);

    // Create a DualView which will be a contiguously stored deep copy of this MV's view.
    dual_view_type dstView ("MV::dual_view", numRows, numColsToCopy);
    Kokkos::View<int*, host_mirror_device_type> whichVecsDst ("whichVecsDst", numColsToCopy);
    Kokkos::View<int*, host_mirror_device_type> whichVecsSrc ("whichVecsSrc", numColsToCopy);

    if (! this->isConstantStride ()) {
      for (size_t j = 0; j < numColsToCopy; ++j) {
        const size_t col = static_cast<size_t> (colRng.lbound ()) + j;
        whichVecsSrc(j) = static_cast<int> (this->whichVectors_[col]);
      }
    }
    else {
      for (size_t j = 0; j < numColsToCopy; ++j) {
        const size_t col = static_cast<size_t> (colRng.lbound ()) + j;
        whichVecsSrc(j) = static_cast<int> (col);
      }
    }
    for (size_t j = 0; j < numColsToCopy; ++j) {
      whichVecsDst(j) = static_cast<int> (j);
    }

    //
    // Do the copy on host first.
    //
    // FIXME (mfh 10 Jul 2014) Exploit contiguity of the desired columns.
    view_.template modify<host_mirror_device_type> ();
    host_view_type srcView =
      Kokkos::subview<host_view_type> (view_.h_view, rowRange, colRange);
    DeepCopySelectedVectors<host_view_type,
      host_view_type, host_mirror_device_type,
      false, false> f (dstView.h_view, srcView, whichVecsDst, whichVecsSrc);
    Kokkos::parallel_for (numRows, f);

    // Sync the output DualView (only) back to device.
    dstView.template sync<device_type> ();

    // Create and return a MultiVector using the new DualView.
    return rcp (new MV (this->getMap (), dstView));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getOrigNumLocalRows () const {
    return origView_.dimension_0 ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getOrigNumLocalCols () const {
    return origView_.dimension_1 ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
              size_t offset) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    const size_t newNumRows = subMap->getNodeNumElements ();
    const bool tooManyElts = newNumRows + offset > this->getOrigNumLocalRows ();
    if (tooManyElts) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > this->getLocalLength (), std::runtime_error,
        "Tpetra::MultiVector::offsetView(NonConst): Invalid input Map.  The "
        "input Map owns " << newNumRows << " entries on process " << myRank <<
        ".  offset = " << offset << ".  Yet, the MultiVector contains only "
        << this->getOrigNumLocalRows () << " rows on this process.");
    }

    const std::pair<size_t, size_t> offsetPair (offset, offset + newNumRows);
    // FIXME (mfh 10 May 2014) Use of origView_ instead of view_ for
    // the second argument may be wrong, if view_ resulted from a
    // previous call to offsetView with offset != 0.
    dual_view_type newView =
      subview<dual_view_type> (origView_, offsetPair, ALL ());
    RCP<const MV> subViewMV;
    if (isConstantStride ()) {
      subViewMV = rcp (new MV (subMap, newView, origView_));
    }
    else {
      subViewMV = rcp (new MV (subMap, newView, origView_, whichVectors_ ()));
    }
    return subViewMV;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  offsetViewNonConst (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
                      size_t offset)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    return Teuchos::rcp_const_cast<MV> (this->offsetView (subMap, offset));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subView (const ArrayView<const size_t>& cols) const
  {
    using Teuchos::Array;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    TEUCHOS_TEST_FOR_EXCEPTION(
      cols.size () == 0, std::runtime_error,
      "Tpetra::MultiVector::subView(ArrayView): "
      "range must include at least one vector.");

    if (isConstantStride ()) {
      return rcp (new MV (this->getMap (), view_, origView_, cols));
    }
    else {
      Array<size_t> newcols (cols.size ());
      for (size_t j = 0; j < static_cast<size_t> (cols.size ()); ++j) {
        newcols[j] = whichVectors_[cols[j]];
      }
      return rcp (new MV (this->getMap (), view_, origView_, newcols ()));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subView (const Teuchos::Range1D& colRng) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::Array;
    using Teuchos::rcp;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    TEUCHOS_TEST_FOR_EXCEPTION(
      colRng.size() == 0, std::runtime_error,
      "Tpetra::MultiVector::subView(Range1D): "
      "range must include at least one vector.");

    // resulting MultiVector is constant stride only if *this is
    if (isConstantStride ()) {
      // view goes from first entry of first vector to last entry of last vector
      std::pair<size_t, size_t> cols (colRng.lbound (), colRng.ubound () + 1);
      return rcp (new MV (this->getMap (),
                          subview<dual_view_type> (view_, ALL (), cols),
                          origView_));
    }
    else {
      // otherwise, use a subset of this whichVectors_ to construct new multivector
      Array<size_t> which (whichVectors_.begin () + colRng.lbound (),
                           whichVectors_.begin () + colRng.ubound () + 1);
      return rcp (new MV (this->getMap (), view_, origView_, which));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subViewNonConst (const ArrayView<const size_t> &cols)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    return Teuchos::rcp_const_cast<MV> (this->subView (cols));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  subViewNonConst (const Teuchos::Range1D &colRng)
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;
    return Teuchos::rcp_const_cast<MV> (this->subView (colRng));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getVector (size_t j) const
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::rcp;
    // using Teuchos::rcp_const_cast;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > V;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( vectorIndexOutOfRange(j), std::runtime_error,
        "Tpetra::MultiVector::getVector(NonConst): index j (== " << j << ") exceeds valid column range for this multivector.");
#endif // HAVE_TPETRA_DEBUG

    // FIXME (mfh 10 May 2014) Why can't Kokkos take size_t instead of
    // unsigned int?
    const unsigned int jj = isConstantStride () ?
      static_cast<unsigned int> (j) :
      static_cast<unsigned int> (whichVectors_[j]);
    // FIXME (mfh 10 May 2014) The same issue shows up here that shows
    // up for offsetView.
    return rcp (new V (this->getMap (),
                       subview<dual_view_type> (view_, ALL (), jj),
                       origView_));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getVectorNonConst (size_t j)
  {
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > V;
    return Teuchos::rcp_const_cast<V> (this->getVector (j));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get1dCopy (Teuchos::ArrayView<Scalar> A, size_t LDA) const
  {
    typedef typename dual_view_type::host_mirror_device_type host_mirror_device_type;
    // The user's array is column major ("LayoutLeft").
    typedef Kokkos::View<Scalar**, Kokkos::LayoutLeft,
      host_mirror_device_type, Kokkos::MemoryUnmanaged> input_view_type;
    typedef typename dual_view_type::t_host host_view_type;

    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();
    const std::pair<size_t, size_t> rowRange (0, numRows);
    const std::pair<size_t, size_t> colRange (0, numCols);

    const char tfecfFuncName[] = "get1dCopy";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      LDA < numRows, std::runtime_error,
      ": LDA = " << LDA << " < numRows = " << numRows << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numRows > static_cast<size_t> (0) &&
      numCols > static_cast<size_t> (0) &&
      static_cast<size_t> (A.size ()) < LDA * (numCols - 1) + numRows,
      std::runtime_error,
      ": A.size() = " << A.size () << ", but its size must be at least "
      << LDA * (numCols - 1) * numRows << " to hold all the entries.");

    // This does a deep copy into A, which is column major with
    // leading dimension LDA.  We assume A is big enough to hold
    // everything.  Copy directly from host mirror of the data, if it
    // exists.

    // Start by sync'ing to host.
    view_.template sync<host_mirror_device_type> ();

    input_view_type dstWholeView (A.getRawPtr (), LDA, numCols);
    input_view_type dstView =
      Kokkos::subview<input_view_type> (dstWholeView, rowRange, Kokkos::ALL ());
    host_view_type srcView =
      Kokkos::subview<host_view_type> (view_.h_view, rowRange, Kokkos::ALL ());

    if (this->isConstantStride ()) {
      Kokkos::deep_copy (dstView, srcView);
    }
    else {
      // FIXME (mfh 10 Jul 2014) Shouldn't we use size_t instead of int here?
      Kokkos::View<int*, host_mirror_device_type> whichVecsDst ("whichVecsDst", numCols);
      Kokkos::View<int*, host_mirror_device_type> whichVecsSrc ("whichVecsSrc", numCols);
      for (size_t j = 0; j < numCols; ++j) {
        whichVecsSrc(j) = static_cast<int> (this->whichVectors_[j]);
        whichVecsDst(j) = static_cast<int> (j);
      }
      DeepCopySelectedVectors<input_view_type, host_view_type,
        host_mirror_device_type, false, false> f (dstView, srcView, whichVecsDst, whichVecsSrc);
      Kokkos::parallel_for (numRows, f);
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get2dCopy (Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const
  {
    typedef typename dual_view_type::host_mirror_device_type
      host_mirror_device_type;
    typedef typename dual_view_type::t_host host_view_type;
    typedef Kokkos::View<Scalar**,
      typename host_view_type::array_layout,
      typename dual_view_type::host_mirror_device_type,
      Kokkos::MemoryUnmanaged> unmanaged_host_view_type;

    const char tfecfFuncName[] = "get2dCopy";
    const size_t numRows = this->getLocalLength ();
    const size_t numCols = this->getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (ArrayOfPtrs.size ()) != numCols,
      std::runtime_error, ": Input array of pointers must contain as many "
      "entries (arrays) as the MultiVector has columns.  ArrayOfPtrs.size() = "
      << ArrayOfPtrs.size () << " != getNumVectors() = " << numCols << ".");

    if (numRows != 0 && numCols != 0) {
      // Start by sync'ing to host.
      view_.template sync<host_mirror_device_type> ();

      // No side effects until we've validated the input.
      for (size_t j = 0; j < numCols; ++j) {
        const size_t dstLen = static_cast<size_t> (ArrayOfPtrs[j].size ());
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          dstLen < numRows, std::invalid_argument, ": Array j = " << j << " of "
          "the input array of arrays is not long enough to fit all entries in "
          "that column of the MultiVector.  ArrayOfPtrs[j].size() = " << dstLen
          << " < getLocalLength() = " << numRows << ".");
      }

      // We've validated the input, so it's safe to start copying.
      for (size_t j = 0; j < numCols; ++j) {
        const size_t col = isConstantStride () ? j : whichVectors_[j];
        host_view_type src =
          Kokkos::subview<host_view_type> (view_.h_view, Kokkos::ALL (), col);
        unmanaged_host_view_type dst (ArrayOfPtrs[j].getRawPtr (),
                                      ArrayOfPtrs[j].size (),
                                      static_cast<size_t> (1));
        Kokkos::deep_copy (dst, src);
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<const Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get1dView () const
  {
    if (getLocalLength () == 0 || getNumVectors () == 0) {
      return Teuchos::null;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isConstantStride (), std::runtime_error, "Tpetra::MultiVector::"
        "get1dView: This MultiVector does not have constant stride, so it is "
        "not possible to view its data as a single array.  You may check "
        "whether a MultiVector has constant stride by calling "
        "isConstantStride().");
      // NOTE (mfh 09 2014) get1dView() and get1dViewNonConst() have
      // always been device->host synchronization points.  We might
      // want to change this in the future.
      typedef typename dual_view_type::host_mirror_device_type host_type;
      view_.template sync<host_type> ();
      // Both get1dView() and get1dViewNonConst() return a host view
      // of the data.
      Teuchos::ArrayRCP<scalar_type> dataAsArcp =
        Kokkos::Compat::persistingView (view_.template view<host_type> ());
      return Teuchos::arcp_const_cast<const scalar_type> (dataAsArcp);

      // Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
      // return node->template viewBuffer<Scalar> (getStride () * (getNumVectors () - 1) +
      //                                           getLocalLength (),
      //                                           MVT::getValues (lclMV_));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Scalar>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get1dViewNonConst ()
  {
    if (getLocalLength () == 0 || getNumVectors () == 0) {
      return Teuchos::null;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isConstantStride (), std::runtime_error, "Tpetra::MultiVector::"
        "get1dViewNonConst: This MultiVector does not have constant stride, so "
        "it is not possible to view its data as a single array.  You may check "
        "whether a MultiVector has constant stride by calling "
        "isConstantStride().");
      // NOTE (mfh 09 May 2014) get1dView() and get1dViewNonConst()
      // have always been device->host synchronization points.  We
      // might want to change this in the future.
      typedef typename dual_view_type::host_mirror_device_type host_type;
      view_.template sync<host_type> ();
      // Both get1dView() and get1dViewNonConst() return a host view
      // of the data.
      Teuchos::ArrayRCP<scalar_type> dataAsArcp =
        Kokkos::Compat::persistingView (view_.template view<host_type> ());
      return dataAsArcp;

      // FIXME (mfh 10 May 2014) In the classic version of
      // MultiVector, once the reference count of the returned result
      // of get1dViewNonConst() reaches zero, any changes in the
      // result are synchronized to the (GPU) device.  We don't do
      // that here, but we could, but adding a custom destructor to
      // the Teuchos::ArrayRCP that synchs to the device.  On the
      // other hand, that would run counter to the semantics of the
      // Kokkos refactor version of Tpetra, where every Tpetra object
      // is supposed to behave like Kokkos::DualView: that is,
      // requiring the user to request synchronization between host
      // and device explicitly.

      // Teuchos::RCP<Node> node = MVT::getNode (lclMV_);
      // return node->template viewBufferNonConst<Scalar> (KokkosClassic::ReadWrite,
      //                                                   getStride () * (getNumVectors () - 1) +
      //                                                   getLocalLength (),
      //                                                   MVT::getValuesNonConst (lclMV_));
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get2dViewNonConst ()
  {
    const size_t numCols = getNumVectors ();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > views (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      const size_t col = isConstantStride () ? j : whichVectors_[j];
      dual_view_type X_col = Kokkos::subview<dual_view_type> (view_, Kokkos::ALL (), col);
      views[j] = Kokkos::Compat::persistingView (X_col.d_view);
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  get2dView () const
  {
    const size_t numCols = getNumVectors ();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > views (numCols);
    for (size_t j = 0; j < numCols; ++j) {
      const size_t col = isConstantStride () ? j : whichVectors_[j];
      dual_view_type X_col = Kokkos::subview<dual_view_type> (view_, Kokkos::ALL (), col);
      views[j] = Kokkos::Compat::persistingView (X_col.d_view);
    }
    return views;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  multiply (Teuchos::ETransp transA,
            Teuchos::ETransp transB,
            const Scalar &alpha,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& B,
            const Scalar &beta)
  {
    using Teuchos::NO_TRANS;      // enums
    using Teuchos::TRANS;
    using Teuchos::CONJ_TRANS;
    using Teuchos::null;
    using Teuchos::ScalarTraits;  // traits
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MV;

    // This routine performs a variety of matrix-matrix multiply operations, interpreting
    // the MultiVector (this-aka C , A and B) as 2D matrices.  Variations are due to
    // the fact that A, B and C can be local replicated or global distributed
    // MultiVectors and that we may or may not operate with the transpose of
    // A and B.  Possible cases are:
    //                                       Num
    //      OPERATIONS                        cases  Notes
    //  1) C(local) = A^X(local) * B^X(local)  4    (X=Trans or Not, No comm needed)
    //  2) C(local) = A^T(distr) * B  (distr)  1    (2D dot product, replicate C)
    //  3) C(distr) = A  (distr) * B^X(local)  2    (2D vector update, no comm needed)
    //
    // The following operations are not meaningful for 1D distributions:
    //
    // u1) C(local) = A^T(distr) * B^T(distr)  1
    // u2) C(local) = A  (distr) * B^X(distr)  2
    // u3) C(distr) = A^X(local) * B^X(local)  4
    // u4) C(distr) = A^X(local) * B^X(distr)  4
    // u5) C(distr) = A^T(distr) * B^X(local)  2
    // u6) C(local) = A^X(distr) * B^X(local)  4
    // u7) C(distr) = A^X(distr) * B^X(local)  4
    // u8) C(local) = A^X(local) * B^X(distr)  4
    //
    // Total of 32 case (2^5).

    const char errPrefix[] = "Tpetra::MultiVector::multiply(transOpA,transOpB,alpha,A,B,beta): ";

    TEUCHOS_TEST_FOR_EXCEPTION( ScalarTraits<Scalar>::isComplex && (transA == TRANS || transB == TRANS), std::invalid_argument,
        errPrefix << "non-conjugate transpose not supported for complex types.");
    transA = (transA == NO_TRANS ? NO_TRANS : CONJ_TRANS);
    transB = (transB == NO_TRANS ? NO_TRANS : CONJ_TRANS);

    // Compute effective dimensions, w.r.t. transpose operations on
    size_t A_nrows = (transA==CONJ_TRANS) ? A.getNumVectors() : A.getLocalLength();
    size_t A_ncols = (transA==CONJ_TRANS) ? A.getLocalLength() : A.getNumVectors();
    size_t B_nrows = (transB==CONJ_TRANS) ? B.getNumVectors() : B.getLocalLength();
    size_t B_ncols = (transB==CONJ_TRANS) ? B.getLocalLength() : B.getNumVectors();

    Scalar beta_local = beta; // local copy of beta; might be reassigned below

    TEUCHOS_TEST_FOR_EXCEPTION(
      getLocalLength() != A_nrows || getNumVectors() != B_ncols || A_ncols != B_nrows,
      std::runtime_error,
      errPrefix << "dimension of *this, op(A) and op(B) must be consistent.  "
      << std::endl << "The local part of *this is "
      << getLocalLength() << " x " << getNumVectors()
      << ", A is " << A_nrows << " x " << A_ncols
      << ", and B is " << B_nrows << " x " << B_ncols << ".");

    bool A_is_local = !A.isDistributed();
    bool B_is_local = !B.isDistributed();
    bool C_is_local = !this->isDistributed();
    bool Case1 = ( C_is_local &&  A_is_local &&  B_is_local);                                           // Case 1: C(local) = A^X(local) * B^X(local)
    bool Case2 = ( C_is_local && !A_is_local && !B_is_local && transA==CONJ_TRANS && transB==NO_TRANS); // Case 2: C(local) = A^T(distr) * B  (distr)
    bool Case3 = (!C_is_local && !A_is_local &&  B_is_local && transA==NO_TRANS  );                     // Case 3: C(distr) = A  (distr) * B^X(local)

    // Test that we are considering a meaningful cases
    TEUCHOS_TEST_FOR_EXCEPTION( !Case1 && !Case2 && !Case3, std::runtime_error,
        errPrefix << "multiplication of op(A) and op(B) into *this is not a supported use case.");

    if (beta != ScalarTraits<Scalar>::zero() && Case2)
    {
      // if Case2, then C is local and contributions must be summed across all nodes
      // however, if beta != 0, then accumulate beta*C into the sum
      // when summing across all nodes, we only want to accumulate this once, so
      // set beta == 0 on all nodes except node 0
      int MyPID = this->getMap()->getComm()->getRank();
      if (MyPID!=0) beta_local = ScalarTraits<Scalar>::zero();
    }

    // Check if A, B, C have constant stride, if not then make temp copy (strided)
    RCP<const MV> Atmp, Btmp;
    RCP<MV>       Ctmp;
    if (isConstantStride() == false) Ctmp = rcp (new MV (*this));
    else Ctmp = rcp(this,false);

    if (A.isConstantStride() == false) Atmp = rcp (new MV (A));
    else Atmp = rcp(&A,false);

    if (B.isConstantStride() == false) Btmp = rcp (new MV (B));
    else Btmp = rcp(&B,false);

#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(!Ctmp->isConstantStride() || !Btmp->isConstantStride() || !Atmp->isConstantStride(), std::logic_error,
        errPrefix << "failed making temporary strided copies of input multivectors.");
#endif

    KMV &C_mv = Ctmp->lclMV_;
    {
      // get local multivectors
      const KMV &A_mv = Atmp->lclMV_;
      const KMV &B_mv = Btmp->lclMV_;
      // do the multiply (GEMM)
      //
      // KR FIXME Need GEMM wrapper
      MVT::GEMM(C_mv,transA,transB,alpha,A_mv,B_mv,beta_local);
    }

    // Dispose of (possibly) extra copies of A, B
    Atmp = null;
    Btmp = null;

    RCP<Node> node = MVT::getNode(lclMV_);
    // If *this was not strided, copy the data from the strided version and then delete it
    if (! isConstantStride ()) {
      // *this is not strided, we must put data from Ctmp into *this
      TEUCHOS_TEST_FOR_EXCEPT(&C_mv != &lclMV_);
      const size_t numVecs = MVT::getNumCols(lclMV_);
      for (size_t j=0; j < numVecs; ++j) {
        node->template copyBuffers<Scalar>(getLocalLength(),MVT::getValues(C_mv,j),MVT::getValuesNonConst(lclMV_,whichVectors_[j]));
      }
    }

    // If Case 2 then sum up *this and distribute it to all processors.
    if (Case2) {
      this->reduce();
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  elementWiseMultiply (Scalar scalarAB,
                       const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& A,
                       const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& B,
                       Scalar scalarThis)
  {
    using Teuchos::arcp_const_cast;
    const char tfecfFuncName[] = "elementWiseMultiply()";

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getLocalLength() != A.getLocalLength() ||
      getLocalLength() != B.getLocalLength(), std::runtime_error,
      ": MultiVectors do not have the same local length.");
#endif // HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      B.getNumVectors() != this->getNumVectors(), std::runtime_error,
      ": MultiVectors 'this' and B must have the same number of vectors.");

    const size_t numVecs = getNumVectors();

    typedef Kokkos::View<Scalar*, DeviceType> view_type;

    if (isConstantStride() && A.isConstantStride()) {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      A.view_.template sync<DeviceType>();
      A.view_.template modify<DeviceType>();
      B.view_.template sync<DeviceType>();
      B.view_.template modify<DeviceType>();
      view_type vector_A = Kokkos::subview<view_type> (A.view_.d_view, Kokkos::ALL (), 0);
      Kokkos::MV_ElementWiseMultiply(scalarThis,view_.d_view,
                                     scalarAB,vector_A,B.view_.d_view);
    }
    else {
      view_.template sync<DeviceType>();
      view_.template modify<DeviceType>();
      A.view_.template sync<DeviceType>();
      A.view_.template modify<DeviceType>();
      B.view_.template sync<DeviceType>();
      B.view_.template modify<DeviceType>();
      view_type vector_A = Kokkos::subview<view_type> (A.view_.d_view, Kokkos::ALL (), 0);
      for (size_t k=0; k < numVecs; ++k) {
        const size_t this_col = isConstantStride () ? k : whichVectors_[k];
        view_type vector_k = Kokkos::subview<view_type> (view_.d_view, Kokkos::ALL (), this_col);
        const size_t B_col = isConstantStride () ? k : B.whichVectors_[k];
        view_type vector_Bk = Kokkos::subview<view_type> (B.view_.d_view, Kokkos::ALL (), B_col);
        Kokkos::V_ElementWiseMultiply(scalarThis,vector_k, scalarAB,vector_A,vector_Bk);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::reduce()
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_SUM;
    typedef typename dual_view_type::t_dev device_view_type;
    typedef typename device_type::host_mirror_device_type host_mirror_device_type;

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->isDistributed (), std::runtime_error,
      "Tpetra::MultiVector::reduce() should only be called with locally "
      "replicated or otherwise not distributed MultiVector objects.");
    const Teuchos::Comm<int>& comm = * (this->getMap ()->getComm ());
    if (comm.getSize () == 1) {
      return;
    }

    const size_t numLclRows = getLocalLength ();
    const size_t numCols = getNumVectors ();

    // FIXME (mfh 16 June 2014) This exception will cause deadlock if
    // it triggers on only some processes.  We don't have a good way
    // to pack this result into the all-reduce below, but this would
    // be a good reason to set a "local error flag" and find other
    // opportunities to let it propagate.
    TEUCHOS_TEST_FOR_EXCEPTION(
      numLclRows > static_cast<size_t> (std::numeric_limits<int>::max ()),
      std::runtime_error, "Tpetra::MultiVector::reduce: On Process " <<
      comm.getRank () << ", the number of local rows " << numLclRows <<
      " does not fit in int.");

    //
    // Use MPI to sum the entries across all local blocks.
    //
    // If this MultiVector's local data are stored contiguously, we
    // can use the local View as the source buffer in the
    // MPI_Allreduce.  Otherwise, we have to allocate a temporary
    // source buffer and pack.
    const bool contig = isConstantStride () && getStride () == numLclRows;
    device_view_type srcBuf;
    if (contig) {
      srcBuf = view_.d_view;
    }
    else {
      srcBuf = device_view_type ("srcBuf", numLclRows, numCols);
      Kokkos::deep_copy (srcBuf, view_.d_view);
    }

    // MPI requires that the send and receive buffers don't alias one
    // another, so we have to copy temporary storage for the result.
    //
    // We expect that MPI implementations will know how to read device
    // pointers.
    device_view_type tgtBuf ("tgtBuf", numLclRows, numCols);

    const int reduceCount = static_cast<int> (numLclRows * numCols);
    reduceAll<int, Scalar> (comm, REDUCE_SUM, reduceCount,
                            srcBuf.ptr_on_device (), tgtBuf.ptr_on_device ());

    // Tell the DualView that we plan to modify the device data.
    view_.template modify<device_type> ();

    const std::pair<size_t, size_t> lclRowRange (0, numLclRows);
    device_view_type d_view =
      subview<device_view_type> (view_.d_view, lclRowRange, ALL ());

    if (contig || isConstantStride ()) {
      Kokkos::deep_copy (d_view, tgtBuf);
    }
    else {
      for (size_t j = 0; j < numCols; ++j) {
        device_view_type d_view_j =
          subview<device_view_type> (d_view, ALL (), j);
        device_view_type tgtBuf_j =
          subview<device_view_type> (tgtBuf, ALL (), j);
        Kokkos::deep_copy (d_view_j, tgtBuf_j);
      }
    }

    // Synchronize the host with changes on the device.
    //
    // FIXME (mfh 16 June 2014) This raises the question of whether we
    // want to synchronize always.  Users will find it reassuring if
    // MultiVector methods always leave the MultiVector in a
    // synchronized state, but it seems silly to synchronize to host
    // if they hardly ever need host data.
    view_.template sync<host_mirror_device_type> ();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  replaceLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    const LocalOrdinal minLocalIndex = this->getMap()->getMinLocalIndex();
    const LocalOrdinal maxLocalIndex = this->getMap()->getMaxLocalIndex();
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow < minLocalIndex || MyRow > maxLocalIndex,
      std::runtime_error,
      "Tpetra::MultiVector::replaceLocalValue: row index " << MyRow
      << " is invalid.  The range of valid row indices on this process "
      << this->getMap()->getComm()->getRank() << " is [" << minLocalIndex
      << ", " << maxLocalIndex << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::replaceLocalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    const size_t colInd = isConstantStride () ?
      VectorIndex : whichVectors_[VectorIndex];
    view_.d_view (static_cast<size_t> (MyRow), colInd) = ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  sumIntoLocalValue (LocalOrdinal MyRow,
                     size_t VectorIndex,
                     const Scalar &ScalarValue)
  {
#ifdef HAVE_TPETRA_DEBUG
    const LocalOrdinal minLocalIndex = this->getMap()->getMinLocalIndex();
    const LocalOrdinal maxLocalIndex = this->getMap()->getMaxLocalIndex();
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow < minLocalIndex || MyRow > maxLocalIndex,
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoLocalValue: row index " << MyRow
      << " is invalid.  The range of valid row indices on this process "
      << this->getMap()->getComm()->getRank() << " is [" << minLocalIndex
      << ", " << maxLocalIndex << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoLocalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    const size_t colInd = isConstantStride () ?
      VectorIndex : whichVectors_[VectorIndex];
    view_.d_view (static_cast<size_t> (MyRow), colInd) += ScalarValue;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  replaceGlobalValue (GlobalOrdinal GlobalRow,
                      size_t VectorIndex,
                      const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap()->getLocalElement(GlobalRow);
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
      std::runtime_error,
      "Tpetra::MultiVector::replaceGlobalValue: global row index " << GlobalRow
      << "is not present on this process " << this->getMap()->getComm()->getRank()
      << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::replaceGlobalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    replaceLocalValue (MyRow, VectorIndex, ScalarValue);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  sumIntoGlobalValue (GlobalOrdinal GlobalRow,
                      size_t VectorIndex,
                      const Scalar &ScalarValue)
  {
    LocalOrdinal MyRow = this->getMap()->getLocalElement(GlobalRow);
#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      MyRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: global row index " << GlobalRow
      << "is not present on this process " << this->getMap()->getComm()->getRank()
      << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      vectorIndexOutOfRange(VectorIndex),
      std::runtime_error,
      "Tpetra::MultiVector::sumIntoGlobalValue: vector index " << VectorIndex
      << " of the multivector is invalid.");
#endif
    sumIntoLocalValue (MyRow, VectorIndex, ScalarValue);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class T>
  Teuchos::ArrayRCP<T>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  getSubArrayRCP (Teuchos::ArrayRCP<T> arr,
                  size_t j) const
  {
    const size_t col = isConstantStride () ? j : whichVectors_[j];
    dual_view_type X_col = Kokkos::subview<dual_view_type> (view_, Kokkos::ALL (), col);
    return Kokkos::Compat::persistingView (X_col.d_view); // ????? host or device???
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  KokkosClassic::MultiVector<Scalar,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getLocalMV () const
  {
    // KR Creates the KokkosClassic thing on the fly.
    // OK to leave this in place as backwards compat.

    KMV kmv (this->getMap ()->getNode ());
    size_t stride[8];
    view_.stride (stride);
    const size_t LDA = view_.dimension_1 () > 1 ? stride[1] : view_.dimension_0 ();
    MVT::initializeValues (kmv, getLocalLength (), getNumVectors (),
                           Kokkos::Compat::persistingView (view_.d_view),
                           LDA,
                           getOrigNumLocalRows (),
                           getOrigNumLocalCols ());
    return kmv;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  typename MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::dual_view_type
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getDualView() const {
    return view_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  TEUCHOS_DEPRECATED
  KokkosClassic::MultiVector<Scalar,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >&
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::getLocalMVNonConst() {
    return lclMV_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::description() const
  {
    using std::endl;
    std::ostringstream oss;
    oss << Teuchos::typeName (*this) << " {"
        << "label: \"" << this->getObjectLabel () << "\""
        << ", numRows: " << getGlobalLength ()
        << ", numCols: " << getNumVectors ()
        << ", isConstantStride: " << isConstantStride ();
    if (isConstantStride ()) {
      oss << ", columnStride: " << getStride ();
    }
    oss << "}";
    return oss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using std::endl;
    using std::setw;

    // Set default verbosity if applicable.
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    RCP<const Teuchos::Comm<int> > comm = this->getMap()->getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();

    if (vl != VERB_NONE) {
      // Don't set the tab level unless we're printing something.
      Teuchos::OSTab tab0 (out);

      if (myImageID == 0) { // >= VERB_LOW prints description()
        out << "Tpetra::MultiVector:" << endl;
        Teuchos::OSTab tab1 (out);
        out << "Template parameters:" << endl;
        {
          Teuchos::OSTab tab2 (out);
          out << "Scalar: " << Teuchos::TypeNameTraits<Scalar>::name () << endl
              << "LocalOrdinal: " << Teuchos::TypeNameTraits<LocalOrdinal>::name () << endl
              << "GlobalOrdinal: " << Teuchos::TypeNameTraits<GlobalOrdinal>::name () << endl
              << "Node: " << Teuchos::TypeNameTraits<Node>::name () << endl;
        }
        out << "label: \"" << this->getObjectLabel () << "\"" << endl
            << "numRows: " << getGlobalLength () << endl
            << "numCols: " << getNumVectors () << endl
            << "isConstantStride: " << isConstantStride () << endl;
        if (isConstantStride ()) {
          out << "columnStride: " << getStride () << endl;
        }
      }
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // At verbosity > VERB_LOW, each process prints something.
            out << "Process " << myImageID << ":" << endl;
            Teuchos::OSTab tab2 (out);

            // >= VERB_MEDIUM: print the local vector length.
            out << "local length: " << getLocalLength();
            if (vl != VERB_MEDIUM) {
              // >= VERB_HIGH: print isConstantStride() and getStride()
              if (isConstantStride()) {
                out << "constant stride: " << getStride() << endl;
              }
              else {
                out << "not constant stride" << endl;
              }
              if (vl == VERB_EXTREME) {
                // VERB_EXTREME: print all the values in the multivector.
                out << "values: " << endl;
                ArrayRCP<ArrayRCP<const Scalar> > X = this->get2dView();
                out << "[";
                for (size_t i = 0; i < getLocalLength(); ++i) {
                  for (size_t j = 0; j < getNumVectors(); ++j) {
                    out << X[j][i];
                    if (j + 1 < getNumVectors()) {
                      out << ", ";
                    }
                  } // for each column
                  if (i + 1 < getLocalLength ()) {
                    out << "; ";
                  } else {
                    out << endl;
                  }
                } // for each row
                out << "]" << endl;
              } // if vl == VERB_EXTREME
            } // if (vl != VERB_MEDIUM)
            else { // vl == VERB_LOW
              out << endl;
            }
          } // if vl != VERB_LOW
        } // if it is my process' turn to print
        comm->barrier();
      } // for each process in the communicator
    } // if vl != VERB_NONE
  }

#if TPETRA_USE_KOKKOS_DISTOBJECT
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViews() const
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViewsNonConst (KokkosClassic::ReadWriteOption rwo)
  {
    // Do nothing in Kokkos::View implementation
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::releaseViews () const
  {
    // Do nothing in Kokkos::View implementation
  }

#else // NOT TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViews() const
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  createViewsNonConst (KokkosClassic::ReadWriteOption /* rwo */ )
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  releaseViews () const
  {}

#endif // TPETRA_USE_KOKKOS_DISTOBJECT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
  {
    replaceMap (newMap);
  }

  template <class Scalar, class LO, class GO, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  createMultiVector (const Teuchos::RCP<const Map<LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >& map,
                     size_t numVectors)
  {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;
    typedef MultiVector<Scalar, LO, GO, node_type> MV;
    return Teuchos::rcp (new MV (map, numVectors));
  }

  /// \brief Nonmember MultiVector constructor with view semantics using user-allocated data.
  /// \relatesalso MultiVector
  /// \relatesalso Vector
  ///
  /// \warning This function is not supported for all Kokkos Node
  ///   types.  Specifically, it is not typically supported for
  ///   GPU accelerator-based nodes like KokkosClassic::ThrustGPUNode.
  ///
  /// \param map [in] The Map describing the distribution of rows of
  ///   the multivector.
  /// \param view [in/out] A pointer to column-major dense matrix
  ///   data.  This will be the multivector's data on the calling
  ///   process.  The multivector will use the pointer directly,
  ///   without copying.
  /// \param LDA [in] The leading dimension (a.k.a. "stride") of the
  ///   column-major input data.
  /// \param numVectors [in] The number of columns in the input data.
  ///   This will be the number of vectors in the returned
  ///   multivector.
  ///
  /// \node To Kokkos and Tpetra developers: If you add a new Kokkos
  ///   Node type that is a host Node type (where memory lives in user
  ///   space, not in a different space as on a GPU), you will need to
  ///   add a specialization of Tpetra::details::ViewAccepter for your
  ///   new Node type.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  createMultiVectorFromView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >& map,
                             const Teuchos::ArrayRCP<Scalar>& view,
                             const size_t LDA,
                             const size_t numVectors)
  {
    (void) map;
    (void) view;
    (void) LDA;
    (void) numVectors;

    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::createMultiVectorFromView: "
      "Not implemented for Node = KokkosDeviceWrapperNode.");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  createCopy (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >& src)
  {
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type> MV;

    MV cpy (src.getMap (), src.getNumVectors ());
    if (src.isConstantStride ()) {
      Kokkos::deep_copy (cpy.getDualView (), src.getDualView ());
    }
    else {
      const char viewName[] = "MultiVector::createCopy::WhichVectors";

      if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
        typedef typename DeviceType::size_type size_type;
        typedef DeepCopySelectedVectors<typename MV::dual_view_type::t_dev,
                                        typename MV::dual_view_type::t_dev,
                                        DeviceType, true, false> functor_type;

        const size_type numWhichVectors = src.whichVectors_.size ();
        Kokkos::View<int*, DeviceType> whichVectors (viewName, numWhichVectors);
        for (size_type i = 0; i < numWhichVectors; ++i) {
          whichVectors(i) = src.whichVectors_[i];
        }
        Kokkos::parallel_for (src.getLocalLength (),
                              functor_type (cpy.getDualView ().template view<DeviceType> (),
                                            src.getDualView ().template view<DeviceType> (),
                                            whichVectors, whichVectors));
      } else {
        typedef typename DeviceType::host_mirror_device_type host_mirror_device_type;
        typedef typename host_mirror_device_type::size_type size_type;
        typedef DeepCopySelectedVectors<typename MV::dual_view_type::t_host,
                                        typename MV::dual_view_type::t_host,
                                        host_mirror_device_type,
                                        true, false> functor_type;

        const size_type numWhichVectors = src.whichVectors_.size ();
        Kokkos::View<int*, host_mirror_device_type> whichVectors (viewName, numWhichVectors);
        for (size_type i = 0; i < numWhichVectors; ++i) {
          whichVectors(i) = src.whichVectors_[i];
        }
        Kokkos::parallel_for (src.getLocalLength (),
                              functor_type (cpy.getDualView ().template view<host_mirror_device_type> (),
                                            src.getDualView ().template view<host_mirror_device_type> (),
                                            whichVectors, whichVectors));
      }
    }
    return cpy;
  }

  template <class DS, class DL, class DG, class DD, class SS, class SL, class SG, class SD>
  void deep_copy( MultiVector<DS,DL,DG,Kokkos::Compat::KokkosDeviceWrapperNode<DD> >& dst,
                  const MultiVector<SS,SL,SG,Kokkos::Compat::KokkosDeviceWrapperNode<SD> >& src) {
    typedef MultiVector<DS,DL,DG,Kokkos::Compat::KokkosDeviceWrapperNode<DD> > MVD;
    typedef const MultiVector<SS,SL,SG,Kokkos::Compat::KokkosDeviceWrapperNode<SD> > MVS;
    if(src.isConstantStride() && dst.isConstantStride()) {
      Kokkos::deep_copy(dst.getDualView(),src.getDualView());
    }
    else {
      if(dst.isConstantStride()) {
        if(src.getDualView().modified_device>=src.getDualView().modified_host) {
          Kokkos::View<int*,DD> whichVectors("MultiVector::createCopy::WhichVectors",src.whichVectors_.size());
          for(int i = 0; i < src.whichVectors_.size(); i++)
            whichVectors(i)=src.whichVectors_[i];
          Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::dual_view_type::t_dev,typename MVS::dual_view_type::t_dev,DD,true,false>
                                                    (dst.getDualView().template view<DD>(),
                                                     src.getDualView().template view<DD>(),
                                                     whichVectors,whichVectors));
        } else {
          Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectors("MultiVector::createCopy::WhichVectors",src.whichVectors_.size());
          for(int i = 0; i < src.whichVectors_.size(); i++)
            whichVectors(i)=src.whichVectors_[i];
          Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::dual_view_type::t_host,typename MVS::dual_view_type::t_host,typename DD::host_mirror_device_type,true,false>
                                                    (dst.getDualView().template view<typename DD::host_mirror_device_type>(),
                                                     src.getDualView().template view<typename DD::host_mirror_device_type>(),
                                                     whichVectors,whichVectors));
        }
      } else {
        if(src.isConstantStride()) {
          if(src.getDualView().modified_device>=src.getDualView().modified_host) {
            Kokkos::View<int*,DD> whichVectors("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectors(i)=dst.whichVectors_[i];
            Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::dual_view_type::t_dev,typename MVS::dual_view_type::t_dev,DD,false,true>
                                                      (dst.getDualView().template view<DD>(),
                                                       src.getDualView().template view<DD>(),
                                                       whichVectors,whichVectors));
          } else {
            Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectors("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectors(i)=dst.whichVectors_[i];
            Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::dual_view_type::t_host,typename MVS::dual_view_type::t_host,typename DD::host_mirror_device_type,false,true>
                                                      (dst.getDualView().template view<typename DD::host_mirror_device_type>(),
                                                       src.getDualView().template view<typename DD::host_mirror_device_type>(),
                                                       whichVectors,whichVectors));
          }
        } else {
          if(src.getDualView().modified_device>=src.getDualView().modified_host) {
            Kokkos::View<int*,DD> whichVectorsDst("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsDst(i)=dst.whichVectors_[i];
            Kokkos::View<int*,DD> whichVectorsSrc("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsSrc(i)=src.whichVectors_[i];
            Kokkos::parallel_for(src.getLocalLength(),DeepCopySelectedVectors<typename MVD::dual_view_type::t_dev,typename MVS::dual_view_type::t_dev,DD,false,false>
                                                      (dst.getDualView().template view<DD>(),
                                                       src.getDualView().template view<DD>(),
                                                       whichVectorsDst,whichVectorsSrc));
          } else {
            Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectorsDst("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsDst(i)=dst.whichVectors_[i];
            Kokkos::View<int*,typename DD::host_mirror_device_type> whichVectorsSrc("MultiVector::createCopy::WhichVectors",dst.whichVectors_.size());
            for(int i = 0; i < dst.whichVectors_.size(); i++)
              whichVectorsSrc(i) = src.whichVectors_[i];
            Kokkos::parallel_for (src.getLocalLength(),
                                  DeepCopySelectedVectors<typename MVD::dual_view_type::t_host,
                                                          typename MVS::dual_view_type::t_host,
                                                          typename DD::host_mirror_device_type, false, false>
                                    (dst.getDualView ().template view<typename DD::host_mirror_device_type> (),
                                     src.getDualView ().template view<typename DD::host_mirror_device_type> (),
                                     whichVectorsDst, whichVectorsSrc));
          }
        }
      }
    }
  }
} // namespace Tpetra


#endif // TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DEF_HPP
