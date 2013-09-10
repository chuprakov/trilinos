/*
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
*/

#ifndef TPETRA_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_DECL_HPP
#define TPETRA_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#if TPETRA_USE_KOKKOS_DISTOBJECT

#include "Kokkos_View.hpp"
#include "Kokkos_Parallel.hpp"
#include "Kokkos_Atomic.hpp"

// Don't include Teuchos_ScalarTraits.hpp here because we want a different
// version for CPU versus GPU

namespace Tpetra {
namespace Details {

  // Functors for implementing packAndPrepare and unpackAndCombine
  // through parallel_for

  template <typename Scalar, typename LocalOrdinal, typename Device>
  struct PackArraySingleColumnConstantStride {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> exportLIDs;
    Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged> src;
    Kokkos::View<Scalar*, device_type> exports;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      exports[k] = src[exportLIDs[k]];
    }

    void pack();
  };

  template <typename Scalar, typename LocalOrdinal, typename Device>
  struct PackArraySingleColumnOffset {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> exportLIDs;
    Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged> src;
    Kokkos::View<Scalar*, device_type> exports;
    size_t offset;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      exports[k] = src[exportLIDs[k] + offset];
    }

    void pack();
  };

  template <typename Scalar, typename LocalOrdinal, typename Device>
  struct PackArrayMultiColumnConstantStride {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> exportLIDs;
    Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged> src;
    Kokkos::View<Scalar*, device_type> exports;
    size_t stride, numCols;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const size_t localRow = exportLIDs[k]; // should use as<size_t>()
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        exports[offset + j] = src[localRow + j*stride];
    }

    void pack();
  };

  template <typename Scalar, typename LocalOrdinal, typename Device>
  struct PackArrayMultiColumnVariableStride {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> exportLIDs;
    Kokkos::View<const size_t*, device_type> srcWhichVectors;
    Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged> src;
    Kokkos::View<Scalar*, device_type> exports;
    size_t stride, numCols;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const size_t localRow = exportLIDs[k]; // should use as<size_t>()
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        exports[offset + j] = src[localRow + srcWhichVectors[j]*stride];
    }

    void pack();
  };

  struct InsertOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      //dest = src;
      Kokkos::atomic_exchange(&dest, src);
    }
  };
  struct AddOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      //dest += src;
      Kokkos::atomic_fetch_add(&dest, src);
    }
  };
  struct AbsMaxOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar max(const Scalar& a, const Scalar& b) const {
      return a > b ? a : b ;
    }

    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      typedef Teuchos::ScalarTraits<Scalar> SCT;
      //dest = max( SCT::magnitude(dest), SCT::magnitude(src) );
      Kokkos::atomic_exchange(
        &dest, max( SCT::magnitude(dest), SCT::magnitude(src) ) );
    }
  };

  template <typename Scalar, typename LocalOrdinal, typename Op, typename Device>
  struct UnpackArrayMultiColumnConstantStride {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> importLIDs;
    Kokkos::View<const Scalar*, device_type> imports;
    Kokkos::View<Scalar*, device_type, Kokkos::MemoryUnmanaged> dest;
    size_t stride, numCols;
    Op op;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const size_t localRow = importLIDs[k]; // should use as<size_t>()
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        op(dest[localRow + j*stride], imports[offset + j]);
    }

    void unpack();
  };

  template <typename Scalar, typename LocalOrdinal, typename Op, typename Device>
  struct UnpackArrayMultiColumnVariableStride {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> importLIDs;
    Kokkos::View<const size_t*, device_type> whichVectors;
    Kokkos::View<const Scalar*, device_type> imports;
    Kokkos::View<Scalar*, device_type, Kokkos::MemoryUnmanaged> dest;
    size_t stride, numCols;
    Op op;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const size_t localRow = importLIDs[k]; // should use as<size_t>()
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        op(dest[localRow + whichVectors[j]*stride], imports[offset + j]);
    }

    void unpack();
  };

  template <typename Scalar, typename LocalOrdinal, typename Device>
  struct PermuteArrayMultiColumnConstantStride {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> permuteToLIDs;
    Kokkos::View<const LocalOrdinal*, device_type> permuteFromLIDs;
    Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged> src;
    Kokkos::View<Scalar*, device_type, Kokkos::MemoryUnmanaged> dest;
    size_t src_stride, dest_stride, numCols;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const size_t toRow = permuteToLIDs[k]; // should use as<size_t>()
      const size_t fromRow = permuteFromLIDs[k]; // should use as<size_t>()
      for (size_t j = 0; j < numCols; ++j)
        dest[toRow + j*dest_stride] = src[fromRow + j*src_stride];
    }

    void permute();
  };

  template <typename Scalar, typename LocalOrdinal, typename Device>
  struct PermuteArrayMultiColumnVariableStride {
    typedef Device device_type;
    typedef typename device_type::size_type size_type;

    Kokkos::View<const LocalOrdinal*, device_type> permuteToLIDs;
    Kokkos::View<const LocalOrdinal*, device_type> permuteFromLIDs;
    Kokkos::View<const size_t*, device_type> src_whichVectors;
    Kokkos::View<const size_t*, device_type> dest_whichVectors;
    Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged> src;
    Kokkos::View<Scalar*, device_type, Kokkos::MemoryUnmanaged> dest;
    size_t src_stride, dest_stride, numCols;

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const size_t toRow = permuteToLIDs[k]; // should use as<size_t>()
      const size_t fromRow = permuteFromLIDs[k]; // should use as<size_t>()
      for (size_t j = 0; j < numCols; ++j)
        dest[toRow + dest_whichVectors[j]*dest_stride] =
          src[fromRow + src_whichVectors[j]*src_stride];
    }

    void permute();
  };

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_USE_KOKKOS_DISTOBJECT

#endif // TPETRA_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_DECL_HPP
