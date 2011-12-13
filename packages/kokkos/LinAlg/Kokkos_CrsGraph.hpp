//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
//@HEADER

#ifndef KOKKOS_CRSGRAPH_HPP
#define KOKKOS_CRSGRAPH_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CompileTimeAssert.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

  //=========================================================================================================================
  // 
  // A host-resident CrsGraph
  // 
  //=========================================================================================================================

  /** \brief A default host-compute compressed-row sparse graph.
      \ingroup kokkos_crs_ops
   */
  template <class Ordinal, 
            class Node,
            class LocalMatOps>
  class CrsGraphHostCompute {
  public:

    typedef Ordinal               OrdinalType;
    typedef Node                  NodeType;
    typedef LocalMatOps           LocalMatOpsType;

    //! @name Constructors/Destructor
    //@{

    //! Default CrsGraphHostCompute constuctor.
    CrsGraphHostCompute(size_t numRows, const RCP<Node> &node);

    //! CrsGraphHostCompute Destructor
    virtual ~CrsGraphHostCompute();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    RCP<Node> getNode() const;

    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Return the number of rows in the graph.
    size_t getNumRows() const;

    //! Return the number of entries in the graph.
    size_t getNumEntries() const;

    //! Indicates that the graph is filled, but empty.
    bool isEmpty() const;

    //! Indicatest that the graph has been finalized.
    bool isFinalized() const;

    //! \brief Indicate that the structure is 1D.
    //! It will never be the case that both is1DStructure() and is2DStructure() return true.
    bool is1DStructure() const;

    //! \brief Indicate that the structure is 2D.
    //! It will never be the case that both is1DStructure() and is2DStructure() return true.
    bool is2DStructure() const;

    //! \brief Indicate that the stucture is optimized.
    bool isOptimized() const;

    //! Submit the indices and offset for 1D storage.
    /** 
          \post is1DStructure() == true
     */
    void set1DStructure(ArrayRCP<Ordinal> inds, 
                        ArrayRCP<size_t>  rowBegs,
                        ArrayRCP<size_t>  rowEnds);
                        
    //! Submit the indices for 2D storage.
    /** 
          \post is2DStructure() == true
     */
    void set2DStructure(ArrayRCP<ArrayRCP<Ordinal> > inds,
                        ArrayRCP<size_t>                      numEntriesPerRow);

    //! Retrieve the structure for 1D storage.
    /** 
          If is1DStructure() == false, then 
          \post inds == rowBegs == rowEnds == null

          Otherwise, 
          \post indices for row \c r are inds[r], where \f$r \in [b,e)\f$, where \f$b = rowBegs[r]\f$ and \f$e = rowEnds[r]\f$
          \post rowBegs has getNumRows()+1 entries; the last entry is inds.size()
          \post rowEnds has getNumRows() entries
     */
    void get1DStructure(ArrayRCP<Ordinal> &inds, 
                        ArrayRCP<size_t>  &rowBegs,
                        ArrayRCP<size_t>  &rowEnds);

    //! Retrieve the structure for 2D storage.
    /** 
          If is2DStructure() == false, then 
          \post inds == numEntriesPerRow == null
     */
    void get2DStructure(ArrayRCP<ArrayRCP<Ordinal> > &inds,
                        ArrayRCP<size_t>                      &numEntriesPerRow);

    //! Instruct the graph to perform any necessary manipulation, including (optionally) optimizing the storage of the graph data.
    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          \post if OptimizeStorage == true, then is2DStructure() == true
     */
    void finalize(bool OptimizeStorage);

    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          @param[in/out] values2D      2D-structured matrix values. Required to be non-null if is2DStructure() is true. Set to null if OptimizeStorage is true.
          @param[in/out] values1D      1D-structured matrix values. Required to be non-null if is1DStructure() is true. Allocated if OptimizeStorage is true.
          \post if OptimizeStorage == true or already is2DStructure(), then is2DStructure() == true.
     */
    template <class Scalar>
    void finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D);

    //! Release data associated with this graph.
    virtual void clear();

    //@}

  protected:
    //! Copy constructor (protected and not implemented) 
    CrsGraphHostCompute(const CrsGraphHostCompute& sources);

    RCP<Node> node_;
    size_t numRows_, numEntries_;
    bool isFinalized_, isEmpty_, is1D_, is2D_, isOpt_;

    // 2D storage
    ArrayRCP<ArrayRCP<Ordinal> >  indices2D_;
    ArrayRCP<size_t>              numEntriesPerRow_;
    // 1D storage
    ArrayRCP<Ordinal>             indices1D_;
    ArrayRCP<size_t>              rowBegs_, rowEnds_;
  };


  //==============================================================================
  template <class Ordinal, class Node, class LocalMatOps>
  CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::CrsGraphHostCompute(size_t numRows, const RCP<Node> &node) 
  : node_(node)
  , numRows_(numRows)
  {
    CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::clear();
  }

  //==============================================================================
  template <class Ordinal, class Node, class LocalMatOps>
  CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::~CrsGraphHostCompute() {
  }

  // ======= clear ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::clear() {
    isFinalized_   = false;
    isEmpty_       = false;
    is1D_          = false;
    is2D_          = false;
    isOpt_         = false;
    numEntries_    = 0;
    indices2D_        = null;
    numEntriesPerRow_ = null;
    rowBegs_          = null;
    rowEnds_          = null;
    indices1D_        = null;
  }

  // ======= node ===========
  template <class Ordinal, class Node, class LocalMatOps>
  RCP<Node> CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::getNode() const {
    return node_;
  }

  // ======= numrows ===========
  template <class Ordinal, class Node, class LocalMatOps>
  size_t CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::getNumRows() const {
    return numRows_;
  }

  // ======= numentries ===========
  template <class Ordinal, class Node, class LocalMatOps>
  size_t CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::getNumEntries() const {
    return numEntries_;
  }

  // ======= isempty ===========
  template <class Ordinal, class Node, class LocalMatOps>
  bool CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::isEmpty() const {
    return isEmpty_;
  }

  // ======= isfinalized ===========
  template <class Ordinal, class Node, class LocalMatOps>
  bool CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::isFinalized() const {
    return isFinalized_;
  }

  // ======= is1d ===========
  template <class Ordinal, class Node, class LocalMatOps>
  bool CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::is1DStructure() const {
    return is1D_;
  }

  // ======= is2d ===========
  template <class Ordinal, class Node, class LocalMatOps>
  bool CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::is2DStructure() const {
    return is2D_;
  }

  // ======= isopt ===========
  template <class Ordinal, class Node, class LocalMatOps>
  bool CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::isOptimized() const {
    return isOpt_;
  }

  // ======= get 1d ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::get1DStructure(ArrayRCP<Ordinal> &inds, 
                                                                     ArrayRCP<size_t>  &rowBegs,
                                                                     ArrayRCP<size_t>  &rowEnds)
  {
    inds = indices1D_;
    rowBegs = rowBegs_;
    rowEnds = rowEnds_;
  }

  // ======= get 2d ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::get2DStructure(ArrayRCP<ArrayRCP<Ordinal> > &inds,
                                                                     ArrayRCP<size_t>                      &numEntriesPerRow) 
  {
    inds = indices2D_;
    numEntriesPerRow = numEntriesPerRow_;
  }

  // ======= set 1d ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::set1DStructure(ArrayRCP<Ordinal> inds, 
                                                                     ArrayRCP<size_t>  rowBegs,
                                                                     ArrayRCP<size_t>  rowEnds)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)rowBegs.size() != numRows_+1 || (size_t)rowEnds.size() != numRows_, std::runtime_error, 
        Teuchos::typeName(*this) << "::set1DStructure(inds,rowBegs,rowEnds): rowBegs and rowEnds are not correctly sized.");
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)rowBegs[numRows_] > (size_t)inds.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::set1DStructure(inds,rowBegs,rowEnds): rowBegs contents to not agree with inds size.");
    this->clear();
    //
    indices1D_ = inds;
    rowBegs_ = rowBegs;
    rowEnds_ = rowEnds;
    if (numRows_ > 0) {
      for (size_t i=0; i < this->getNumRows(); ++i) {
        numEntries_ += (this->rowEnds_[i] - this->rowBegs_[i]);
#ifdef HAVE_KOKKOS_DEBUG
        // row i goes like [ begs[i] , ends[i] )
        // sanity        : begs[i] <= ends[i]
        // ordering      : begs[i] <= begs[i+1]
        // no overlapping: ends[i] <= begs[i+1]
        TEUCHOS_TEST_FOR_EXCEPTION( rowBegs_[i+1] < rowBegs_[i] || rowEnds_[i] < rowBegs_[i] || rowEnds_[i] > rowBegs_[i+1], std::runtime_error,
            Teuchos::typeName(*this) << "::set1DStructure(inds,rowBegs,rowEnds): ends and begs are not consistent.");
#endif
      }
    }
    is1D_ = true;
  }

  // ======= set 2d ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::set2DStructure(ArrayRCP<ArrayRCP<Ordinal> > inds,
                                                                     ArrayRCP<size_t>                      numEntriesPerRow)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)inds.size() != numRows_ || (size_t)numEntriesPerRow.size() != numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::set2DStructure(inds,numEntriesPerRow): numEntriesPerRow and inds must have as many entries as the number of rows specified to the constructor.");
    this->clear();
    //
    indices2D_  = inds;
    if (indices2D_ != null) {
      numEntriesPerRow_ = numEntriesPerRow;
      numEntries_ = std::accumulate(this->numEntriesPerRow_.begin(), this->numEntriesPerRow_.end(), 0);
#ifdef HAVE_KOKKOS_DEBUG
      for (size_t i=0; i<numRows_; ++i) {
        TEUCHOS_TEST_FOR_EXCEPTION( (size_t)inds[i].size() < numEntriesPerRow[i], std::runtime_error,
            Teuchos::typeName(*this) << "::set2DStructure(): inds[" << i << "] == " << inds[i] 
            << " is not large enough for the specified number of entries, "
            << " numEntriesPerRow[" << i << "] == " << numEntriesPerRow[i]);
      }
#endif
    }
    is2D_ = true;
  }

  // ======= finalize ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage)
  {
    // allocations not done using the Node. no current need for host-based nodes, and 
    // this leads to incorrect behavior when we try to reuse this code from child CrsGraphDeviceCompute
    if (isFinalized() && !(OptimizeStorage == true && isOptimized() == false)) return;
    if ((indices1D_ == null && indices2D_ == null) || (this->getNumEntries() == 0)) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
      if (OptimizeStorage) {
        // move into packed 1D storage
        if (is1DStructure() == false) {
          // allocate 1D storage
          // we these are for host use, so we'll forgo the view
          indices1D_ = arcp<Ordinal>(this->getNumEntries());
        }
        ArrayRCP<size_t> offsets = arcp<size_t>(numRows_+1);
        // copy/pack data
        size_t curoffset = 0;
        size_t curnuminds;
        typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
        newinds = indices1D_.begin();
        for (size_t i=0; i < numRows_; ++i) {
          offsets[i] = curoffset;
          if (is1DStructure()) {
            curnuminds = rowEnds_[i] - rowBegs_[i];
            oldinds = indices1D_.begin() + rowBegs_[i];
          }
          else {
            curnuminds = numEntriesPerRow_[i];
            oldinds = indices2D_[i].begin();
          }
          std::copy(oldinds, oldinds+curnuminds, newinds);
          newinds += curnuminds;
          curoffset += curnuminds;
        }
        offsets[numRows_] = curoffset;
        TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
            Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        // done with the original row beg/end offsets, can point to the new overlapping one
        rowBegs_   = offsets;
        rowEnds_   = offsets.persistingView(1,numRows_);
        isOpt_     = true;
        is1D_      = true;
        // delete 2D storage (if there was any)
        is2D_      = false;
        numEntriesPerRow_ = null;
        indices2D_        = null;
      }
    }
    isFinalized_ = true;
  }


  // ======= finalize ===========
  // finalize() storage for the graph with associated matrix values
  // this is called from a CrsMatrix, and we're doing the finalize the for the graph and matrix at the same time, so the matrix doesn't have to.
  template <class Ordinal, class Node, class LocalMatOps>
  template <class Scalar>
  void CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D)
  {
    if (isFinalized() && !(OptimizeStorage == true && isOptimized() == false)) return;
    if ((indices1D_ == null && indices2D_ == null) || (this->getNumEntries() == 0)) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
      // move into packed 1D storage
      if (OptimizeStorage) {
        if (is1DStructure() == false) {
          // allocate 1D storage
          // we know this is a host-base node, so we'll forgo the view of rowBegs_,rowEnds_
          indices1D_ = arcp<Ordinal>(this->getNumEntries());
          values1D   = arcp<Scalar >(this->getNumEntries());
        }
        ArrayRCP<size_t> offsets = arcp<size_t>(numRows_+1);
        // copy/pack data
        size_t curoffset = 0;
        size_t curnuminds;
        typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
        typename ArrayRCP<Scalar >::iterator oldvals, newvals;
        newinds = indices1D_.begin();
        newvals = values1D.begin();
        for (size_t i=0; i < numRows_; ++i) {
          offsets[i] = curoffset;
          if (is1DStructure()) {
            curnuminds = rowEnds_[i] - rowBegs_[i];
            oldinds = indices1D_.begin() + rowBegs_[i];
            oldvals = values1D.begin() + rowBegs_[i];
          }
          else {
            curnuminds = numEntriesPerRow_[i];
            oldinds = indices2D_[i].begin();
            oldvals = values2D[i].begin();
          }
          std::copy(oldinds, oldinds+curnuminds, newinds);
          std::copy(oldvals, oldvals+curnuminds, newvals);
          newinds += curnuminds;
          newvals += curnuminds;
          curoffset += curnuminds;
        }
        offsets[numRows_] = curoffset;
        TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
            Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        // done with the original row beg/end offsets, can point to the new overlapping one
        rowBegs_   = offsets;
        rowEnds_   = offsets.persistingView(1,numRows_);
        is1D_      = true;
        isOpt_     = true;
        // delete 2D storage (if there was any)
        is2D_      = false;
        numEntriesPerRow_ = null;
        indices2D_        = null;
        values2D          = null;
      }
    }
    isFinalized_ = true;
  }


  //=========================================================================================================================
  // 
  // A device-resident CrsGraph
  // 
  //=========================================================================================================================


  /** \brief A default device-compute compressed-row sparse graph.
      \ingroup kokkos_crs_ops

      This is externally identical to the host-based graph; in fact, it
      derives from CrsGraphHostCompute. The difference is that that it
      contains additional storage and logic for device-bound compute buffers, and it over-rides finalize() to fill these.
   */
  template <class Ordinal, 
            class Node,
            class LocalMatOps>
  class CrsGraphDeviceCompute : public CrsGraphHostCompute<Ordinal,Node,LocalMatOps> {
  public:

    //! @name Constructors/Destructor
    //@{

    //! Default CrsGraphDeviceCompute constuctor.
    CrsGraphDeviceCompute(size_t numRows, const RCP<Node> &node);

    //! CrsGraphDeviceCompute Destructor
    ~CrsGraphDeviceCompute();

    //@}

    //! @name Methods over-riding CrsGraphDeviceCompute.
    //@{

    //! Instruct the graph to perform any necessary manipulation, including (optionally) optimizing the storage of the graph data.
    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          \post if OptimizeStorage == true, then is2DStructure() == true
     */
    void finalize(bool OptimizeStorage);

    //! Instruct the graph to perform any necessary manipulation, including (optionally) optimizing the storage of the graph data, performing identical transformation on matrix values as well.
    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          @param[in/out] values2D      2D-structured matrix values. Required to be non-null if is2DStructure() is true. Set to null if OptimizeStorage is true.
          @param[in/out] values1D      1D-structured matrix values. Required to be non-null if is1DStructure() is true. Allocated if OptimizeStorage is true.
          @param[out]    d_values1D    1D-structured matrix values, resident on the device. Allocated and filled, regardless of OptimizeStorage.
          \post if OptimizeStorage == true or already is2DStructure(), then is2DStructure() == true.
     */
    template <class Scalar>
    void finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D, ArrayRCP<Scalar> &d_values1D);

    //! Return the device-bound buffers.
    void getDeviceBuffers(ArrayRCP<Ordinal> &d_inds, ArrayRCP<size_t> &d_offs) const;

    //! Release data associated with this graph.
    virtual void clear();

    //@}

  protected:
    //! Copy constructor (protected and not implemented) 
    CrsGraphDeviceCompute(const CrsGraphDeviceCompute& sources);

    // device storage (always 1D packed)
    ArrayRCP<Ordinal> pbuf_indices_;
    ArrayRCP<size_t > pbuf_offsets_;
  };

  //==============================================================================
  template <class Ordinal, class Node, class LocalMatOps>
  CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps>::CrsGraphDeviceCompute(size_t numRows, const RCP<Node> &node) 
  : CrsGraphHostCompute<Ordinal,Node,LocalMatOps>(numRows,node) 
  {}

  //===== destructor =====
  template <class Ordinal, class Node, class LocalMatOps>
  CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps>::~CrsGraphDeviceCompute() 
  {}

  //===== clear =====
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps>::clear() { 
    CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::clear();
    pbuf_indices_ = null;
    pbuf_offsets_ = null;
  }

  //==============================================================================
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage)
  {
    if (this->isFinalized() && !(OptimizeStorage == true && this->isOptimized() == false)) return;
    // call "normal" finalize(). this handles the re-structuring of data. below, we will handle the movement to device.
    CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::finalize(OptimizeStorage);
    // all we're doing here now is copying data to device.
    // copy into a 1D structure on the device, regardless of host format
    if (this->isEmpty()) {
      pbuf_indices_ = null;
      pbuf_offsets_ = null;
    }
    else {
      // allocate space on the device and copy data there, in a packed format
      pbuf_offsets_ = this->getNode()->template allocBuffer<size_t>(this->getNumRows()+1);
      pbuf_indices_ = this->getNode()->template allocBuffer<Ordinal>(this->getNumEntries());
      if (this->isOptimized()) {
        // should be packed now; single copy should do, and offsets are rowBegs_
        this->getNode()->template copyToBuffer<size_t >(this->getNumRows()+1, this->rowBegs_(),   pbuf_offsets_);
        this->getNode()->template copyToBuffer<Ordinal>(this->getNumEntries(),this->indices1D_(), pbuf_indices_);
      }
      else {
        ArrayRCP<size_t > view_offsets = this->getNode()->template viewBufferNonConst<size_t >(WriteOnly,pbuf_offsets_.size(),pbuf_offsets_);
        ArrayRCP<Ordinal> view_indices = this->getNode()->template viewBufferNonConst<Ordinal>(WriteOnly,pbuf_indices_.size(),pbuf_indices_);
        typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
        newinds = view_indices.begin();
        size_t curnuminds, curoffset = 0;
        for (size_t i=0; i < this->getNumRows(); ++i) {
          view_offsets[i] = curoffset;
          if (this->is1DStructure()) {
            curnuminds = this->rowEnds_[i] - this->rowBegs_[i];
            oldinds = this->indices1D_.begin() + this->rowBegs_[i];
          }
          else {
            curnuminds = this->numEntriesPerRow_[i];
            oldinds = this->indices2D_[i].begin();
          }
          std::copy(oldinds, oldinds+curnuminds, newinds);
          newinds += curnuminds;
          curoffset += curnuminds;
        }
        view_offsets[this->getNumRows()] = curoffset;
        TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
            Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        view_offsets = null;
        view_indices = null;
      }
    }
  }

  // ======= get device ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps>::getDeviceBuffers(ArrayRCP<Ordinal> &d_inds, ArrayRCP<size_t> &d_offs) const
  {
    d_inds = pbuf_indices_;
    d_offs = pbuf_offsets_;
  }


  //==============================================================================
  template <class Ordinal, class Node, class LocalMatOps>
  template <class Scalar>
  void CrsGraphDeviceCompute<Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &h_vals2D, ArrayRCP<Scalar> &h_vals1D, ArrayRCP<Scalar> &d_valsPacked) 
  {
    if (this->isFinalized() && !(OptimizeStorage == true && this->isOptimized() == false)) return;
    // call "normal" finalize(). this handles the re-structuring of data. below, we will handle the movement to device.
    CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::finalize(OptimizeStorage,h_vals2D,h_vals1D);
    if (this->isEmpty()) {
      pbuf_indices_ = null;
      pbuf_offsets_ = null;
      d_valsPacked  = null;
    }
    else {
      // allocate space on the device and copy data there, in a packed format
      pbuf_offsets_ = this->getNode()->template allocBuffer<size_t>(this->getNumRows()+1);
      pbuf_indices_ = this->getNode()->template allocBuffer<Ordinal>(this->getNumEntries());
      d_valsPacked  = this->getNode()->template allocBuffer<Scalar >(this->getNumEntries());
      if (this->isOptimized()) {
        // should be packed now; single copy should do, and offsets are rowBegs_
        this->getNode()->template copyToBuffer<size_t >(this->getNumRows()+1 ,  this->rowBegs_(), pbuf_offsets_);
        this->getNode()->template copyToBuffer<Ordinal>(this->getNumEntries(),this->indices1D_(), pbuf_indices_);
        this->getNode()->template copyToBuffer<Scalar >(this->getNumEntries(),        h_vals1D(), d_valsPacked);
      }
      else {
        ArrayRCP<size_t > view_offsets = this->getNode()->template viewBufferNonConst<size_t >(WriteOnly,pbuf_offsets_.size(),pbuf_offsets_);
        ArrayRCP<Ordinal> view_indices = this->getNode()->template viewBufferNonConst<Ordinal>(WriteOnly,pbuf_indices_.size(),pbuf_indices_);
        ArrayRCP<Scalar >  view_values = this->getNode()->template viewBufferNonConst<Scalar >(WriteOnly, d_valsPacked.size(),d_valsPacked);
        typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
        typename ArrayRCP<Scalar >::iterator oldvals, newvals;
        newinds = view_indices.begin();
        newvals = view_values.begin();
        size_t curnuminds, curoffset = 0;
        for (size_t i=0; i < this->getNumRows(); ++i) {
          view_offsets[i] = curoffset;
          if (this->is1DStructure()) {
            curnuminds = this->rowEnds_[i] - this->rowBegs_[i];
            oldinds = this->indices1D_.begin() + this->rowBegs_[i];
            oldvals = h_vals1D.begin() + this->rowBegs_[i];
          }
          else {
            curnuminds = this->numEntriesPerRow_[i];
            oldinds = this->indices2D_[i].begin();
            oldvals = h_vals2D[i].begin();
          }
          std::copy(oldinds, oldinds+curnuminds, newinds);
          std::copy(oldvals, oldvals+curnuminds, newvals);
          newinds += curnuminds;
          newvals += curnuminds;
          curoffset += curnuminds;
        }
        view_offsets[this->getNumRows()] = curoffset;
        TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
            Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        view_offsets = null;
        view_indices = null;
        view_values  = null;
      }
    }
  }



  //=========================================================================================================================
  // 
  // A first-touch allocation host-resident CrsGraph
  // 
  //=========================================================================================================================

  /** \brief A host-compute compressed-row sparse graph with first-touch allocation.
      \ingroup kokkos_crs_ops
   */
  template <class Ordinal, 
            class Node,
            class LocalMatOps>
  class FirstTouchHostCrsGraph : public CrsGraphHostCompute<Ordinal,Node,LocalMatOps> {
  public:

    typedef typename CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::OrdinalType      OrdinalType;
    typedef typename CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::NodeType         NodeType;
    typedef typename CrsGraphHostCompute<Ordinal,Node,LocalMatOps>::LocalMatOpsType  LocalMatOpsType;

    //! @name Constructors/Destructor
    //@{

    //! Default FirstTouchHostCrsGraph constuctor.
    FirstTouchHostCrsGraph(size_t numRows, const RCP<Node> &node);

    //! FirstTouchHostCrsGraph Destructor
    virtual ~FirstTouchHostCrsGraph();

    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Instruct the graph to perform any necessary manipulation, including (optionally) optimizing the storage of the graph data.
    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          \post if OptimizeStorage == true, then is2DStructure() == true
     */
    void finalize(bool OptimizeStorage);

    /** 
          @param[in] OptimizeStorage   Permit the graph to reallocate storage on the host in order to provide optimal storage and/or performance.
          @param[in/out] values2D      2D-structured matrix values. Required to be non-null if is2DStructure() is true. Set to null if OptimizeStorage is true.
          @param[in/out] values1D      1D-structured matrix values. Required to be non-null if is1DStructure() is true. Allocated if OptimizeStorage is true.
          \post if OptimizeStorage == true or already is2DStructure(), then is2DStructure() == true.
     */
    template <class Scalar>
    void finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D);

    //@}

  protected:
    //! Copy constructor (protected and not implemented) 
    FirstTouchHostCrsGraph(const FirstTouchHostCrsGraph& sources);
  };

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif
  /// \class FirstTouchCopyIndicesKernel
  ///
  /// Kokkos kernel for copying array indices using a first-touch
  /// initialization strategy for CPU memory.
  ///
  /// \note We have to store member data as raw pointers, rather than
  /// ArrayRCPs, because ArrayRCP are not thread safe, and the arrays
  /// get accessed inside a Kokkos kernel.
  template <class T>
  struct FirstTouchCopyIndicesKernel {
    const size_t * numEntriesPerRow;
    const size_t * offsets1D;
    T * data1D;
    const ArrayRCP<T> * data2D;

    inline KERNEL_PREFIX void execute(size_t row) {
      const size_t rowNumInds = numEntriesPerRow[row];
      const T* const oldinds = data2D[row].getRawPtr();
      T* const newinds = data1D + offsets1D[row];
      std::copy(oldinds, oldinds+rowNumInds, newinds);
    }
  };

  //==============================================================================
  template <class Ordinal, class Node, class LocalMatOps>
  FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps>::FirstTouchHostCrsGraph(size_t numRows, const RCP<Node> &node) 
  : CrsGraphHostCompute<Ordinal,Node,LocalMatOps>(numRows,node) 
  {
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  //==============================================================================
  template <class Ordinal, class Node, class LocalMatOps>
  FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps>::~FirstTouchHostCrsGraph() {}

  // ======= finalize ===========
  template <class Ordinal, class Node, class LocalMatOps>
  void FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage)
  {
    // allocations not done using the Node. no current need for host-based nodes, and 
    // this leads to incorrect behavior when we try to reuse this code from child CrsGraphDeviceCompute
    if (this->isFinalized() && !(OptimizeStorage == true && this->isOptimized() == false)) return;
    if ((this->indices1D_ == null && this->indices2D_ == null) || (this->getNumEntries() == 0)) {
      this->isEmpty_ = true;
    }
    else {
      this->isEmpty_ = false;
      if (OptimizeStorage) {
        // move into packed 1D storage
        ArrayRCP<size_t> offsets = arcp<size_t>(this->numRows_+1);
        if (this->is2DStructure() == true) {
          // 2D to 1D packed: first-touch allocation
          // allocate 1D storage
          this->indices1D_ = arcp<Ordinal>(this->getNumEntries());
          // compute offset array on host thread
          {
            size_t curoffset = 0;
            for (size_t i=0; i < this->numRows_; ++i) {
              offsets[i] = curoffset;
              curoffset += this->numEntriesPerRow_[i];
            }
            offsets[this->numRows_] = curoffset;
            TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
                Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
          }
          FirstTouchCopyIndicesKernel<Ordinal> kern;
          kern.offsets1D = offsets.getRawPtr();
          kern.numEntriesPerRow = this->numEntriesPerRow_.getRawPtr();
          kern.data1D = this->indices1D_.getRawPtr();
          kern.data2D = this->indices2D_.getRawPtr();
          this->getNode()->template parallel_for<FirstTouchCopyIndicesKernel<Ordinal> >(0,this->numRows_,kern);
        }
        else {
          // 1D to 1D packed: no first-touch
          // copy/pack data
          size_t curoffset = 0;
          size_t curnuminds;
          typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
          newinds = this->indices1D_.begin();
          for (size_t i=0; i < this->numRows_; ++i) {
            offsets[i] = curoffset;
            curnuminds = this->rowEnds_[i] - this->rowBegs_[i];
            oldinds = this->indices1D_.begin() + this->rowBegs_[i];
            std::copy(oldinds, oldinds+curnuminds, newinds);
            newinds += curnuminds;
            curoffset += curnuminds;
          }
          offsets[this->numRows_] = curoffset;
          TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
              Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
        }
        // done with the original row beg/end offsets, can point to the new overlapping one
        this->rowBegs_   = offsets;
        this->rowEnds_   = offsets.persistingView(1,this->numRows_);
        this->isOpt_     = true;
        this->is1D_      = true;
        // delete 2D storage (if any)
        this->is2D_      = false;
        this->numEntriesPerRow_ = null;
        this->indices2D_        = null;
      }
    }
    this->isFinalized_ = true;
  }


  // ======= finalize ===========
  // finalize() storage for the graph with associated matrix values
  // this is called from a CrsMatrix, and we're doing the finalize the for the graph and matrix at the same time, so the matrix doesn't have to.
  template <class Ordinal, class Node, class LocalMatOps>
  template <class Scalar>
  void FirstTouchHostCrsGraph<Ordinal,Node,LocalMatOps>::finalize(bool OptimizeStorage, ArrayRCP<ArrayRCP<Scalar> > &values2D, ArrayRCP<Scalar> &values1D)
  {
    if (this->isFinalized() && !(OptimizeStorage == true && this->isOptimized() == false)) return;
    if ((this->indices1D_ == null && this->indices2D_ == null) || (this->getNumEntries() == 0)) {
      this->isEmpty_ = true;
    }
    else {
      this->isEmpty_ = false;
      // move into packed 1D storage
      if (OptimizeStorage) {
        ArrayRCP<size_t> offsets = arcp<size_t>(this->numRows_+1);
        if (this->is2DStructure() == true) {
          // 2D to 1D packed: first-touch allocation
          // allocate 1D storage
          this->indices1D_ = arcp<Ordinal>(this->getNumEntries());
          values1D         = arcp<Scalar >(this->getNumEntries());
          // compute offset array on host thread
          {
            size_t curoffset = 0;
            for (size_t i=0; i < this->numRows_; ++i) {
              offsets[i] = curoffset;
              curoffset += this->numEntriesPerRow_[i];
            }
            offsets[this->numRows_] = curoffset;
            TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
                Teuchos::typeName(*this) << "::finalize(): Internal logic error. Please contact Kokkos team.");
          }
          {
            FirstTouchCopyIndicesKernel<Ordinal> indskern;
            indskern.offsets1D = offsets.getRawPtr();
            indskern.numEntriesPerRow = this->numEntriesPerRow_.getRawPtr();
            indskern.data1D = this->indices1D_.getRawPtr();
            indskern.data2D = this->indices2D_.getRawPtr();
            this->getNode()->template parallel_for<FirstTouchCopyIndicesKernel<Ordinal> >(0,this->numRows_,indskern);
          }
          {
            FirstTouchCopyIndicesKernel<Scalar> valskern;
            valskern.offsets1D = offsets.getRawPtr();
            valskern.numEntriesPerRow = this->numEntriesPerRow_.getRawPtr();
            valskern.data1D = values1D.getRawPtr();
            valskern.data2D = values2D.getRawPtr();
            this->getNode()->template parallel_for<FirstTouchCopyIndicesKernel<Scalar> >(0,this->numRows_,valskern);
          }
        }
        else {
          // copy/pack data
          // 1D to 1D packed: no first-touch
          size_t curoffset = 0;
          size_t curnuminds;
          typename ArrayRCP<Ordinal>::iterator oldinds, newinds;
          typename ArrayRCP<Scalar >::iterator oldvals, newvals;
          newinds = this->indices1D_.begin();
          newvals = values1D.begin();
          for (size_t i=0; i < this->numRows_; ++i) {
            offsets[i] = curoffset;
            curnuminds = this->rowEnds_[i] - this->rowBegs_[i];
            oldinds = this->indices1D_.begin() + this->rowBegs_[i];
            oldvals = values1D.begin() + this->rowBegs_[i];
            std::copy(oldinds, oldinds+curnuminds, newinds);
            std::copy(oldvals, oldvals+curnuminds, newvals);
            newinds += curnuminds;
            newvals += curnuminds;
            curoffset += curnuminds;
          }
          offsets[this->numRows_] = curoffset;
          TEUCHOS_TEST_FOR_EXCEPTION( curoffset != this->getNumEntries(), std::logic_error, 
			      Teuchos::typeName(*this) << "::finalize(): "
			      "Internal logic error: curoffset (= " 
			      << curoffset << ") != this->getNumEntries() (= "
			      << this->getNumEntries() 
			      << ").  Please contact Kokkos team.");
        }
        // done with the original row beg/end offsets, can point to the new overlapping one
        this->rowBegs_   = offsets;
        this->rowEnds_   = offsets.persistingView(1,this->numRows_);
        this->is1D_      = true;
        this->isOpt_     = true;
        // delete 2D storage (if there was any)
        this->is2D_      = false;
        this->numEntriesPerRow_ = null;
        this->indices2D_        = null;
        values2D          = null;
      }
    }
    this->isFinalized_ = true;
  }



  //=========================================================================================================================
  // 
  // Specializations
  // 
  //=========================================================================================================================

  /** \brief Kokkos compressed-row sparse graph class.
      \ingroup kokkos_crs_ops
      
      Default specialization is a host-bound CrsGraphHostCompute object.
    */
  template <class Ordinal, 
            class Node,
            class LocalMatOps>
  class CrsGraph : public CrsGraphHostCompute<Ordinal,Node,LocalMatOps> {
  public:
    CrsGraph(size_t numRows, const RCP<Node> &node) : CrsGraphHostCompute<Ordinal,Node,LocalMatOps>(numRows,node) {}
  private:
    CrsGraph(const CrsGraph<Ordinal,Node,LocalMatOps> &graph); // not implemented
  };

#ifndef HAVE_KOKKOS_NO_FIRST_TOUCH_MATVEC_ALLOCATION
#ifdef HAVE_KOKKOS_TBB
  /** \brief Kokkos compressed-row sparse graph class.
      \ingroup kokkos_crs_ops
      
      Specialization is a first-touch host-bound FirstTouchHostCrsGraph object.
    */
  class TBBNode;
  template <class Ordinal, 
            class LocalMatOps>
  class CrsGraph<Ordinal,TBBNode,LocalMatOps> : public FirstTouchHostCrsGraph<Ordinal,TBBNode,LocalMatOps> {
  public:
    CrsGraph(size_t numRows, const RCP<TBBNode> &node) : FirstTouchHostCrsGraph<Ordinal,TBBNode,LocalMatOps>(numRows,node) {}
  private:
    CrsGraph(const CrsGraph<Ordinal,TBBNode,LocalMatOps> &graph); // not implemented
  };
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
  class TPINode;
  /** \brief Kokkos compressed-row sparse graph class.
      \ingroup kokkos_crs_ops
      
      Specialization is a first-touch host-bound FirstTouchHostCrsGraph object.
    */
  template <class Ordinal, 
            class LocalMatOps>
  class CrsGraph<Ordinal,TPINode,LocalMatOps> : public FirstTouchHostCrsGraph<Ordinal,TPINode,LocalMatOps> {
  public:
    CrsGraph(size_t numRows, const RCP<TPINode> &node) : FirstTouchHostCrsGraph<Ordinal,TPINode,LocalMatOps>(numRows,node) {}
  private:
    CrsGraph(const CrsGraph<Ordinal,TPINode,LocalMatOps> &graph); // not implemented
  };
#endif
#endif

  /** \brief Kokkos compressed-row sparse graph class.
      \ingroup kokkos_crs_ops

      Specialization for device-based graph operation is a CrsGraphDeviceCompute.
    */
  template <class S, class O, class N> class DefaultDeviceSparseOps;
  template <class S,
            class Ordinal,
            class Node>
  class CrsGraph<Ordinal,Node,DefaultDeviceSparseOps<S,Ordinal,Node> > : public CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<S,Ordinal,Node> > {
  public:
    CrsGraph(size_t numRows, const RCP<Node> &node) : CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<S,Ordinal,Node> >(numRows,node) {}
  private:
    CrsGraph(const CrsGraph<Ordinal,Node,DefaultDeviceSparseOps<S,Ordinal,Node> > &graph); // not implemented
  };

} // namespace Kokkos

#endif /* KOKKOS_CRSGRAPH_HPP */
