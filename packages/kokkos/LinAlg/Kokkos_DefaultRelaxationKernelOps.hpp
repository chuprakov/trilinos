//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

#ifndef KOKKOS_DEFAULTRELAXATION_KERNELOPS_HPP
#define KOKKOS_DEFAULTRELAXATION_KERNELOPS_HPP

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

#ifdef __CUDACC__
#include <Teuchos_ScalarTraitsCUDA.hpp>
#else
#include <Teuchos_ScalarTraits.hpp>
#endif

// NTS: Need to remove Scalar divisions and replace with appropriate inverse listed in Traits

namespace Kokkos {

  // Extract Matrix Diagonal for Type 1 storage
  template <class Scalar, class Ordinal>
  struct ExtractDiagonalOp1 {

    // mat data
    size_t numRows;
    const size_t  *begs;
    const size_t  *ends;
    const Ordinal *inds;
    const Scalar  *vals;
    Scalar * diag;


    inline KERNEL_PREFIX void execute(size_t row) {
      for (size_t c=begs[row]; c<ends[row]; ++c) {
        if(row==(size_t)inds[c]) {
          diag[row]=vals[c];
          break;
        }
      }
    }
  };

  // Extract Matrix Diagonal for Type 2 storage
  template <class Scalar, class Ordinal>
  struct ExtractDiagonalOp2 {

    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    Scalar * diag;
    size_t numRows;

    inline KERNEL_PREFIX void execute(size_t row) {
      const Scalar  *curval = vals_beg[row];
      const Ordinal *curind = inds_beg[row];
      for (size_t j=0;j<numEntries[row];j++) {
        if(row==(size_t)curind[j]){
          diag[row]=curval[j];
          break;
        }
      }
    }
  };


  /************************************************************************************/
  /********************************* Jacobi Kernels ***********************************/
  /************************************************************************************/
  // Jacobi for Type 1 storage.
  template <class Scalar, class Ordinal>
  struct DefaultJacobiOp1 {
    size_t numRows;
    const size_t  *begs;
    const size_t  *ends;
    const Ordinal *inds;
    const Scalar  *vals;
    const Scalar  *diag;
    // vector data (including multiple rhs)
    Scalar       *x;
    const Scalar *x0;
    const Scalar *b;
    Scalar damping_factor;
    size_t xstride, bstride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row  = i % numRows;
      const size_t rhs  = (i - row) / numRows;
      Scalar       *xj  = x + rhs * xstride;
      const Scalar *x0j = x0 + rhs * xstride;
      const Scalar *bj  = b + rhs * bstride;

      Scalar tmp = bj[row];
      for (size_t c=begs[row]; c<ends[row]; ++c) {
        tmp -= vals[c] * x0j[inds[c]];
      }
      xj[row]=x0j[row]+damping_factor*tmp/diag[row];
    }
  };


  // Jacobi for Type 2 storage.
  template <class Scalar, class Ordinal>
  struct DefaultJacobiOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    const Scalar  *diag;
    size_t numRows;
    // vector data (including multiple rhs)    
    Scalar        *x;
    const Scalar *x0;
    const Scalar  *b;
    Scalar damping_factor;
    size_t xstride, bstride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      const size_t rhs = (i - row) / numRows;
      Scalar       *xj = x + rhs * xstride;
      const Scalar *x0j = x0 + rhs * xstride;
      const Scalar *bj = b + rhs * bstride;
      Scalar tmp = bj[row];
      const Scalar  *curval = vals_beg[row];
      const Ordinal *curind = inds_beg[row];
      for (size_t j=0; j != numEntries[row]; ++j) {
        tmp -= (curval[j]) * x0j[curind[j]];
      }
      xj[row]=x0j[row]+damping_factor*tmp/diag[row];
    }
  };

  /************************************************************************************/
  /************************ Fine-Grain Gauss-Seidel Kernels ***************************/
  /************************************************************************************/

  // Fine-grain "hybrid" Gauss-Seidel for Type 1 storage.
  // Note: This is actually real Gauss-Seidel for a serial node, and hybrid for almost any other kind of node.
  template <class Scalar, class Ordinal>
  struct DefaultFineGrainHybridGaussSeidelOp1 {
    const size_t  *offsets;
    const Ordinal *inds;
    const Scalar  *vals;
    const Scalar  *diag;
    size_t numRows;
    // vector data (including multiple rhs)
    Scalar       *x;
    const Scalar *b;
    Scalar damping_factor;
    size_t xstride, bstride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      const size_t rhs = (i - row) / numRows;
      Scalar       *xj = x + rhs * xstride;
      const Scalar *bj = b + rhs * bstride;
      Scalar tmp = bj[row];
      for (size_t c=offsets[row];c<offsets[row+1];c++) {
        tmp -= vals[c] * xj[inds[c]];
      }
      xj[row]+=damping_factor*tmp/diag[row];
    }
  };


  // Fine-grain "hybrid" Gauss-Seidel for Type 2 storage.
  // Note: This is actually real Gauss-Seidel for a serial node, and hybrid for almost any other kind of node.
  template <class Scalar, class Ordinal>
  struct DefaultFineGrainHybridGaussSeidelOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    const Scalar  *diag;
    size_t numRows;
    // vector data (including multiple rhs)    
    Scalar        *x;
    const Scalar  *b;
    Scalar damping_factor;
    size_t xstride, bstride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      const size_t rhs = (i - row) / numRows;
      Scalar       *xj = x + rhs * xstride;
      const Scalar *bj = b + rhs * bstride;
      Scalar tmp = bj[row];
      const Scalar  *curval = vals_beg[row];
      const Ordinal *curind = inds_beg[row];
      for (size_t j=0; j != numEntries[row]; ++j) {
        tmp -= (curval[j]) * xj[curind[j]];
      }
      xj[row]+=damping_factor*tmp/diag[row];
    }
  };


  /************************************************************************************/
  /************************ Coarse-Grain Gauss-Seidel Kernels *************************/
  /************************************************************************************/

  // Coarse-grain "hybrid" Gauss-Seidel for Type 1 storage.
  // Note: This is actually real Gauss-Seidel for a serial node, and hybrid for almost any other kind of node.
  template <class Scalar, class Ordinal>
  struct DefaultCoarseGrainHybridGaussSeidelOp1 {
    const size_t  *offsets;
    const Ordinal *inds;
    const Scalar  *vals;
    const Scalar  *diag;
    size_t numRows;
    size_t numChunks;
    // vector data (including multiple rhs)
    Scalar       *x;
    const Scalar *b;
    Scalar damping_factor;
    size_t xstride, bstride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t chunk = i % numChunks;
      const size_t rhs = (i - chunk) / numChunks;
      const size_t start_r = chunk * numRows / numChunks;
      const size_t stop_r  = TEUCHOS_MIN((chunk+1)*numRows/numChunks,numRows);
      Scalar       *xj = x + rhs * xstride;
      const Scalar *bj = b + rhs * bstride;
      for (size_t row=start_r;row<stop_r;row++){
        Scalar tmp = bj[row];
        for (size_t c=offsets[row];c<offsets[row+1];c++) {
          tmp -= vals[c] * xj[inds[c]];
        }
        xj[row]+=damping_factor*tmp/diag[row];
      }
    }
  };


  // Coarse-grain "hybrid" Gauss-Seidel for Type 2 storage.
  // Note: This is actually real Gauss-Seidel for a serial node, and hybrid for almost any other kind of node.
  template <class Scalar, class Ordinal>
  struct DefaultCoarseGrainHybridGaussSeidelOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    const Scalar  *diag;
    size_t numRows;
    size_t numChunks;
    // vector data (including multiple rhs)    
    Scalar        *x;
    const Scalar  *b;
    Scalar damping_factor;
    size_t xstride, bstride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t chunk = i % numChunks;
      const size_t rhs = (i - chunk) / numChunks;
      const size_t start_r = chunk * numRows / numChunks;
      const size_t stop_r  = TEUCHOS_MIN((chunk+1)*numRows/numChunks,numRows);
      Scalar       *xj = x + rhs * xstride;
      const Scalar *bj = b + rhs * bstride;
      for (size_t row=start_r;row<stop_r;row++){
        Scalar tmp = bj[row];
        const Scalar  *curval = vals_beg[row];
        const Ordinal *curind = inds_beg[row];
        for (size_t j=0; j!=numEntries[row];j++) {
          tmp -= (curval[j]) * xj[curind[j]];
        }
        xj[row]+=damping_factor*tmp/diag[row];
      }
    }
  };
 
  /************************************************************************************/
  /******************************** Chebyshev Kernels *********************************/
  /************************************************************************************/

  template <class Scalar, class Ordinal>
  struct DefaultChebyshevOp1 {
    const size_t  *offsets;
    const Ordinal *inds;
    const Scalar  *vals;
    const Scalar  *diag;
    size_t numRows;
    // vector data (including multiple rhs)
    Scalar       *x,*w;
    const Scalar *x0,*b;
    Scalar oneOverTheta,dtemp1,dtemp2;
    size_t stride;
    bool first_step;
    bool zero_initial_guess;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row  = i % numRows;
      const size_t rhs  = (i - row) / numRows;
      Scalar       *xj  = x + rhs * stride;
      const Scalar *x0j = x0 + rhs * stride;
      Scalar       *wj  = w + rhs * stride;      
      const Scalar *bj  = b + rhs * stride;
      Scalar        vj;

      if(first_step){
	if(zero_initial_guess)
	  // x= theta^{-1} D^{-1} b
	  xj[row]=wj[row]=bj[row] / diag[row] *oneOverTheta;
	else{
	  // v=Ax
	  vj=Teuchos::ScalarTraits<Scalar>::zero();
	  for (size_t c=offsets[row];c<offsets[row+1];c++) {
	    vj += vals[c] * x0j[inds[c]];
	  }
	  // w=theta^{-1} D^{-1} (b -Ax)
	  wj[row]=(bj[row]-vj)/diag[row]*oneOverTheta;
	  // x+=w
	  xj[row]+=wj[row];
	}
      }
      else{
	//v=Ax
	vj=Teuchos::ScalarTraits<Scalar>::zero();
	for (size_t c=offsets[row];c<offsets[row+1];c++) {
	  vj += vals[c] * x0j[inds[c]];
	}
	// w=dtemp1*w +  D^{-1}*dtemp2*(b-Ax)
	wj[row]=dtemp1*wj[row]+dtemp2*(bj[row]-vj)/diag[row];
	// x+=w
	xj[row]+=wj[row];
      }

      //      printf("[%3d-%d] x=%11.4e v=%11.4e w=%11.4e x0=%11.4e\n",row,first_step,xj[row],vj,wj[row],x0j[row]);

    }

  };

}// namespace Kokkos

#endif /* KOKKOS_DEFAULTRELAXATION_KERNELOPS_HPP */
