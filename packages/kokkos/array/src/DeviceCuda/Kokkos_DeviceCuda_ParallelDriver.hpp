/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_DEVICECUDA_PARALLELDRIVER_HPP
#define KOKKOS_DEVICECUDA_PARALLELDRIVER_HPP

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

struct DeviceCudaTraits {
  enum { WarpSize       = 32      /* 0x0020 */ };
  enum { WarpIndexMask  = 0x001f  /* Mask for warpindex */ };
  enum { WarpIndexShift = 5       /* WarpSize == 1 << WarpShift */ };
  enum { SharedMemoryBanks_13 = 16 /* Compute device 1.3 */ };
  enum { SharedMemoryBanks_20 = 32 /* Compute device 2.0 */ };
 

  enum { ConstantMemoryCapacity = 0x010000 /* 64k bytes */ };
  enum { ConstantMemoryCache    = 0x002000 /*  8k bytes */ };

  typedef unsigned long
    ConstantGlobalBufferType[ ConstantMemoryCapacity / sizeof(unsigned long) ];
};

}
}

//----------------------------------------------------------------------------

#include <impl/Kokkos_DeviceCuda_macros.hpp>

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

__device__ __constant__
Kokkos::Impl::DeviceCudaTraits::ConstantGlobalBufferType
kokkos_device_cuda_constant_memory_buffer ;

template< class ParallelDriver >
__global__
void kokkos_device_cuda_parallel_driver()
{
  ((const ParallelDriver *) kokkos_device_cuda_constant_memory_buffer )->run_on_device();
}

namespace Kokkos {
namespace Impl {

template< class ParallelDriver >
void device_cuda_run( const ParallelDriver & driver ,
                      const dim3           & block ,
                      const dim3           & grid ,
                      const DeviceCuda::size_type shmem = 0 )
{
  cudaMemcpyToSymbol( kokkos_device_cuda_constant_memory_buffer , & driver , sizeof(ParallelDriver) );

  kokkos_device_cuda_parallel_driver< ParallelDriver ><<< block , grid , shmem >>>();
}

}
}

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

#include <impl/Kokkos_DeviceClear_macros.hpp>

#endif /* KOKKOS_DEVICECUDA_PARALLELDRIVER_HPP */

