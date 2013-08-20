// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_Storage_i_h)
#define Intrepid_MiniTensor_Storage_i_h

namespace Intrepid {

//
// Raw pointer storage
//
template<typename T>
inline
Storage<T, DYNAMIC>::Storage() :
  storage_(NULL),
  size_(0)
{
  return;
}

template<typename T>
inline
Storage<T, DYNAMIC>::Storage(Index const number_entries) :
  storage_(NULL),
  size_(0)
{
  resize(number_entries);
  return;
}

template<typename T>
inline
Storage<T, DYNAMIC>::~Storage()
{
  clear();
  return;
}

template<typename T>
inline
T const &
Storage<T, DYNAMIC>::operator[](Index const i) const
{
  assert(i < size());
  return storage_[i];
}

template<typename T>
inline
T &
Storage<T, DYNAMIC>::operator[](Index const i)
{
  assert(i < size());
  return storage_[i];
}

template<typename T>
inline
Index
Storage<T, DYNAMIC>::size() const
{
  return size_;
}

template<typename T>
inline
void
Storage<T, DYNAMIC>::resize(Index const number_entries)
{
  if (number_entries == size()) return;

  clear();
  storage_ = new T[number_entries];
  size_ = number_entries;

  return;
}

template<typename T>
inline
void
Storage<T, DYNAMIC>::clear()
{
  if (storage_ == NULL) return;

  delete [] storage_;
  storage_ = NULL;
  size_ = 0;

  return;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Storage_i_h
