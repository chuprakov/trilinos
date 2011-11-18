// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_AlltoAll.hpp

    \brief AlltoAll communication methods
*/

#ifndef _ZOLTAN2_ALLTOALL_HPP_
#define _ZOLTAN2_ALLTOALL_HPP_

/*! \file Zoltan2_AlltoAll.hpp
*/

#include <vector>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Standards.hpp>

//
// TODO: doxygen comments and error handling and timing
//

namespace Zoltan2
{

/*! \brief AlltoAll sends/receives a fixed number of objects to/from all processes.
 *
 * The data type T of the objects must be a type for which 
 * Teuchos::SerializationTraits are defined.  This is most likely every 
 * fundamental data type plus std::pair<T1,T2>. It does not
 * include std::vector<T2>.
 */

template <typename T, typename LNO>
void AlltoAll(const Comm<int> &comm,
              Zoltan2::Environment &env,
              const ArrayView<const T> &sendBuf,  // input
              LNO count,                          // input
              ArrayRCP<T> &recvBuf)         // output - allocated here
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (count == 0) return;   // count is the same on all procs

  LNO n = nprocs * count;
  T *ptr = new T [n]; 
  Z2_GLOBAL_MEMORY_ASSERTION(env, n, ptr);
  recvBuf = Teuchos::arcp<T>(ptr, 0, n, true);

  // Do self messages

  for (LNO i=0, offset = rank*count; i < count; i++, offset++){
    recvBuf.get()[offset] = sendBuf.getRawPtr()[offset];
  }

#ifdef HAVE_MPI
  // Post receives

  RCP<CommRequest> r;
  Array<RCP<CommRequest> > req(nprocs-1);

  for (int p=0; p < nprocs; p++){
    if (p != rank){
      ArrayRCP<T> recvBufPtr(recvBuf.get() + p*count, 0, count, false);
      try{
        r  = Teuchos::ireceive<int, T>(comm, recvBufPtr, p);
      }
      catch (const std::exception &e){
        Z2_THROW_OUTSIDE_ERROR(env, e);
      }
    
      req.push_back(r);
    }
  }

  // Wait until all are posted

  Teuchos::barrier(comm);

  // Do ready sends.

  for (int p=0; p < nprocs; p++){
    if (p != rank){
      try {
        Teuchos::readySend<int, T>(comm, sendBuf.view(p*count, count), p);
      } 
      catch (std::exception &e){
        Z2_THROW_OUTSIDE_ERROR(env, e);
      }
    }
  }

  if (req.size() > 0){
    try {
      Teuchos::waitAll<int>(comm, req);
    }
    catch (std::exception &e){
      Z2_THROW_OUTSIDE_ERROR(env, e);
    }
  }
#endif
}

/*! \brief AlltoAllv sends/receives a variable number of objects to/from all processes.
 *
 * The data type T of the objects must be a type for which 
 * Teuchos::SerializationTraits are defined.  This is most likely every 
 * fundamental data type plus std::pair<T1,T2>. It does not
 * include std::vector<T2>.
 */

template <typename T, typename LNO>
void AlltoAllv(const Comm<int> &comm,
              Zoltan2::Environment &env,  
              const ArrayView<const T> &sendBuf,      // input
              const ArrayView<const LNO> &sendCount,  // input
              ArrayRCP<T> &recvBuf,      // output, allocated here
              ArrayRCP<LNO> &recvCount)  // output, allocated here
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  try{
    AlltoAll<LNO, LNO>(comm, env, sendCount, 1, recvCount);
  }
  catch (const std::exception &e){
    Z2_THROW_ZOLTAN2_ERROR(env, e);
  }

  size_t totalIn=0, offsetIn=0, offsetOut=0;

  for (int i=0; i < nprocs; i++){
    totalIn += recvCount[i];
    if (i < rank){
      offsetIn += recvCount[i];
      offsetOut += sendCount[i];
    }
  }

  T *ptr = NULL;
  if (totalIn){
    ptr = new T [totalIn]; 
  }
  Z2_GLOBAL_MEMORY_ASSERTION(env, totalIn, !totalIn||ptr);
  recvBuf = Teuchos::arcp<T>(ptr, 0, totalIn, true);

  T *in = recvBuf.get() + offsetIn;           // Copy self messages
  const T *out = sendBuf.getRawPtr() + offsetOut;

  for (LNO i=0; i < recvCount[rank]; i++){
    in[i] = out[i];
  }

#ifdef HAVE_MPI
  // Post receives

  RCP<CommRequest> r;
  Array<RCP<CommRequest> > req(nprocs-1);

  offsetIn = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && recvCount[p] > 0){
      ArrayRCP<T> recvBufPtr(recvBuf.get() + offsetIn, 0, recvCount[p], false);

      try{
        r  = Teuchos::ireceive<int, T>(comm, recvBufPtr, p);
      }
      catch (const std::exception &e){
        Z2_THROW_OUTSIDE_ERROR(env, e);
      }
    
      req.push_back(r);
    }
    offsetIn += recvCount[p];
  }

  // Wait till all are posted

  Teuchos::barrier(comm);

  // Do all ready sends

  offsetOut = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && sendCount[p] > 0){
      try{
        Teuchos::readySend<int, T>(comm, sendBuf.view(offsetOut, sendCount[p]), p);
      }
      catch(const std::exception &e){
        Z2_THROW_OUTSIDE_ERROR(env, e);
      }
    }
    offsetOut += sendCount[p];
  }

  if (req.size() > 0){
    try{
      Teuchos::waitAll<int>(comm, req);
    }
    catch(const std::exception &e){
      Z2_THROW_OUTSIDE_ERROR(env, e);
    }
  }
#endif
}

/*! \brief Serialization for std::vector<T>
 *
 * Teuchos::SerializationTraits exist for types that can be
 * sent in a Teuchos message, such as std::pair<T1, T2>.  It
 * does not exist for std::vector<T>, and it seems the serialization
 * interface does not make this possible.
 * 
 * These four methods are what the SerializationTraits interface
 * might look like if it could support type std::vector<T>.
 *
 * Buffer layout for std::vector<T> of size N, variable length vectors -
 *      LNO numberOfVectors      
 *      LNO offsetToStartsOfVectorElements[N] 
 *      T first element of first vector
 *        ...
 *      T last element of last vector
 *
 * Buffer layout for std::vector<T> of size N, identical length vectors -
 *      LNO numberOfVectors      
 *      T first element of first vector
 *        ...
 *      T last element of last vector
 *
 * Important: number of bytes returned is always a multiple of sizeof(T)
 */

template <typename T, typename LNO>
  LNO fromObjectsToIndirectBytes(const LNO count, 
    std::vector<T> const v[], 
    LNO vLen=0)   // set vLen to vector length if all are the same length
{
  LNO nelements=0, nbytes=0, preamble=0;

  if (vLen == 0){
    for (LNO i=0; i < count; i++)
      nelements += v[i].size();
    preamble = sizeof(LNO) * (1 + count);
  }
  else{
    nelements = vLen * count;
    preamble = sizeof(LNO);
  }

  nbytes = preamble + (preamble % sizeof(T));  // T alignment

  nbytes += nelements * sizeof(T);

  return nbytes;
}

template <typename T, typename LNO>
  void serialize(
    const LNO count, const std::vector<T> v[], const LNO bytes, char buf[], 
      LNO vLen=0)
{
  LNO preamble = sizeof(LNO);

  if (vLen == 0){
    preamble *= (1 + count);
  }

  LNO nbytes = preamble + (preamble % sizeof(T));  // T alignment

  LNO offset = nbytes / sizeof(T);

  LNO *info = reinterpret_cast<LNO *>(buf);
  T* elements = reinterpret_cast<T *>(buf) + offset;

  *info++ = count;

  if (vLen == 0){
    for (LNO i=0; i < count; i++){
      int nelements = v[i].size();
      *info++ = offset;
  
      for (LNO j=0; j < nelements; j++)
        *elements++ = v[i][j];
  
      offset += nelements;
    }
  }
  else{
    for (LNO i=0; i < count; i++){
      for (LNO j=0; j < vLen; j++){
        *elements++ = v[i][j];
      }
    }
  }
}

template <typename T, typename LNO>
LNO fromIndirectBytesToObjectCount(const LNO bytes, char buf[]) 
{
  LNO *count = reinterpret_cast<LNO *>(buf);
  return count[0];
}

template <typename T, typename LNO>
  void deserialize(const LNO bytes, const char buf[], 
     const LNO count, std::vector<T> v[], LNO vLen=0)
{
  LNO preamble = sizeof(LNO);

  if (vLen == 0){
    preamble *= (1 + count);
  }

  LNO nbytes = preamble + (preamble % sizeof(T));  // T alignment
  LNO offset = nbytes / sizeof(T);

  const T* elements = reinterpret_cast<const T *>(buf) + offset;

  if (vLen > 0){
    for (LNO i=0; i < count; i++){
      v[i].resize(vLen);
      v[i].clear();
      for (LNO j=0; j < vLen; j++){
        v[i].push_back(*elements++);
      }
    }
  }
  else{
    const LNO *info = reinterpret_cast<const LNO *>(buf) + 1;
    LNO lastOffset = LNO(bytes/sizeof(T));
  
    for (LNO i=0; i < count; i++){
  
      LNO length = ((i == count-1) ? lastOffset : info[i+1]) - info[i];
  
      v[i].resize(length);
      v[i].clear();
      for (LNO j=0; j < length; j++){
        v[i].push_back(*elements++);
      }
    }
  }

}

/*! \brief AlltoAllv sends/receives a std::vector<T> to/from all processes.
 *
 * The vectors need not be the same length. The data type T must be a type 
 * for which Teuchos::SerializationTraits are defined.  
 */

template <typename T, typename LNO>
void AlltoAllv(const Comm<int>     &comm,
  Zoltan2::Environment &env,
  const ArrayView<const std::vector<T> > &sendBuf,
  const ArrayView<const LNO>             &sendCount,
  ArrayRCP<std::vector<T> >        &recvBuf,
  ArrayRCP<LNO>                    &recvCount,
  LNO            vLen=0)      // set if all vectors are the same length
{
  int nprocs = comm.getSize();
  size_t totalSendSize = 0;
  LNO offset = 0;
  using Teuchos::is_null;

  LNO *sendSize = new LNO [nprocs];
  Z2_GLOBAL_MEMORY_ASSERTION(env, nprocs, sendSize);

  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0){
      sendSize[p] = 
        fromObjectsToIndirectBytes<T, LNO>(sendCount[p], 
          sendBuf.getRawPtr() + offset, vLen);

      offset += sendCount[p];
      totalSendSize += sendSize[p];
    }
    else{
      sendSize[p] = 0;
    }
  }

  size_t bufSize = totalSendSize/sizeof(T);
  
  T *buf = NULL;
  if (bufSize)
    buf = new T [bufSize];
  
  Z2_GLOBAL_MEMORY_ASSERTION(env, bufSize, !bufSize || buf);

  const std::vector<T> *vptr = sendBuf.getRawPtr();

  char *charBuf = reinterpret_cast<char *>(buf);

  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0){
      serialize<T, LNO>(sendCount[p], vptr, sendSize[p], charBuf, vLen);
      vptr += sendCount[p];
      charBuf += sendSize[p];
      sendSize[p] /= sizeof(T);
    }
  }

  ArrayRCP<T> recvT;
  ArrayRCP<LNO> recvSize;
  ArrayView<const T> bufView(buf, bufSize);
  ArrayView<const LNO> sendSizeView(sendSize, nprocs);

  try{
    AlltoAllv<T, LNO>(comm, env, bufView, sendSizeView, recvT, recvSize);
  }
  catch (const std::exception &e){
    Z2_THROW_ZOLTAN2_ERROR(env, e);
  }

  delete [] sendSize;
  if (bufSize)
    delete [] buf;

  LNO *vectorCount = new LNO [nprocs];
  Z2_GLOBAL_MEMORY_ASSERTION(env, nprocs, vectorCount);

  LNO totalCount = 0;

  charBuf = reinterpret_cast<char *>(recvT.get());

  for (int p=0; p < nprocs; p++){
    if (recvSize[p] > 0){
      LNO bytecount = recvSize[p] * sizeof(T);
      vectorCount[p] = 
        fromIndirectBytesToObjectCount<T, LNO>(bytecount, charBuf);

      charBuf += bytecount;
      totalCount += vectorCount[p];
    }
    else{
      vectorCount[p] = 0;
    }
  }

  std::vector<T> *inVectors = NULL;
  if (totalCount)
    inVectors = new std::vector<T> [totalCount];
  Z2_GLOBAL_MEMORY_ASSERTION(env, nprocs, !totalCount || inVectors);

  charBuf = reinterpret_cast<char *>(recvT.get());
  std::vector<T> *inv = inVectors;

  for (int p=0; p < nprocs; p++){
    if (recvSize[p] > 0){
      LNO bytecount = recvSize[p] * sizeof(T);
      deserialize<T, LNO>(bytecount, charBuf, vectorCount[p], inv, vLen);

      charBuf += bytecount;
      inv += vectorCount[p];
    }
  }

  recvBuf = Teuchos::arcp(inVectors, 0, totalCount);
  recvCount = Teuchos::arcp(vectorCount, 0, nprocs);
}

}                   // namespace Z2
#endif
