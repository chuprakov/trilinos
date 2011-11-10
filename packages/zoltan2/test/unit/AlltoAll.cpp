// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER

// TODO: doxygen comments
//     make this a real unit test that gives helpful information if it fails
//     and uses different template values

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <Teuchos_RCP.hpp>   
#include <Teuchos_ArrayRCP.hpp>   
#include <Teuchos_Comm.hpp>   
#include <Teuchos_ParameterList.hpp>   
#include <Teuchos_DefaultComm.hpp>   
#include <Zoltan2_Environment.hpp>   
#include <Zoltan2_AlltoAll.hpp>   

using namespace std;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  Teuchos::ParameterList params;
  params.set(std::string("ERROR_CHECK_LEVEL"), 1);
  params.set(std::string("DEBUG_LEVEL"), 0);
        
  Teuchos::RCP<Zoltan2::Environment> envPtr = 
    Teuchos::rcp(new Zoltan2::Environment(params, comm));

  // In this test, our local IDs are ints and our global IDs are longs.

  int errcode = 0;
  if (!errcode){

    // test of Zoltan2::AlltoAll using a Scalar type (ints)
  
    int *sendBuf = new int [2*nprocs];
    
    Teuchos::ArrayView<const int> sendBufView(sendBuf, 2*nprocs);
  
    for (int i=0, j=1; i < 2*nprocs ; i+=2,j++){
      sendBuf[i] = j*10;
      sendBuf[i+1] = j*10 + 1;
    } 
  
    Teuchos::ArrayRCP<int> recvBuf;
    
    Zoltan2::AlltoAll<int, int>(*comm, *envPtr,
        sendBufView,    // ints to send from this process to all the others
        2,              // two ints per process
        recvBuf);       // will be allocated and filled in AlltoAll

    delete [] sendBuf;
  
    int *inBuf = recvBuf.get();
  
    int myvals[2] = {(rank+1) * 10, (rank+1) * 10 + 1};
  
    for (int i=0; i < 2*nprocs; i+=2){
      if (inBuf[i] != myvals[0] && inBuf[i+1] != myvals[1]){
        errcode = 1;
        break;
      }
    }
  }

  if (!errcode){

    // test of Zoltan2::AlltoAll using a non-Scalar type for which
    //  SerializationTraits are defined in Teuchos (std::pair).
  
    std::pair<int, int> *outBufPair = new std::pair<int, int> [nprocs];

    Teuchos::ArrayView< const std::pair<int, int> > sendBufPair(
      outBufPair, nprocs);
  
    for (int i=0,j=1; i < nprocs ; i++,j++){
      outBufPair[i].first = j*10;
      outBufPair[i].second = j*10 + 1;
    } 
  
    Teuchos::ArrayRCP< std::pair<int, int> > recvBufPair;
    
    Zoltan2::AlltoAll<std::pair<int, int> , int>(*comm, *envPtr,
        sendBufPair,    // ints to send from this process to all the others
        1,              // one pair per process
        recvBufPair);   // will be allocated and filled in AlltoAll

    delete [] outBufPair;
  
    std::pair<int, int > *inBufPair = recvBufPair.get();

    int myvals[2] = {(rank+1) * 10, (rank+1) * 10 + 1};
  
    for (int i=0; i < nprocs; i++){
      if (inBufPair[i].first != myvals[0] && inBufPair[i].second != myvals[1]){
        errcode = 1;
        break;
      }
    }
  }

  if (!errcode){

    // test of Zoltan2::AlltoAllv (different sized messages) using a Scalar type

    int myMsgSizeBase=rank*nprocs + 1;
    int *outMsgSizes = new int [nprocs];
    long totalOut = 0;

    for (int p=0; p < nprocs; p++){
      outMsgSizes[p] = myMsgSizeBase + p;
      totalOut += outMsgSizes[p];
    }
    Teuchos::ArrayView<const int> sendCount(outMsgSizes, nprocs);
  
    int *outBuf = new int [totalOut];
  
    int *out = outBuf;
    for (int p=0; p < nprocs; p++){
      for (int i=0; i < outMsgSizes[p]; i++){
        *out++ = p+rank;
      }
    }
    Teuchos::ArrayView<const int> sendBuf(outBuf, totalOut);
  
    Teuchos::ArrayRCP<int> recvCount;
    Teuchos::ArrayRCP<int> recvBuf;
    
    Zoltan2::AlltoAllv<int, int>(*comm, *envPtr,
                  sendBuf,    
                  sendCount,   
                  recvBuf,
                  recvCount);

    delete [] outBuf;
    delete [] outMsgSizes;
  
    int *inMsgSizes = recvCount.get();
    int *inBuf = recvBuf.get();

    for (int p=0; p < nprocs; p++){
      for (int i=0; i < inMsgSizes[p]; i++){
        if (*inBuf++ != rank+p){
          errcode = 1;
          break;
        }
      }
    }
  }

  if (!errcode){

    // test of Zoltan2::AlltoAllv using vectors - which can not
    //    be serialized by Teuchos.

    int myMsgSizeBase=rank*nprocs + 1;
    int *outMsgSizes = new int [nprocs];
    long totalOut = 0;

    for (int p=0; p < nprocs; p++){
      outMsgSizes[p] = myMsgSizeBase + p;
      totalOut += outMsgSizes[p];
    }
    Teuchos::ArrayView<const int> sendCount(outMsgSizes, nprocs);
  
    std::vector<float> *outBuf = new std::vector<float> [totalOut];
  
    std::vector<float> *out = outBuf;
    for (int p=0; p < nprocs; p++){
      for (int i=0; i < outMsgSizes[p]; i++, out++){
        (*out).reserve(2);
        (*out).push_back(p+rank);
        (*out).push_back(p+rank + 0.5);
      }
    }
    Teuchos::ArrayView<const std::vector<float> > sendBuf(outBuf, totalOut);
  
    Teuchos::ArrayRCP<int> recvCount;
    Teuchos::ArrayRCP<std::vector<float> > recvBuf;
    
    Zoltan2::AlltoAllv<float , int>(*comm, *envPtr,
        sendBuf,    
        sendCount,   
        recvBuf,
        recvCount,
        2);        // optimization: all vectors are size 2

    delete [] outMsgSizes;
    delete [] outBuf;
  
    int *inMsgSizes = recvCount.get();
    std::vector<float> *inBuf = recvBuf.get();

    for (int p=0; p < nprocs; p++){
      for (int i=0; i < inMsgSizes[p]; i++, inBuf++){
        std::vector<float> &v = *inBuf;
        if ( (v[0] != rank + p) || (v[1] != rank +p +0.5)){
          errcode = 1;
          break;
        }
      }
    }
  }

  if (rank == 0){
    if (errcode)
      std::cout << "FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return errcode;
}
