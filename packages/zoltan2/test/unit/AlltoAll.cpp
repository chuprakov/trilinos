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
#include <Zoltan2_AlltoAll.hpp>
#include <Teuchos_GlobalMPISession.hpp>   
#include <Teuchos_RCP.hpp>   
#include <Teuchos_ArrayRCP.hpp>   
#include <Teuchos_Comm.hpp>   
#include <Teuchos_DefaultComm.hpp>   

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int rank = session.getRank();
  int nprocs = session.getNProc();
  bool pass = true;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  if (pass){

    // test of Z2::AlltoAll using a ordinal type
  
    int *sendBuf = new int [2*nprocs];
    
    Teuchos::ArrayView<int> sendBufView(sendBuf, 2*nprocs);
  
    for (int i=0, j=1; i < 2*nprocs ; i+=2,j++){
      sendBuf[i] = j*10;
      sendBuf[i+1] = j*10 + 1;
    } 
  
    Teuchos::ArrayRCP<int> recvBuf;
    
    Z2::AlltoAll<int, long, int>(*comm, 
        sendBufView,    // ints to send from this process to all the others
        2,              // two ints per process
        recvBuf);       // will be allocated and filled in AlltoAll

    delete [] sendBuf;
  
    int *inBuf = recvBuf.get();
  
    int myvals[2] = {(rank+1) * 10, (rank+1) * 10 + 1};
  
    for (int i=0; i < 2*nprocs; i+=2){
      if (inBuf[i] != myvals[0] && inBuf[i+1] != myvals[1]){
        pass = false;
        break;
      }
    }
  }

  if (pass){

    // test of Z2::AlltoAll using a non-ordinal type
  
    std::pair<int, int> *outBufPair = new std::pair<int, int> [nprocs];

    Teuchos::ArrayView< std::pair<int, int> > sendBufPair(outBufPair, nprocs);
  
    for (int i=0,j=1; i < nprocs ; i++,j++){
      outBufPair[i].first = j*10;
      outBufPair[i].second = j*10 + 1;
    } 
  
    Teuchos::ArrayRCP< std::pair<int, int> > recvBufPair;
    
    Z2::AlltoAll<std::pair<int, int> , long, int>(*comm, 
        sendBufPair,    // ints to send from this process to all the others
        1,              // one pair per process
        recvBufPair);   // will be allocated and filled in AlltoAll

    delete [] outBufPair;
  
    std::pair<int, int > *inBufPair = recvBufPair.get();

    int myvals[2] = {(rank+1) * 10, (rank+1) * 10 + 1};
  
    for (int i=0; i < nprocs; i++){
      if (inBufPair[i].first != myvals[0] && inBufPair[i].second != myvals[1]){
        pass = false;
        break;
      }
    }
  }

  if (pass){

    // test of Z2::AlltoAllv (different sized messages) using a ordinal type

    long myMsgSizeBase=rank*nprocs + 1;
    long *outMsgSizes = new long [nprocs];
    long totalOut = 0;

    for (int p=0; p < nprocs; p++){
      outMsgSizes[p] = myMsgSizeBase + p;
      totalOut += outMsgSizes[p];
    }

    Teuchos::ArrayView<long> sendCount(outMsgSizes, nprocs);
  
    int *outBuf = new int [totalOut];
    Teuchos::ArrayView<int> sendBuf(outBuf, totalOut);
  
    int *out = outBuf;
    for (int p=0; p < nprocs; p++){
      for (int i=0; i < outMsgSizes[p]; i++){
        *out++ = p+rank;
      }
    }
  
    Teuchos::ArrayRCP<long> recvCount;
    Teuchos::ArrayRCP<int> recvBuf;
    
    Z2::AlltoAllv<int, long, int>(*comm, 
                  sendBuf,    
                  sendCount,   
                  recvBuf,
                  recvCount);

    delete [] outBuf;
    delete [] outMsgSizes;
  
    long *inMsgSizes = recvCount.get();
    int *inBuf = recvBuf.get();

    for (int p=0; p < nprocs; p++){
      for (int i=0; i < inMsgSizes[p]; i++){
        if (*inBuf++ != rank+p){
          pass = false;
          break;
        }
      }
    }
  }

  if (pass){

    // test of Z2::AlltoAllv using vectors - which can not
    //    be serialized by Teuchos (so they can't be Scalars).

    long myMsgSizeBase=rank*nprocs + 1;
    long *outMsgSizes = new long [nprocs];
    long totalOut = 0;

    for (int p=0; p < nprocs; p++){
      outMsgSizes[p] = myMsgSizeBase + p;
      totalOut += outMsgSizes[p];
    }

    Teuchos::ArrayView<long> sendCount(outMsgSizes, nprocs);
  
    std::vector<float> *outBuf = new std::vector<float> [totalOut];
    Teuchos::ArrayView<std::vector<float> > sendBuf(outBuf, totalOut);
  
    std::vector<float> *out = outBuf;
    for (int p=0; p < nprocs; p++){
      for (int i=0; i < outMsgSizes[p]; i++){
        (*out).push_back(p+rank);
        (*out).push_back(p+rank + 0.5);
      }
    }
  
    Teuchos::ArrayRCP<long> recvCount;
    Teuchos::ArrayRCP<std::vector<float> > recvBuf;
    
    Z2::AlltoAllv<std::vector<float> , long, int>(*comm, 
                  sendBuf,    
                  sendCount,   
                  recvBuf,
                  recvCount);

    delete [] outMsgSizes;
    delete [] outBuf;
  
    long *inMsgSizes = recvCount.get();
    std::vector<float> *inBuf = recvBuf.get();

    for (int p=0; p < nprocs; p++){
      for (int i=0; i < inMsgSizes[p]; i++){
        std::vector<float> &v = *inBuf++;
        if ( (v[0] != rank + p) || (v[1] != rank +p +0.5)){
          pass = false;
          break;
        }
      }
    }
  }
#endif

  if (!pass)
    std::cout << "FAIL" << std::endl;
  else
    std::cout << "PASS" << std::endl;


  return 0;
}
