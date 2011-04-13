//@HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
*/

#include <stdexcept>
#include <vector>
#include <map>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>
#include <Zoltan2_IdentifierTraits.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_AlltoAll.hpp>

namespace Z2
{

template<typename AppLID, typename AppGID, typename LNO, typename GNO> 
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::ArrayRCP<AppGID> &gids, 
    Teuchos::ArrayRCP<AppLID> &lids) 
         : _comm(in_comm), _myGids(gids), _myLids(lids),
           _globalNumberOfIds(0), _localNumberOfIds(0), _haveLocalIds(false),
           _myRank(0), _numProcs(0)
{
  _numProcs = _comm->getSize(); 
  _myRank = _comm->getRank(); 

  _localNumberOfIds = _myGids.size();

  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<std::string, LNO> id2index_hash_t;

  Teuchos::Tuple<teuchos_size_t, 4> counts;

  counts[0] = _myLids.size();
  counts[1] = _localNumberOfIds;
  counts[2] = counts[3] = 0;

  Teuchos::reduceAll<int, teuchos_size_t>(*_comm, Teuchos::REDUCE_SUM, 
    2, counts.getRawPtr(), counts.getRawPtr()+2);

  _haveLocalIds = (counts[2] > 0);
  _globalNumberOfIds = counts[3];

  TEST_FOR_EXCEPTION(_haveLocalIds && (counts[0] != _localNumberOfIds),
             std::runtime_error,
          "local IDs are provided but number of global IDs"
          " does not equal number of local IDs");

  if (_haveLocalIds){   // hash LID to index in LID vector
    id2index_hash_t *p = new id2index_hash_t(_localNumberOfIds);

    AppLID *lidPtr = _myLids.get();  // for performance

    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      p->put(IdentifierTraits<AppLID>::key(lidPtr[i]), i);
    }

    _lidHash = Teuchos::RCP<id2index_hash_t>(p);
  }

  // If the application's global ID data type (AppGID) is a Teuchos Ordinal,
  // we will be able to use it as our internal global numbers (GNO).  
  // Otherwise we will have to map their global IDs to valid global numbers.

  if (IdentifierTraits<AppGID>::isGlobalOrdinalType()){

    // Are the AppGIDs consecutive and increasing with process rank? 
    // If so GID/proc lookups can be optimized.

    AppGID min(0), max(0), globalMin(0), globalMax(0);
    AppGID *gidPtr = _myGids.get();  // for performance
    min = max = gidPtr[0];
    AppGID checkVal = min;
    bool consecutive = true;
  
    for (teuchos_size_t i=1; i < _localNumberOfIds; i++){
      if (consecutive && (gidPtr[i] != ++checkVal)){
        consecutive=false;
        break;
      }
      if (gidPtr[i] < min)
        min = gidPtr[i];
      else if (gidPtr[i] > max)
        max = gidPtr[i];
    }

    Teuchos::Tuple<AppGID, 4> results;

    results[0] = static_cast<AppGID>(consecutive ? 1 : 0);
    results[1] = min;
    results[2] = results[3] = 0;

    Teuchos::reduceAll<int, AppGID>(*_comm, Teuchos::REDUCE_MIN, 2, 
      results.getRawPtr(), results.getRawPtr()+2);

    if (results[2] != 1)       // min of consecutive flags
      consecutive=false;

    if (consecutive){
      globalMin = results[3];
      Teuchos::reduceAll<int, AppGID>(*_comm, Teuchos::REDUCE_MAX, 1, &max, 
        &globalMax);
      if (globalMax - globalMin + 1 != static_cast<AppGID>(_globalNumberOfIds))
        consecutive = false;   // there are gaps in the gids
  
      if (consecutive){
        AppGID myStart = _myGids[0];

        _gnoDist = Teuchos::ArrayRCP<GNO>(_numProcs + 1);

        AppGID *startGID = static_cast<AppGID *>(_gnoDist.getRawPtr());
      
        Teuchos::gatherAll<int, AppGID>(*_comm, 1, &myStart, _numProcs, 
          startGID);
      
        for (int p=1; p < _numProcs; p++){
          if (startGID[p] < startGID[p-1]){
            consecutive = false;  // gids do not increase with process rank
            break;
          }
        }
        if (consecutive){
          startGID[_numProcs] = globalMax + 1;
        }
      }
    }
  }
  else{

    // AppGIDs are not Ordinals.  We map them to consecutive 
    // global numbers starting with 0. 

    _gnoDist = Teuchos::ArrayRCP<GNO>(_numProcs + 1, 0);
    GNO myNum = static_cast<GNO>(_localNumberOfIds);

    Teuchos::scan<int, GNO>(*_comm, Teuchos::REDUCE_SUM, 1, &myNum, 
      _gnoDist.getRawPtr() + 1);
  }

  if (_gnoDist.size() == 0){

    // We need a hash table mapping the application global ID
    // to its index in _myGids.

    id2index_hash_t *p =  new id2index_hash_t(_localNumberOfIds);

    AppGID *gidPtr = _myGids.get();  // for performance

    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      p->put(IdentifierTraits<AppGID>::key(gidPtr[i]), i);
    }

    _gidHash = Teuchos::RCP<id2index_hash_t>(p);
  }
}

  /*! Constructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap()  
         : _comm(), _myGids(), _myLids(),
           _globalNumberOfIds(0), _localNumberOfIds(0), _haveLocalIds(false),
           _myRank(0), _numProcs(0)
{
}

  /*! Destructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::~IdentifierMap() 
  {
  }

  /*! Copy Constructor */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO>::IdentifierMap(const IdentifierMap &id)
{
    //TODO
}

  /*! Assignment operator */
template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  IdentifierMap<AppLID,AppGID,LNO,GNO> &
    IdentifierMap<AppLID,AppGID,LNO,GNO>::operator=(const IdentifierMap<AppLID,
                  AppGID,LNO,GNO> &id)
{
    //TODO
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID,AppGID,LNO,GNO>::initialize(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::ArrayRCP<AppGID> &gids, 
    Teuchos::ArrayRCP<AppLID> &lids) 
{
  _gnoDist.relese();
  _gidHash.release();
  _lidHash.release();

  this->IdentifierMap(in_comm, gids, lids);
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  bool IdentifierMap<AppLID, AppGID, LNO, GNO>::gnosAreGids()
{
  return IdentifierTraits<AppGID>::isGlobalOrdinalType();
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidTranslate(
    Teuchos::ArrayView<AppGID> &gid, 
    Teuchos::ArrayView<GNO > &gno,
    TranslationType tt)
{
  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  teuchos_size_t len=gid.size();

  if (len == 0){
    return;
  }

  TEST_FOR_EXCEPTION((tt!=TRANSLATE_GNO_TO_GID) && (tt!=TRANSLATE_GID_TO_GNO),
                  std::runtime_error,
                 "TranslationType is invalid");

  TEST_FOR_EXCEPTION(
    ((tt==TRANSLATE_GNO_TO_GID) && (gid.size() < gno.size())) ||
    ((tt==TRANSLATE_GID_TO_GNO) && (gno.size() < gid.size())),
                  std::runtime_error,
                 "Destination array is too small");


  if (IdentifierTraits<AppGID>::isGlobalOrdinalType()){   
                              // our gnos are the app gids
    if (tt == TRANSLATE_GNO_TO_GID)
      for (teuchos_size_t i=0; i < len; i++)
        gid[i] = static_cast<AppGID>(gno[i]);
    else
      for (teuchos_size_t i=0; i < len; i++)
        gno[i] = static_cast<GNO>(gid[i]);
  }
  else{              // we mapped gids to consecutive gnos
  
    GNO firstGNO = _gnoDist[_myRank];
    if (tt == TRANSLATE_GNO_TO_GID)
      for (teuchos_size_t i=0; i < len; i++)
        gid[i] = _myGids[gno[i] - firstGNO];
    else
      for (teuchos_size_t i=0; i < len; i++){
        const LNO &idx = _gidHash->get(
          Z2::IdentifierTraits<AppGID>::key(gid[i]));
        gno[i] = firstGNO + idx;
      }
  }
  return;
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::lidTranslate(
    Teuchos::ArrayView<AppLID> &lid, 
    Teuchos::ArrayView<GNO> &gno, 
    TranslationType tt)
{
  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  teuchos_size_t len=lid.size();

  if (len == 0){
    return;
  }

  TEST_FOR_EXCEPTION((tt!=TRANSLATE_GNO_TO_LID) && (tt!=TRANSLATE_LID_TO_GNO),
                  std::runtime_error,
                 "TranslationType is invalid");

  TEST_FOR_EXCEPTION(
    ((tt==TRANSLATE_GNO_TO_LID) && (lid.size() < gno.size())) ||
    ((tt==TRANSLATE_LID_TO_GNO) && (gno.size() < lid.size())),
                  std::runtime_error,
                 "Destination array is too small");

  TEST_FOR_EXCEPTION(!_haveLocalIds,
                  std::runtime_error,
                 "local ID translation is requested but none were provided");

  GNO firstGNO;
  if (_gnoDist.size() > 0)
    firstGNO = _gnoDist[_myRank];
  
  if (TRANSLATE_GNO_TO_LID){
    for (teuchos_size_t i=0; i < len; i++){
      LNO idx = 0;
      if (_gnoDist.size() > 0)  // gnos are consecutive
        idx = gno[i] - firstGNO;
      else                      // gnos must be the app gids
        idx = _gidHash->get(
          IdentifierTraits<AppGID>::key(static_cast<AppGID>(gno[i])));
      
      lid[i] = _myLids[idx];
    }
  }
  else{
    for (teuchos_size_t i=0; i < len; i++){
      LNO idx = _lidHash->get(IdentifierTraits<AppLID>::key(lid[i]));

      if (_gnoDist.size() > 0)  // gnos are consecutive
        gno[i] = firstGNO + idx;
      else                     // gnos must be the app gids
        gno[i] = static_cast<GNO>(_myGids[idx]);
    }
  }
}

template<typename AppLID, typename AppGID, typename LNO, typename GNO>
  void IdentifierMap<AppLID, AppGID, LNO, GNO>::gidGlobalTranslate(
    Teuchos::ArrayView<const AppGID> &in_gid,
    Teuchos::ArrayView<GNO> &out_gno,
    Teuchos::ArrayView<int> &out_proc)
{
  typedef typename Teuchos::Array<GNO>::size_type teuchos_size_t;
  typedef typename Teuchos::Hashtable<std::string, LNO> id2index_hash_t;
  typedef typename Teuchos::Hashtable<std::string, Teuchos::Array<LNO> > 
    id2index_array_hash_t;

  teuchos_size_t len=in_gid.size();

  if (len == 0){
    return;
  }

  TEST_FOR_EXCEPTION( (out_gno.size() < len) || (out_proc.size() < len),
                  std::runtime_error,
                 "Destination array is too small");


  if (IdentifierTraits<AppGID>::isGlobalOrdinalType() && (_gnoDist.size() > 0)){

    // Easy case - communication is not needed.
    // Global numbers are the application global IDs and
    // they are increasing consecutively with rank.
 
    typename std::map<GNO, int> firstGnoToProc;
    typename std::map<GNO, int>::iterator pos;

    for (int p=0; p <= _numProcs; p++){
      firstGnoToProc[_gnoDist[p]] = p;
    }

    for (teuchos_size_t i=0; i < len; i++){
      GNO globalNumber = static_cast<GNO>(in_gid[i]);
      out_gno[i] = globalNumber;
      pos = firstGnoToProc.upper_bound(globalNumber);
      out_proc[i] = pos->first - 1;
    }

    return;
  }

  bool needGnoInfo = !IdentifierTraits<AppGID>::isGlobalOrdinalType();

  ///////////////////////////////////////////////////////////////////////
  // First: Hash each of my AppGIDs to a process that will answer
  // for it.  Send my AppGIDs (and the Gnos if they are different)
  // to their assigned processes.  Build a search structure for
  // the AppGIDs that were assigned to me, so I can reply with
  // with the process owning them (and their Gnos if they are different).
  ///////////////////////////////////////////////////////////////////////

  Teuchos::Array<int> hashProc(0);
  Teuchos::Array<AppGID> gidOutBuf(0);
  Teuchos::Array<GNO> gnoOutBuf(0);
  Teuchos::Array<LNO> countOutBuf(_numProcs, 0);
  Teuchos::Array<LNO> offsetBuf(_numProcs + 1, 0);

  Teuchos::ArrayRCP<AppGID> gidInBuf();
  Teuchos::ArrayRCP<GNO> gnoInBuf();
  Teuchos::ArrayRCP<LNO> countInBuf();

  if (_localNumberOfIds > 0){

    hashProc.reserve(_localNumberOfIds);
    gidOutBuf.reserve(_localNumberOfIds);

    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      hashProc[i] = IdentifierTraits<AppGID>::hashCode(_myGids[i]) % _numProcs;
      countOutBuf[hashProc[i]]++;
    }
  
    for (int p=1; p <= _numProcs; p++){
      offsetBuf[p] = offsetBuf[p-1] + countOutBuf[p-1];
    }
  
    if (needGnoInfo){   
      // The gnos are not the gids, which also implies that
      // gnos are consecutive numbers given by _gnoDist.
      gnoOutBuf.resize(_localNumberOfIds);
    }
  
    for (teuchos_size_t i=0; i < _localNumberOfIds; i++){
      LNO offset = offsetBuf[hashProc[i]];
      gidOutBuf[offset] = _myGids[i];
      if (needGnoInfo)
        gnoOutBuf[offset] = _gnoDist[_myRank] + i;
      offsetBuf[hashProc[i]] = offset + 1;
    }
    hashProc.clear();
  }

  // Teuchos comment #1: An Array can be passed for an ArrayView parameter.
  // Teuchos comment #2: AppGID need not be a Teuchos Packet type,
  //                     so we wrote our own AlltoAllv.
  // Z2::AlltoAllv comment: Buffers are in process rank order.

  AlltoAllv(*_comm, gidOutBuf, countOutBuf, gidInBuf, countInBuf);

  gidOutBuf.clear();
  
  if (needGnoInfo){
    countInBuf.release();
    AlltoAllv(*_comm, gnoOutBuf, countOutBuf, gnoInBuf, countInBuf);
  }

  gnoOutBuf.clear();
  countOutBuf.clear();

  //
  // Save the information that was hashed to me so I can do lookups.
  //

  std::map<LNO, int> firstIndexToProc;
  LNO total = 0;

  for (int p=0; p < _numProcs; p++){
    firstIndexToProc[total] = p;
    total += countInBuf[p];
  }

  firstIndexToProc[total] = _numProcs;

  id2index_hash_t gidToIndex(total);

  total = 0;
  for (int p=0; p < _numProcs; p++){
    for (LNO i=countInBuf[p]; i < countInBuf[p+1]; i++, total++){
      gidToIndex.put(IdentifierTraits<AppGID>::key(gidInBuf[total]), total);
    }
  }

  // Keep gnoInBuf.  We're done with the others.

  gidInBuf.release();
  countInBuf.release();

  ///////////////////////////////////////////////////////////////////////
  // Send a request for information to the "answer process" for each 
  // of the GIDs in in_gid.  First remove duplicate GIDs from the list.
  ///////////////////////////////////////////////////////////////////////

  
  Teuchos::Array<std::string> uniqueGidQueries(0);
  Teuchos::Array<Teuchos::Array<LNO> > uniqueGidQueryIndices(0);
  teuchos_size_t numberOfUniqueGids = 0;
  Teuchos::Array<LNO> gidLocation(0);

  countOutBuf.resize(_numProcs, 0);

  if (len > 0){
  
    // For efficiency, guess a reasonable size for the array.
    // (In an input adapter, how many objects will have the same neighbor?)

    teuchos_size_t sizeChunk = 4;

    teuchos_size_t tableSize = len / sizeChunk;

    tableSize =  (tableSize < 1) ? 1 : tableSize;

    id2index_array_hash_t *gidIndices = new id2index_array_hash_t(tableSize);
  
    for (LNO i=0; i < len; i++){

      std::string uniqueKey(IdentifierTraits<AppGID>::key(in_gid[i]));

      if (gidIndices->containsKey(uniqueKey)){
        Teuchos::Array<LNO> &v = gidIndices->get(uniqueKey);
        teuchos_size_t n = v.size();
        if (n % sizeChunk == 0){
          v.reserve(n + sizeChunk);
        }
        v.push_back(i);
      }
      else{
        Teuchos::Array<LNO> v(sizeChunk);
        v[0] = i;
        gidIndices->put(uniqueKey, v);
      }
    }
  
    numberOfUniqueGids = gidIndices->size();

    gidIndices->arrayify(uniqueGidQueries, uniqueGidQueryIndices);
  
    delete gidIndices;
  
    gidOutBuf.reserve(numberOfUniqueGids);
    hashProc.reserve(numberOfUniqueGids);
  
    for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){
      hashProc[i] = Teuchos::hashCode(uniqueGidQueries[i]) % _numProcs;
      countOutBuf[hashProc[i]]++;
    }
  
    offsetBuf[0] = 0;
  
    for (int p=0; p < _numProcs; p++){
      offsetBuf[p+1] = offsetBuf[p] + countOutBuf[p];
    }
  
    gidLocation.reserve(numberOfUniqueGids);
  
    for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){
      AppGID gid = IdentifierTraits<AppGID>::keyToGid(uniqueGidQueries[i]);
      gidLocation[i] = offsetBuf[hashProc[i]];
      gidOutBuf[gidLocation[i]] = gid;
      offsetBuf[hashProc[i]] = gidLocation[i] + 1;
    }

    hashProc.clear();
  }

  AlltoAllv(*_comm, gidOutBuf, countOutBuf, gidInBuf, countInBuf);

  gidOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Create and send answers to the processes that made requests of me.
  ///////////////////////////////////////////////////////////////////////

  total = 0;

  for (int p=0; p < _numProcs; p++){
    countOutBuf[p] = countInBuf[p];
    total += countOutBuf[p];
  }

  Teuchos::Array<int> procOutBuf(total);
  Teuchos::ArrayRCP<int> procInBuf();

  if (needGnoInfo){
    gnoOutBuf.reserve(total);
  }

  if (total > 0){
  
    total=0;
  
    for (int p=0; p < _numProcs; p++){
      for (LNO i=0; i < countInBuf[p]; i++, total++){
        std::string s(IdentifierTraits<AppGID>::key(gidInBuf[total]));
        LNO index = gidToIndex.get(s);
        int proc = firstIndexToProc.upper_bound(index);
        procOutBuf[total] = proc-1;
  
        if (needGnoInfo){
          gnoOutBuf[total] = gnoInBuf[index];
        }
      }
    }
  }

  gidInBuf.release();
  if (needGnoInfo){
    gnoInBuf.release();
  }

  AlltoAllv(*_comm, procOutBuf, countOutBuf, procInBuf, countInBuf);

  procOutBuf.clear();

  if (needGnoInfo){
    AlltoAllv(*_comm, gnoOutBuf, countOutBuf, gnoInBuf, countInBuf);
    gnoOutBuf.clear();
  }

  countOutBuf.clear();

  ///////////////////////////////////////////////////////////////////////
  // Done.  Process the replies to my queries
  ///////////////////////////////////////////////////////////////////////

  for (teuchos_size_t i=0; i < numberOfUniqueGids; i++){

    std::string s(uniqueGidQueries[i]);
    Teuchos::Array<LNO> v(uniqueGidQueryIndices[i]);

    int gidProc = procInBuf[gidLocation[i]];

    LNO gno;
    if (needGnoInfo){
      gno = gnoInBuf[gidLocation[i]];
    }
    else{
      gno = static_cast<GNO>(IdentifierTraits<AppGID>::keyToGid(s));
    }

    for (teuchos_size_t j=0; j < v.size(); j++){
      out_gno[v[j]] = gno;
      out_proc[v[j]] = gidProc;
    }
  }
}

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_ */
