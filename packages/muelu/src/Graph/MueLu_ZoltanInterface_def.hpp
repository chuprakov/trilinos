#ifndef MUELU_ZOLTANINTERFACE_DEF_HPP
#define MUELU_ZOLTANINTERFACE_DEF_HPP

#include "MueLu_ZoltanInterface_decl.hpp"
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#include <Teuchos_Utils.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {


  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  ZoltanInterface(RCP<const FactoryBase> AFact, RCP<const FactoryBase> TransferFact)
    : AFact_(AFact), TransferFact_(TransferFact)
  {}

  //-------------------------------------------------------------------------------------------------------------
  // DeclareInput
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  DeclareInput(Level & level) const
  {
    level.DeclareInput("A", AFact_.get());
    level.DeclareInput("Coordinates", NoFactory::get()); //FIXME JJH
    //level.DeclareInput("Coordinates", TransferFact_.get()); //FIXME JJH
  } //DeclareInput()

  //-------------------------------------------------------------------------------------------------------------
  // Build
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  Build(Level &level) const
  {

    FactoryMonitor m(*this, "ZoltanInterface", level);

    RCP<Operator> A = level.Get< RCP<Operator> >("A",AFact_.get());
    // Tell Zoltan what kind of local/global IDs we will use.
    // In our case, each GID is two ints and there are no local ids.
    // One can skip this step if the IDs are just single ints.
    RCP<const Teuchos::MpiComm<int> > mpiComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A->getRowMap()->getComm());
    float zoltanVersion_;
    Zoltan_Initialize(0, NULL, &zoltanVersion_);
    //TODO define zoltanComm_ as a subcommunicator?!;
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > zoltanComm_ = mpiComm->getRawMpiComm();
    RCP<Zoltan> zoltanObj_ = rcp( new Zoltan( (*zoltanComm_)() ) );  //extract the underlying MPI_Comm handle and create a Zoltan object
    if (zoltanObj_==Teuchos::null) throw(Exceptions::RuntimeError("MueLu::Zoltan : Unable to create Zoltan data structure"));
    int rv;
    if ((rv=zoltanObj_->Set_Param("num_gid_entries", "1")) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_gid_entries' returned error code " + Teuchos::toString(rv)));
    if ( (rv=zoltanObj_->Set_Param("num_lid_entries", "0") ) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'num_lid_entries' returned error code " + Teuchos::toString(rv)));
    if ( (rv=zoltanObj_->Set_Param("obj_weight_dim", "1") ) != ZOLTAN_OK )
      throw(Exceptions::RuntimeError("MueLu::Zoltan::Setup : setting parameter 'obj_weight_dim' returned error code " + Teuchos::toString(rv)));

    GO numPartitions_ = level.Get<GO>("number of partitions");
    std::stringstream ss;
    ss << numPartitions_;
    zoltanObj_->Set_Param("num_global_partitions",ss.str());

    //if (level.IsAvailable("Coordinates",TransferFact_.get()) == false) //FIXME JJH
    //~~ if (level.IsAvailable("Coordinates",NoFactory::get()) == false) //FIXME JJH
    //~~  throw(Exceptions::HaltRepartitioning("MueLu::ZoltanInterface : no coordinates available"));
    //RCP<MultiVector> XYZ = level.Get< RCP<MultiVector> >("Coordinates",TransferFact_.get()); //FIXME JJH
    //~~    RCP<MultiVector> XYZ = level.Get< RCP<MultiVector> >("Coordinates");

    //TODO: coordinates should be const

    Array<ArrayRCP<SC> > XYZ; // Using this format because no communications needed here. No need for a map and a Xpetra::MultiVector

    // Build XYZ from XCoordinates, YCoordinates and ZoltanInterface
    if (level.IsAvailable("XCoordinates")) {
      
      { 
        XYZ.push_back(level.Get< ArrayRCP<SC> >("XCoordinates"));
      }
      
      if (level.IsAvailable("YCoordinates")) {
        XYZ.push_back(level.Get< ArrayRCP<SC> >("YCoordinates"));
      }
      
      if (level.IsAvailable("ZCoordinates")) {
        XYZ.push_back(level.Get< ArrayRCP<SC> >("ZCoordinates"));
      }

    } else if (level.IsAvailable("Coordinates")) {

      RCP<MultiVector> multiVectorXYZ = level.Get< RCP<MultiVector> >("Coordinates");
      for (int i=0; i< (int)multiVectorXYZ->getNumVectors(); i++) { //FIXME cast
        XYZ.push_back(multiVectorXYZ->getDataNonConst(i)); // It's OK to leave 'open' the MultiVector until the destruction of XYZ because no communications using Xpetra
      }

    } else {
      throw(Exceptions::HaltRepartitioning("MueLu::ZoltanInterface : no coordinates available"));
    }

    //~~ size_t problemDimension_ = XYZ->getNumVectors();
    size_t problemDimension_ = XYZ.size();

    zoltanObj_->Set_Num_Obj_Fn(GetLocalNumberOfRows,(void *) &*A);
    zoltanObj_->Set_Obj_List_Fn(GetLocalNumberOfNonzeros,(void *) &*A);
    zoltanObj_->Set_Num_Geom_Fn(GetProblemDimension, (void *) &problemDimension_);
    //~~ zoltanObj_->Set_Geom_Multi_Fn(GetProblemGeometry, (void *) &*XYZ);
    zoltanObj_->Set_Geom_Multi_Fn(GetProblemGeometry, (void *) &XYZ);

    // Data pointers that Zoltan requires.
    ZOLTAN_ID_PTR import_gids = NULL;  // Global nums of objs to be imported   
    ZOLTAN_ID_PTR import_lids = NULL;  // Local indices to objs to be imported 
    int   *import_procs = NULL;        // Proc IDs of procs owning objs to be imported.
    int   *import_to_part = NULL;      // Partition #s to which imported objs should be assigned.
    ZOLTAN_ID_PTR export_gids = NULL;  // Global nums of objs to be exported   
    ZOLTAN_ID_PTR export_lids = NULL;  // local indices to objs to be exported 
    int   *export_procs = NULL;        // Proc IDs of destination procs for objs to be exported.
    int   *export_to_part = NULL;      // Partition #s for objs to be exported.
    int   num_imported;                // Number of objs to be imported.          
    int   num_exported;                // Number of objs to be exported.          
    int   newDecomp;                   // Flag indicating whether the decomposition has changed 
    int   num_gid_entries;             // Number of array entries in a global ID.  
    int   num_lid_entries;

    rv = zoltanObj_->LB_Partition(newDecomp, num_gid_entries, num_lid_entries,
                                  num_imported, import_gids, import_lids, import_procs, import_to_part,
                                  num_exported, export_gids, export_lids, export_procs, export_to_part);
    if (rv == ZOLTAN_FATAL) {
      throw(Exceptions::RuntimeError("Zoltan::LB_Partition() returned error code"));
    }

    //TODO check that A's row map is 1-1.  Zoltan requires this.
    RCP<const Map> rowMap = A->getRowMap();
    int mypid = rowMap->getComm()->getRank();
    RCP<Xpetra::Vector<GO,LO,GO,NO> > decomposition;
    if (newDecomp) {
      decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(rowMap,false); //Don't bother initializing, as this will just be overwritten.
      ArrayRCP<GO> decompEntries = decomposition->getDataNonConst(0);
      for (typename ArrayRCP<GO>::iterator i = decompEntries.begin(); i != decompEntries.end(); ++i)
        *i = mypid;
      LocalOrdinal blockSize = A->GetFixedBlockSize();
      for (int i=0; i< num_exported; ++i) {
        LO localEl = rowMap->getLocalElement(export_gids[i]);
        int partNum = export_to_part[i];
        for (LO j=0; j<blockSize; ++j)
          decompEntries[ localEl + j ] = partNum;
          //decompEntries[ rowMap->getLocalElement(export_gids[i]) + j ] = export_to_part[i];
      }
    } else {
      //Running on one processor, so decomposition is the trivial one, all zeros.
      decomposition = Xpetra::VectorFactory<GO,LO,GO,NO>::Build(rowMap,true);
    } // if (newDecomp) ... else
    level.Set<RCP<Xpetra::Vector<GO,LO,GO,NO> > >("Partition",decomposition,this);

    zoltanObj_->LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_to_part);
    zoltanObj_->LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_to_part);

  } //Build()

  //-------------------------------------------------------------------------------------------------------------
  // GetLocalNumberOfRows
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  int ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetLocalNumberOfRows(void *data, int *ierr)
  {
    if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return -1;
    }
    *ierr = ZOLTAN_OK;
    //TODO is there a safer way to cast?
    //Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> *A = (Operator*) data;
    Operator *A = (Operator*) data;
    LocalOrdinal blockSize = A->GetFixedBlockSize(); //FIXME
    if (blockSize==0) throw(Exceptions::RuntimeError("MueLu::Zoltan : Operator has block size 0."));
    return (A->getRowMap()->getNodeNumElements() / blockSize); //FIXME
  } //GetLocalNumberOfRows()

  //-------------------------------------------------------------------------------------------------------------
  // GetLocalNumberOfNonzeros
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetLocalNumberOfNonzeros(void *data, int NumGidEntries, int NumLidEntries, ZOLTAN_ID_PTR gids,
                           ZOLTAN_ID_PTR lids, int wgtDim, float *weights, int *ierr)
  {
    if (data == NULL || NumGidEntries < 1) {
      *ierr = ZOLTAN_FATAL;
      return;
    } else {
      *ierr = ZOLTAN_OK;
    }

    //TODO is there a safer way to cast?
    Operator *A = (Operator*) data;
    RCP<const Map> map = A->getRowMap();
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    LocalOrdinal blockSize = A->GetFixedBlockSize(); //FIXME
    if (blockSize==0) throw(Exceptions::RuntimeError("MueLu::Zoltan : Operator has block size 0."));
    if (blockSize == 1) {
      for (size_t i=0; i<map->getNodeNumElements(); ++i) {
        gids[i] = (ZOLTAN_ID_TYPE) map->getGlobalElement(i);
        A->getLocalRowView(i,cols,vals);
        weights[i] = cols.size();
      }
    } else {
      LocalOrdinal numBlocks = A->getRowMap()->getNodeNumElements() / blockSize;
      std::set<LocalOrdinal> uniqueColsInBlockRow;
      Teuchos::ArrayView<LO> nonconstCols =Teuchos::av_const_cast<LO>(cols);
      for (LocalOrdinal i=0; i<numBlocks; ++i) {
        gids[i] = (ZOLTAN_ID_TYPE) map->getGlobalElement(i*blockSize);
        for (LocalOrdinal j=i*blockSize; j<(i+1)*blockSize; ++j) {
          A->getLocalRowView(j,cols,vals);
          LO *tt = nonconstCols.getRawPtr();
          uniqueColsInBlockRow.insert(tt,tt+nonconstCols.size());  //yes, cols.size() is correct, one past last entry
        }
        weights[i] = uniqueColsInBlockRow.size();
        uniqueColsInBlockRow.clear();
      } //for (size_t i=0; i<numBlocks; ++i)
    }

  } //GetLocalNumberOfNonzeros()

  //-------------------------------------------------------------------------------------------------------------
  // GetProblemDimension
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  int ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetProblemDimension(void *data, int *ierr)
  {
    //TODO is there a safer way to cast?
    int dim = *((int*)data);
    *ierr = ZOLTAN_OK; /* set error flag */
    return(dim);
  } //GetProblemDimension

  //-------------------------------------------------------------------------------------------------------------
  // GetProblemGeometry
  //-------------------------------------------------------------------------------------------------------------

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>     
  void ZoltanInterface<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  GetProblemGeometry(void *data, int numGIDEntries, int numLIDEntries, int numObjectIDs, 
                     ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int dim, double *coordinates, int *ierr)
  {
    if (data == NULL) {
      *ierr = ZOLTAN_FATAL;
      return;
    }

    //TODO is there a safer way to cast?
    //~~ MultiVector *XYZ = (MultiVector*) data;
    Array<ArrayRCP<SC> > * XYZpt = (Array<ArrayRCP<SC> > *)data;
    Array<ArrayRCP<SC> > & XYZ   = *XYZpt;

    //~~ if (dim != (int) XYZ->getNumVectors()) {
    if (dim != (int) XYZ.size()) { //FIXME: cast to size_t instead?
      //FIXME I'm assuming dim should be 1,2,or 3 coming in?!
      *ierr = ZOLTAN_FATAL;
      return;
    }

    //~ assert(numObjectIDs == XYZ->getLocalLength());
    for(int j=0; j<dim; j++) {
      assert(numObjectIDs == XYZ[j].size()); //FIXME: TEST_FOR_EXCEPTION instead?      
    }

    /*~~
    ArrayRCP<ArrayRCP<const SC> > XYZdata(dim);
    for (int j=0; j<dim; ++j) XYZdata[j] = XYZ->getData(j);
    for (size_t i=0; i<XYZ->getLocalLength(); ++i) {
      for (int j=0; j<dim; ++j) {
        coordinates[i*dim+j] = (double) XYZdata[j][i];
      }
    }
    */

    for (size_t i=0; i<(size_t)XYZ[0].size(); ++i) { //FIXME cast OK?
      for (int j=0; j<dim; ++j) {
        coordinates[i*dim+j] = (double) XYZ[j][i];
      }
    }

    *ierr = ZOLTAN_OK;

  } //GetProblemGeometry

} //namespace MueLu

#endif //if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

#endif // MUELU_ZOLTANINTERFACE_DEF_HPP
