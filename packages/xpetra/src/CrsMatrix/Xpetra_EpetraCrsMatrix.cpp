#include <Teuchos_Array.hpp>
#include "Xpetra_EpetraCrsMatrix.hpp"

namespace Xpetra {

  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }
  
  // TODO: convert array size_t to int
  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist) {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }
  
  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist)
    : mtx_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), toEpetra(colMap), maxNumEntriesPerRow, toEpetra(pftype)))) { }
  
  EpetraCrsMatrix::EpetraCrsMatrix(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype, const Teuchos::RCP< Teuchos::ParameterList > &plist) {
    Teuchos::Array<int> numEntriesPerRowToAlloc(NumEntriesPerRowToAlloc.begin(), NumEntriesPerRowToAlloc.end()); // convert array of "size_t" to array of "int"
    mtx_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, toEpetra(rowMap), toEpetra(colMap), numEntriesPerRowToAlloc.getRawPtr(), toEpetra(pftype)));
  }

  void EpetraCrsMatrix::insertGlobalValues(int globalRow, const ArrayView<const int> &cols, const ArrayView<const double> &vals) { 
    XPETRA_ERR_CHECK(mtx_->InsertGlobalValues(globalRow, vals.size(), vals.getRawPtr(), cols.getRawPtr())); 
  }

  //TODO: throw same exception as Tpetra
  void EpetraCrsMatrix::getLocalRowCopy(int LocalRow, const ArrayView<int> &Indices, const ArrayView<double> &Values, size_t &NumEntries) const { 
    int numEntries=-1;
    XPETRA_ERR_CHECK(mtx_->ExtractMyRowCopy(LocalRow, Indices.size(), numEntries, Values.getRawPtr(), Indices.getRawPtr()));
    NumEntries = numEntries;
  }

  void EpetraCrsMatrix::getGlobalRowView(int GlobalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const { 
    int      numEntries;
    double * eValues;
    int    * eIndices;
      
    XPETRA_ERR_CHECK(mtx_->ExtractGlobalRowView(GlobalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  void EpetraCrsMatrix::getLocalRowView(int LocalRow, ArrayView<const int> &indices, ArrayView<const double> &values) const { 
    int      numEntries;
    double * eValues;
    int    * eIndices;
      
    XPETRA_ERR_CHECK(mtx_->ExtractMyRowView(LocalRow, numEntries, eValues, eIndices));
    if (numEntries == 0) { eValues = NULL; eIndices = NULL; } // Cf. TEUCHOS_TEST_FOR_EXCEPT( p == 0 && size_in != 0 ) in Teuchos ArrayView constructor.

    indices = ArrayView<const int>(eIndices, numEntries);
    values  = ArrayView<const double>(eValues, numEntries);
  }

  void EpetraCrsMatrix::apply(const MultiVector<double,int,int> & X, MultiVector<double,int,int> &Y, Teuchos::ETransp mode, double alpha, double beta) const { 
    TEUCHOS_TEST_FOR_EXCEPTION((alpha != 1) || (beta != 0), Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraCrsMatrix.multiply() only accept alpha==1 and beta==0");
      
    XPETRA_DYNAMIC_CAST(const EpetraMultiVector, X, eX, "Xpetra::EpetraCrsMatrix->apply() only accept Xpetra::EpetraMultiVector as input arguments.");
    XPETRA_DYNAMIC_CAST(      EpetraMultiVector, Y, eY, "Xpetra::EpetraCrsMatrix->apply() only accept Xpetra::EpetraMultiVector as input arguments.");

    TEUCHOS_TEST_FOR_EXCEPTION((mode != Teuchos::NO_TRANS) && (mode != Teuchos::TRANS), Xpetra::Exceptions::NotImplemented, "Xpetra::EpetraCrsMatrix->apply() only accept mode == NO_TRANS or mode == TRANS");
    bool eTrans = toEpetra(mode);

    // /!\ UseTranspose value
    TEUCHOS_TEST_FOR_EXCEPTION(mtx_->UseTranspose(), Xpetra::Exceptions::NotImplemented, "An exception is throw to let you know that Xpetra::EpetraCrsMatrix->apply() do not take into account the UseTranspose() parameter of Epetra_CrsMatrix.");
      
    XPETRA_ERR_CHECK(mtx_->Multiply(eTrans, *eX.getEpetra_MultiVector(), *eY.getEpetra_MultiVector()));
  }

  std::string EpetraCrsMatrix::description() const { 
    // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
    std::ostringstream oss;
    //TODO: oss << DistObject<char, LocalOrdinal,GlobalOrdinal>::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
      
  } 
    
  void EpetraCrsMatrix::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { 

    // This implementation come from Tpetra_CrsMatrix_def.hpp (without modification)
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
      numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,11) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    // 
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl; 
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << getGlobalNumDiags() << std::endl;
        out << "Global max number of entries = " << getGlobalMaxNumRowEntries() << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        getRowMap()->describe(out,vl);
        //
        if (getColMap() != null) {
          if (getColMap() == getRowMap()) {
            if (myImageID == 0) out << "\nColumn map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nColumn map: " << std::endl;
            getColMap()->describe(out,vl);
          }
        }
        if (getDomainMap() != null) {
          if (getDomainMap() == getRowMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else if (getDomainMap() == getColMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nDomain map: " << std::endl;
            getDomainMap()->describe(out,vl);
          }
        }
        if (getRangeMap() != null) {
          if (getRangeMap() == getDomainMap()) {
            if (myImageID == 0) out << "\nRange map is domain map." << std::endl;
          }
          else if (getRangeMap() == getRowMap()) {
            if (myImageID == 0) out << "\nRange map is row map." << std::endl;
          }
          else {
            if (myImageID == 0) out << "\nRange map: " << std::endl;
            getRangeMap()->describe(out,vl);
          }
        }
        if (myImageID == 0) out << std::endl;
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl;
            // TODO: need a graph
            //               if (staticGraph_->indicesAreAllocated() == false) {
            //                 out << "Node not allocated" << std::endl;
            //               }
            //               else {
            //                 out << "Node number of allocated entries = " << staticGraph_->getNodeAllocationSize() << std::endl;
            //               }

            // TMP:
            //            const Epetra_CrsGraph & staticGraph_ = mtx_->Graph();
            // End of TMP

            out << "Node number of entries = " << getNodeNumEntries() << std::endl;
            if (isFillComplete()) {
              out << "Node number of diagonals = " << getNodeNumDiags() << std::endl;
            }
            out << "Node max number of entries = " << getNodeMaxNumRowEntries() << std::endl;
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row" 
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << std::setw(width) << "(Index,Value)";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              const size_t nE = getNumEntriesInLocalRow(r);
              GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
              out << std::setw(width) << myImageID 
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getGlobalRowView(gid,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << rowinds[j]
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getLocalRowView(r,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << getColMap()->getGlobalElement(rowinds[j]) 
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
    
  }

  // TODO: use toEpetra()    
  void EpetraCrsMatrix::doImport(const DistObject<char, int, int> &source, 
                                 const Import<int, int> &importer, CombineMode CM) {

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, source, tSource, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraCrsMatrix::doExport(const DistObject<char, int, int> &dest,
                                 const Import<int, int>& importer, CombineMode CM) {
      
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, dest, tDest, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM)); 
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

  void EpetraCrsMatrix::doImport(const DistObject<char, int, int> &source,
                                 const Export<int, int>& exporter, CombineMode CM) {

    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, source, tSource, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tSource.getEpetra_CrsMatrix();
    int err = mtx_->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");

  }

  void EpetraCrsMatrix::doExport(const DistObject<char, int, int> &dest,
                                 const Export<int, int>& exporter, CombineMode CM) {
      
    XPETRA_DYNAMIC_CAST(const EpetraCrsMatrix, dest, tDest, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraCrsMatrix as input arguments.");
    XPETRA_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Xpetra::EpetraCrsMatrix::doImport only accept Xpetra::EpetraImport as input arguments.");

    RCP<const Epetra_CrsMatrix> v = tDest.getEpetra_CrsMatrix();
    int err = mtx_->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM)); 
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
  }

}
