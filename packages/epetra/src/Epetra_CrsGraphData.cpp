
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_CrsGraphData.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
//#include "Epetra_ConfigDefs.h" //DATA_DEBUG

//=============================================================================
Epetra_CrsGraphData::Epetra_CrsGraphData(Epetra_DataAccess CV, const Epetra_BlockMap& RowMap, bool StaticProfile)
  // maps
  : RowMap_(RowMap),
    ColMap_(RowMap),
    DomainMap_(RowMap),
    RangeMap_(RowMap),
    // importer & exporter
    Importer_(0),
    Exporter_(0),
    // booleans
    HaveColMap_(false),
    Filled_(false),
    Allocated_(false),
    Sorted_(false),
    StorageOptimized_(false),
    NoRedundancies_(false),
    IndicesAreGlobal_(false),
    IndicesAreLocal_(false),
    IndicesAreContiguous_(false),
    LowerTriangular_(true),
    UpperTriangular_(true),
    NoDiagonal_(true),
    GlobalConstantsComputed_(false),
    StaticProfile_(StaticProfile),
    // ints
    IndexBase_(RowMap.IndexBase()),
    NumGlobalEntries_(0),
    NumGlobalBlockRows_(RowMap.NumGlobalElements()),
    NumGlobalBlockCols_(NumGlobalBlockRows_),
    NumGlobalBlockDiagonals_(0),
    NumMyEntries_(0),
    NumMyBlockRows_(RowMap.NumMyElements()),
    NumMyBlockCols_(NumMyBlockRows_),
    NumMyBlockDiagonals_(0),
    MaxRowDim_(RowMap.MaxElementSize()),
    MaxColDim_(MaxRowDim_),
    GlobalMaxRowDim_(RowMap.MaxElementSize()),
    GlobalMaxColDim_(GlobalMaxRowDim_),
    MaxNumNonzeros_(0),
    GlobalMaxNumNonzeros_(0),
    NumGlobalNonzeros_(0),
    NumGlobalRows_(RowMap.NumGlobalPoints()),
    NumGlobalCols_(NumGlobalRows_),
    NumGlobalDiagonals_(0),
    NumMyNonzeros_(0),
    NumMyRows_(RowMap.NumMyPoints()),
    NumMyCols_(NumMyRows_),
    NumMyDiagonals_(0),	
    MaxNumIndices_(0),
    GlobalMaxNumIndices_(0),
    Indices_(new int *[NumMyBlockRows_]),
    NumAllocatedIndicesPerRow_(0),
    NumIndicesPerRow_(0),
    IndexOffset_(0),
    All_Indices_(0),
    CV_(CV)
{
  //cout << "--CRSGD created(rowmap ctr), addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_CrsGraphData::Epetra_CrsGraphData(Epetra_DataAccess CV, 
					 const Epetra_BlockMap& RowMap, 
					 const Epetra_BlockMap& ColMap, bool StaticProfile)
  // maps
  : RowMap_(RowMap),
    ColMap_(ColMap),
    DomainMap_(ColMap),
    RangeMap_(RowMap),
    // importer & exporter
    Importer_(0),
    Exporter_(0),
    // booleans
    HaveColMap_(true),
    Filled_(false),
    Allocated_(false),
    Sorted_(false),
    StorageOptimized_(false),
    NoRedundancies_(false),
    IndicesAreGlobal_(false),
    IndicesAreLocal_(false),
    IndicesAreContiguous_(false),
    LowerTriangular_(true),
    UpperTriangular_(true),
    NoDiagonal_(true),
    GlobalConstantsComputed_(false),
    StaticProfile_(StaticProfile),
    // ints
    IndexBase_(RowMap.IndexBase()),
    NumGlobalEntries_(0),
    NumGlobalBlockRows_(RowMap.NumGlobalElements()),
    NumGlobalBlockCols_(ColMap.NumGlobalElements()),
    NumGlobalBlockDiagonals_(0),
    NumMyEntries_(0),
    NumMyBlockRows_(RowMap.NumMyElements()),
    NumMyBlockCols_(ColMap.NumMyElements()),
    NumMyBlockDiagonals_(0),
    MaxRowDim_(RowMap.MaxElementSize()),
    MaxColDim_(ColMap.MaxElementSize()),
    GlobalMaxRowDim_(RowMap.MaxElementSize()),
    GlobalMaxColDim_(ColMap.MaxElementSize()),
    MaxNumNonzeros_(0),
    GlobalMaxNumNonzeros_(0),
    NumGlobalNonzeros_(0),
    NumGlobalRows_(RowMap.NumGlobalPoints()),
    NumGlobalCols_(ColMap.NumGlobalPoints()),
    NumGlobalDiagonals_(0),
    NumMyNonzeros_(0),
    NumMyRows_(RowMap.NumMyPoints()),
    NumMyCols_(ColMap.NumMyPoints()),
    NumMyDiagonals_(0),	
    MaxNumIndices_(0),
    GlobalMaxNumIndices_(0),
    Indices_(new int *[NumMyBlockRows_]),
    NumAllocatedIndicesPerRow_(0),
    NumIndicesPerRow_(0),
    IndexOffset_(0),
    All_Indices_(0),
    CV_(CV)
{
  //cout << "--CRSGD created(rowmap&colmap ctr), addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_CrsGraphData::~Epetra_CrsGraphData() {

  if(Indices_ != 0 && !StorageOptimized_) {
    for (int i=0; i<NumMyBlockRows_; i++) {
      if (Indices_[i]!=0 && CV_==Copy && !StaticProfile_) 
	delete [] Indices_[i]; 
      Indices_[i] = 0;
    } 
    delete[] Indices_;
    Indices_ = 0;
  }

  if(Importer_ != 0) {
    delete Importer_;
    Importer_ = 0;
  }
  if(Exporter_ != 0) {
    delete Exporter_;
    Importer_ = 0;
  }

  NumMyBlockRows_ = 0;	// are these needed?
  Filled_ = false;      // they're about to go out of scope, after all
  Allocated_ = false;

  //cout << "--CRSGD destroyed, addr: " << this << endl; //DATA_DEBUG
}

//==========================================================================
int Epetra_CrsGraphData::MakeImportExport() {
  // Create Import object for use by matrix classes.    This is only needed if ColMap and DomainMap are different
  if (!ColMap_.SameAs(DomainMap_)) {
    if (Importer_ != 0) {
      delete Importer_;
      Importer_ = 0;
    }
    Importer_ = new Epetra_Import(ColMap_, DomainMap_);
  }
  
  // Now see if we need to define an export map.  This is only needed if RowMap and RangeMap are different
  if (!RowMap_.SameAs(RangeMap_)) {
    if (Exporter_ != 0) {
      delete Exporter_;
      Exporter_ = 0;
    }
    Exporter_ = new Epetra_Export(RowMap_, RangeMap_); // Create Export object. 
  }
   
  return(0);
}

//==========================================================================
int Epetra_CrsGraphData::ReAllocateAndCast(char*& UserPtr, int& Length, const int IntPacketSizeTimesNumTrans) {
  if(IntPacketSizeTimesNumTrans > Length) {
    if(Length > 0) 
      delete[] UserPtr;
    Length = IntPacketSizeTimesNumTrans;
    int* newPtr = new int[Length];
    UserPtr = reinterpret_cast<char*> (newPtr);
  }
  return(0);
}

//==========================================================================
void Epetra_CrsGraphData::Print(ostream& os, int level) const {
  bool four_bit = (level >= 4);      // 4-bit = BlockMaps
  bool two_bit = ((level % 4) >= 2); // 2-bit = Indices
  bool one_bit = ((level % 2) == 1); // 1-bit = Everything else

  os << "\n***** CrsGraphData (output level " << level << ") *****" << endl;

  if(four_bit) {
    os << "RowMap_:\n" << RowMap_ << endl;
    os << "ColMap_:\n" << ColMap_ << endl;
    os << "DomainMap_:\n" << DomainMap_ << endl;
    os << "RangeMap_:\n" << RangeMap_ << endl;
  }
	
  if(one_bit) {
    os.width(26); os << "HaveColMap_: "              << HaveColMap_;
    os.width(25); os << "Filled_: "                  << Filled_;
    os.width(25); os << "Allocated_: "               << Allocated_;
    os.width(25); os << "Sorted_: "                  << Sorted_ << endl;
    os.width(26); os << "StorageOptimized_: "        << StorageOptimized_;
    os.width(25); os << "NoRedundancies_: "          << NoRedundancies_;
    os.width(25); os << "IndicesAreGlobal_: "        << IndicesAreGlobal_;
    os.width(25); os << "IndicesAreLocal_: "         << IndicesAreLocal_ << endl;
    os.width(26); os << "IndicesAreContiguous_: "    << IndicesAreContiguous_;
    os.width(25); os << "LowerTriangular_: "         << LowerTriangular_;
    os.width(25); os << "UpperTriangular_: "         << UpperTriangular_;
    os.width(25); os << "NoDiagonal_: "              << NoDiagonal_ << endl;
    os.width(25); os << "GlobalConstantsComputed_: " << GlobalConstantsComputed_ << endl;
    os.width(25); os << "StaticProfile_: " << StaticProfile_ << endl << endl;
		
    os.width(10); os << "NGBR_: " << NumGlobalBlockRows_;
    os.width(10); os << "NGBC_: " << NumGlobalBlockCols_;
    os.width(10); os << "NGBD_: " << NumGlobalBlockDiagonals_;
    os.width(10); os << "NGE_: "  << NumGlobalEntries_;
    os.width(10); os << "NGR_: "  << NumGlobalRows_;
    os.width(10); os << "NGC_: "  << NumGlobalCols_;
    os.width(10); os << "NGD_: "  << NumGlobalDiagonals_;
    os.width(10); os << "NGN_: "  << NumGlobalNonzeros_;
    os.width(10); os << "IB_: "   << IndexBase_ << endl;
    os.width(10); os << "GMRD_: " << GlobalMaxRowDim_;
    os.width(11); os << "GMCD_: " << GlobalMaxColDim_;
    os.width(11); os << "GMNI_: " << GlobalMaxNumIndices_;
    os.width(11); os << "NMBR_: " << NumMyBlockRows_;
    os.width(10); os << "NMBC_: " << NumMyBlockCols_;
    os.width(10); os << "NMBD_: " << NumMyBlockDiagonals_;
    os.width(10); os << "NME_: "  << NumMyEntries_;
    os.width(10); os << "NMR_: "  << NumMyRows_;
    os.width(10); os << "CV_: " << CV_ << endl;
    os.width(10); os << "NMC_: "  << NumMyCols_;
    os.width(10); os << "NMD_: "  << NumMyDiagonals_;
    os.width(10); os << "NMN_: "  << NumMyNonzeros_;
    os.width(10); os << "MRD_: "  << MaxRowDim_;
    os.width(11); os << "MCD_: "  << MaxColDim_;
    os.width(11); os << "MNI_: "  << MaxNumIndices_;
    os.width(11); os << "MNN_: "  << MaxNumNonzeros_;
    os.width(11); os << "GMNN_: " << GlobalMaxNumNonzeros_;
    os.width(11); os << "RC: " << ReferenceCount() << endl << endl;
		
    os << "NIPR_: " << NumIndicesPerRow_ << endl;
    os << "NAIPR_: " << NumAllocatedIndicesPerRow_ << endl;
    os << "IndexOffset_: " << IndexOffset_ << endl;
    os << "All_Indices_: " << All_Indices_ << endl;
  }
		
  if(two_bit) {
    os << "Indices_: " << Indices_ << endl;
    if(Indices_ != 0) {
      for(int i = 0; i < NumMyBlockRows_; i++) {
	os << "Indices_[" << i << "]: (" << Indices_[i] << ") ";
	if(Indices_[i] != 0) {
	  for(int j = 0; j < NumAllocatedIndicesPerRow_[i]; j++)
	    os << Indices_[i][j] << " ";
	}
	os << endl;
      }
    }
  }
	
  os << "***** End CrsGraphData *****" << endl;
}
