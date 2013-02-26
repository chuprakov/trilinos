//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include <EpetraExt_ConfigDefs.h>
#include <EpetraExt_MatrixMatrix.h>

#include <EpetraExt_MMHelpers.h>

#include <EpetraExt_Transpose_RowMatrix.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Directory.h>
#include <Epetra_HashTable.h>
#include <Epetra_Distributor.h>
#include <Epetra_IntSerialDenseVector.h>

#ifdef HAVE_VECTOR
#include <vector>
#endif

#ifdef HAVE_MPI
#include <Epetra_MpiDistributor.h>
#endif


#include <Teuchos_TimeMonitor.hpp>

namespace EpetraExt {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static inline int C_estimate_nnz(const Epetra_CrsMatrix & A, const Epetra_CrsMatrix &B){
  // Follows the NZ estimate in ML's ml_matmatmult.c
  int Aest=(A.NumMyRows()>0)? A.MaxNumEntries()/A.NumMyRows():100;
  int Best=(B.NumMyRows()>0)? B.MaxNumEntries()/B.NumMyRows():100;

  int nnzperrow=(int)(sqrt(Aest) + sqrt(Best) - 1);
  nnzperrow*=nnzperrow;
 
  return (int)(A.NumMyRows()*nnzperrow*0.75 + 100);
 }

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
// Template params are <PID,GID>
static inline bool lessthan12(std::pair<int,int> i, std::pair<int,int> j){
  return ((i.first<j.first) || (i.first==j.first && i.second <j.second));
}

#ifdef LIGHTWEIGHT_MATRIX
int aztecoo_and_ml_compatible_map_union(const Epetra_CrsMatrix &B, const LightweightCrsMatrix &Bimport, Epetra_Map*& unionmap)
#else
int aztecoo_and_ml_compatible_map_union(const Epetra_CrsMatrix &B, const Epetra_CrsMatrix &Bimport, Epetra_Map*& unionmap)
#endif
{
#ifdef HAVE_MPI
  // As per Mike Heroux, AztecOO is very picky about how matrix column maps are ordered.  Specifically:
  // 1) Local entries (row is owned by the processor) must come first.
  // 2) Off-processor entries are ordered by processor
  //
  // ML and IFPACK add the following constraint for blocking reasons:
  // 3) Unknown ordering within a processor must be preserved. That is, the local numbering of b is less than c if the global numbering of b 
  //    is less than c and if both b and c are owned by the same processor. 

  // This routine will merge two colmaps into "unionmap" that will meet those criteria ordering criteria


  // NTS: Should this be re-implemented in the same way as MakeColMap?
  
  Epetra_Util util;
  int MyPID = B.Comm().MyPID();

  // Since we're getting called, we know we have to be using an MPI implementation of Epetra.  Which means we should have an MpiDistributor for both B and Bimport
#ifdef LIGHTWEIGHT_MATRIX
  if(!B.Importer()) EPETRA_CHK_ERR(-1);
#else
  if(!B.Importer() || !Bimport.Importer()) EPETRA_CHK_ERR(-1);
#endif

  // Get the (PID,GID) pairs
  std::vector< std::pair<int,int> > Bpgids, Bimportpgids,Mpgids;
  Epetra_MpiDistributor *Bdist=dynamic_cast<Epetra_MpiDistributor*>(&B.Importer()->Distributor());
  if(!Bdist) EPETRA_CHK_ERR(-2);

  util.GetPidGidPairs(*B.Importer(),Bpgids,true);
#ifdef LIGHTWEIGHT_MATRIX
  int Nimport=Bimport.ColMap_.NumMyElements();
  if(Nimport != (int) Bimport.ColMapOwningPIDs_.size()) {
    //    printf("[%d] Nimport = %d  Bimport.ColMapOwningPIDs_.size() = %d\n",B.Comm().MyPID(),Nimport, Bimport.ColMapOwningPIDs_.size());
    EPETRA_CHK_ERR(-21);

  }
  Bimportpgids.resize(Nimport);
  for(int i=0;i<Nimport;i++){
    Bimportpgids[i].first = (Bimport.ColMapOwningPIDs_[i]==MyPID)?(-1):(Bimport.ColMapOwningPIDs_[i]);
    Bimportpgids[i].second= Bimport.ColMap_.GID(i);
  }
#else
  Epetra_MpiDistributor *Bidist=dynamic_cast<Epetra_MpiDistributor*>(&Bimport.Importer()->Distributor());
  if(!Bidist) EPETRA_CHK_ERR(-2);
  util.GetPidGidPairs(*Bimport.Importer(),Bimportpgids,true);
#endif
 
  // Sort the (PID,GID) pairs by PID then GID (pre-merge)
  std::sort(Bpgids.begin(),Bpgids.end(),lessthan12);
  std::sort(Bimportpgids.begin(),Bimportpgids.end(),lessthan12);

  // Make sure the merge list is big enough
  Mpgids.resize(Bpgids.size()+Bimportpgids.size());

  // Merge, then unique the sorted lists.
  // The default == for std::pair is OK, the default < is not
  std::merge(Bpgids.begin(),Bpgids.end(),Bimportpgids.begin(),Bimportpgids.end(),Mpgids.begin(),lessthan12);
  std::vector< std::pair<int,int> >::iterator last_el=std::unique(Mpgids.begin(),Mpgids.end());
  int NumUniqueCols=last_el-Mpgids.begin();

  // Flag all the colids actually in my domain map, copy all into vector for Map constructor
  const Epetra_Map& DomainMap=B.DomainMap();
  std::vector<bool> LocalGCIDs(DomainMap.NumMyElements(),false);
  std::vector<int> MyGCIDs(NumUniqueCols);

  for(int i=0;i<NumUniqueCols; i++){
    MyGCIDs[i]=Mpgids[i].second;
    if(Mpgids[i].first==-1){
      int LID=DomainMap.LID(Mpgids[i].second);
      if(LID!=-1) LocalGCIDs[LID]=true;
      else EPETRA_CHK_ERR(-3);
    }
  }

  // Re-order the local unknowns to be in the same order as the domain map
  // This code is more or less what FillComplete does.
  for(int i=0;i<DomainMap.NumMyElements();i++){
    if(LocalGCIDs[i]){
      MyGCIDs[i]=DomainMap.GID(i);
    }
  }

  // Make the map
  unionmap=new Epetra_Map(-1,NumUniqueCols,&MyGCIDs[0],B.ColMap().IndexBase(),B.Comm());
  return 0;
#else
  return -1;
#endif
}



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
// Provide a "resize" operation for double*'s. 
inline void resize_doubles(int nold,int nnew,double*& d){
  if(nnew > nold){
    double *tmp = new double[nnew];
    for(int i=0; i<nold; i++)
      tmp[i]=d[i];
    delete [] d;
    d=tmp;
  }
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int  MatrixMatrix::mult_A_B_newmatrix(const Epetra_CrsMatrix & A,
		       const Epetra_CrsMatrix & B,
		       CrsMatrixStruct& Bview,
		       std::vector<int> & Bcol2Ccol,
		       std::vector<int> & Bimportcol2Ccol,
		       Epetra_CrsMatrix& C){

  // *****************************
  // Improved Parallel Gustavson in Local IDs
  // *****************************
  const Epetra_Map * colmap_C = &(C.ColMap());

  int m=A.NumMyRows();
  int n=colmap_C->NumMyElements();
  int i,j,k;

  // DataPointers for A
  int *Arowptr, *Acolind;
  double *Avals;
  EPETRA_CHK_ERR(A.ExtractCrsDataPointers(Arowptr,Acolind,Avals));

  // MemorySetup: If somebody else is sharing this C's graphdata, make a new one.
  // This is needed because I'm about to walk all over the CrsGrapData...
  C.ExpertMakeUniqueCrsGraphData();

  // The status array will contain the index into colind where this entry was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with -1's indicating that there are no entries yet.
  std::vector<int> c_status(n, -1);

  // Classic csr assembly (low memory edition)
  int CSR_alloc=C_estimate_nnz(A,B);
  int CSR_ip=0,OLD_ip=0;
  Epetra_IntSerialDenseVector & CSR_rowptr = C.ExpertExtractIndexOffset();
  Epetra_IntSerialDenseVector & CSR_colind = C.ExpertExtractIndices();  
  double *&                     CSR_vals   = C.ExpertExtractValues();

  CSR_rowptr.Resize(m+1);
  CSR_colind.Resize(CSR_alloc);
  resize_doubles(0,CSR_alloc,CSR_vals);

  // Static Profile stuff
  std::vector<int> NumEntriesPerRow(m);

  // For each row of A/C
  for(i=0; i<m; i++){			       
    CSR_rowptr[i]=CSR_ip;

    for(k=Arowptr[i]; k<Arowptr[i+1]; k++){
      int Ak=Acolind[k];
      double Aval = Avals[k];
      if(Aval==0) continue;

      if(!Bview.remote[Ak]){
	// Local matrix
	int* Bcol_inds = Bview.indices[Ak];
	double* Bvals_k= Bview.values[Ak];
	for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
	  int Cj=Bcol2Ccol[Bcol_inds[j]];

	  if(c_status[Cj]<OLD_ip){
	    // New entry
	    c_status[Cj]=CSR_ip;
	    CSR_colind[CSR_ip]=Cj;
	    CSR_vals[CSR_ip]=Aval*Bvals_k[j];
	    CSR_ip++;
	  }
	  else
	    CSR_vals[c_status[Cj]]+=Aval*Bvals_k[j];
	}
      }
      else{
	// Remote matrix
	int* Bcol_inds = Bview.indices[Ak];
	double* Bvals_k= Bview.values[Ak];
	
	for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
	  int Cj=Bimportcol2Ccol[Bcol_inds[j]];

	  if(c_status[Cj]<OLD_ip){
	    // New entry
	    c_status[Cj]=CSR_ip;
	    CSR_colind[CSR_ip]=Cj;
	    CSR_vals[CSR_ip]=Aval*Bvals_k[j];
	    CSR_ip++;
	  }
	  else
	    CSR_vals[c_status[Cj]]+=Aval*Bvals_k[j];
	}
      }
    }
    NumEntriesPerRow[i]=CSR_ip-CSR_rowptr[i];

    // Resize for next pass if needed
    if(CSR_ip + n > CSR_alloc){
      resize_doubles(CSR_alloc,2*CSR_alloc,CSR_vals);
      CSR_alloc*=2;
      CSR_colind.Resize(CSR_alloc);
    }
    OLD_ip=CSR_ip;
  }
  CSR_rowptr[m]=CSR_ip;

#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=M.getNewTimer("Final Sort");
  mtime->start();
#endif

  // Sort the entries
  sort_crs_entries(m, &CSR_rowptr[0], &CSR_colind[0], &CSR_vals[0]);

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=M.getNewTimer("ESFC");
  mtime->start();
#endif

  // Update the CrsGraphData
  C.ExpertStaticFillComplete(B.DomainMap(),A.RangeMap());

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif

  return 0;
}








int MatrixMatrix::mult_A_B_reuse(const Epetra_CrsMatrix & A,
		   const Epetra_CrsMatrix & B,
		   CrsMatrixStruct& Bview,
		   std::vector<int> & Bcol2Ccol,
		   std::vector<int> & Bimportcol2Ccol,
		   Epetra_CrsMatrix& C){

  // *****************************
  // Improved Parallel Gustavson in Local IDs
  // *****************************
  const Epetra_Map * colmap_C = &(C.ColMap());

  int m=A.NumMyRows();
  int n=colmap_C->NumMyElements();
  int i,j,k;

  // DataPointers for A
  int *Arowptr, *Acolind;
  double *Avals;
  EPETRA_CHK_ERR(A.ExtractCrsDataPointers(Arowptr,Acolind,Avals));


  // DataPointers for C
  int *CSR_rowptr, *CSR_colind;
  double *CSR_vals;
  EPETRA_CHK_ERR(C.ExtractCrsDataPointers(CSR_rowptr,CSR_colind,CSR_vals));


  // The status array will contain the index into colind where this dude was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with -1's indicating that there are no entries yet.
  std::vector<int> c_status(n, -1);

  // Classic csr assembly
  int CSR_alloc=CSR_rowptr[m] - CSR_rowptr[0];
  int CSR_ip=0,OLD_ip=0;

  // For each row of A/C
  for(i=0; i<m; i++){			       
    for(k=Arowptr[i]; k<Arowptr[i+1]; k++){
      int Ak=Acolind[k];

      double Aval = Avals[k];
      if(Aval==0) continue;

      if(!Bview.remote[Ak]){
	// Local matrix
	int* Bcol_inds = Bview.indices[Ak];
	double* Bvals_k= Bview.values[Ak];
	for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
	  int Cj=Bcol2Ccol[Bcol_inds[j]];

	  if(c_status[Cj]<OLD_ip){
	    // New entry
	    if(CSR_ip >= CSR_alloc) EPETRA_CHK_ERR(-13);
	    c_status[Cj]=CSR_ip;
	    CSR_colind[CSR_ip]=Cj;
	    CSR_vals[CSR_ip]=Aval*Bvals_k[j];
	    CSR_ip++;
	  }
	  else
	    CSR_vals[c_status[Cj]]+=Aval*Bvals_k[j];
	}
      }
      else{
	// Remote matrix
	int* Bcol_inds = Bview.indices[Ak];
	double* Bvals_k= Bview.values[Ak];
	
	for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
	  int Cj=Bimportcol2Ccol[Bcol_inds[j]];

	  if(c_status[Cj]<OLD_ip){
	    // New entry
	    if(CSR_ip >= CSR_alloc) EPETRA_CHK_ERR(-14);
	    c_status[Cj]=CSR_ip;
	    CSR_colind[CSR_ip]=Cj;
	    CSR_vals[CSR_ip]=Aval*Bvals_k[j];
	    CSR_ip++;
	  }
	  else
	    CSR_vals[c_status[Cj]]+=Aval*Bvals_k[j];
	}
      }
    }
    OLD_ip=CSR_ip;
  }

  // Sort the entries
  sort_crs_entries(m, &CSR_rowptr[0], &CSR_colind[0], &CSR_vals[0]);

  return 0;
}



//kernel method for computing the local portion of C = A*B
  int mult_A_B_general(const Epetra_CrsMatrix & A,
		       CrsMatrixStruct & Aview,
		       const Epetra_CrsMatrix & B,
		       CrsMatrixStruct& Bview,
		       Epetra_CrsMatrix& C,
		       bool call_FillComplete_on_result)
{
 printf("A*B general\n");

  int C_firstCol = Bview.colMap->MinLID();
  int C_lastCol = Bview.colMap->MaxLID();

  int C_firstCol_import = 0;
  int C_lastCol_import = -1;

  int* bcols = Bview.colMap->MyGlobalElements();
  int* bcols_import = NULL;
  if (Bview.importColMap != NULL) {
    C_firstCol_import = Bview.importColMap->MinLID();
    C_lastCol_import = Bview.importColMap->MaxLID();

    bcols_import = Bview.importColMap->MyGlobalElements();
  }

  int C_numCols = C_lastCol - C_firstCol + 1;
  int C_numCols_import = C_lastCol_import - C_firstCol_import + 1;

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  // Allocate workspace memory
  double* dwork = new double[C_numCols];
  int* iwork = new int[C_numCols];
  int *c_cols=iwork;
  double *c_vals=dwork;
  int *c_index=new int[C_numCols];

  int  i, j, k;
  bool C_filled=C.Filled();


  //To form C = A*B we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc.
  
  // Mark indices as empty w/ -1
  for(k=0;k<C_numCols;k++) c_index[k]=-1;

  //loop over the rows of A.
  for(i=0; i<Aview.numRows; ++i) {
    //only navigate the local portion of Aview... (It's probable that we
    //imported more of A than we need for A*B, because other cases like A^T*B 
    //need the extra rows.)
    if (Aview.remote[i]) {
      continue;
    }

    int* Aindices_i = Aview.indices[i];
    double* Aval_i  = Aview.values[i];
    int global_row = Aview.rowMap->GID(i);

    //loop across the i-th row of A and for each corresponding row
    //in B, loop across colums and accumulate product
    //A(i,k)*B(k,j) into our partial sum quantities C_row_i. In other words,
    //as we stride across B(k,:) we're calculating updates for row i of the
    //result matrix C.

    /* Outline of the revised, ML-inspired algorithm
 
    C_{i,j} = \sum_k A_{i,k} B_{k,j}

    This algorithm uses a "middle product" formulation, with the loop ordering of 
    i, k, j.  This means we compute a row of C at a time, but compute partial sums of 
    each entry in row i until we finish the k loop.

    This algorithm also has a few twists worth documenting.  

    1) The first major twist involves the c_index, c_cols and c_vals arrays.  The arrays c_cols 
    and c_vals store the *local* column index and values accumulator respectively.  These 
    arrays are allocated to a size equal to the max number of local columns in C, namely C_numcols.
    The value c_current tells us how many non-zeros we currently have in this row.
    
    So how do we take a LCID and find the right accumulator?  This is where the c_index array
    comes in.  At the start (and stop) and the i loop, c_index is filled with -1's.  Now 
    whenever we find a LCID in the k loop, we first loop at c_index[lcid].  If this value is
    -1 we haven't seen this entry yet.  In which case we add the appropriate stuff to c_cols
    and c_vals and then set c_index[lcid] to the location of the accumulator (c_current before
    we increment it).  If the value is NOT -1, this tells us the location in the c_vals/c_cols
    arrays (namely c_index[lcid]) where our accumulator lives.

    This trick works because we're working with local ids.  We can then loop from 0 to c_current
    and reset c_index to -1's when we're done, only touching the arrays that have changed.
    While we're at it, we can switch to the global ids so we can call [Insert|SumInto]GlobalValues.
    Thus, the effect of this trick is to avoid passes over the index array.

    2) The second major twist involves handling the remote and local components of B separately.
    (ML doesn't need to do this, because its local ordering scheme is consistent between the local
    and imported components of B.)  Since they have different column maps, they have inconsistent 
    local column ids.  This means the "second twist" won't work as stated on both matrices at the 
    same time.  While this could be handled any number of ways, I have chosen to do the two parts 
    of B separately to make the code easier to read (and reduce the memory footprint of the MMM).
    */

    // Local matrix: Zero Current counts for matrix
    int c_current=0;    

    // Local matrix: Do the "middle product"
    for(k=0; k<Aview.numEntriesPerRow[i]; ++k) {
      int Ak=Aindices_i[k];
      double Aval = Aval_i[k];
      // We're skipping remote entries on this pass.
      if(Bview.remote[Ak] || Aval==0) continue;
      
      int* Bcol_inds = Bview.indices[Ak];
      double* Bvals_k = Bview.values[Ak];

      for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
	int col=Bcol_inds[j];
	if(c_index[col]<0){
	  // We haven't seen this entry before; add it.  (In ML, on
	  // the first pass, you haven't seen any of the entries
	  // before, so they are added without the check.  Not sure
	  // how much more efficient that would be; depends on branch
	  // prediction.  We've favored code readability here.)
	  c_cols[c_current]=col;	      
	  c_vals[c_current]=Aval*Bvals_k[j];
	  c_index[col]=c_current;
	  c_current++;
	}
	else{ 
	  // We've already seen this entry; accumulate it.
	  c_vals[c_index[col]]+=Aval*Bvals_k[j];	    
	}
      }
    }    
    // Local matrix: Reset c_index and switch c_cols to GIDs
    for(k=0; k<c_current; k++){
      c_index[c_cols[k]]=-1;
      c_cols[k]=bcols[c_cols[k]]; // Switch from local to global IDs.     
    }
    // Local matrix: Insert.
    //
    // We should check global error results after the algorithm is
    // through.  It's probably safer just to let the algorithm run all
    // the way through before doing this, since otherwise we have to
    // remember to free all allocations carefully.
    int err = C_filled ?
      C.SumIntoGlobalValues(global_row,c_current,c_vals,c_cols)
      :
      C.InsertGlobalValues(global_row,c_current,c_vals,c_cols);   

    // Remote matrix: Zero current counts again for matrix
    c_current=0;    

    // Remote matrix: Do the "middle product"
    for(k=0; k<Aview.numEntriesPerRow[i]; ++k) {
      int Ak=Aindices_i[k];
      double Aval = Aval_i[k];
      // We're skipping local entries on this pass.
      if(!Bview.remote[Ak] || Aval==0) continue;
      
      int* Bcol_inds = Bview.indices[Ak];
      double* Bvals_k = Bview.values[Ak];

      for(j=0; j<Bview.numEntriesPerRow[Ak]; ++j) {
	int col=Bcol_inds[j];
	if(c_index[col]<0){
	  c_cols[c_current]=col;	      
	  c_vals[c_current]=Aval*Bvals_k[j];
	  c_index[col]=c_current;
	  c_current++;
	}
	else{
	  c_vals[c_index[col]]+=Aval*Bvals_k[j];	    
	}
      }
    }    
    // Remote matrix: Reset c_index and switch c_cols to GIDs
    for(k=0; k<c_current; k++){
      c_index[c_cols[k]]=-1;
      c_cols[k]=bcols_import[c_cols[k]];      
    }
    // Remove matrix: Insert
    //
    // See above (on error handling).
    err = C_filled ?
      C.SumIntoGlobalValues(global_row,c_current,c_vals,c_cols)
      :
      C.InsertGlobalValues(global_row,c_current,c_vals,c_cols);    
  }
  
  // Since Multiply won't do this
  if(call_FillComplete_on_result)
    C.FillComplete(B.DomainMap(),A.RangeMap());

  delete [] dwork;
  delete [] iwork;
  delete [] c_index;
  return(0);
}







int MatrixMatrix::mult_A_B(const Epetra_CrsMatrix & A,
			   CrsMatrixStruct & Aview,
			   const Epetra_CrsMatrix & B,
			   CrsMatrixStruct& Bview,
			   Epetra_CrsMatrix& C,
			   bool call_FillComplete_on_result){

  int i,rv;
  Epetra_Map* mapunion1 = 0;

  // DEBUG
#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
#endif  



  // If the user doesn't want us to call FillComplete, use the general routine
  if(!call_FillComplete_on_result) {
    rv=mult_A_B_general(A,Aview,B,Bview,C,false);
    return rv;
  }

  // Is this a "clean" matrix
  bool NewFlag=!C.IndicesAreLocal() && !C.IndicesAreGlobal();


  // Does ExtractCrsDataPointers work?
  int *C_rowptr, *C_colind;
  double * C_vals;
  C.ExtractCrsDataPointers(C_rowptr,C_colind,C_vals);
  bool ExtractFailFlag=!C_rowptr || !C_colind || !C_vals;

  // It's a new matrix that hasn't been fill-completed, use the general routine
  if(!NewFlag && ExtractFailFlag){
    rv=mult_A_B_general(A,Aview,B,Bview,C,call_FillComplete_on_result);
    return rv;
  }

#ifdef ENABLE_MMM_TIMINGS
  if(NewFlag) mtime=M.getNewTimer("M5 CMap");
  else mtime=M.getNewTimer("M5r CMap");
  mtime->start();
#endif
  // If new, build & clobber a colmap for C
  if(NewFlag){
    if(Bview.importMatrix) {
      EPETRA_CHK_ERR( aztecoo_and_ml_compatible_map_union(B,*Bview.importMatrix,mapunion1) );
      EPETRA_CHK_ERR( C.ReplaceColMap(*mapunion1) );
    }
    else
      EPETRA_CHK_ERR( C.ReplaceColMap(B.ColMap()) );
  }

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  if(NewFlag) mtime=M.getNewTimer("M5 Lookups");
  else mtime=M.getNewTimer("M5r Lookups");
  mtime->start();
#endif

  // ********************************************
  // Setup Bcol2Ccol / Bimportcol2Ccol lookups
  // ********************************************
  const Epetra_Map * colmap_B = &(B.ColMap());
  const Epetra_Map * colmap_C = &(C.ColMap());

  std::vector<int> Bcol2Ccol(colmap_B->NumMyElements());

  if(colmap_B->SameAs(*colmap_C)){
    // Maps are the same: Use local IDs as the hash
    for(i=0;i<colmap_B->NumMyElements();i++)
      Bcol2Ccol[i]=i;				
  }
  else {
    // Maps are not the same:  Use the map's hash
    for(i=0;i<colmap_B->NumMyElements();i++){
      Bcol2Ccol[i]=colmap_C->LID(colmap_B->GID(i));
      if(Bcol2Ccol[i]==-1) EPETRA_CHK_ERR(-11);
    }
  }

  std::vector<int> Bimportcol2Ccol;
  if(Bview.importMatrix){
#ifdef LIGHTWEIGHT_MATRIX
    Bimportcol2Ccol.resize(Bview.importMatrix->ColMap_.NumMyElements());
    for(i=0;i<Bview.importMatrix->ColMap_.NumMyElements();i++){
      Bimportcol2Ccol[i]=colmap_C->LID(Bview.importMatrix->ColMap_.GID(i));
      if(Bimportcol2Ccol[i]==-1) EPETRA_CHK_ERR(-12);
    }
#else
    Bimportcol2Ccol.resize(Bview.importMatrix->ColMap().NumMyElements());
    for(i=0;i<Bview.importMatrix->ColMap().NumMyElements();i++){
      Bimportcol2Ccol[i]=colmap_C->LID(Bview.importMatrix->ColMap().GID(i));
      if(Bimportcol2Ccol[i]==-1) EPETRA_CHK_ERR(-12);
    }
#endif
  }		

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  if(NewFlag) mtime=M.getNewTimer("M5 Core");
  else mtime=M.getNewTimer("M5r Core");
  mtime->start();
#endif

  // Call the appropriate core routine
  if(NewFlag) {
    EPETRA_CHK_ERR(mult_A_B_newmatrix(A,B,Bview,Bcol2Ccol,Bimportcol2Ccol,C));
  }
  else {
    EPETRA_CHK_ERR(mult_A_B_reuse(A,B,Bview,Bcol2Ccol,Bimportcol2Ccol,C));
  }
  
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif

  // Cleanup      
  delete mapunion1;
  return 0;
}




}//namespace EpetraExt
