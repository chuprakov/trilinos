// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Aztec2Petra.h"
#include "Vbr2Petra.h"
#include "ml_include.h"

#include "ml_common.h"

#include "ml_epetra_utils.h"


#include <unistd.h>
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int VbrMatrix2PetraMatrix(int blk_size, 
		AZ_MATRIX * Amat, 
		Epetra_Comm *  comm,
		Epetra_CrsMatrix * &F,
	        Epetra_CrsMatrix * &B,
	        Epetra_CrsMatrix * &Bt,
	        Epetra_CrsMatrix * &C,
	        Epetra_Map  **VBR_map, 
		Epetra_Map **petra_maps) {

  int * MyGlobalBlocks, * global_bindx, *update;
  int proc_config[AZ_PROC_SIZE];
   proc_config[AZ_node   ] = comm->MyPID();
   proc_config[AZ_N_procs] = comm->NumProc();

#ifdef EPETRA_MPI
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
  AZ_set_comm(proc_config, MPI_COMM_WORLD);
#endif

 if (!Amat->has_global_indices) {
   //   create a global bindx
   AZ_revert_to_global(proc_config, Amat, &global_bindx, &update);
   MyGlobalBlocks = update;
 }
 else // Already have global ordering
   {
     global_bindx = Amat->bindx;
     MyGlobalBlocks = Amat->update;
     if (MyGlobalBlocks==0) EPETRA_CHK_ERR(-1);
   }

  // Get matrix information
  int NumMyBlocks = 0;

  NumMyBlocks = (Amat->data_org[AZ_N_int_blk] + Amat->data_org[AZ_N_bord_blk]);

  int * bpntr = Amat->bpntr;
  int * rpntr = Amat->rpntr;
  int * indx = Amat->indx;
  double * val = Amat->val;

  int NumGlobalBlocks;
  comm->SumAll(&NumMyBlocks, &NumGlobalBlocks, 1);

  Epetra_Map  *F_map, *Brow_map;

  F_map = petra_maps[0];
  if (F_map==0) EPETRA_CHK_ERR(-2); // Ran out of memory

  Brow_map = petra_maps[1];
  if (Brow_map==0) EPETRA_CHK_ERR(-2); // Ran out of memory

  *VBR_map = new Epetra_Map(NumGlobalBlocks*3, NumMyBlocks*3, 0, *comm);
  if (*VBR_map==0) EPETRA_CHK_ERR(-2); // Ran out of memory

  Epetra_CrsMatrix *CC   = new Epetra_CrsMatrix(Copy, *Brow_map, 0);
  if (CC==0) EPETRA_CHK_ERR(-3); // Ran out of memory

  Epetra_CrsMatrix *FF = new Epetra_CrsMatrix(Copy, *F_map, 0);
  if (FF==0) EPETRA_CHK_ERR(-3); // Ran out of memory

  Epetra_CrsMatrix *BB = new Epetra_CrsMatrix(Copy, *Brow_map, 0);
  if (BB==0) EPETRA_CHK_ERR(-3); // Ran out of memory

  Epetra_CrsMatrix *BtBt = new Epetra_CrsMatrix(Copy, *F_map, 0);
  if (BtBt==0) EPETRA_CHK_ERR(-3); // Ran out of memory
  
  int count = 0;
  int   C_empty = 1;
  for (int i=0; i<NumMyBlocks; i++) {
    int BlockRow = MyGlobalBlocks[i];
    int NumBlockEntries = bpntr[i+1] - bpntr[i];
    int *BlockIndices = global_bindx + bpntr[i];
    int therow, thecolumn, globalcol;

    for (int ii = 0; ii < NumBlockEntries; ii++) {
    for (int jj=0; jj < blk_size; jj++) {
      for (int kk=0; kk < blk_size; kk++) {
	if (val[count] != 0.0) {
	if (jj == blk_size-1) {
	  if (kk == blk_size-1) {
	    /* C block */

	    C_empty = 0;
	    therow = BlockRow;
	    thecolumn = *BlockIndices;
     int ierr = CC->InsertGlobalValues(therow, 1, &(val[count]), &thecolumn);
	  }
	  else {
	    /* Bt block */

	    therow = (blk_size-1)*(BlockRow)+kk;
	    thecolumn = *BlockIndices;
	    int ierr = BtBt->InsertGlobalValues(therow, 1, &(val[count]), &thecolumn);
	  }
	}
	else {
	  if (kk == blk_size-1) {
	    /* B block */

	    thecolumn = (blk_size-1)*(*BlockIndices)+jj;
	    int ierr = BB->InsertGlobalValues(BlockRow, 1, &(val[count]), &thecolumn);
	  }
	  else {
	    /* F block */

	    therow = (blk_size-1)*(BlockRow)+kk;
	    thecolumn = (blk_size-1)*(*BlockIndices)+jj;
	    int ierr = FF->InsertGlobalValues(therow, 1, &(val[count]), &thecolumn);
	  }
	}
      }
	count++;
      }
    }
    BlockIndices++;
    }
  }

  int ierr=BB->FillComplete(*F_map, *Brow_map);
  if (ierr!=0) {
    cerr <<"Error in Epetra_VbrMatrix TransformToLocal" << ierr << endl;
    EPETRA_CHK_ERR(ierr);
  }
  ierr=FF->FillComplete();    
  if (ierr!=0) {
    cerr <<"Error in Epetra_VbrMatrix TransformToLocal" << ierr << endl;
    EPETRA_CHK_ERR(ierr);
  }
    ierr=BtBt->FillComplete(*Brow_map, *F_map);
  if (ierr!=0) {
    cerr <<"Error in Epetra_VbrMatrix TransformToLocal" << ierr << endl;
    EPETRA_CHK_ERR(ierr);
  }
  if (!C_empty) {
    ierr=CC->FillComplete();        
    if (ierr!=0) {
      cerr <<"Error in Epetra_VbrMatrix TransformToLocal" << ierr << endl;
      EPETRA_CHK_ERR(ierr);
    }
  }
  else {
    delete CC;
    CC = NULL;
  }

    B  = BB;
    F  = FF; 
    Bt = BtBt;
    C  = CC; 

  // RPP: Can not use the OperatorRangeMap in the ctor of the "b" vector 
  // below.  In MPSalsa, we delete the VbrMatrix yet still use the vector "b".
  // Deleting the matrix deletes the OperatorRangeMap that the b vector is 
  // based on.  Losing the map means "b" and all vectors that are created 
  // with the copy constructor of "b" break.  Mike has suggested 
  // using reference counting (Boost smart pointers) so the map is not 
  // deleted.  For now we will use the "map" variable as the base map for "b". 
  //b = new Epetra_Vector (View, A->OperatorRangeMap(), az_b);

  if (!Amat->has_global_indices) {
    AZ_free((void *) global_bindx);
    AZ_free((void *) update);
  }

  return 0;
}

int VbrVector2BlockVector(Epetra_Vector *vbr, 
			  Epetra_Vector *velocity, 
			  Epetra_Vector *pressure, int blk_size)
{
  double *vbr_data, *velocity_data, *pressure_data;

  vbr->ExtractView(&vbr_data);
  velocity->ExtractView(&velocity_data);
  pressure->ExtractView(&pressure_data);
  int Npressure = pressure->Map().NumMyElements();

  for (int i = 0; i < Npressure; i++) {
    pressure_data[i] = vbr_data[(i+1)*blk_size-1];
    for (int j = 0; j < blk_size-1; j++) {
      velocity_data[ (blk_size-1)*i + j] = vbr_data[i*blk_size+j];
    }
  }
  return 0;
}

int BlockVector2VbrVector(Epetra_Vector *vbr, 
			  Epetra_Vector *velocity, 
			  Epetra_Vector *pressure, int blk_size)
{
  double *vbr_data, *velocity_data, *pressure_data;

  vbr->ExtractView(&vbr_data);
  velocity->ExtractView(&velocity_data);
  pressure->ExtractView(&pressure_data);
  int Npressure = pressure->Map().NumMyElements();

  for (int i = 0; i < Npressure; i++) {
    vbr_data[(i+1)*blk_size-1] = pressure_data[i];
    for (int j = 0; j < blk_size-1; j++) {
      vbr_data[i*blk_size+j] = velocity_data[ (blk_size-1)*i + j];
    }
  }
  return 0;
}
