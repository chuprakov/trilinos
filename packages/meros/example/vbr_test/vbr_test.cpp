#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Vbr2Petra.h"
#include "ReadPetra.h"
#include "AztecOO.h"

main() {

  AZ_MATRIX *Amat;
  Epetra_SerialComm comm;

  Epetra_CrsMatrix *F, *B, *Bt, *C;
  Epetra_Map       *FMap, *BrowMap, *VbrMap;


  /* Read file and shove data into an Aztec Matrix */
  /* This Aztec matrix should be very similar to   */
  /* what MPSalsa provides.                        */

  int Nnz, Nblks;
  printf("Enter the total number of matrix nonzeros (4228 for mac-vbr and 327519 for salsa)\n");
  scanf("%d",&Nnz);
  printf("Enter the total number of block rows (256 for mac-vbr and 4131 for salsa)\n");
  scanf("%d",&Nblks);

  ReadAztecVbr( 3, Nnz, Nblks, &Amat, "K_reordered");

  /* If the file 'AZ_write_matrix_now' is present, print the */
  /* Aztec matrix to a set of data files                     */

  AZ_capt_matrix(Amat->val, Amat->indx, Amat->bindx, Amat->rpntr, 
		 Amat->cpntr, Amat->bpntr, Amat->data_org);

  int NumMyBlocks = Amat->data_org[AZ_N_int_blk] + Amat->data_org[AZ_N_bord_blk];

  /* Convert the Aztec VBR file to a 2x2 set of block Petra matrices */

  Epetra_Map **petra_maps;
  petra_maps = (Epetra_Map **) malloc(sizeof(Epetra_Map *)*2);
  petra_maps[0] = new Epetra_Map(NumMyBlocks*(3-1), 0, comm);
  petra_maps[1] = new Epetra_Map(NumMyBlocks, 0, comm);


  VbrMatrix2PetraMatrix(3, /* proc_config, */ Amat, &comm, F, B, Bt, C, &VbrMap
			, petra_maps
			);
  FMap    = (Epetra_Map *) &(F->OperatorDomainMap());
  BrowMap = (Epetra_Map *) &(B->OperatorRangeMap());


  /* Read a petra vector from a file */

  Epetra_Vector *guess_vbr, *velocity, *pressure, *vel_rhs, *pres_rhs, *vel_tmp, *pres_tmp;
  guess_vbr = new Epetra_Vector ( (Epetra_BlockMap &) *VbrMap);

  ReadPetraVector(guess_vbr, "guess_reordered");

  velocity = new Epetra_Vector ( (Epetra_BlockMap &) *FMap);
  pressure = new Epetra_Vector ( (Epetra_BlockMap &) *BrowMap);

  /* Convert the Petra-VBR ordered vector to 2 smaller vectors */
  /* corresponding to a 2x2 block partitioning of the matrix   */

  VbrVector2BlockVector(guess_vbr, velocity, pressure, 3);

  /* Try out the matrix-vector product */

  vel_rhs  = new Epetra_Vector ( (Epetra_BlockMap &) *FMap);
  pres_rhs = new Epetra_Vector ( (Epetra_BlockMap &) *BrowMap);
  vel_tmp  = new Epetra_Vector ( (Epetra_BlockMap &) *FMap);
  pres_tmp = new Epetra_Vector ( (Epetra_BlockMap &) *BrowMap);


  F->Multiply(0, *velocity, *vel_rhs); 
  B->Multiply(0, *velocity, *pres_rhs); 
  Bt->Multiply(0, *pressure, *vel_tmp); 
  vel_rhs->Update(1., *vel_tmp, 1.);
  if (C != NULL) {
    C->Multiply(0, *pressure, *pres_tmp);
    pres_rhs->Update(1., *pres_tmp, 1.);
  }



  BlockVector2VbrVector(guess_vbr, vel_rhs, pres_rhs,3);
  double *junk;
  guess_vbr->ExtractView(&junk);

  for (int i = 0; i < Nblks*3; i++) printf("%20.13e\n",junk[i]);



}
