
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Aztec2TSF.h"
#include "Vbr2Petra.h"
#include "TSFMPI.h"
#include "TSFOut.h"
#include "TSFVectorType.h"
#include "TSFProductSpace.h"
#include "TSFBlockLinearOperator.h"
#include "PetraVectorType.h"
#include "PetraVectorSpace.h"
#include "PetraVector.h"
#include "PetraMatrix.h"
#include "GMRESSolver.h"
#include "AZTECSolver.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"
#include "TSFPreconditioner.h"
#include "TSFMatrixOperator.h"
#include "TSFLinearOperator2EpetraRowMatrix.h"

Epetra_RowMatrix *Aztec2TSF(   AZ_MATRIX * Amat, 
			Epetra_Comm * & junkcomm,
			Epetra_BlockMap * & VbrMap,
			      Epetra_Map **petra_maps) {
  Epetra_RowMatrix *tmp;
  Aztec2TSF(Amat, junkcomm, VbrMap, tmp, NULL, NULL, NULL, NULL, petra_maps);
  return tmp;
}



int Aztec2TSF(	AZ_MATRIX * Amat, 
		Epetra_Comm * & junkcomm,
		Epetra_BlockMap * & map,
		Epetra_RowMatrix * &jacobian,
		double * x,            Epetra_Vector ** tmpSolution,
		double * resid_vector, Epetra_Vector ** residual,
		Epetra_Map **petra_maps)
{

  Epetra_CrsMatrix *F_crs, *B_crs, *Bt, *C;
  int blk_size = 3;

  Epetra_Map *VbrMap;

  // 1) Make a series of petra matrices corresponding to blocks.

  VbrMatrix2PetraMatrix(blk_size, Amat, junkcomm, F_crs, B_crs, Bt, C,
                        &VbrMap, petra_maps);
  map = VbrMap;

  // Build TSF vector spaces

  Epetra_Map *vel1 = (Epetra_Map *) &(F_crs->OperatorDomainMap());
  Epetra_Map *vel2 = (Epetra_Map *) &(F_crs->RowMatrixColMap());
  Epetra_Import *vel3  = (Epetra_Import *) F_crs->RowMatrixImporter();
  TSFSmartPtr<Epetra_Map> tv1 = TSFSmartPtr<Epetra_Map>(vel1,false);
  TSFSmartPtr<Epetra_Map> tv2 = TSFSmartPtr<Epetra_Map>(vel2,false);
  TSFSmartPtr<Epetra_Import> tv3 = TSFSmartPtr<Epetra_Import>(vel3,false);
  TSFVectorSpace velocitySpace = new PetraVectorSpace( tv1, tv2, tv3);

  Epetra_Map *p1 = (Epetra_Map *) &(Bt->OperatorDomainMap());
  Epetra_Map *p2 = (Epetra_Map *) &(Bt->RowMatrixColMap());
  Epetra_Import *p3 = (Epetra_Import *) Bt->RowMatrixImporter();
  TSFSmartPtr<Epetra_Map> tp1 = TSFSmartPtr<Epetra_Map>(p1,false);
  TSFSmartPtr<Epetra_Map> tp2 = TSFSmartPtr<Epetra_Map>(p2,false);
  TSFSmartPtr<Epetra_Import> tp3 = TSFSmartPtr<Epetra_Import>(p3,false);
  TSFVectorSpace pressureSpace = new PetraVectorSpace( tp1, tp2, tp3);

  // Convert block matrices to epetra/TSF matrices

  PetraMatrix* F_petra = new PetraMatrix(velocitySpace, velocitySpace);
  F_petra->setPetraMatrix(F_crs,true);    // insert the epetra matrix into our TSF
  TSFLinearOperator F_tsf = F_petra; // matrix and make it a TSF linear op

  PetraMatrix* B_petra = new PetraMatrix(velocitySpace, pressureSpace);
  B_petra->setPetraMatrix(B_crs,true);
  TSFLinearOperator B_tsf = B_petra;

  PetraMatrix* Bt_petra = new PetraMatrix(pressureSpace, velocitySpace);
  Bt_petra->setPetraMatrix(Bt,true);
  TSFLinearOperator Bt_tsf = Bt_petra;
  TSFLinearOperator C_tsf;

  if (C != NULL) {
    PetraMatrix* C_petra = new PetraMatrix(pressureSpace, pressureSpace);
    C_petra->setPetraMatrix(C,true);
    C_tsf = C_petra;
  }

  // Set up 2x2 block operator

  TSFVectorSpace rangeBlockSpace = new TSFProductSpace(velocitySpace, pressureSpace);
  TSFVectorSpace domainBlockSpace = new TSFProductSpace(velocitySpace, pressureSpace);

  TSFLinearOperator saddleA_tsf =  new TSFBlockLinearOperator(rangeBlockSpace, domainBlockSpace);

  saddleA_tsf.setBlock(0,0,F_tsf);
  saddleA_tsf.setBlock(0,1,Bt_tsf);
  saddleA_tsf.setBlock(1,0,B_tsf);
  if (C != NULL)  saddleA_tsf.setBlock(1,1,C_tsf);

  // Build a map to go between vbr-style vectors to TSF-style vectors
  // Specifically, map[i] != 1 indicates that the ith variable is in the
  //                           first block.
  //               map[i] == 1 indicates that the ith variable is in the
  //                           second block.

 int Npressure = B_crs->NumMyRows();
 int Nvelocity = F_crs->NumMyRows();

 int *imap = (int *) malloc(sizeof(int)*(Nvelocity + Npressure));
 for (int i = 0; i < Nvelocity + Npressure; i++) imap[i] = 0;
 for (int i = 0; i < Npressure; i++)  imap[(i+1)*blk_size - 1] = 1;


 // build the epetra matrix wrapper around the TSF block matrix

 TSFLinearOperator2EpetraRowMatrix *tmp = 
   new TSFLinearOperator2EpetraRowMatrix(saddleA_tsf, junkcomm, VbrMap, imap,
					 EPETRA_MATRIX);
  jacobian = tmp;

  // 2) Make a big petra vector out of 'x', 'resid_vector', 'tmpSolution'
  //    and 'residual'. This vector should correspond to the new ordering.

  if (resid_vector != NULL) 
    *residual     = new Epetra_Vector(View, *VbrMap, resid_vector);
  if (x != NULL )
    *tmpSolution =  new Epetra_Vector(View, *VbrMap, x);





 return 0;
}
#include "ml_epetra_utils.h"

int TSF_MatrixMult(const TSFLinearOperator& B, const TSFLinearOperator& Bt,
		   TSFLinearOperator& result)
{
  if (B.isMatrixOperator()) {
    const TSFSmartPtr<const TSFMatrixOperator> M = B.getMatrix();
    const PetraMatrix* pm = dynamic_cast<const PetraMatrix*>(M.get());
    if (pm==0) {
      printf("TSF_MatrixMult: first argument is not a Petra_Matrix\n");
      exit(1);
    }
  }
  else {
      printf("TSF_MatrixMult: first argument is not a Matrix\n");
      exit(1);
  }
  Epetra_CrsMatrix  *B_crs  = PetraMatrix::getConcrete(B);

  if (Bt.isMatrixOperator()) {
    const TSFSmartPtr<const TSFMatrixOperator> M = Bt.getMatrix();
    const PetraMatrix* pm = dynamic_cast<const PetraMatrix*>(M.get());
    if (pm==0) {
      printf("TSF_MatrixMult: second argument is not a Petra_Matrix\n");
      exit(1);
    }
  }
  else {
      printf("TSF_MatrixMult: second argument is not a Matrix\n");
      exit(1);
  }



  Epetra_CrsMatrix  *Bt_crs = PetraMatrix::getConcrete(Bt);
  Epetra_CrsMatrix *result_crs = Epetra_MatrixMult(B_crs,Bt_crs);

  PetraMatrix* result_petra = new PetraMatrix(B.range(), Bt.domain());
  result_petra->setPetraMatrix(result_crs,true);  // insert the epetra matrix into our TSF
  result = result_petra;                          // matrix and make it a TSF linear op

  return 0;
}

int TSF_MatrixAdd(const TSFLinearOperator& B, const TSFLinearOperator& Bt,
		  double scalar,  TSFLinearOperator& result)
{
  if (B.isMatrixOperator()) {
    const TSFSmartPtr<const TSFMatrixOperator> M = B.getMatrix();
    const PetraMatrix* pm = dynamic_cast<const PetraMatrix*>(M.get());
    if (pm==0) {
      printf("TSF_MatrixMult: first argument is not a Petra_Matrix\n");
      exit(1);
    }
  }
  else {
      printf("TSF_MatrixMult: first argument is not a Matrix\n");
      exit(1);
  }
  Epetra_CrsMatrix  *B_crs  = PetraMatrix::getConcrete(B);

  if (Bt.isMatrixOperator()) {
    const TSFSmartPtr<const TSFMatrixOperator> M = Bt.getMatrix();
    const PetraMatrix* pm = dynamic_cast<const PetraMatrix*>(M.get());
    if (pm==0) {
      printf("TSF_MatrixMult: second argument is not a Petra_Matrix\n");
      exit(1);
    }
  }
  else {
      printf("TSF_MatrixMult: second argument is not a Matrix\n");
      exit(1);
  }



  Epetra_CrsMatrix  *Bt_crs = PetraMatrix::getConcrete(Bt);
  Epetra_CrsMatrix *result_crs = Epetra_MatrixAdd(B_crs,Bt_crs,scalar);

  PetraMatrix* result_petra = new PetraMatrix(B.range(), Bt.domain());
  result_petra->setPetraMatrix(result_crs,true);  // insert the epetra matrix into our TSF
  result = result_petra;                          // matrix and make it a TSF linear op

  return 0;
}
