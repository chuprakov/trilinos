
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
#include "Epetra2TSFutils.h"
#include "Aztec2Petra.h"

#include "Epetra_Vector.h" //temp!
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Map.h"

Epetra_RowMatrix *Aztec2TSF(   AZ_MATRIX * Amat, 
			Epetra_Comm * & junkcomm,
			Epetra_BlockMap * & VbrMap,
			      Epetra_Map **petra_maps) {
  Epetra_RowMatrix *tmp;
  Aztec2TSF(Amat, junkcomm, VbrMap, tmp, NULL, NULL, NULL, NULL, petra_maps);
  return tmp;
}

TSFLinearOperator Aztec1x1VBR_2_TSF(AZ_MATRIX *Fp, 
				    //TSFLinearOperator Fp_tsf,
		      TSFVectorSpace pressureSpace,
		      Epetra_Comm *comm, int proc_config[])
{
  AZ_MATRIX *Fpmat;
  Epetra_CrsMatrix *Fp_crs;
  Epetra_RowMatrix *Fp_row;
  int **junk_indices;
  Epetra_BlockMap *junkmap;
 
  int *bindx;
  double *val;

  int Nrows, nz_ptr, ival, iblk_row, j, jblk, NnzTotal, oldnz_ptr;
  Nrows  = Fp->data_org[AZ_N_border]+Fp->data_org[AZ_N_internal];
  nz_ptr = Nrows + 1;

  NnzTotal = ival = Fp->indx[Fp->bpntr[Nrows]];
  bindx = (int    *) malloc(sizeof(int)*(NnzTotal+2));
  val   = (double *) malloc(sizeof(double)*(NnzTotal+2));
    
  bindx[0] = nz_ptr;

  /* loop over block rows */

  for (iblk_row = 0; iblk_row < Nrows; iblk_row++) {
    oldnz_ptr = nz_ptr;
    ival = Fp->indx[Fp->bpntr[iblk_row]];
    for (j = Fp->bpntr[iblk_row]; j < Fp->bpntr[iblk_row+1]; j++) {
      jblk = Fp->bindx[j];
      if (jblk == iblk_row) {
   //      	if(jblk=iblk_row) 
//         	if(iblk_row == 0) //Pinned 	
//         	val[iblk_row] = 1.0; //Pinned Fp
//         	else
	val[iblk_row] = Fp->val[ival];
      }
      else {
	bindx[nz_ptr] = jblk;
  //	if(iblk_row == 0) //Pinned
  //	val[nz_ptr++] = 0.0; //Pinned
  //	else
	val[nz_ptr++] = Fp->val[ival];
      }
      ival += 1;
    }
    bindx[iblk_row+1] = bindx[iblk_row] + nz_ptr - oldnz_ptr;
  }

  Fpmat = AZ_matrix_create(Fp->data_org[AZ_N_internal] +
			   Fp->data_org[AZ_N_border]);

  int old_matrix_type = Fp->data_org[AZ_matrix_type];
  Fp->data_org[AZ_matrix_type] = AZ_MSR_MATRIX;
  AZ_set_MSR(Fpmat, bindx, val, Fp->data_org, 0, NULL, AZ_LOCAL);

  Epetra_Vector *tmpSolution = NULL, *residual = NULL;

  Aztec2Petra (proc_config, Fpmat, NULL, NULL, comm, *& junkmap, 
	       Fp_row, *& tmpSolution, *& residual, junk_indices);

  free( val);
  free( bindx);
  AZ_matrix_destroy(&Fpmat);
  Fp->data_org[AZ_matrix_type] = old_matrix_type;

  Fp_crs = dynamic_cast<Epetra_CrsMatrix*>(Fp_row);

  TSFLinearOperator Fp_tsf  = EpetraCRS2TSF(pressureSpace,pressureSpace,
					    Fp_crs);

  return Fp_tsf;
}

TSFLinearOperator Aztec2TSFpin(TSFLinearOperator Fp_tsf, TSFLinearOperator C_tsf, Epetra_Comm *comm)
{
  cerr << "\n About to pin Fp";

  Epetra_CrsMatrix *C_crs = PetraMatrix::getConcrete(C_tsf);
  Epetra_CrsMatrix *Fp_crs = PetraMatrix::getConcrete(Fp_tsf);

  // Pinning Stuff
  Epetra_Vector x_ran(C_crs->RowMap());
  Epetra_Vector y_ran(C_crs->RowMap());

  x_ran.Random();

  C_crs->Multiply(false, x_ran, y_ran);

  cerr << "\n C multiplied";

  double *before, *after;
  int nzeros = 0;
  int tmp;
  int ind[x_ran.MyLength()];
  int i;

  x_ran.ExtractView(&before);
  y_ran.ExtractView(&after);

  cerr << "\n Length of x vec is: " << x_ran.MyLength(); 

  for(i=0; i<x_ran.MyLength(); i++)
    {
      if(before[i] == after[i])
        {
          ind[nzeros] = i;
          tmp = i;
          cerr << "\n Pinning: " << ind[nzeros];
          nzeros = nzeros + 1;
        }
    }

  if(nzeros>1) {
cerr << "\n WARNING! Number of pinned are: " << nzeros; 
 exit(1);
  }

  cerr << "\n row to pin is: " << ind[nzeros-1];
  cerr << "\n row to pin is (tmp): " << tmp;

  if(nzeros == 0) return Fp_tsf;
  else 
    {
      int MyRowNum = tmp; //need to get GID of this with comm!!!
      int MyColNum = Fp_crs->NumMyCols();
      double *val = new double[MyColNum];
      int *ind = new int[MyColNum];
      int NumIndices;

      cerr << "\n Rownum is: " << MyRowNum << " Col num is: " << MyColNum;

      Fp_crs->ExtractMyRowView(MyRowNum, NumIndices, val, ind);

      cerr << "\n Extraction complete! ";
      cerr << "\n numentries is: " << NumIndices;

      for(int kkk=0; kkk<NumIndices; kkk++)
        {
          if(kkk==MyRowNum) val[MyRowNum] = 1.0;
          else         val[kkk] = 0.0;

          ind[kkk] = kkk;
          cerr << "\n ind is: " << ind[kkk] << " val is: " << val[kkk]; 
        }
      cerr << "\n About to insert ";

      //      Fp_crs->InsertMyValues(MyRowNum, MyColNum, val, ind); 

      //  TSFVectorSpace pressureSpace = C_tsf.range();
  
      //  Fp_tsf  = EpetraCRS2TSF(pressureSpace,pressureSpace,
      //	    Fp_crs);
   
  return Fp_tsf;
    }
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

  // Note: TSF constructor is PetraMatrix(domain, range)
  PetraMatrix* F_petra = new PetraMatrix(velocitySpace, velocitySpace);
  F_petra->setPetraMatrix(F_crs,true);    // insert the epetra matrix into our TSF
  TSFLinearOperator F_tsf = F_petra; // matrix and make it a TSF linear op

  // Note: TSF constructor is PetraMatrix(domain, range)
  PetraMatrix* B_petra = new PetraMatrix(velocitySpace, pressureSpace);
  B_petra->setPetraMatrix(B_crs,true);
  TSFLinearOperator B_tsf = B_petra;

  // Note: TSF constructor is PetraMatrix(domain, range)
  PetraMatrix* Bt_petra = new PetraMatrix(pressureSpace, velocitySpace);
  Bt_petra->setPetraMatrix(Bt,true);
  TSFLinearOperator Bt_tsf = Bt_petra;
  TSFLinearOperator C_tsf;

  if (C != NULL) {
    // Note: TSF constructor is PetraMatrix(domain, range)
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

  // Note: TSF constructor is PetraMatrix(domain, range)
  PetraMatrix* result_petra = new PetraMatrix(Bt.domain(), B.range());
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
  //  cerr << "got here in matrix add " << endl;
  Epetra_CrsMatrix *result_crs = Epetra_MatrixAdd(B_crs,Bt_crs,scalar);


  // Note: TSF constructor is PetraMatrix(domain, range)
  PetraMatrix* result_petra = new PetraMatrix(Bt.domain(), B.range());
  result_petra->setPetraMatrix(result_crs,true);  // insert the epetra matrix into our TSF
  result = result_petra;                          // matrix and make it a TSF linear op

  return 0;
}


// #include "ml_epetra_operator.h"
// #include "ml_aztec_utils.h"
// int ML_TSF_defaults(TSF::TSFLinearSolver &FSolver, 
// 		    ML_solverData *solver_data,
// 		    bool symmetric, Epetra_RowMatrix *F)
// {
//   ML *ml_handle;
//   ML_Aggregate *agg_object;



//   int N_levels = 10;
//    ML_Set_PrintLevel(10);
//    ML_Create(&ml_handle, N_levels);
//    solver_data->ml = ml_handle;
//    EpetraMatrix2MLMatrix(ml_handle, 0, F);
//    ML_Aggregate_Create(&agg_object);
//    ML_Aggregate_Set_MaxCoarseSize(agg_object,30);
//    if (symmetric != true) {
//      ML_Set_Symmetrize(ml_handle, ML_TRUE);
//      ML_Aggregate_Set_DampingFactor(agg_object,0.25);
//    }
//    N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
//                                                   ML_INCREASING, agg_object);
//    if (symmetric != true) {
//      if (solver_data->aztec_status.get() == 0) {
//        double *dtemp = new double[AZ_STATUS_SIZE];
//        solver_data->aztec_status = TSFSmartPtr<double>(dtemp, true);
//      }
//      if (solver_data->aztec_proc_config.get() == 0) {
//        int *itemp = new int[AZ_PROC_SIZE];
//        solver_data->aztec_proc_config = TSFSmartPtr<int>(itemp, true);
//      }


//      int options[AZ_OPTIONS_SIZE];
//      double params[AZ_PARAMS_SIZE];
// #ifdef EPETRA_MPI
//      AZ_set_proc_config(solver_data->aztec_proc_config.get(), MPI_COMM_WORLD);
// #else
//      AZ_set_proc_config(solver_data->aztec_proc_config.get(), AZ_NOT_MPI);
// #endif
//      AZ_defaults(options, params);
//      options[AZ_precond] = AZ_dom_decomp;
//      options[AZ_subdomain_solve] = AZ_ilut;
//      options[AZ_overlap] = 1;
//      params[AZ_ilut_fill] = 2.;


//      for (int level = 0; level < N_levels; level++) {
//        ML_Gen_SmootherAztec(ml_handle, level, options, params, 
// 			    solver_data->aztec_proc_config.get(), solver_data->aztec_status.get(), 
// 			    AZ_ONLY_PRECONDITIONER,  ML_BOTH, NULL);
//      }
//    }

//    else 
//      ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS, ML_BOTH, 1, 1);

//    ML_Gen_Solver    (ml_handle, ML_MGV, 0, N_levels-1);

//    Epetra_ML_Operator  *MLop = new Epetra_ML_Operator(ml_handle,
// 				  (F->OperatorDomainMap().Comm()),
// 				   (F->OperatorDomainMap()),
// 			           (F->OperatorDomainMap()));
//    MLop->SetOwnership(true);
//    ML_Aggregate_Destroy(&agg_object);

//    if (symmetric == true)
//      solver_data->azOptions.put(AZ_solver, AZ_cg);
//    else
//      solver_data->azOptions.put(AZ_solver, AZ_gmres);

//    solver_data->azOptions.put(AZ_kspace, 100);
//    solver_data->azOptions.put(AZ_conv, AZ_r0);
//    solver_data->azParams.put(AZ_tol, 1e-8);
//    solver_data->azOptions.put(AZ_max_iter, 200);
//    solver_data->azOptions.put(AZ_output, 1);
//    TSFSmartPtr<Epetra_Operator> Smart_MLprec = TSFSmartPtr<Epetra_Operator>(MLop, true);

//    printf("commenting out FSolver due to compilation problems\n");
//    exit(1);
//    //   FSolver = new AZTECSolver(solver_data->azOptions, solver_data->azParams, Smart_MLprec);

//    return N_levels;
// }
