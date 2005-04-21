
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

#include "ReadPetra.h"

int ReadPetraMatrix(Epetra_Map *row_map, Epetra_Map *col_map, 
		    Epetra_CrsMatrix **A, char *fname) { 



    Epetra_CrsMatrix *AA = new Epetra_CrsMatrix(Copy, *row_map, *col_map, 0);
    if (AA==0) EPETRA_CHK_ERR(-3); // Ran out of memory

    FILE *fp;
    int row, col, count = 0;
    double value;
    int ret_val;

    fp = fopen(fname,"r");
    if (fp == NULL) {
      printf("File %s not found\n",fname);
      exit(1);
    }
    while (   (ret_val = fscanf(fp,"%d%d%lf",&row, &col, &value)) == 3) {
      row--;
      col--;
      AA->InsertGlobalValues(row, 1, &value, &col);
    }
    fclose(fp);

    int ierr=AA->FillComplete(*col_map, *row_map);
    if (ierr!=0) {
      cerr <<"Error in Epetra_VbrMatrix FillComplete" << ierr << endl;
      EPETRA_CHK_ERR(ierr);
    }
    *A = AA;   

    return 0;
}
int ReadPetraVector(Epetra_Vector *AA, char *fname) { 

  FILE *fp;
  int location;
  double value;

  fp = fopen(fname,"r");
  if (fp == NULL) {
    printf("File %s not found\n",fname);
    exit(1);
  }
  int ret_val;
  while (   (ret_val = fscanf(fp,"%d%lf",&location, &value)) == 2) {
      location--;
      AA->ReplaceGlobalValues(1, &value, &location);
  }
  fclose(fp);

  return 0;
}

#include "az_aztec_defs.h"
#define BLOCK_SIZE 3
#define MAX_BLOCKS_PER_ROW 30
void ReadAztecVbr(int bsize, int Nnz, int Nbrows,
		  AZ_MATRIX **Amat, char *fname)
{
  int    *rpntr, *cpntr, *indx, *bindx, *bpntr;
  double *aaa;
  int Nnz_per_block = 2;
  int TotalBlocks;
  int row_in_next_block;
  int Nread;
  FILE *fp;
  int row,col;
  double val;
  int Nblocks_in_row = 0;
  int block_row = 0;
  int block_col[MAX_BLOCKS_PER_ROW];
  double blocks[MAX_BLOCKS_PER_ROW][BLOCK_SIZE][BLOCK_SIZE];
  int i,j,k, bcol;
  int all_done = 0;
  int count = 0, bcount = 0;
  int total_size;
  int *data_org;

  if (bsize != BLOCK_SIZE) {
   printf("block size in input argument does not match #define BLOCK_SIZE.\n");
   printf("Recompile with correct block size and hopefully this will work.\n");
   exit(1);
  }

  data_org = (int *) malloc(sizeof(int)*AZ_COMMLESS_DATA_ORG_SIZE); 
  rpntr = (int *) malloc(sizeof(int)*(Nbrows+1));  
  cpntr = (int *) malloc(sizeof(int)*(Nbrows+1));  
  TotalBlocks = Nnz/Nnz_per_block;
  indx = (int *) malloc(sizeof(int)*(TotalBlocks+1));  
  bindx = (int *) malloc(sizeof(int)*TotalBlocks); 
  total_size = TotalBlocks*BLOCK_SIZE*BLOCK_SIZE;
  aaa = (double *) malloc(sizeof(double)*total_size);
  bpntr = (int *) malloc(sizeof(int)*(Nbrows+1));    
  bpntr[0] = 0;
  
  for (i = 0; i <= Nbrows; i++) {
    rpntr[i] = BLOCK_SIZE*i;
    cpntr[i] = BLOCK_SIZE*i;
  }

  row_in_next_block = BLOCK_SIZE+1;

  fp = fopen(fname,"r");
  if (fp == NULL) {
    printf("File %s not found\n",fname);
    exit(1);
  }

  Nread = fscanf(fp,"%d%d%lf",&row,&col,&val);
  while ( all_done == 0) {
    if (Nread < 3) {
      row = row_in_next_block;
      all_done = 1;
    }
    if (row < row_in_next_block) {
      /* record the entry */
      bcol = (int) floor( (((double)(col)) - .5)/((double) BLOCK_SIZE));

      /* see if bcolumn already exists */
      for (i = 0; i < Nblocks_in_row; i++) {
	if (block_col[i] == bcol) break;
      }
      if (i == Nblocks_in_row) {
	if (i == MAX_BLOCKS_PER_ROW) {
	  printf("Exceeded the total number of allowable blocks per row: ");
	  printf("%d.\nChange MAX_BLOCKS_PER_ROW and recompile again.\n",
		 MAX_BLOCKS_PER_ROW);
	  exit(1);
	}

	block_col[i] = bcol;
	Nblocks_in_row++;
	for (j = 0; j < BLOCK_SIZE; j++) {
	  for (k = 0; k < BLOCK_SIZE; k++) {
	    blocks[i][j][k] = 0.;
	  }
	}
      }
      blocks[i][row-block_row*BLOCK_SIZE - 1][col- bcol*BLOCK_SIZE-1] = val;
      Nread = fscanf(fp,"%d%d%lf",&row,&col,&val);
    }
    else {

      /* print all the blocks in this block row */

      for (i = 0; i < Nblocks_in_row; i++) {
	if (bcount >= TotalBlocks) {
	  printf("We've run out of space. Space was allocated assuming\n");
	  printf("that there are %d nonzero per block.  ",Nnz_per_block);
	  printf("If this is not\ncorrect, lower Nnz_per_block, recompile ");
	  printf("and try again.\n");
	  exit(1);
	}
		
	indx[bcount++] = count;
	bindx[bpntr[block_row] + i] = block_col[i];
	/* printf("block row and column = %d %d\n",block_row, block_col[i]);*/
	for (k = 0; k < BLOCK_SIZE; k++) {
	  for (j = 0; j < BLOCK_SIZE; j++) {
	    /* printf("   %20.13e ",blocks[i][j][k]); */
	    aaa[count++] = blocks[i][j][k];
	  }
	  /* printf("\n"); */
	}
      }
      row_in_next_block += BLOCK_SIZE;
      block_row++;
      Nblocks_in_row = 0;
      bpntr[block_row] = bcount;
    }
  }
  indx[bcount++] = count;
  data_org[AZ_matrix_type] = AZ_VBR_MATRIX;
  data_org[AZ_N_int_blk] = block_row;
  data_org[AZ_N_bord_blk] = 0;
  data_org[AZ_N_internal] = block_row*BLOCK_SIZE;
  data_org[AZ_N_border] = 0;
  data_org[AZ_N_external] = 0;
  data_org[AZ_N_ext_blk] = 0;
  data_org[AZ_N_neigh] = 0;
  data_org[AZ_total_send] = 0;
  data_org[AZ_name] = 666;

  *Amat = AZ_matrix_create(data_org[AZ_N_internal]+data_org[AZ_N_border]);
  AZ_set_VBR(*Amat, rpntr, cpntr, bpntr, indx, bindx, aaa, data_org,
  	     0, NULL, AZ_LOCAL);



}



void AZ_capt_matrix(double val[], int indx[], int bindx[], int rpntr[],
		       int cpntr[], int bpntr[], 
		       int data_org[])

     /*******************************************************************************

										     Routine to capture matrix in a form acceptable to MATLAB.  Prints matrix out
										     in i,j,a(i,j) format.  Currently implemented only for one processor.  Based
										     on the routine AZ_output_matrix.
										     Test to see if we should capture matrix, rhs and partitioning info
										     in an ASCII data file.  If the file "AZ_write_matrix_now" exists in
										     the current working directory, then the files 

										     - AZ_capture_matrix.dat
										     - AZ_capture_rhs.dat
										     - AZ_capture_partition.dat (VBR only)
  
										     will be appended with the current matrix in (i,j,val) format, the
										     current RHS and the current partition information.  The existence
										     of "AZ_write_matrix_now" is check each time.  Thus, capturing can
										     be turned on and off at will during the run of a simulation.


										     Author:          Michael A. Heroux, 9222, SNL.
										     =======

										     Return code:     void
										     ============

										     Parameter list:
										     ===============

										     val:             Array containing the nonzero entries of the matrix.

										     indx,
										     bindx,
										     rpntr,
										     cpntr,
										     bpntr:           Arrays used for DMSR and DVBR sparse matrix storage.

										     data_org:        Array containing information on the distribution of the
										     matrix to this processor as well as communication parameters
										     (see file params.txt).
										     b:               Right hand side values.

     *******************************************************************************/

{

  /* local variables */

  int  iblk_row, i, j, ib1, ib2, n1, jblk, m1, ipoint, jpoint;
  int  ival = 0;
  int  num_total_nonzeros;
  int  num_total_nodes, num_total_equations = 0;
  int  Proc, Num_Proc;

  /********** execution begins **********/
  { FILE *AZ_capture_flag;
  AZ_capture_flag = fopen("AZ_write_matrix_now","r");
  if(AZ_capture_flag)    {
    FILE *AZ_capt_matrix, *AZ_capture_rhs, *AZ_capture_partition;

    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
      num_total_nodes    = data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk];
      num_total_equations = rpntr[num_total_nodes];
      num_total_nonzeros = indx[bpntr[num_total_nodes]];

      /***** Print out the VBR partitioning information for the matrix *****/

      AZ_capture_partition = fopen("AZ_capture_partition.dat","a");
	  
      fprintf(AZ_capture_partition, "Start of partition\n");

      for (i = 0; i < num_total_nodes + 1; i++)
	fprintf(AZ_capture_partition,"%d\n", rpntr[i]);

      fclose(AZ_capture_partition);

      /***** Print out the VBR i,j,a(i,j) information for the matrix *****/

      AZ_capt_matrix = fopen("AZ_capture_matrix.dat","a");
      fprintf(AZ_capt_matrix, "Start of VBR matrix\n");
      fprintf(AZ_capt_matrix, "%d %d\n", 
	      num_total_equations, num_total_nonzeros);

      /* loop over block rows */

      for (iblk_row = 0; iblk_row < num_total_nodes; iblk_row++) {

	/* the starting point row index of the current block */

	ib1 = rpntr[iblk_row];      

	/* number of rows in the current row block */

	m1 = rpntr[iblk_row+1] - ib1;

	/* starting index of current row block */

	ival = indx[bpntr[iblk_row]];

	/* loop over all the blocks in the current block-row */

	for (j = bpntr[iblk_row]; j < bpntr[iblk_row+1]; j++) {
	  jblk = bindx[j];

	  /* the starting point column index of the current block */

	  ib2 = cpntr[jblk];

	  /* number of columns in the current block */

	  n1 = cpntr[jblk+1] - ib2;

	  for (jpoint = 0; jpoint < n1; jpoint++)
	    for (ipoint = 0; ipoint < m1; ipoint++) {
	      fprintf(AZ_capt_matrix,"%d %d %22.16e\n", 
		      ib1+ipoint+1, 
		      ib2+jpoint+1, 
		      val[ival+jpoint*m1+ipoint]);
	  
	    }

	  ival += m1*n1;

	}
      }
      fclose(AZ_capt_matrix);
    }

    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
      num_total_equations = data_org[AZ_N_internal]+data_org[AZ_N_border];

      num_total_nonzeros = bindx[num_total_equations]-1;

      /***** Print out the MSR i,j,a(i,j) information for the matrix *****/

      AZ_capt_matrix = fopen("AZ_capture_matrix.dat","a");
      fprintf(AZ_capt_matrix, "Start of MSR matrix\n");
      fprintf(AZ_capt_matrix, "%d %d\n", 
	      num_total_equations, num_total_nonzeros);
      for (i = 0; i < num_total_equations; i++) {
	fprintf(AZ_capt_matrix,"%d %d %22.16e\n", i+1, i+1, val[i]);
	for (j = bindx[i]; j < bindx[i+1]; j++ )
	  fprintf(AZ_capt_matrix,"%d %d %22.16e\n", 
		  i+1, bindx[j]+1, val[j]);
      }
      fclose(AZ_capt_matrix);
    }
    fclose(AZ_capture_flag);
  }
  }
} /* AZ_capture_matrix */
       

  
