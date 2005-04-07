
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

#ifndef AZTEC
#define AZTEC
#endif
#ifndef _AZTEC2TSF_H_
#define _AZTEC2TSF_H_

#ifndef __cplusplus
#define __cplusplus
#endif

#include "az_aztec.h"

#ifdef AZTEC_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "TSFLinearOperator.h"
#include "TSFHashtable.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"

extern int Aztec2TSF(   AZ_MATRIX * Amat,
                        Epetra_Comm * & junkcomm,
                        Epetra_BlockMap * & VbrMap,
                        Epetra_RowMatrix * &jacobian,
                        double * x,            Epetra_Vector ** tmpSolution,
                        double * resid_vector, Epetra_Vector ** residual,
                        Epetra_Map **);

extern TSF::TSFLinearOperator Aztec1x1VBR_2_TSF(AZ_MATRIX *Fp,
                                                TSF::TSFVectorSpace pressureSpace,
                                                Epetra_Comm *comm, int proc_config[]);

//extern TSF::TSFLinearOperator Aztec1x1VBR_2_TSF(AZ_MATRIX *Fp,
//                                                Epetra_CrsMatrix *C_crs,
//                                                TSF::TSFVectorSpace pressureSpace,
//                                                Epetra_Comm *comm, int proc_config[]);

extern TSF::TSFLinearOperator Aztec2TSFpin(TSF::TSFLinearOperator Fp_tsf,
                                            TSF::TSFLinearOperator C_tsf,
                                            Epetra_Comm *comm);

extern TSF::TSFLinearOperator TSFDirichletpin(TSF::TSFLinearOperator C_tsf,
					      int *coordinate2,
					      int numbc,				    
                                            Epetra_Comm *comm);

extern Epetra_RowMatrix *Aztec2TSF(   AZ_MATRIX * Amat,
                                      Epetra_Comm * & junkcomm,
                                      Epetra_BlockMap * & VbrMap,
                                      Epetra_Map **);

extern int TSF_MatrixMult(const TSF::TSFLinearOperator& B,const TSF::TSFLinearOperator& Bt,
                          TSF::TSFLinearOperator& result);
extern int TSF_MatrixAdd(const TSF::TSFLinearOperator& B,const TSF::TSFLinearOperator& Bt, double scalar,
                         TSF::TSFLinearOperator& result);

struct ML_TSF_data {
  TSF::TSFHashtable<int, int>  azOptions;
  TSF::TSFHashtable<int, double>  azParams;
  ML *ml;
  TSF::TSFSmartPtr<double> aztec_status;
  TSF::TSFSmartPtr<int> aztec_proc_config;
};
typedef struct ML_TSF_data ML_solverData;



extern int ML_TSF_defaults(TSF::TSFLinearSolver &FSolver,
                           ML_solverData *,
                           bool symmetric,
                           Epetra_RowMatrix *);




#endif /* _VBR2PETRA_H_ */
