
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
