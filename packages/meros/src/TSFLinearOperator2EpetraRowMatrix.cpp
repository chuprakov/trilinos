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

#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "TSFMPI.h"
#include "TSFOut.h"
#include "TSFVectorType.h"
#include "TSFProductSpace.h"
#include "TSFLinearOperator.h"
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


TSFLinearOperator2EpetraRowMatrix::TSFLinearOperator2EpetraRowMatrix(TSFLinearOperator A_tsf, const Epetra_Comm *Comm_pet,
	     Epetra_Map *Map_pet, int *map, int matrix_type) {

  OwnMap_      = true;

    // This is some garbage code to prevent both the matrix and its inverse
    // delete the map. Ugh. It might be better to use a smart pointer for
    // the map .... so then it will go away automatically at the right time

    if (matrix_type == EPETRA_INVERSE)      OwnMap_ = false;

    matrix_type_ = matrix_type;
    map_         = map;
    Matrix_      = A_tsf;
    Comm_        = Comm_pet;
    OperatorDomainMap_ = Map_pet;

    rangeBlockSpace_  = A_tsf.range();
    domainBlockSpace_ = A_tsf.domain();

    numRangeBlocks_   = rangeBlockSpace_.numBlocks();
    numDomainBlocks_  = domainBlockSpace_.numBlocks();

    rangeLength_ = rangeBlockSpace_.dim();
    domainLength_ = domainBlockSpace_.dim();

    if (!A_tsf.isZeroOperator())
      nnz_ = 9999;
    

    myname = (char *) malloc(sizeof(char)*(1+strlen("TSFLinearOperator2EpetraRowMatrix: '")));
    strcpy(myname,"TSFLinearOperator2EpetraRowMatrix: ");
    nop = (char *) malloc(sizeof(char)*strlen("' Not Implemented!"));
    strcpy(nop," Not Implemented!");
  }


  TSFLinearOperator2EpetraRowMatrix::TSFLinearOperator2EpetraRowMatrix(const Epetra_Comm *Comm_pet, int matrix_type) {

    matrix_type_ = matrix_type;
    Comm_        = Comm_pet;

    myname = (char *) malloc(sizeof(char)*(1+strlen("TSFLinearOperator2EpetraRowMatrix: '")));
    strcpy(myname,"TSFLinearOperator2EpetraRowMatrix: ");
    nop = (char *) malloc(sizeof(char)*strlen("' Not Implemented!"));
    strcpy(nop," Not Implemented!");
  }
  TSFLinearOperator2EpetraRowMatrix::~TSFLinearOperator2EpetraRowMatrix() {
    if ((OwnMap_ == true) && (map_ != NULL)) free(map_);
    free(nop);
    free(myname);
  }

  int TSFLinearOperator2EpetraRowMatrix::NumMyRows()                        const { return(rangeLength_); }
  int TSFLinearOperator2EpetraRowMatrix::NumMyCols()                        const { return(domainLength_);} 
  double TSFLinearOperator2EpetraRowMatrix::NormInf()                       const { return(-1);} 
  bool TSFLinearOperator2EpetraRowMatrix::HasNormInf()                      const { return false;}
  const Epetra_Map & TSFLinearOperator2EpetraRowMatrix::OperatorDomainMap() const {return(*OperatorDomainMap_); } 
  const Epetra_Map & TSFLinearOperator2EpetraRowMatrix::OperatorRangeMap()  const {return(*OperatorDomainMap_); }
  const Epetra_Comm & TSFLinearOperator2EpetraRowMatrix::Comm()             const { return(*Comm_);             }
  int TSFLinearOperator2EpetraRowMatrix::SetUseTranspose(bool UseTranspose)       { return(-1); }
  char  *TSFLinearOperator2EpetraRowMatrix::Label()                         const { return("TSFLinearOperator2EpetraRowMatrix");}

  // NumGlobalNonzeros is checked by AztecOO.cpp. I don't see a simple way to get nnz from the tsf operator, 
  // so since they only seem to check that it is not < 1, if operator is not a ZeroOperator, we return 9999. VEH
  int TSFLinearOperator2EpetraRowMatrix::NumGlobalNonzeros()  const {return(nnz_);}

  int TSFLinearOperator2EpetraRowMatrix::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const { 
    if (matrix_type_ == EPETRA_INVERSE) {
      printf("%sOnly ApplyInverse() (and not Apply()) can be\n",myname);
      printf("invoked on an EPETRA_INVERSE matrix.\n");
      exit(1);
    }
    else return(applyoperator(X, Y)); 
  }
  int TSFLinearOperator2EpetraRowMatrix::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

    if (matrix_type_ != EPETRA_INVERSE) {
      printf("%sOnly Apply() (and not ApplyInverse()) can be\n",myname);
      printf("invoked on an EPETRA_MATRIX matrix.\n");
      exit(1);
    }
    else 
    return(applyoperator(X, Y)); 
  }



  //***************************************************************************
  // NOT CALLED BY US!!!!
  //***************************************************************************
  int TSFLinearOperator2EpetraRowMatrix::MaxNumEntries () const {printf("%sMaxNumEntries%s\n",myname,nop);        exit(1); return 0;}

  const Epetra_BlockMap & TSFLinearOperator2EpetraRowMatrix::Map () const {
                             printf("%sMap%s\n",myname,nop);        exit(1); 
			     return(*OperatorDomainMap_);}


  int TSFLinearOperator2EpetraRowMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const {
                             printf("%sNumMyRowEntries%s\n",myname,nop);        exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const {
                             printf("%sExtractMyRowCopy%s\n",myname,nop);        exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {
                             printf("%sExtractDiagonalCopy%s\n",myname,nop);        exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
                             printf("%sMultiply%s\n",myname,nop);        exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
	    Epetra_MultiVector& Y) const {
                             printf("%sSolve%s\n",myname,nop);        exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::InvRowSums(Epetra_Vector& x) const {
                             printf("%sInvRowSums%s\n",myname,nop);        exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::LeftScale(const Epetra_Vector& x) {
                             printf("%sLeftScale%s\n",myname,nop);        exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::InvColSums(Epetra_Vector& x) const {
                                  printf("%sInvColSums%s\n",myname,nop);        exit(1); return 0;}
  bool    TSFLinearOperator2EpetraRowMatrix::Filled()         const {printf("%sFilled%s\n",myname,nop);            exit(1); return 0;}
  double TSFLinearOperator2EpetraRowMatrix::NormOne()         const {printf("%sNormOne%s\n",myname,nop);           exit(1); return 0;}
// int TSFLinearOperator2EpetraRowMatrix::NumGlobalNonzeros()  const {printf("%sNumGlobalNonzeros%s\n",myname,nop); exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::NumGlobalRows()      const {printf("%sNumGlobalRows%s\n",myname,nop);     exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::NumGlobalCols()      const {printf("%sNumGlobalCols%s\n",myname,nop);     exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::NumGlobalDiagonals() const {printf("%sNumGlobalDiagonals%s\n",myname,nop);exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::NumMyNonzeros()      const {printf("%sNumMyNonzeros%s\n",myname,nop);     exit(1); return 0;}
  int TSFLinearOperator2EpetraRowMatrix::NumMyDiagonals()     const {printf("%sNumMyDiagonals%s\n",myname,nop);    exit(1); return 0;}
  bool TSFLinearOperator2EpetraRowMatrix::LowerTriangular()   const {printf("%sLowerTriangular%s\n",myname,nop);   exit(1); return 0;}
  bool TSFLinearOperator2EpetraRowMatrix::UpperTriangular()   const {printf("%sUpperTriangular%s\n",myname,nop);   exit(1); return 0;}
  bool TSFLinearOperator2EpetraRowMatrix::UseTranspose()      const {printf("%sUseTransposer%s\n",myname,nop);     exit(1); return 0;}

  int TSFLinearOperator2EpetraRowMatrix::RightScale(const Epetra_Vector& x) {
                                    printf("%sRightScale%s\n",myname,nop);        exit(1); return(0); }
  const Epetra_Map & TSFLinearOperator2EpetraRowMatrix::RowMatrixRowMap() const {
                                  printf("%sRowMatrixRowMap%s\n",myname,nop);   exit(1); return(*OperatorDomainMap_);}
  const Epetra_Map & TSFLinearOperator2EpetraRowMatrix::RowMatrixColMap() const {
                                  printf("%sRowMatrixColMap%s\n",myname,nop);   exit(1); return(*OperatorDomainMap_);}
  const Epetra_Import * TSFLinearOperator2EpetraRowMatrix::RowMatrixImporter() const {
   	    return(0); }	 

  // The current applyoperator() is very specific to the needs of MPSalsa. In particular,
  // it is assumed that the TSF matrix that is being wrapped is a block matrix that
  // operates on a series of petra type vector spaces. It is also assumed that the 
  // data needs to be reshuffled in the vectors. For example, the 'X' and 'Y' vectors
  // correspond to something like
  //                ( u0, v0, p0, u1, v1, p1, u2, v2, p2, .... )
  // and we want to reorder this like
  //                ( u0, v0, u1,  v1, u2, v2, .... )  and (p0, p1, p2, ...)
  // Edited so applyoperator() can use composed TSF matrices that are compositions
  // of block matrices that operate on a series of petra type vector spaces. VEH


  int TSFLinearOperator2EpetraRowMatrix::applyoperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    Epetra_Vector *in_Epet  = (Epetra_Vector *) &X;
    Epetra_Vector *out_Epet = (Epetra_Vector *) &Y;

    // make two tsf vectors corresponding to the subblocks comprising 'in'

    TSFVector in_TsfVelocity = domainBlockSpace_.getBlock(0).createMember();
    TSFVector in_TsfPressure = domainBlockSpace_.getBlock(1).createMember();
    //    TSFVector in_TsfVelocity = Matrix_.getBlock(0,0).domain().createMember();
    //    TSFVector in_TsfPressure = Matrix_.getBlock(0,1).domain().createMember();
    TSFSmartPtr<Epetra_Vector> in_EpetVelocity = PetraVector::getPetraPtr(in_TsfVelocity);
    TSFSmartPtr<Epetra_Vector> in_EpetPressure = PetraVector::getPetraPtr(in_TsfPressure);   

  // make two tsf vectors corresponding to the subblocks comprising 'out'

    TSFVector out_TsfVelocity = rangeBlockSpace_.getBlock(0).createMember();
    TSFVector out_TsfPressure = rangeBlockSpace_.getBlock(1).createMember();
    //    TSFVector out_TsfVelocity = Matrix_.getBlock(0,0).domain().createMember();
    //    TSFVector out_TsfPressure = Matrix_.getBlock(0,1).domain().createMember();
    TSFSmartPtr<Epetra_Vector> out_EpetVelocity = PetraVector::getPetraPtr(out_TsfVelocity);
    TSFSmartPtr<Epetra_Vector> out_EpetPressure = PetraVector::getPetraPtr(out_TsfPressure);   

  // make two large tsf vectors for 'in' and 'out'

    TSFVector in_Tsf = domainBlockSpace_.createMember();
    //    TSFVector in_Tsf = rangeBlockSpace_.createMember();
    in_Tsf.setBlock(0, in_TsfVelocity);
    in_Tsf.setBlock(1, in_TsfPressure);
    TSFVector out_Tsf = rangeBlockSpace_.createMember();
    out_Tsf.setBlock(0, out_TsfVelocity);
    out_Tsf.setBlock(1, out_TsfPressure);

  // put 'in_Epet' into subblocks used by TSF 

    int VelocityCount = 0;
    int PressureCount = 0;
    double *DataVelocity, *DataPressure, *Data;
    int errcode;

    errcode = in_Epet->ExtractView(&Data);
    errcode = in_EpetVelocity.get()->ExtractView(&DataVelocity);
    errcode = in_EpetPressure.get()->ExtractView(&DataPressure);

    for (int i= 0 ; i < in_Epet->MyLength(); i++) {

      if (map_[i] == 1) {
	DataPressure[PressureCount++] = Data[i];
      }
      else {              
	DataVelocity[VelocityCount++] = Data[i];
      }
    }


    // Do the multiply

    Matrix_.apply(in_Tsf, out_Tsf);  // Operator overloading doesn't work
                                     // as a new copy of 'out_Tsf' is created.


  // put 'out_Tsf' data in 'out_Epet'

    VelocityCount = 0;
    PressureCount = 0;

    out_Epet->ExtractView(&Data);
    out_EpetVelocity.get()->ExtractView(&DataVelocity);
    out_EpetPressure.get()->ExtractView(&DataPressure);
    for (int i= 0 ; i < out_Epet->MyLength(); i++) {
      if (map_[i] == 1)  Data[i] = DataPressure[PressureCount++];
      else               Data[i] = DataVelocity[VelocityCount++];
    }

    return 0;
  }
  
  TSFLinearOperator TSFLinearOperator2EpetraRowMatrix::getTSF( ) { return(Matrix_); }
  int *  TSFLinearOperator2EpetraRowMatrix::getBlockAssignments() { return(map_); }


  int TSFLinearOperator2EpetraRowMatrix::initialize(TSFLinearOperator A_tsf,  Epetra_Map *Map_pet, int *map) {



    OwnMap_      = false;
    map_         = map;
    Matrix_      = A_tsf;
    OperatorDomainMap_ = Map_pet;

    rangeBlockSpace_  = A_tsf.range();
    domainBlockSpace_ = A_tsf.domain();

    numRangeBlocks_   = rangeBlockSpace_.numBlocks();
    numDomainBlocks_  = domainBlockSpace_.numBlocks();

    rangeLength_ = rangeBlockSpace_.dim();
    domainLength_ = domainBlockSpace_.dim();
    return(0);

  }

