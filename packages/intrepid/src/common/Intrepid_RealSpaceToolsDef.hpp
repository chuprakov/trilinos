// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_RealSpaceToolsDef.hpp
    \brief  Definition file for utility classes providing basic linear algebra functionality.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/


namespace Intrepid {



template<class Scalar>
void RealSpaceTools<Scalar>::absval(Scalar* absArray, const Scalar* inArray, const int size) {
  for (int i=0; i<size; i++) {
    absArray[i] = std::abs(inArray[i]);
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::absval(Scalar* inoutAbsArray, const int size) {
  for (int i=0; i<size; i++) {
    inoutAbsArray[i] = std::abs(inoutAbsArray[i]);
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::absval(ArrayScalar & absArray, const ArrayScalar & inArray) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.getRank() != absArray.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::absval): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.getRank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.getDimension(i) != absArray.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::absval): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.getSize(); i++) {
    absArray[i] = std::abs(inArray[i]);
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::absval(ArrayScalar & inoutAbsArray) {
  for (int i=0; i<inoutAbsArray.getSize(); i++) {
    inoutAbsArray[i] = std::abs(inoutAbsArray[i]);
  }
}



template<class Scalar>
Scalar RealSpaceTools<Scalar>::vectorNorm(const Scalar* inVec, const int dim, const ENorm normType) {
  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:
      for(int i = 0; i < dim; i++){
        temp += inVec[i]*inVec[i];
      }
      temp = std::sqrt(temp);
      break;
    case NORM_INF:
      temp = std::abs(inVec[0]);
      for(int i = 1; i < dim; i++){
        Scalar absData = std::abs(inVec[i]);
        if (temp < absData) temp = absData;
      }
      break;
    case NORM_ONE:
      for(int i = 0; i < dim; i++){
        temp += std::abs(inVec[i]);
      }
      break;
    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
  return temp;
}



template<class Scalar>
template<class VecArray>
Scalar RealSpaceTools<Scalar>::vectorNorm(const VecArray & inVec, const ENorm normType) {

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inVec.getRank() != 1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Vector argument must have rank 1!");
#endif

  int size = inVec.getSize();

  Scalar temp = (Scalar)0;
  switch(normType) {
    case NORM_TWO:
      for(int i = 0; i < size; i++){
        temp += inVec[i]*inVec[i];
      }
      temp = std::sqrt(temp);
      break;
    case NORM_INF:
      temp = std::abs(inVec[0]);
      for(int i = 1; i < size; i++){
        Scalar absData = std::abs(inVec[i]);
        if (temp < absData) temp = absData;
      }
      break;
    case NORM_ONE:
      for(int i = 0; i < size; i++){
        temp += std::abs(inVec[i]);
      }
      break;
    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
  return temp;
}



template<class Scalar>
template<class NormArray, class VecArray>
void RealSpaceTools<Scalar>::vectorNorm(NormArray & normArray, const VecArray & inVecs, const ENorm normType) {

  int arrayRank = inVecs.getRank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != normArray.getRank()+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Ranks of norm and vector array arguments are incompatible!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::vectorNorm): Rank of vector array must be 2 or 3!");
    for (int i=0; i<arrayRank-1; i++) {
      TEST_FOR_EXCEPTION( ( inVecs.getDimension(i) != normArray.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Dimensions of norm and vector arguments do not agree!");
    }
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inVecs.getDimension(arrayRank-1); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 3:
      dim_i0 = inVecs.getDimension(0);
      dim_i1 = inVecs.getDimension(1);
      break;
    case 2:
      dim_i1 = inVecs.getDimension(0);
      break;
  }

  switch(normType) {
    case NORM_TWO: {
      int offset_i0, offset, normOffset;
      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          for(int i = 0; i < dim; i++){
            temp += inVecs[i]*inVecs[i];
          }
          normArray[normOffset] = std::sqrt(temp);
        }
      }
      break;
    } // case NORM_TWO

    case NORM_INF: {
      int offset_i0, offset, normOffset;
      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          temp = std::abs(inVecs[0]);
          for(int i = 1; i < dim; i++){
            Scalar absData = std::abs(inVecs[i]);
            if (temp < absData) temp = absData;
          }
          normArray[normOffset] = temp;
        }
      }
      break;
    } // case NORM_INF

    case NORM_ONE: {
      int offset_i0, offset, normOffset;
      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset      = offset_i0 + i1;
          normOffset  = offset;
          offset     *= dim;
          Scalar temp = (Scalar)0;
          for(int i = 0; i < dim; i++){
            temp += std::abs(inVecs[i]);
          }
          normArray[normOffset] = temp;
        }
      }
      break;
    } // case NORM_ONE

    default:
      TEST_FOR_EXCEPTION( ( (normType != NORM_TWO) && (normType != NORM_INF) && (normType != NORM_ONE) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::transpose(Scalar* transposeMat, const Scalar* inMat, const int dim) {
  for(int i=0; i < dim; i++){
    transposeMat[i*dim+i]=inMat[i*dim+i];    // Set diagonal elements
    for(int j=i+1; j < dim; j++){
      transposeMat[i*dim+j]=inMat[j*dim+i];  // Set off-diagonal elements
      transposeMat[j*dim+i]=inMat[i*dim+j];
    }
  }
}



template<class Scalar>
template<class MatArray>
void RealSpaceTools<Scalar>::transpose(MatArray & transposeMats, const MatArray & inMats) {
  int arrayRank = inMats.getRank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != transposeMats.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::transpose): Matrix array arguments do not have identical ranks!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::transpose): Rank of matrix array must be 2, 3, or 4!");
    for (int i=0; i<arrayRank; i++) {
      TEST_FOR_EXCEPTION( ( inMats.getDimension(i) != transposeMats.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::transpose): Dimensions of matrix arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.getDimension(arrayRank-2) != inMats.getDimension(arrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::transpose): Matrices are not square!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.getDimension(arrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 4:
      dim_i0 = inMats.getDimension(0);
      dim_i1 = inMats.getDimension(1);
      break;
    case 3:
      dim_i1 = inMats.getDimension(0);
      break;
  }

  int offset_i0, offset;

  for (int i0=0; i0<dim_i0; i0++) {
    offset_i0 = i0*dim_i1;
    for (int i1=0; i1<dim_i1; i1++) {
      offset  = offset_i0 + i1;
      offset *= (dim*dim);

      for(int i=0; i < dim; i++){
        transposeMats[offset+i*dim+i]=inMats[offset+i*dim+i];    // Set diagonal elements
        for(int j=i+1; j < dim; j++){
          transposeMats[offset+i*dim+j]=inMats[offset+j*dim+i];  // Set off-diagonal elements
          transposeMats[offset+j*dim+i]=inMats[offset+i*dim+j];
        }
      }

    } // i1
  } // i0

}



template<class Scalar>
void RealSpaceTools<Scalar>::inverse(Scalar* inverseMat, const Scalar* inMat, const int dim) {

  switch(dim) {
    case 3: {
      int i, j, rowID = 0, colID = 0;
      int rowperm[3]={0,1,2};
      int colperm[3]={0,1,2}; // Complete pivoting
      Scalar emax(0);

      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( inMat[i*3+j] ) >  emax){
            rowID = i;  colID = j; emax = std::abs( inMat[i*3+j] );
          }
        }
      }
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
      if( rowID ){
        rowperm[0] = rowID;
        rowperm[rowID] = 0;
      }
      if( colID ){
        colperm[0] = colID;
        colperm[colID] = 0;
      }
      Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          B[i][j] = inMat[rowperm[i]*3+colperm[j]];
        }
      }
      B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
      for(i=0; i < 2; ++i){
        for(j=0; j < 2; ++j){
          S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
        }
      }
      Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif

      Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
      Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

      for(j=0; j<2;j++)
        Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
      for(i=0; i<2;i++)
        Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

      Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
      Bi[1][1] =  Si[0][0];
      Bi[1][2] =  Si[0][1];
      Bi[2][1] =  Si[1][0];
      Bi[2][2] =  Si[1][1];
      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          inverseMat[i*3+j] = Bi[colperm[i]][rowperm[j]]; // set inverse
        }
      }
      break;
    } // case 3

    case 2: {

      Scalar determinant    = inMat[0]*inMat[3]-inMat[1]*inMat[2];;
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (inMat[0]==(Scalar)0) && (inMat[1]==(Scalar)0) &&
                            (inMat[2]==(Scalar)0) && (inMat[3]==(Scalar)0) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
      TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
      inverseMat[0] =   inMat[3] / determinant;
      inverseMat[1] = - inMat[1] / determinant;
      //
      inverseMat[2] = - inMat[2] / determinant;
      inverseMat[3] =   inMat[0] / determinant;
      break;
    } // case 2

    case 1: {
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( inMat[0] == (Scalar)0 ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
      inverseMat[0] = (Scalar)1 / inMat[0];
      break;
    } // case 1

  } // switch (dim)
}



template<class Scalar>
template<class MatArray>
void RealSpaceTools<Scalar>::inverse(MatArray & inverseMats, const MatArray & inMats) {

  int arrayRank = inMats.getRank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != inverseMats.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Matrix array arguments do not have identical ranks!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Rank of matrix array must be 2, 3, or 4!");
    for (int i=0; i<arrayRank; i++) {
      TEST_FOR_EXCEPTION( ( inMats.getDimension(i) != inverseMats.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::inverse): Dimensions of matrix arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.getDimension(arrayRank-2) != inMats.getDimension(arrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Matrices are not square!");
    TEST_FOR_EXCEPTION( ( (inMats.getDimension(arrayRank-2) < 1) || (inMats.getDimension(arrayRank-2) > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::inverse): Spatial dimension must be 1, 2, or 3!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.getDimension(arrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 4:
      dim_i0 = inMats.getDimension(0);
      dim_i1 = inMats.getDimension(1);
      break;
    case 3:
      dim_i1 = inMats.getDimension(0);
      break;
  }

  switch(dim) {
    case 3: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;
          offset *= 9;

          int i, j, rowID = 0, colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMats[offset+i*3+j] ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMats[offset+i*3+j] );
              }
            }
          }
#ifdef HAVE_INTREPID_DEBUG
          TEST_FOR_EXCEPTION( ( emax == (Scalar)0 ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
          if( rowID ){
            rowperm[0] = rowID;
            rowperm[rowID] = 0;
          }
          if( colID ){
            colperm[0] = colID;
            colperm[colID] = 0;
          }
          Scalar B[3][3], S[2][2], Bi[3][3]; // B=rowperm inMat colperm, S=Schur complement(Boo)
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              B[i][j] = inMats[offset+rowperm[i]*3+colperm[j]];
            }
          }
          B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
          for(i=0; i < 2; ++i){
            for(j=0; j < 2; ++j){
              S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
            }
          }
          Scalar detS = S[0][0]*S[1][1]- S[0][1]*S[1][0], Si[2][2];
#ifdef HAVE_INTREPID_DEBUG
          TEST_FOR_EXCEPTION( ( detS == (Scalar)0 ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif

          Si[0][0] =  S[1][1]/detS;                  Si[0][1] = -S[0][1]/detS;
          Si[1][0] = -S[1][0]/detS;                  Si[1][1] =  S[0][0]/detS;

          for(j=0; j<2;j++)
            Bi[0][j+1] = -( B[0][1]*Si[0][j] + B[0][2]* Si[1][j])/B[0][0];
          for(i=0; i<2;i++)
            Bi[i+1][0] = -(Si[i][0]*B[1][0] + Si[i][1]*B[2][0]);

          Bi[0][0] =  ((Scalar)1/B[0][0])-Bi[0][1]*B[1][0]-Bi[0][2]*B[2][0];
          Bi[1][1] =  Si[0][0];
          Bi[1][2] =  Si[0][1];
          Bi[2][1] =  Si[1][0];
          Bi[2][2] =  Si[1][1];
          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              inverseMats[offset+i*3+j] = Bi[colperm[i]][rowperm[j]]; // set inverse
            }
          }
        } // for i1
      } // for i0
      break;
    } // case 3

    case 2: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;;
          offset *= 4;

          Scalar determinant    = inMats[offset]*inMats[offset+3]-inMats[offset+1]*inMats[offset+2];;
#ifdef HAVE_INTREPID_DEBUG
          TEST_FOR_EXCEPTION( ( (inMats[offset]==(Scalar)0)   && (inMats[offset+1]==(Scalar)0) &&
                                (inMats[offset+2]==(Scalar)0) && (inMats[offset+3]==(Scalar)0) ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
          TEST_FOR_EXCEPTION( ( determinant == (Scalar)0 ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
          inverseMats[offset]   = inMats[offset+3] / determinant;
          inverseMats[offset+1] = - inMats[offset+1] / determinant;
          //
          inverseMats[offset+2] = - inMats[offset+2] / determinant;
          inverseMats[offset+3] =   inMats[offset] / determinant;
        } // for i1
      } // for i0
      break;
    } // case 2

    case 1: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;;
#ifdef HAVE_INTREPID_DEBUG
          TEST_FOR_EXCEPTION( ( inMats[offset] == (Scalar)0 ),
                              std::invalid_argument,
                              ">>> ERROR (Matrix): Inverse of a zero matrix is undefined!");
#endif
          inverseMats[offset] = (Scalar)1 / inMats[offset];
        } // for i1
      } // for i2
      break;
    } // case 1

  } // switch (dim)
}



template<class Scalar>
Scalar RealSpaceTools<Scalar>::det(const Scalar* inMat, const int dim) {
  Scalar determinant(0);

  switch (dim) {
    case 3: {
      int i,j,rowID = 0;
      int colID = 0;
      int rowperm[3]={0,1,2};
      int colperm[3]={0,1,2}; // Complete pivoting
      Scalar emax(0);

      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( inMat[i*dim+j] ) >  emax){
            rowID = i;  colID = j; emax = std::abs( inMat[i*dim+j] );
          }
        }
      }
      if( emax > 0 ){
        if( rowID ){
          rowperm[0] = rowID;
          rowperm[rowID] = 0;
        }
        if( colID ){
          colperm[0] = colID;
          colperm[colID] = 0;
        }
        Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
        for(i=0; i < 3; ++i){
          for(j=0; j < 3; ++j){
            B[i][j] = inMat[rowperm[i]*dim+colperm[j]];
          }
        }
        B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
        for(i=0; i < 2; ++i){
          for(j=0; j < 2; ++j){
            S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
          }
        }
        determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
        if( rowID ) determinant = -determinant;
        if( colID ) determinant = -determinant;
      }
      break;
    } // case 3

    case 2:
      determinant = inMat[0]*inMat[3]-
                    inMat[1]*inMat[2];
      break;

    case 1:
      determinant = inMat[0];
      break;

    default:
      TEST_FOR_EXCEPTION( ( (dim != 1) && (dim != 2) && (dim != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  } // switch (dim)

  return determinant;
}



template<class Scalar>
template<class ArrayScalar>
Scalar RealSpaceTools<Scalar>::det(const ArrayScalar & inMat) {

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( (inMat.getRank() != 2),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Rank of matrix argument must be 2!");
    TEST_FOR_EXCEPTION( ( inMat.getDimension(0) != inMat.getDimension(1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Matrix is not square!");
    TEST_FOR_EXCEPTION( ( (inMat.getDimension(0) < 1) || (inMat.getDimension(0) > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif

  int dim = inMat.getDimension(0);
  Scalar determinant(0);

  switch (dim) {
    case 3: {
      int i,j,rowID = 0;
      int colID = 0;
      int rowperm[3]={0,1,2};
      int colperm[3]={0,1,2}; // Complete pivoting
      Scalar emax(0);

      for(i=0; i < 3; ++i){
        for(j=0; j < 3; ++j){
          if( std::abs( inMat[i*dim+j] ) >  emax){
            rowID = i;  colID = j; emax = std::abs( inMat[i*dim+j] );
          }
        }
      }
      if( emax > 0 ){
        if( rowID ){
          rowperm[0] = rowID;
          rowperm[rowID] = 0;
        }
        if( colID ){
          colperm[0] = colID;
          colperm[colID] = 0;
        }
        Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
        for(i=0; i < 3; ++i){
          for(j=0; j < 3; ++j){
            B[i][j] = inMat[rowperm[i]*dim+colperm[j]];
          }
        }
        B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
        for(i=0; i < 2; ++i){
          for(j=0; j < 2; ++j){
            S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
          }
        }
        determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
        if( rowID ) determinant = -determinant;
        if( colID ) determinant = -determinant;
      }
      break;
    } // case 3

    case 2:
      determinant = inMat[0]*inMat[3]-
                    inMat[1]*inMat[2];
      break;

    case 1:
      determinant = inMat[0];
      break;

    default:
      TEST_FOR_EXCEPTION( ( (dim != 1) && (dim != 2) && (dim != 3) ),
                          std::invalid_argument,
                          ">>> ERROR (Matrix): Invalid matrix dimension.");
  } // switch (dim)

  return determinant;
}




template<class Scalar>
template<class DetArray, class MatArray>
void RealSpaceTools<Scalar>::det(DetArray & detArray, const MatArray & inMats) {

  int matArrayRank = inMats.getRank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( matArrayRank != detArray.getRank()+2 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Determinant and matrix array arguments do not have compatible ranks!");
    TEST_FOR_EXCEPTION( ( (matArrayRank < 3) || (matArrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Rank of matrix array must be 3 or 4!");
    for (int i=0; i<matArrayRank-2; i++) {
      TEST_FOR_EXCEPTION( ( inMats.getDimension(i) != detArray.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::det): Dimensions of determinant and matrix array arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.getDimension(matArrayRank-2) != inMats.getDimension(matArrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Matrices are not square!");
    TEST_FOR_EXCEPTION( ( (inMats.getDimension(matArrayRank-2) < 1) || (inMats.getDimension(matArrayRank-2) > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.getDimension(matArrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(matArrayRank) {
    case 4:
      dim_i0 = inMats.getDimension(0);
      dim_i1 = inMats.getDimension(1);
      break;
    case 3:
      dim_i1 = inMats.getDimension(0);
      break;
  }

  switch(dim) {
    case 3: {
      int offset_i0, offset, detOffset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset     = offset_i0 + i1;
          detOffset  = offset;
          offset    *= 9;

          int i,j,rowID = 0;
          int colID = 0;
          int rowperm[3]={0,1,2};
          int colperm[3]={0,1,2}; // Complete pivoting
          Scalar emax(0), determinant(0);

          for(i=0; i < 3; ++i){
            for(j=0; j < 3; ++j){
              if( std::abs( inMats[offset+i*3+j] ) >  emax){
                rowID = i;  colID = j; emax = std::abs( inMats[offset+i*3+j] );
              }
            }
          }
          if( emax > 0 ){
            if( rowID ){
              rowperm[0] = rowID;
              rowperm[rowID] = 0;
            }
            if( colID ){
              colperm[0] = colID;
              colperm[colID] = 0;
            }
            Scalar B[3][3], S[2][2]; // B=rowperm inMat colperm, S=Schur complement(Boo)
            for(i=0; i < 3; ++i){
              for(j=0; j < 3; ++j){
                B[i][j] = inMats[offset+rowperm[i]*3+colperm[j]];
              }
            }
            B[1][0] /= B[0][0]; B[2][0] /= B[0][0];// B(:,0)/=pivot
            for(i=0; i < 2; ++i){
              for(j=0; j < 2; ++j){
                S[i][j] = B[i+1][j+1] - B[i+1][0] * B[0][j+1]; // S = B -z*y'
              }
            }
            determinant = B[0][0] * (S[0][0] * S[1][1] - S[0][1] * S[1][0]); // det(B)
            if( rowID ) determinant = -determinant;
            if( colID ) determinant = -determinant;
          }
          detArray[detOffset] = determinant;
        } // for i1
      } // for i0
      break;
    } // case 3

    case 2: {
      int offset_i0, offset, detOffset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset     = offset_i0 + i1;
          detOffset  = offset;
          offset    *= 4;

          detArray[detOffset] = inMats[offset]*inMats[offset+3]-inMats[offset+1]*inMats[offset+2];;
        } // for i1
      } // for i0
      break;
    } // case 2

    case 1: {
      int offset_i0, offset;

      for (int i0=0; i0<dim_i0; i0++) {
        offset_i0 = i0*dim_i1;
        for (int i1=0; i1<dim_i1; i1++) {
          offset  = offset_i0 + i1;;
          detArray[offset] = inMats[offset];
        } // for i1
      } // for i2
      break;
    } // case 1

  } // switch (dim)
}



template<class Scalar>
void RealSpaceTools<Scalar>::add(Scalar* sumArray, const Scalar* inArray1, const Scalar* inArray2, const int size) {
  for (int i=0; i<size; i++) {
    sumArray[i] = inArray1[i] + inArray2[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::add(Scalar* inoutSumArray, const Scalar* inArray, const int size) {
  for (int i=0; i<size; i++) {
    inoutSumArray[i] += inArray[i];
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::add(ArrayScalar & sumArray, const ArrayScalar & inArray1, const ArrayScalar & inArray2) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( (inArray1.getRank() != inArray2.getRank()) || (inArray1.getRank() != sumArray.getRank()) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
    for (int i=0; i<inArray1.getRank(); i++) {
      TEST_FOR_EXCEPTION( ( (inArray1.getDimension(i) != inArray2.getDimension(i)) || (inArray1.getDimension(i) != sumArray.getDimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray1.getSize(); i++) {
    sumArray[i] = inArray1[i] + inArray2[i];
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::add(ArrayScalar & inoutSumArray, const ArrayScalar & inArray) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.getRank() != inoutSumArray.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.getRank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.getDimension(i) != inoutSumArray.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.getSize(); i++) {
    inoutSumArray[i] += inArray[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::subtract(Scalar* diffArray, const Scalar* inArray1, const Scalar* inArray2, const int size) {
  for (int i=0; i<size; i++) {
    diffArray[i] = inArray1[i] - inArray2[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::subtract(Scalar* inoutDiffArray, const Scalar* inArray, const int size) {
  for (int i=0; i<size; i++) {
    inoutDiffArray[i] -= inArray[i];
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::subtract(ArrayScalar & diffArray, const ArrayScalar & inArray1, const ArrayScalar & inArray2) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( (inArray1.getRank() != inArray2.getRank()) || (inArray1.getRank() != diffArray.getRank()) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
    for (int i=0; i<inArray1.getRank(); i++) {
      TEST_FOR_EXCEPTION( ( (inArray1.getDimension(i) != inArray2.getDimension(i)) || (inArray1.getDimension(i) != diffArray.getDimension(i)) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray1.getSize(); i++) {
    diffArray[i] = inArray1[i] - inArray2[i];
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::subtract(ArrayScalar & inoutDiffArray, const ArrayScalar & inArray) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.getRank() != inoutDiffArray.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.getRank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.getDimension(i) != inoutDiffArray.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.getSize(); i++) {
    inoutDiffArray[i] -= inArray[i];
  }
}




template<class Scalar>
void RealSpaceTools<Scalar>::scale(Scalar* scaledArray, const Scalar* inArray, const int size, const Scalar scalar) {
  for (int i=0; i<size; i++) {
    scaledArray[i] = scalar*inArray[i];
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::scale(Scalar* inoutScaledArray, const int size, const Scalar scalar) {
  for (int i=0; i<size; i++) {
    inoutScaledArray[i] *= scalar;
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::scale(ArrayScalar & scaledArray, const ArrayScalar & inArray, const Scalar scalar) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( inArray.getRank() != scaledArray.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::scale): Array arguments must have identical ranks!");
    for (int i=0; i<inArray.getRank(); i++) {
      TEST_FOR_EXCEPTION( ( inArray.getDimension(i) != scaledArray.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::scale): Dimensions of array arguments do not agree!");
    }
#endif

  for (int i=0; i<inArray.getSize(); i++) {
    scaledArray[i] = scalar*inArray[i];
  }
}



template<class Scalar>
template<class ArrayScalar>
void RealSpaceTools<Scalar>::scale(ArrayScalar & inoutScaledArray, const Scalar scalar) {
  for (int i=0; i<inoutScaledArray.getSize(); i++) {
    inoutScaledArray[i] *= scalar;
  }
}




template<class Scalar>
Scalar RealSpaceTools<Scalar>::dot(const Scalar* inArray1, const Scalar* inArray2, const int size) {
  Scalar dot(0);
  for (int i=0; i<size; i++) {
    dot += inArray1[i]*inArray2[i];
  }
  return dot;  
}



template<class Scalar>
template<class VecArray>
Scalar RealSpaceTools<Scalar>::dot(const VecArray & inVec1, const VecArray & inVec2) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( (inVec1.getRank() != 1) || (inVec2.getRank() != 1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Vector arguments must have rank 1!");
    TEST_FOR_EXCEPTION( ( inVec1.getDimension(0) != inVec2.getDimension(0) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Dimensions of vector arguments must agree!");
#endif

  Scalar dot(0);
  for (int i=0; i<inVec1.getSize(); i++) {
    dot += inVec1[i]*inVec2[i];
  }
  return dot;  

}



template<class Scalar>
template<class DotArray, class VecArray>
void RealSpaceTools<Scalar>::dot(DotArray & dotArray, const VecArray & inVecs1, const VecArray & inVecs2) {

  int arrayRank = inVecs1.getRank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( arrayRank != dotArray.getRank()+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Ranks of norm and vector array arguments are incompatible!");
    TEST_FOR_EXCEPTION( ( arrayRank != inVecs2.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Ranks of input vector arguments must be identical!");
    TEST_FOR_EXCEPTION( ( (arrayRank < 2) || (arrayRank > 3) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::dot): Rank of input vector arguments must be 2 or 3!");
    for (int i=0; i<arrayRank; i++) {
      TEST_FOR_EXCEPTION( ( inVecs1.getDimension(i) != inVecs2.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::dot): Dimensions of input vector arguments do not agree!");
    }
    for (int i=0; i<arrayRank-1; i++) {
      TEST_FOR_EXCEPTION( ( inVecs1.getDimension(i) != dotArray.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::dot): Dimensions of dot-product and vector arrays do not agree!");
    }
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inVecs1.getDimension(arrayRank-1); // spatial dimension

  // determine i0 and i1 dimensions
  switch(arrayRank) {
    case 3:
      dim_i0 = inVecs1.getDimension(0);
      dim_i1 = inVecs1.getDimension(1);
      break;
    case 2:
      dim_i1 = inVecs1.getDimension(0);
      break;
  }

  int offset_i0, offset, dotOffset;
  for (int i0=0; i0<dim_i0; i0++) {
    offset_i0 = i0*dim_i1;
    for (int i1=0; i1<dim_i1; i1++) {
      offset      = offset_i0 + i1;
      dotOffset   = offset;
      offset     *= dim;
      Scalar dot(0);
      for (int i=0; i<dim; i++) {
        dot += inVecs1[offset+i]*inVecs2[offset+i];
      }
      dotArray[dotOffset] = dot;
    }
  }
}



template<class Scalar>
void RealSpaceTools<Scalar>::matvec(Scalar* matVec, const Scalar* inMat, const Scalar* inVec, const int dim) {
  for (int i=0; i<dim; i++) {
    Scalar sumdot(0);
    for (int j=0; j<dim; j++) {
      sumdot += inMat[i*dim+j]*inVec[j];
    }
    matVec[i] = sumdot; 
  }
}



template<class Scalar>
template<class MatArray, class VecArray>
void RealSpaceTools<Scalar>::matvec(VecArray & matVecs, const MatArray & inMats, const VecArray & inVecs) {
  int matArrayRank = inMats.getRank();

#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( matArrayRank != inVecs.getRank()+1 ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Vector and matrix array arguments do not have compatible ranks!");
    TEST_FOR_EXCEPTION( ( (matArrayRank < 3) || (matArrayRank > 4) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Rank of matrix array must be 3 or 4!");
    TEST_FOR_EXCEPTION( ( matVecs.getRank() != inVecs.getRank() ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Vector arrays must be have the same rank!");
    for (int i=0; i<matArrayRank-1; i++) {
      TEST_FOR_EXCEPTION( ( inMats.getDimension(i) != inVecs.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector and matrix array arguments do not agree!");
    }
    for (int i=0; i<inVecs.getRank(); i++) {
      TEST_FOR_EXCEPTION( ( matVecs.getDimension(i) != inVecs.getDimension(i) ),
                          std::invalid_argument,
                          ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector array arguments do not agree!");
    }
    TEST_FOR_EXCEPTION( ( inMats.getDimension(matArrayRank-2) != inMats.getDimension(matArrayRank-1) ),
                        std::invalid_argument,
                        ">>> ERROR (RealSpaceTools::matvec): Matrices are not square!");
#endif

  int dim_i0 = 1; // first  index dimension (e.g. cell index)
  int dim_i1 = 1; // second index dimension (e.g. point index)
  int dim    = inMats.getDimension(matArrayRank-2); // spatial dimension

  // determine i0 and i1 dimensions
  switch(matArrayRank) {
    case 4:
      dim_i0 = inMats.getDimension(0);
      dim_i1 = inMats.getDimension(1);
      break;
    case 3:
      dim_i1 = inMats.getDimension(0);
      break;
  }

  int offset_i0, offset, vecOffset;

  for (int i0=0; i0<dim_i0; i0++) {
    offset_i0 = i0*dim_i1;
    for (int i1=0; i1<dim_i1; i1++) {
      offset     = offset_i0 + i1;
      vecOffset  = offset*dim;
      offset     = vecOffset*dim;

      for (int i=0; i<dim; i++) {
        Scalar sumdot(0);
        for (int j=0; j<dim; j++) {
          sumdot += inMats[offset+i*dim+j]*inVecs[vecOffset+j];
        }
        matVecs[vecOffset+i] = sumdot;
      }
    }
  }
}



} // namespace Intrepid
