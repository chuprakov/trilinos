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
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert C. Kirby (robert.c.kirby@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit tests for the Intrepid::Basis_F0_TRI_C1_FEM_FIAT class.
\author Created by R. Kirby via FIAT.
*/

#include "Intrepid_DefaultBasisFactory.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                 Unit Test (Basis_F0_TRI_C1_FEM_FIAT)                                   |\n" \
    << "|                                                                             |\n" \
    << "|     1) Basis creation, conversion of Dof tags into enumeration and back     |\n" \
    << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
    << "|                      Robert C. Kirby (robert.c.kirby@ttu.edu)               |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|  FIAT website:       http://www.fenics.org/fiat                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Basis creation, exception testing                                   |\n"\
    << "===============================================================================\n";

  int errorFlag = 0;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 1;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  endThrowNumber += 5;
#endif

  // Reference element points are the standard equispaced lattice of degree 1
  Teuchos::Array< Point<double> > elNodes;
  Point<double> tempPt(2, FRAME_REFERENCE);
  elNodes.assign(3, tempPt);  
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  elNodes[0] = Point<double>( 0.0000000000000000e+00 , 0.0000000000000000e+00 , FRAME_REFERENCE);
  elNodes[1] = Point<double>( 1.0000000000000000e+00 , 0.0000000000000000e+00 , FRAME_REFERENCE);
  elNodes[2] = Point<double>( 0.0000000000000000e+00 , 1.0000000000000000e+00 , FRAME_REFERENCE);
#endif


  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_0, 
                                                            CELL_TRI, 
                                                            RECONSTRUCTION_SPACE_COMPLETE, 
                                                            1, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);
    // Exception 1: DIV cannot be applied to scalar functions
    try {
      triBasis->getValues(vals, elNodes, OPERATOR_DIV);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 1----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    // Exceptions 2-6: all tags/bd Ids below are wrong and should cause getLocalDofEnumeration() and 
    // getLocalDofTag() to access invalid array elements thereby causing Teuchos bounds check exception
    try {      
      LocalDofTag myTag = {{3,0,0,1}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 2----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    try {
      LocalDofTag myTag = {{1,1,1,0}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 3----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      LocalDofTag myTag = {{0,3,0,0}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 4----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    }; 

    try {
      int bfId = 3;
      triBasis -> getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 5----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    try {
      int bfId = -1;
      triBasis -> getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 6----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  // Check if number of thrown exceptions matches the one we expect (1 + 5)
  if (TestForException_getThrowNumber() != endThrowNumber) {
    errorFlag++;
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"\
    << "===============================================================================\n";

  try{
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_0, 
                                                            CELL_TRI, 
                                                            RECONSTRUCTION_SPACE_COMPLETE, 
                                                            1, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);
    Teuchos::Array<LocalDofTag> allTags;
    allTags = triBasis -> getAllLocalDofTags();

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int        bfId  = triBasis -> getLocalDofEnumeration(allTags[i]);
      LocalDofTag myTag = triBasis -> getLocalDofTag(bfId);
       if( !( (myTag.tag_[0] == allTags[i].tag_[0]) &&
             (myTag.tag_[1] == allTags[i].tag_[1]) &&
             (myTag.tag_[2] == allTags[i].tag_[2]) &&
             (myTag.tag_[3] == allTags[i].tag_[3]) ) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getLocalDofEnumeration( {" 
          << allTags[i].tag_[0] << ", " 
          << allTags[i].tag_[1] << ", " 
          << allTags[i].tag_[2] << ", " 
          << allTags[i].tag_[3] << "}) = " << bfId <<" but \n";   
        *outStream << " getLocalDofTag(" << bfId << ") = { "
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "}\n";        
      }
    }

    // Now do the same but loop over basis functions
    for( int bfId = 0; bfId < triBasis -> getNumLocalDof(); bfId++) {
      LocalDofTag myTag  = triBasis -> getLocalDofTag(bfId);
      int myBfId = triBasis -> getLocalDofEnumeration(myTag);
      if( bfId != myBfId) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getLocalDofTag(" << bfId << ") = { "
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "} but getLocalDofEnumeration({" 
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "} ) = " << myBfId << "\n";
      }
    }
  }
  catch (std::logic_error err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 3: correctness of basis function values                                |\n"\
    << "===============================================================================\n";

  outStream -> precision(20);

  // VALUE: Each correct basis function at each correct point (point increasing fastest)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  double basisValues[] = { 1.0000000000000000e+00 , 2.7755575615628914e-17 , 0.0000000000000000e+00 , 2.7755575615628914e-17 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , 5.5511151231257827e-17 , 5.5511151231257827e-17 , 1.0000000000000000e+00 };
#endif

  // DERIVS: Components of derivatives increasing fastest, followed by points, then basis functions
#ifndef DOXYGEN_SHOULD_SKIP_THIS
double basisD1s[] = { -1.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 1.0000000000000000e+00 , -1.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 1.0000000000000000e+00 , -1.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 1.0000000000000000e+00 };
#endif


  // CURL: each correct values of the curls of the basis functions at the points (point increasing fastest
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  double basisCurls[] = { -1.0000000000000000e+00 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , -0.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , -0.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , -1.0000000000000000e+00 , 1.0000000000000000e+00 , -0.0000000000000000e+00 };
#endif


  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_0, 
                                                            CELL_TRI, 
                                                            RECONSTRUCTION_SPACE_COMPLETE, 
                                                            1, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);

    // Check VALUE of basis functions
    triBasis -> getValues(vals, elNodes, OPERATOR_VALUE);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisValues[i]) < INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisValues[i]) > INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisValues[i]) / std::abs(basisValues[i] ) > INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed value: " << vals[i] 
          << " but reference value: " << basisValues[i] << "\n";
      }
    }

    // Check GRAD of basis function
    triBasis -> getValues(vals, elNodes, OPERATOR_GRAD);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisD1s[i]) < 10.0 * INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisD1s[i]) > 10.0 * INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisD1s[i]) / std::abs(basisD1s[i] ) > 10.0 * INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed grad component: " << vals[i] 
          << " but reference grad component: " << basisD1s[i] << "\n";
      }
    }

    // Check CURL of basis function
    triBasis -> getValues(vals, elNodes, OPERATOR_CURL);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisCurls[i]) < 10.0 * INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisCurls[i]) > 10.0 * INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisCurls[i]) / std::abs(basisCurls[i] ) > 10.0 * INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed curl component: " << vals[i] 
          << " but reference curl component: " << basisCurls[i] << "\n";
      }
    }


    // Check D1 of basis function
    triBasis -> getValues(vals, elNodes, OPERATOR_D1);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisD1s[i]) < pow(15.0,1) * INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisD1s[i]) > pow(15.0,1) * INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisD1s[i]) / std::abs(basisD1s[i] ) > pow(15.0,1) * INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed D1 component: " << vals[i] 
          << " but reference D1 component: " << basisD1s[i] << "\n";
      }
    }

      // Check D2 of basis function: should be zero 
      triBasis -> getValues(vals, elNodes, OPERATOR_D2);
      for (int i=0; i < vals.getSize(); i++) {
        if ( std::abs(vals[i])  > pow( 15.0 , 2 ) * INTREPID_FIAT_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // Get the multi-index of the value where the error is and the operator order
          Teuchos::Array<int> myIndex;
          vals.getMultiIndex(myIndex,i);
          *outStream << " At multi-index { ";
          for(int j = 0; j < vals.getRank(); j++) {
            *outStream << myIndex[j] << " ";
          }
          *outStream << "}  computed D2 component: " << vals[i] 
            << " but reference D2 component:  0 \n";
        }
      }

      // Check D3 of basis function: should be zero 
      triBasis -> getValues(vals, elNodes, OPERATOR_D3);
      for (int i=0; i < vals.getSize(); i++) {
        if ( std::abs(vals[i])  > pow( 15.0 , 3 ) * INTREPID_FIAT_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // Get the multi-index of the value where the error is and the operator order
          Teuchos::Array<int> myIndex;
          vals.getMultiIndex(myIndex,i);
          *outStream << " At multi-index { ";
          for(int j = 0; j < vals.getRank(); j++) {
            *outStream << myIndex[j] << " ";
          }
          *outStream << "}  computed D3 component: " << vals[i] 
            << " but reference D3 component:  0 \n";
        }
      }


  }
  
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}

