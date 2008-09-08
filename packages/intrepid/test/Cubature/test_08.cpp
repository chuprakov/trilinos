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


/** \file
\brief  Unit test (CubatureDirect,CubatureTensor): correctness of
        integration of monomials for 3D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_CubatureSparse.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

using namespace Intrepid;


/*
  Monomial evaluation.
    in 1D, for point p(x)    : x^xDeg
    in 2D, for point p(x,y)  : x^xDeg * y^yDeg
    in 3D, for point p(x,y,z): x^xDeg * y^yDeg * z^zDeg
*/
double computeMonomial(Point<double> p, int xDeg, int yDeg=0, int zDeg=0) {
  double val = 1.0;
  int polydeg[3];
  polydeg[0] = xDeg; polydeg[1] = yDeg; polydeg[2] = zDeg;
  for (int i=0; i<p.getDim(); i++) {
    val *= std::pow(p.getCoordinates()[i],polydeg[i]);
  }
  return val;
}


/*
  Computes integrals of monomials over a given reference cell.
*/
void computeIntegral(Teuchos::Array<double>& testIntFixDeg, ECell cellType, int cubDegree) {

  Teuchos::RCP< Cubature<double> > myCub;  

  int ambientDim =  MultiCell<double>::getCellDim(cellType);

  switch (cellType) {

    case CELL_HEX:
        myCub = Teuchos::rcp(new CubatureSparse<double,3>(cubDegree));
      break;

    default:
      TEST_FOR_EXCEPTION((cellType != CELL_HEX),
                          std::invalid_argument,
                          ">>> ERROR (Unit Test -- Cubature -- 3D Monomial): Invalid cell type.");
  } // end switch

  int numCubPoints = myCub->getNumPoints();
  int numPolys     = (cubDegree+1)*(cubDegree+2)*(cubDegree+3)/6;

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;
  Teuchos::SerialDenseMatrix<int, double> functValues(numPolys, numCubPoints);

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  myCub->getCubature(cubPoints, cubWeights);

  int polyCt = 0;
  for (int xDeg=0; xDeg <= cubDegree; xDeg++) {
    for (int yDeg=0; yDeg <= cubDegree-xDeg; yDeg++) {
      for (int zDeg=0; zDeg <= cubDegree-xDeg-yDeg; zDeg++) {
        for (int i=0; i<numCubPoints; i++) {
          functValues(polyCt,i) = computeMonomial(cubPoints[i], xDeg, yDeg, zDeg);
        }
        polyCt++;
      }
    }
  }

  Teuchos::BLAS<int, double> myblas;
  int inc = 1;
  double alpha = 1.0;
  double beta  = 0.0;
  myblas.GEMV(Teuchos::NO_TRANS, numPolys, numCubPoints, alpha, functValues.values(), numPolys,
              &cubWeights[0], inc, beta, &testIntFixDeg[0], inc);
}


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
  << "|                 Unit Test (CubatureDirect,CubatureTensor)                   |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing integrals of monomials on reference cells in 3D            |\n" \
  << "|                         - using Level 2 BLAS -                              |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov),                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Matthew Keegan (mskeega@sandia.gov).                   |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: integrals of monomials in 3D (Level 2 BLAS version) using Sparse    |\n"\
  << "|           Grid Construction                                                 |\n"\
  << "===============================================================================\n";

  // >>> ASSUMPTION: max polynomial degree integrated exactly is the same for
  // >>>             edges and triangles !!!
  // internal variables:
  int                                      errorFlag = 0;
  int                                      polyCt = 0;
  int                                      offset = 0;
  Teuchos::Array< Teuchos::Array<double> > testInt;
  Teuchos::Array< Teuchos::Array<double> > analyticInt;
  Teuchos::Array<double>                   tmparray(1);
  double                                   reltol = 1.0e+04 * INTREPID_TOL;
  int                                      numPoly = (INTREPID_MAX_CUBATURE_DEGREE_EDGE+1)*
                                                     (INTREPID_MAX_CUBATURE_DEGREE_EDGE+2)*
                                                     (INTREPID_MAX_CUBATURE_DEGREE_EDGE+3)/6;
  testInt.assign(numPoly, tmparray);
  analyticInt.assign(numPoly, tmparray);

  // get names of files with analytic values
  std::string basedir = "./data";
  std::stringstream namestream;
  std::string filename;
  namestream << basedir << "/HEX_integrals" << ".dat";
  namestream >> filename;

  // reference cells tested
  ECell             testType = CELL_HEX;
  // format of data files with analytic values
  TypeOfExactData dataFormat = INTREPID_UTILS_FRACTION;

  // compute and compare integrals
  try {
    //for (int cellCt=0; cellCt < 3; cellCt++) {
      *outStream << "\nIntegrals of monomials on a reference " << MultiCell<double>::getCellName(testType) << ":\n";
      std::ifstream filecompare(&filename[0]);
      // compute integrals
      for (int cubDeg=0; cubDeg <= INTREPID_MAX_CUBATURE_DEGREE_SPARSE3D; cubDeg++) {
        int numMonomials = (cubDeg+1)*(cubDeg+2)*(cubDeg+3)/6; 
        testInt[cubDeg].resize(numMonomials);
        computeIntegral(testInt[cubDeg], testType, cubDeg);
      }
      // get analytic values
      if (filecompare.is_open()) {
        getAnalytic(analyticInt, filecompare, dataFormat);
        // close file
        filecompare.close();
      }
      // perform comparison
      for (int cubDeg=0; cubDeg <= INTREPID_MAX_CUBATURE_DEGREE_SPARSE3D; cubDeg++) {
        polyCt = 0;
        offset = 0;
        int oldErrorFlag = errorFlag;
        for (int xDeg=0; xDeg <= cubDeg; xDeg++) {
          for (int yDeg=0; yDeg <= cubDeg-xDeg; yDeg++) {
            for (int zDeg=0; zDeg <= cubDeg-xDeg-yDeg; zDeg++) {
              double abstol = ( analyticInt[polyCt+offset][0] == 0.0 ? reltol : std::fabs(reltol*analyticInt[polyCt+offset][0]) );
              double absdiff = std::fabs(analyticInt[polyCt+offset][0] - testInt[cubDeg][polyCt]);
              if (absdiff > abstol) {
                *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                           << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg
                           << " * z^" << std::setw(2) << zDeg << ":" << "   "
                           << std::scientific << std::setprecision(16)
                           << testInt[cubDeg][polyCt] << "   " << analyticInt[polyCt+offset][0] << "   "
                           << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
                errorFlag++;
                *outStream << std::right << std::setw(118) << "^^^^---FAILURE!\n";
              }
              polyCt++;
            }
            offset = offset + INTREPID_MAX_CUBATURE_DEGREE_EDGE - cubDeg;
          }
          offset = offset + (INTREPID_MAX_CUBATURE_DEGREE_EDGE - cubDeg)*(INTREPID_MAX_CUBATURE_DEGREE_EDGE - cubDeg + 1)/2;
        }
        *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg;
        if (errorFlag == oldErrorFlag)
         *outStream << " passed.\n";
        else
         *outStream << " failed.\n";
      }
      *outStream << "\n";
    //}  // end for cellCt
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  return errorFlag;
}
