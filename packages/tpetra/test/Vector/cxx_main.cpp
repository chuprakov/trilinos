// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#include "../tpetra_test_util.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI

// function prototype
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);
void checkOutputs(); 
// checkOutputs does not need to be passed verbose/debug 
// because it is only called if debug is true

int main(int argc, char* argv[]) {
	int myImageID = 0; // assume we are on serial
	int numImages = 1; // if MPI, will be reset later
  
	// initialize MPI if needed
#ifdef TPETRA_MPI
	numImages = -1;
	myImageID = -1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numImages);
	MPI_Comm_rank(MPI_COMM_WORLD, &myImageID);
#endif // TPETRA_MPI

	// initialize verbose & debug flags
	bool verbose = false;
	bool debug = false;
	if(argc > 1) {
		if(argv[1][0] == '-' && argv[1][1] == 'v')
			verbose = true;
		if(argv[1][0] == '-' && argv[1][1] == 'd') {
			debug = true;
			verbose = true;
		}
	}
  
	// change verbose to only be true on Image 0
	verbose = (verbose && (myImageID == 0));
  
	// start the testing
	if(verbose) outputStartMessage("Vector");
	int ierr = 0;
  
	// call the actual test routines
	ierr += unitTests<int, float>(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);
    
	if(debug) 
  		checkOutputs();

	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Vector", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
	int ierr = 0;
	int returnierr = 0;
#ifdef TPETRA_MPI
	const Tpetra::MpiPlatform<OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
	const Tpetra::MpiPlatform<OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
	const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
	const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif // TPETRA_MPI
	if(verbose) cout << "Starting unit tests for Vector<" 
					 << Teuchos::OrdinalTraits<OrdinalType>::name() << "," 
					 << Teuchos::ScalarTraits<ScalarType>::name() << ">." << endl;
	
	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
	if(verbose) cout << "Starting code coverage section..." << endl;
	
	// constructors
	if(debug) cout << "Constructors..." << endl;
	// taking a VectorSpace
	OrdinalType const ESlength = 10;
	OrdinalType const ESindexBase = 2;
	Tpetra::ElementSpace<OrdinalType> elementspace(ESlength, ESindexBase, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
	Tpetra::Vector<OrdinalType, ScalarType> vector(vectorspace);
	// taking a VectorSpace and a user array of entries
	OrdinalType myLength = vectorspace.getNumMyEntries();
	std::vector<ScalarType> scalarArray(myLength); // allocate to size myLength
	for(OrdinalType i = 0; i < myLength; i++)
		scalarArray[i] = Teuchos::ScalarTraits<ScalarType>::random();

	ScalarType* ptr = 0;
	if(!scalarArray.empty())
		ptr = &scalarArray[0];
	Tpetra::Vector<OrdinalType, ScalarType> vector1a(ptr, myLength, vectorspace);
	// cpy ctr
	Tpetra::Vector<OrdinalType, ScalarType> v2(vector);
	
	// print
	if(debug) {
		cout << "Overloaded << operator..." << endl;
		cout << vector << endl;
	}
	
	// attribute access
	if(debug) cout << "Attribute access methods..." << endl;
	OrdinalType temp = 0;
	temp = vector.getNumGlobalEntries();
	temp = vector.getNumMyEntries();
	
	// element access - [] operator
	if(debug) cout << "Element access methods..." << endl;
	if(vector.getNumMyEntries() > 1) {
		ScalarType const temp1 = vector[1];
		vector[0] = temp1;
	}
	
	// set all to scalar
	if(debug) cout << "setAllToScalar..." << endl;
	ScalarType scalar3 = Teuchos::ScalarTraits<ScalarType>::one();
	scalar3 += (scalar3 + scalar3); // 1 += (1+1) , should equal 3
	vector.setAllToScalar(scalar3);
	
	// set all to random
	if(debug) cout << "setAllToRandom..." << endl;
	vector.setAllToRandom();
	
	// math functions
	if(debug) cout << "Mathematical operations..." << endl;
	
	// dot product
	if(debug) cout << "dot product..." << endl;
	vector.dotProduct(v2); // throw away return value
	// absolute value
	if(debug) cout << "absolute value..." << endl;
	vector.absoluteValue(v2);
	// reciprocal
	if(debug) cout << "reciprocal..." << endl;
	vector.reciprocal(v2);
	// scale
	if(debug) cout << "scale..." << endl;
	vector.scale(Teuchos::ScalarTraits<ScalarType>::random());
	// scale #2
	if(debug) cout << "scale #2..." << endl;
	vector.scale(Teuchos::ScalarTraits<ScalarType>::random(), v2);
	// update
	if(debug) cout << "update..." << endl;
	vector.update(Teuchos::ScalarTraits<ScalarType>::random(), v2, Teuchos::ScalarTraits<ScalarType>::random());
	// update #2
	if(debug) cout << "update #2..." << endl;
	vector.update(Teuchos::ScalarTraits<ScalarType>::random(), v2, 
				  Teuchos::ScalarTraits<ScalarType>::random(), vector1a, Teuchos::ScalarTraits<ScalarType>::random());
	// 1-norm
	if(debug) cout << "1-norm..." << endl;
	vector.norm1(); // throw away return value
	// 2-norm
	if(debug) cout << "2-norm..." << endl;
	vector.norm2(); // throw away return value
	// Infinity-norm
	if(debug) cout << "Infinity-norm..." << endl;
	vector.normInf(); // throw away return value
	// Weighted 2-norm
	if(debug) cout << "Weighted 2-norm (RMS norm)..." << endl;
	vector.normWeighted(v2); // throw away return value
	// min value
	if(debug) cout << "minValue..." << endl;
	vector.minValue();
	// max value
	if(debug) cout << "maxValue..." << endl;
	vector.maxValue();
	// mean value
	if(debug) cout << "meanValue..." << endl;
	vector.meanValue();
	// elementwiseMultiply
	if(debug) cout << "elementwiseMultiply..." << endl;
	ScalarType const temp1 = Teuchos::ScalarTraits<ScalarType>::random();
	vector.elementwiseMultiply(temp1, vector1a, v2, scalar3);
	// elementwiseReciprocalMultiply
	if(debug) cout << "elementwiseReciprocalMultiply..." << endl;
	vector.elementwiseMultiply(temp1, vector1a, v2, scalar3);
	
	// assignment operator
	if(debug) cout << "assignment operator..." << endl;
	v2 = vector;
	
	if(verbose) cout << "Code coverage section finished." << endl;
	
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================
	
	if(verbose) cout << "Starting actual testing section..." << endl;
	
	// default ctr initializing to zero
	if(verbose) cout << "Checking to see that default constructor initializes to zeros... ";
	Tpetra::Vector<OrdinalType, ScalarType> testVector1(vectorspace);
	OrdinalType length = testVector1.getNumMyEntries();
  for(OrdinalType i = 0; i < length; i++) {
    if (testVector1[i] != 0) {
      if(debug) {
        cout << "element " << i << " = " << testVector1[i] << ", should be zero" << endl;
      }
      ierr++;
    }
  }
	if(verbose) {
		if(ierr == 0) {
			cout << "Passed" << endl;
    }
		else {
			cout << "Failed" << endl;
    }
  }
	returnierr += ierr;
	ierr = 0;
	
	// user array ctr initializing correctly
	if(verbose) cout << "Checking to see that user array constructor initializes values correctly... ";
	length = vector1a.getNumMyEntries();
	for(OrdinalType i = 0; i < length; i++)
		if(vector1a[i] != scalarArray[i])
			ierr++;
	if(verbose) {
		if(ierr == 0) {
			cout << "Passed" << endl;
    }
		else {
			cout << "Failed" << endl;
    }
  }
	returnierr += ierr;
	ierr = 0;
	
	// changing data
	if(verbose) cout << "Changing data... ";
	OrdinalType nME = testVector1.getNumMyEntries();
	
	if(nME > 2) {
		ScalarType value2 = Teuchos::ScalarTraits<ScalarType>::random();
		testVector1[2] = value2;
		if(testVector1[2] != value2) {
			if(debug) cout << "element 2 = " << testVector1[2] << ", should be " << value2 << endl;
			ierr++;
		}
	}

	if(nME > 5) {
		ScalarType value5 = Teuchos::ScalarTraits<ScalarType>::random();
		testVector1[5] = value5;
		if(testVector1[5] != value5) {
			if(debug) cout << "element 5 = " << testVector1[5] << ", should be " << value5 << endl;
			ierr++;
		}
	}

	if(nME > 0) {
		ScalarType value0 = Teuchos::ScalarTraits<ScalarType>::random();
		testVector1[0] = value0;
		if(testVector1[0] != value0) {
			if(debug) cout << "element 0 = " << testVector1[0] << ", should be " << value0 << endl;
			ierr++;
		}
	}
	
  if (verbose) {
    if (ierr == 0) {
      cout << "Passed" << endl;
    }
    else {
      cout << "Failed" << endl;
    }
  }
	returnierr += ierr;
	ierr = 0;
	
	// check updates
	Tpetra::Vector<OrdinalType, ScalarType> u1(vectorspace);
	Tpetra::Vector<OrdinalType, ScalarType> u2(vectorspace);
	Tpetra::Vector<OrdinalType, ScalarType> u3(vectorspace);
	int const uLength = 10;

	std::vector<ScalarType> u1Vals(uLength);
	u1Vals[0] = -1;
	u1Vals[1] = 2;
	u1Vals[2] = -1;
	u1Vals[3] = -1;
	u1Vals[4] = 2;
	u1Vals[5] = -1;
	u1Vals[6] = -1;
	u1Vals[7] = 2;
	u1Vals[8] = -1;
	u1Vals[9] = 2;
	std::vector<ScalarType> u2Vals(uLength);
	u2Vals[0] = 2;
	u2Vals[1] = 3;
	u2Vals[2] = 4;
	u2Vals[3] = 5;
	u2Vals[4] = 4;
	u2Vals[5] = 3;
	u2Vals[6] = 2;
	u2Vals[7] = 1;
	u2Vals[8] = 8;
	u2Vals[9] = 9;

	// u1 and u2 have the same distribution, so we can combine their calls
	// into the same loopset
	for(int i = 0; i < uLength; i++) {
		if(vectorspace.isMyGlobalIndex(i)) {
			u1[vectorspace.getLocalIndex(i)] = u1Vals[i];
			u2[vectorspace.getLocalIndex(i)] = u2Vals[i];
		}
	}
	
	if(debug) {
		cout << "before update:" << endl;
		cout << "u1:" << endl << u1 << endl;
		cout << "u2:" << endl << u2 << endl;
		cout << "u3:" << endl << u3 << endl << endl;
	}
	
	u3.update(1.0, u1, 2.0, u2, 0.0);
	
	if(debug) {
		cout << "after update:" << endl;
		cout << "u1:" << endl << u1 << endl;
		cout << "u2:" << endl << u2 << endl;
		cout << "u3:" << endl << u3 << endl << endl;
	}
	
	// finish up
	if (verbose) {
		if(returnierr == 0) {
			cout << "VectorTest <" 
			     << Teuchos::OrdinalTraits<OrdinalType>::name() << ", " 
			     << Teuchos::ScalarTraits<ScalarType>::name() << "> passed." << endl;
    }
		else {
			cout << "VectorTest <" 
			     << Teuchos::OrdinalTraits<OrdinalType>::name() << ", " 
			     << Teuchos::ScalarTraits<ScalarType>::name() << "> failed." << endl;
    }
  }
	return(returnierr);
}

//======================================================================
void checkOutputs() {
	cout << "Doing checkOutput..." << endl;
	int const length = 5;
	int const indexBase = 0;
#ifdef TPETRA_MPI
	const Tpetra::MpiPlatform<int, int> platformE(MPI_COMM_WORLD);
	const Tpetra::MpiPlatform<int, double> platformV(MPI_COMM_WORLD);
#else
	const Tpetra::SerialPlatform <int, int> platformE;
	const Tpetra::SerialPlatform <int, double> platformV;
#endif // TPETRA_MPI
	Tpetra::ElementSpace<int> elementspace(length, indexBase, platformE);
	Tpetra::VectorSpace<int, double> vectorspace(elementspace, platformV);
	Tpetra::Vector<int, double> v1(vectorspace);
	cout << "Created v1, default constructor" << endl;
	v1.printValues(cout);
	cout << "Setting all to 3.0" << endl;
	v1.setAllToScalar(3.0);
	v1.printValues(cout);
	cout << "Setting all to random" << endl;
	v1.setAllToRandom();
	v1.printValues(cout);
	cout << "Min value: " << v1.minValue() << endl;
	cout << "Max value: " << v1.maxValue() << endl;
	cout << "Mean value: " << v1.meanValue() << endl;
	cout << "1-norm: " << v1.norm1() << endl;
	cout << "2-norm: " << v1.norm2() << endl;
	cout << "Inf-norm: " << v1.normInf() << endl;
	cout << "Creating v2, setting all to random" << endl;
	Tpetra::Vector<int, double> v2(vectorspace);
	v2.setAllToRandom();
	v2.printValues(cout);
	cout << "dot product of v1.v2 : " << v1.dotProduct(v2) << endl;
	cout << "Setting v1 by 2.0 (scale):" << endl;
	cout << "v1(before): "; v1.printValues(cout);
	v1.scale(2.0);
	cout << "v1(after): "; v1.printValues(cout);

	// extracting values should only be done if we own some
	int const nME = v1.getNumMyEntries();
	if(nME > 0) {
		cout << "copying v1 to user array" << endl;
		double* dblArray = new double[nME];
		v1.extractCopy(dblArray);
		cout << "values of dblArray: ";
		for(int i = 0; i < nME; i++)
			cout << dblArray[i] << " ";
		cout << endl;
		delete[] dblArray;
	}

	cout << "Setting v1 to reciprocal of v2" << endl;
	// only set value of an entry on the Image that owns it
	if(vectorspace.isMyGlobalIndex(0))
		v2[vectorspace.getLocalIndex(0)] = 0.25;
	if(vectorspace.isMyGlobalIndex(1))
		v2[vectorspace.getLocalIndex(1)] = 5.0;
	if(vectorspace.isMyGlobalIndex(2))
		v2[vectorspace.getLocalIndex(2)] = 0.125;
	if(vectorspace.isMyGlobalIndex(3))
		v2[vectorspace.getLocalIndex(3)] = 2.0;
	if(vectorspace.isMyGlobalIndex(4))
		v2[vectorspace.getLocalIndex(4)] = 0.75;

	cout << "v2: "; v2.printValues(cout);
	v1.reciprocal(v2);
	cout << "v1: "; v1.printValues(cout);
	
	cout << "Elementwise multiply:" << endl;
	// only set value of an entry on the Image that owns it
	if(vectorspace.isMyGlobalIndex(0))
		v1[vectorspace.getLocalIndex(0)] = 6.25; 
	if(vectorspace.isMyGlobalIndex(1))
		v1[vectorspace.getLocalIndex(1)] = 18.0; 
	if(vectorspace.isMyGlobalIndex(2))
		v1[vectorspace.getLocalIndex(2)] = 0.0; 
	if(vectorspace.isMyGlobalIndex(3))
		v1[vectorspace.getLocalIndex(3)] = 3.0; 
	if(vectorspace.isMyGlobalIndex(4))
		v1[vectorspace.getLocalIndex(4)] = 1.0;

	cout << "v3 = v1 @ v2" << endl;
	Tpetra::Vector<int, double> v3(vectorspace);
	v3.elementwiseMultiply(1.0, v1, v2, 1.0);
	cout << "v1: "; v1.printValues(cout);
	cout << "v2: "; v2.printValues(cout);
	cout << "v3: "; v3.printValues(cout);
	cout << "v3 = v3 @ v1" << endl;
	v3.elementwiseMultiply(1.0, v1, v3, 0.0);
	cout << "v3: "; v3.printValues(cout);
	
}
