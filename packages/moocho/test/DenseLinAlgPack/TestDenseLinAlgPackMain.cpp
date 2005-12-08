// Run all the test software for DenseLinAlgPack

//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
#include <iostream>

#include "DenseLinAlgPack_TestDenseLinAlgPack.hpp"

int main() {

//	DenseLinAlgPack::TestingPack::TestVectorBasicOp(std::cout);
//	DenseLinAlgPack::TestingPack::TestGenMatrixBasicOp(std::cout);
//	DenseLinAlgPack::TestingPack::TestDenseLinAlgPackIO(std::cin,std::cout);
//	DenseLinAlgPack::TestingPack::TestPivotVecMat(std::cout);
	bool result = DenseLinAlgPack::TestingPack::TestDenseLinAlgPack( &std::cout );
	if(result)
		std::cerr << "DenseLinAlgPack checks out!\n";
	else
		std::cerr << "Oops!  At least one of the tests for DenseLinAlgPack failed!\n";
	return ( result == true ? 0 : -1 );
}
