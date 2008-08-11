// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "ELRDFadUnitTests.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION( ELRDFadOpsUnitTest );

ELRDFadOpsUnitTest::ELRDFadOpsUnitTest() :
  urand(0.0, 1.0), n(5), tol_a(1.0e-15), tol_r(1.0e-14) {}

ELRDFadOpsUnitTest::ELRDFadOpsUnitTest(int numComponents, 
					   double absolute_tolerance, 
					   double relative_tolerance) :
  urand(0.0, 1.0), 
  n(numComponents), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) {}

void ELRDFadOpsUnitTest::setUp() {
  double val;

  val = urand.number();
  a_dfad = DFadType(n,val);
  a_fad = FAD::Fad<double>(n,val);
  
  val = urand.number();
  b_dfad = DFadType(n,val);
  b_fad = FAD::Fad<double>(n,val);

  for (int i=0; i<n; i++) {
    val = urand.number();
    a_dfad.fastAccessDx(i) = val;
    a_fad.fastAccessDx(i) = val;

    val = urand.number();
    b_dfad.fastAccessDx(i) = val;
    b_fad.fastAccessDx(i) = val;
  }
}

void ELRDFadOpsUnitTest::tearDown() {}

void ELRDFadOpsUnitTest::compareFads(const DFadType& x_dfad,
				       const FAD::Fad<double>& x_fad) {

  if (x_dfad.size() != x_fad.size()) {
    std::cout << "x_dfad.size() = " << x_dfad.size() << std::endl
	      << "x_fad.size() = " << x_fad.size() << std::endl
	      << "x_dfad = " << x_dfad << std::endl
	      << "x_fad = " << x_fad << std::endl;
  }

  // Compare sizes
  CPPUNIT_ASSERT(x_dfad.size() == x_fad.size());
  
  // Compare hasFastAccess
  CPPUNIT_ASSERT(x_dfad.hasFastAccess() == x_fad.hasFastAccess());
  
  // Compare values
  compareDoubles(x_dfad.val(), x_fad.val());
  
  for (int i=0; i<x_fad.size(); i++) {
    
    // Compare dx
    compareDoubles(x_dfad.dx(i), x_fad.dx(i));
    
    // Compare fastAccessDx
    compareDoubles(x_dfad.fastAccessDx(i), x_fad.fastAccessDx(i));
  }
}

void ELRDFadOpsUnitTest::compareDoubles(double a, double b) {
  if ( fabs(a-b) >= tol_a + tol_r*fabs(a) )
    std::cout << "a = " << a << std::endl
	      << "b = " << b << std::endl;
  CPPUNIT_ASSERT( fabs(a-b) < tol_a + tol_r*fabs(a) );
}

void ELRDFadOpsUnitTest::testMax() {
  double val;

  DFadType aa_dfad = a_dfad + 1.0;
  c_dfad = max(aa_dfad, a_dfad);
  compareDoubles(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), aa_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }
  
  c_dfad = max(a_dfad, aa_dfad);
  compareDoubles(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), aa_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = max(a_dfad+1.0, a_dfad);
  compareDoubles(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), aa_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }
  
  c_dfad = max(a_dfad, a_dfad+1.0);
  compareDoubles(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), aa_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }
  
  val = a_dfad.val() + 1;
  c_dfad = max(a_dfad, val);
  compareDoubles(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_dfad.dx(i), 0.0);
  
  val = a_dfad.val() - 1;
  c_dfad = max(a_dfad, val);
  compareDoubles(c_dfad.val(), a_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), a_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() + 1;
  c_dfad = max(val, b_dfad);
  compareDoubles(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_dfad.dx(i), 0.0);
  
  val = b_dfad.val() - 1;
  c_dfad = max(val, b_dfad);
  compareDoubles(c_dfad.val(), b_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), b_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}

void ELRDFadOpsUnitTest::testMin() {
  double val;

  DFadType aa_dfad = a_dfad - 1.0;
  c_dfad = min(aa_dfad, a_dfad);
  compareDoubles(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), aa_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  c_dfad = min(a_dfad, aa_dfad);
  compareDoubles(c_dfad.val(), aa_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), aa_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), aa_dfad.fastAccessDx(i));
  }

  val = a_dfad.val() - 1;
  c_dfad = min(a_dfad, val);
  compareDoubles(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_dfad.dx(i), 0.0);
  
  val = a_dfad.val() + 1;
  c_dfad = min(a_dfad, val);
  compareDoubles(c_dfad.val(), a_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), a_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), a_dfad.fastAccessDx(i));
  }

  val = b_dfad.val() - 1;
  c_dfad = min(val, b_dfad);
  compareDoubles(c_dfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_dfad.dx(i), 0.0);
  
  val = b_dfad.val() + 1;
  c_dfad = min(val, b_dfad);
  compareDoubles(c_dfad.val(), b_dfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_dfad.dx(i), b_dfad.dx(i));
    compareDoubles(c_dfad.fastAccessDx(i), b_dfad.fastAccessDx(i));
  }
}
