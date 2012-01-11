// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Teuchos {


/*
TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentNames ) {
  // Lists with different names should not be equal
  ParameterList A("Tom");
  ParameterList B("Bob");
  TEST_ASSERT( A != B );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( A != B );
}

TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentNames ) {
  ParameterList A("Julie");
  ParameterList B("Shannon");
  TEST_ASSERT( !haveSameValues(A,B) );
  A.set("Hello","World");
  B.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
}
*/


TEUCHOS_UNIT_TEST( createParameterList, empty )
{
  RCP<ParameterList> pl = createParameterList();
  TEST_ASSERT(nonnull(pl));
}


TEUCHOS_UNIT_TEST( createParameterList, withName )
{
  RCP<ParameterList> pl = createParameterList("dummyName");
  TEST_ASSERT(nonnull(pl));
  TEST_EQUALITY_CONST(pl->name(), "dummyName");
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityWithEmpty )
{
  // An empty list should not be equal to a full list
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( A == B );
  A.set("Hello","World");
  TEST_ASSERT( A != B );
  B.set("Hello","World");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentSublistNames )
{
  // Sublists with different names should not be equal
  ParameterList A;
  ParameterList B;
  A.sublist("Bob");
  B.sublist("Tom");
  TEST_ASSERT( A != B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, operatorEqualityDifferentLengths )
{
  ParameterList A;
  ParameterList B;
  A.set("A","a");
  A.set("B","b");
  A.set("C","c");
  A.print();

  B.set("A","a");
  B.set("B","b");
  B.print();

  TEST_ASSERT( A != B );

  B.set("C","c");
  TEST_ASSERT( A == B );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesWithEmpty )
{
  ParameterList A;
  ParameterList B;
  TEST_ASSERT( haveSameValues(A,B) );
  A.set("Hello","World");
  TEST_ASSERT( !haveSameValues(A,B) );
  B.set("Hello","World");
  TEST_ASSERT( haveSameValues(A,B) );
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, haveSameValuesDifferentSublistNames )
{
  ParameterList A;
  ParameterList B;
  A.sublist("Smith").set("People",4);
  B.sublist("Jones").set("People",4);
  TEST_ASSERT( !haveSameValues(A,B) ); // sublist names matter
}


} // namespace Teuchos



