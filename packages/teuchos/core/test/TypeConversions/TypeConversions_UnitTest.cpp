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

#include "Teuchos_as.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <limits>

// Putting the unit tests in an anonymous namespace avoids name collisions.
namespace {

// Test for overflow when converting an std::string (containing a
// positive integer too large to fit in int) to int.
TEUCHOS_UNIT_TEST( asSafe, stringToIntPositiveOverflow ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  const int maxInt = std::numeric_limits<int>::max ();

  // Write a long int to the given string that is guaranteed to
  // overflow int, if long is strictly bigger than int.
  if (sizeof (int) < sizeof (long)) {
    const long maxIntPlusOne = static_cast<long> (maxInt) + static_cast<long> (1);
    os << maxIntPlusOne;

    // Now attempt to convert the string to int.  The conversion
    // should fail, but leave the string unaffected.
    int intVal = 0;
    TEST_THROW(intVal = asSafe<int> (os.str ()), std::range_error);
    (void) intVal; // Silence compiler warning.

    // Since the string is unaffected, conversion to long should work
    // just fine, and return the correct result.
    long longVal = 0;
    TEST_NOTHROW(longVal = asSafe<long> (os.str ()));
    TEST_EQUALITY_CONST(longVal, maxIntPlusOne);
  }
  else { // C++ standard requires then that sizeof(int) == sizeof(long).
    os << maxInt;
    // Converting the string to int should not throw and should return
    // the correct result.
    int intVal = 0;
    TEST_NOTHROW(intVal = asSafe<int> (os.str ()));
    TEST_EQUALITY_CONST(intVal, maxInt);
  }
}

// Test for overflow when converting an std::string (containing a
// negative integer too negative to fit in int) to int.
TEUCHOS_UNIT_TEST( asSafe, stringToIntNegativeOverflow ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  const int minInt = std::numeric_limits<int>::min ();

  // Write a long int to the given string that is guaranteed to
  // overflow int, if long is strictly bigger than int.
  if (sizeof (int) < sizeof (long)) {
    const long minIntMinusOne = static_cast<long> (minInt) - static_cast<long> (1);
    os << minIntMinusOne;

    // Now attempt to convert the string to int.  The conversion
    // should fail, but leave the string unaffected.
    int intVal = 0;
    TEST_THROW(intVal = asSafe<int> (os.str ()), std::range_error);
    (void) intVal; // Silence compiler warning

    // Since the string is unaffected, conversion to long should work
    // just fine, and return the correct result.
    long longVal = 0;
    TEST_NOTHROW(longVal = asSafe<long> (os.str ()));
    TEST_EQUALITY_CONST(longVal, minIntMinusOne);
  }
  else { // C++ standard requires then that sizeof(int) == sizeof(long).
    os << minInt;
    // Converting the string to int should not throw and should return
    // the correct result.
    int intVal = 0;
    TEST_NOTHROW(intVal = asSafe<int> (os.str ()));
    TEST_EQUALITY_CONST(intVal, minInt);
  }
}

// Unit test for conversion from std::string (containing a positive
// integer) to built-in integer types (may be signed or unsigned).
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, stringToIntegerPositive, IntegerType ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  os << static_cast<IntegerType> (42);
  IntegerType val = 0;
  TEST_NOTHROW(val = asSafe<IntegerType> (os.str ()));
  TEST_EQUALITY_CONST(val, static_cast<IntegerType> (42));
}

// Unit test for conversion from std::string (containing a negative
// integer) to built-in integer types (must be signed).
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, stringToIntegerNegative, IntegerType ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  os << static_cast<IntegerType> (-42);
  IntegerType val = 0;
  TEST_NOTHROW(val = asSafe<IntegerType> (os.str ()));
  TEST_EQUALITY_CONST(val, static_cast<IntegerType> (-42));
}

// Unit test for conversion from std::string (NOT containing an
// integer) to built-in integer types (may be signed or unsigned).
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( asSafe, stringToIntegerShouldThrow, IntegerType ) {
  using Teuchos::asSafe;

  std::ostringstream os;
  os << "This string definitely does not contain an integer.";
  IntegerType val = 0;
  TEST_THROW(val = asSafe<IntegerType> (os.str ()), std::invalid_argument);
  (void) val; // Silence compiler warning
}

// Macros to instantiate templated unit tests for conversion from
// std::string to built-in integer types.  AnyIntegerType may be
// signed or unsigned; SignedIntegerType must be signed.

#define UNIT_TEST_GROUP_ANY_INTEGER( AnyIntegerType ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, stringToIntegerPositive, AnyIntegerType ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, stringToIntegerShouldThrow, AnyIntegerType )

#define UNIT_TEST_GROUP_SIGNED_INTEGER( SignedIntegerType ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( asSafe, stringToIntegerNegative, SignedIntegerType )

// Instantiations of templated unit tests for conversion from
// std::string to built-in integer types.

UNIT_TEST_GROUP_ANY_INTEGER( int )
UNIT_TEST_GROUP_ANY_INTEGER( long )

//
// Hack to work around Bug 5757
//
typedef unsigned int unsigned_int_type;
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned int )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_int_type )

//
// Hack to work around Bug 5757
//
typedef unsigned long unsigned_long_type;
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned long )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_long_type )


#ifdef HAVE_TEUCHOS_LONG_LONG_INT
//
// Hack to work around Bug 5757
//
typedef long long long_long_type;
//UNIT_TEST_GROUP_ANY_INTEGER( long long )
UNIT_TEST_GROUP_ANY_INTEGER( long_long_type )
//UNIT_TEST_GROUP_SIGNED_INTEGER( long long )
UNIT_TEST_GROUP_SIGNED_INTEGER( long_long_type )

//
// Hack to work around Bug 5757
//
typedef unsigned long long unsigned_long_long_type;
//UNIT_TEST_GROUP_ANY_INTEGER( unsigned long long )
UNIT_TEST_GROUP_ANY_INTEGER( unsigned_long_long_type )

#endif // HAVE_TEUCHOS_LONG_LONG_INT

} // namespace (anonymous)



