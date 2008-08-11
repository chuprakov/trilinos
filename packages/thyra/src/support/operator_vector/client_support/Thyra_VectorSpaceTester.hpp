// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_VECTOR_SPACE_TESTER_HPP
#define THYRA_VECTOR_SPACE_TESTER_HPP

#include "Thyra_VectorSpaceTesterDecl.hpp"
#include "Thyra_TestingTools.hpp"

namespace Thyra {

template <class Scalar>
VectorSpaceTester<Scalar>::VectorSpaceTester(
  const ScalarMag     warning_tol
  ,const ScalarMag    error_tol
  ,const int          num_random_vectors
  ,const int          num_mv_cols
  ,const bool         show_all_tests
  ,const bool         dump_all
  )
  :num_mv_cols_(num_mv_cols) // declared first!
  ,warning_tol_(warning_tol)
  ,error_tol_(error_tol)
  ,num_random_vectors_(num_random_vectors)
  ,show_all_tests_(show_all_tests)
  ,dump_all_(dump_all)
{
  vectorTester_.warning_tol(warning_tol);
  vectorTester_.error_tol(error_tol);
  vectorTester_.num_random_vectors(num_random_vectors);
  vectorTester_.show_all_tests(show_all_tests);
  vectorTester_.dump_all(dump_all);
}

template <class Scalar>
bool VectorSpaceTester<Scalar>::check(
  const VectorSpaceBase<Scalar>  &vs
  ,Teuchos::FancyOStream         *out_arg
  ) const
{
  using std::endl;
  using Teuchos::describe;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;

  Teuchos::RCP<FancyOStream> out = Teuchos::rcp(out_arg,false);
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  OSTab tab(out,1,"THYRA");

  bool result, success = true;

  if(out.get()) *out <<endl<< "*** Entering Thyra::VectorSpaceTester<"<<ST::name()<<">::check(vs,...) ...\n";

  if(out.get()) *out <<endl<< "Testing a vector space vs described as:\n" << describe(vs,verbLevel);

  if(out.get())
    *out
      <<endl<< "A) Calling basic query functions ...\n"
      <<endl<< "vs.dim() = " << vs.dim()
      <<endl<< "vs.hasInCoreView() = " << vs.hasInCoreView() << std::endl;

  if(out.get()) *out <<endl<< "B) Checking that vs is compatible with itself ...\n";

  if(out.get()) *out <<endl<< "vs.isCompatible(vs)=";
  result = vs.isCompatible(vs);
  if(!result) success = false;
  if(out.get()) *out << result << " == true : " << passfail(result) << std::endl;

  if(out.get()) *out <<endl<< "C) Creating a randomized vector member v ...\n";
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    v = createMember(vs);
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),&*v);

  if(out.get()) *out <<endl<< "D) Testing the VectorBase interface of v ...\n";

  result = vectorTester_.check(*v,out.get());
  if(!result) success = false;

  if(out.get()) *out <<endl<< "C) Creating a randomized MultiVector member mv ...\n";
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
    mv = createMembers(vs,num_mv_cols());
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),&*mv);

  if(out.get()) *out <<endl<< "D) Testing the MultiVectorBase interface of mv ...\n";

  result = vectorTester_.multiVectorTester().check(*mv,out.get());
  if(!result) success = false;

  if(out.get()) {
    if(success)
      *out << endl <<"Congratulations, this VectorSpaceBase object seems to check out!\n";
    else
      *out << endl <<"Oh no, at least one of the tests performed with this VectorSpaceBase object failed (see above failures)!\n";
    *out << endl << "*** Leaving VectorSpaceTester<"<<ST::name()<<">::check(vs,...)\n";
  }
  
  return success;

}

} // namespace Thyra

#endif // THYRA_VECTOR_SPACE_TESTER_HPP
