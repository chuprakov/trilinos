/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;

SymbolicTransformation::SymbolicTransformation()
{}

RCP<ScalarExpr> SymbolicTransformation::chooseSign(int sign, 
                                                           const RCP<ScalarExpr>& expr) 
{
  /* return expr if sign == 1, -expr if sign == -1. No other
   * cases should happen. */
  switch(sign)
    {
    case 1:
      return expr;
    case -1:
      {
        Expr e = -Expr::handle(expr);
        RCP<ScalarExpr> rtn = rcp_dynamic_cast<ScalarExpr>(e.ptr());
        TEUCHOS_TEST_FOR_EXCEPTION(rtn.get() == NULL, std::logic_error,
                           "Non-scalar expr "
                           << e.toString() 
                           << " detected in SymbolicTransformation::chooseSign");
        return rtn;
      }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
                         "sign != +/- 1 in Expr::transformSign()");
    }
  return expr;
}

Expr SymbolicTransformation::chooseSign(int sign, 
                                        const Expr& expr) 
{
  /* return expr if sign == 1, -expr if sign == -1. No other
   * cases should happen. */
  switch(sign)
    {
    case 1:
      return expr;
    case -1:
      return -expr;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
                         "sign != +/- 1 in Expr::transformSign()");
    }
  return expr;
}

RCP<ScalarExpr> SymbolicTransformation::getScalar(const Expr& expr)
{
  RCP<ScalarExpr> s = rcp_dynamic_cast<ScalarExpr>(expr.ptr());

  TEUCHOS_TEST_FOR_EXCEPTION(s.get()==NULL, std::logic_error,
                     "non-scalar detected in SymbolicTransformation::getScalar");

  return s;
}
