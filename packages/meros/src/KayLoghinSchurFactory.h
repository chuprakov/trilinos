// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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

#ifndef KAYLOGHINSCHURFACTORY_H
#define KAYLOGHINSCHURFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolver.h"
#include "TSFOperatorSourceBase.h"
#include "RightBlockNSOperatorSource.h"
#include "KayLoghinRightOperatorSource.h"
#include "SchurFactoryBase.h"
#include "SchurFactory.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace Meros
{

  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * KayLoghinSchurFactory
   *
   */

  class KayLoghinSchurFactory : public SchurFactoryBase
    {
    public:
      /** ctor */
      KayLoghinSchurFactory(TSFLinearSolver& ApSolver);

      //      /** ctor (with Mp solver) */
      //      KayLoghinSchurFactory(TSFLinearSolver& ApSolver,
      //                            TSFLinearSolver& MpSolver);

      // maybe have a default that takes no solvers?

      /** virtual destructor */
      virtual ~KayLoghinSchurFactory(){;}

      /** get Schur complement inverse approximation */
      virtual TSFLinearOperator getSchurInvApprox(const TSFOperatorSource& klfac) const;

      /** write to a string */
      virtual string toString() const ;

    private:
      TSFLinearSolver ApSolver_;
    };

  /** \relates KayLoghinSchurFactory
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const KayLoghinSchurFactory& x);

  /** \relates KayLoghinSchurFactory
   * write to a string
   */
  string toString(const KayLoghinSchurFactory& x);



}


#endif
