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

#ifndef KAYLOGHINRIGHTOPERATORSOURCE_H
#define KAYLOGHINRIGHTOPERATORSOURCE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
//#include "KayLoghinSchurFactory.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace Meros
{
  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * KayLoghinRightOperatorSource
   *
   */

  // Need a better name for this class.

  class KayLoghinRightOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      // KayLoghinRightOperatorSource(TSFLinearOperator& S);

      /** ctor (no Mp) */
      KayLoghinRightOperatorSource(TSFLinearOperator& S,
                                   TSFLinearOperator& Fp,
                                   TSFLinearOperator& Ap);

      // /** ctor (with Mp) */
      // KayLoghinRightOperatorSource(TSFLinearOperator& S,
      //                             TSFLinearOperator& Fp,
      //                             TSFLinearOperator& Ap,
      //                             TSFLinearOperator& Mp);

      // /** ctor, we build Fp (no Mp) */
      // KayLoghinRightOperatorSource(TSFLinearOperator& S,
      //                             TSFLinearOperator& Ap);

      // I don't see how we can have an option where we build Fp but they supply Mp
      // No way to distinguish the constructors.

      /** virtual destructor */
      virtual ~KayLoghinRightOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** get pressure Poisson operator */
      TSFLinearOperator getAp() const;

      /** get pressure convection diffusion operator */
      TSFLinearOperator getFp() const;

      /** get mass matrix */
      TSFLinearOperator getMp() const;

      /** write to a string */
      string toString() const ;

    private:
      TSFLinearOperator S_;
      TSFLinearOperator Fp_;
      TSFLinearOperator Ap_;
      TSFLinearOperator Mp_;
    };


}


#endif
