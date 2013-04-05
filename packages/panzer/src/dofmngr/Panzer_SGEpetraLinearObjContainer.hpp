// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_SGEpetraLinearObjContainer_hpp__
#define __Panzer_SGEpetraLinearObjContainer_hpp__

#include "Panzer_config.hpp"
#ifdef HAVE_STOKHOS

#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Teuchos_RCP.hpp"

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"

#include <vector>

namespace panzer {

/** Linear object container for SG-Epetra objects.
  */
class SGEpetraLinearObjContainer : public LinearObjContainer {
public:
   typedef std::vector<Teuchos::RCP<EpetraLinearObjContainer> > CoeffVector;
   typedef CoeffVector::iterator iterator;
   typedef CoeffVector::const_iterator const_iterator;

   SGEpetraLinearObjContainer(const CoeffVector & coeffs,
                              const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & basis);

   virtual void initialize();

   virtual std::size_t size() const
   { return coeffs_.size(); }

   CoeffVector::iterator begin() { return coeffs_.begin(); }
   CoeffVector::iterator end() { return coeffs_.end(); }

   CoeffVector::const_iterator begin() const { return coeffs_.begin(); }
   CoeffVector::const_iterator end() const { return coeffs_.end(); }

   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > getExpansion() const
   { return expansion_; }

   /** Used to get around a minor flaw in the Stokhos::EpetraVectorOrthogPoly where
     * it contains a view of the data it uses and distributes through an RCP mechansism.
     * When it is deallocated then the eventual content of those RCPs (while still
     * technically valid pointers to Epetra_Vector objects that are views of another 
     * Epetra_Vector that have been deallocated) is no longer valid. Yuck!
     */
   void storeForSafeKeeping(const Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> & v)
   { sg_vecs_.push_back(v); }

private:
   CoeffVector coeffs_;
   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion_;
   std::vector<Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> > sg_vecs_;
};

}

#endif
#endif
