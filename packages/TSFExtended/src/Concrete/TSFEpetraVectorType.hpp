/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/

#ifndef TSFEPETRAVECTORTYPE_HPP
#define TSFEPETRAVECTORTYPE_HPP

#include "TSFCoreEpetraVectorSpaceFactory.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "TSFVectorTypeExtensions.hpp"
#include "TSFLinearOperator.hpp"


namespace TSFExtended
{
  using namespace Teuchos;
  /**
   * TSF extension of TSFCore::EpetraVectorSpaceFactory, allowing 
   * use in handles and more extensive capability for creating
   * distributed spaces.
   * This class derives
   * from TSFCore::EpetraVectorSpaceFactory, so it can be used 
   * seamlessly in any 
   * TSFCore-based code.
   */
  class EpetraVectorType : public VectorTypeExtensions<double>,
                           public TSFCore::EpetraVectorSpaceFactory,
                           public Handleable<VectorTypeExtensions<double> >,
                           public Printable,
                           public Describable
  {
  public:
    /** Ctor needs no arguments */
    EpetraVectorType();
      
    /** virtual dtor */
    virtual ~EpetraVectorType() {;}

    /** create a distributed vector space.
     * @param dimension the dimension of the space 
     * @param nLocal number of indices owned by the local processor
     * @param locallyOwnedIndices array of indices owned by this processor  
     */
    RefCountPtr<const TSFCore::VectorSpace<double> > 
    createSpace(int dimension, 
                int nLocal,
                const int* locallyOwnedIndices) const ;

    /**  
     * Create an importer for accessing ghost elements.
     * @param space the distributed vector space on which ghost elements
     * are to be shared
     * @param nGhost number of ghost elements needed by this processor
     * @param ghostIndices read-only C array of off-processor indices needed
     * by this processor.
     * @return A RCP to a GhostImporter object.
     */
    RefCountPtr<GhostImporter<double> > 
    createGhostImporter(const VectorSpace<double>& space,
                        int nGhost,
                        const int* ghostIndices) const ;


    /**
     * Create an empty matrix of type compatible with this vector type,
     * sized according to the given domain and range spaces.
     */
    LinearOperator<double>
    createMatrix(const VectorSpace<double>& domain,
                 const VectorSpace<double>& range) const ;
      
    

    /** \name Describable interface */
    //@{
    /** Return a short description  */
    string describe() const {return "EpetraVectorType";}
    //@}

    /** \name Printable interface */
    //@{
    /** Print to stream */
    void print(ostream& os) const {os << describe();}
    //@}

    GET_RCP(VectorTypeExtensions<double>);

  };
  
}

#endif
