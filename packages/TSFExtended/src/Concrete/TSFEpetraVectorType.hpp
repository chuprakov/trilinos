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

    /** \name Handleable interface */
    //@{
    /** Return a ref count pointer to a newly created object */
    virtual RefCountPtr<VectorTypeExtensions<double> > getRcp() 
    {return rcp(this);}
    //@}
  };
  
}

#endif
