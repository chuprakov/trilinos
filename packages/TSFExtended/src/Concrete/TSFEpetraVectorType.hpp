#ifndef TSFEPETRAVECTORTYPE_HPP
#define TSFEPETRAVECTORTYPE_HPP

#include "TSFCoreEpetraVectorSpaceFactory.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFHandleable.hpp"


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
  class EpetraVectorType : public TSFCore::EpetraVectorSpaceFactory,
                           public Handleable<TSFCore::VectorSpaceFactory<double> >
  {
  public:
    /** Ctor needs no arguments */
    EpetraVectorType();
      
    /** virtual dtor */
    virtual ~EpetraVectorType() {;}

    /** create a vector space of dimension dim */
    virtual RefCountPtr<const TSFCore::VectorSpace<double> > 
    createVecSpc(int dimension) const ;

    /** create a vector space in which the local processor owns
     * indices \f$[firstLocal, firstLocal+nLocal]f$. */
    virtual RefCountPtr<const TSFCore::VectorSpace<double> > 
    createVecSpc(int dimension, 
                 int nLocal,
                 int firstLocal) const  ;
      
    /** create a vector space in which the given local indices are owned by 
     * this processor */
    virtual RefCountPtr<const TSFCore::VectorSpace<double> > 
    createVecSpc(int dimension, 
                 int nLocal,
                 const int* localIndices) const  ;
      
    

    /** \name Describable interface */
    //@{
    /** Return a short description  */
    string describe() const {return "EpetraVectorType";}
    //@}

    /** \name Handleable interface */
    virtual RefCountPtr<TSFCore::VectorSpaceFactory<double> > getRcp() 
    {return rcp(this);}
  };
  
}

#endif
