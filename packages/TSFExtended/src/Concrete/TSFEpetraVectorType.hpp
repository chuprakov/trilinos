#ifndef TSFEPETRAVECTORTYPE_HPP
#define TSFEPETRAVECTORTYPE_HPP

#include "TSFCoreEpetraVectorSpace.hpp"


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
                           public VectorTypeExtensions<double>
  {
  public:
    /** Ctor needs no arguments */
    EpetraVectorType();
      
    /** virtual dtor */
    virtual ~EpetraVectorType() {;}

    /** create a vector space in which the local processor owns
     * indices \f$[firstLocal, firstLocal+nLocal]f$. */
    virtual RefCountPtr<TSFCore::VectorSpace<Scalar> > 
    createVecSpc(int dimension, 
                 int nLocal,
                 int firstLocal) const  ;
      
    /** create a vector space in which the given local indices are owned by 
     * this processor */
    virtual RefCountPtr<TSFCore::VectorSpace<Scalar> > 
    createVecSpc(int dimension, 
                 int nLocal,
                 const int* localIndices) const  ;
      
    /** create a vector space in which the given local indices are owned
     * by this processor, and the given ghost indices are available but
     * not owned. */
    virtual RefCountPtr<TSFCore::VectorSpace<Scalar> > 
    createVecSpc(int dimension, 
                 int nLocal,
                 const int* localIndices,
                 int nGhost,
                 const int* ghostIndices) const  ;

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
