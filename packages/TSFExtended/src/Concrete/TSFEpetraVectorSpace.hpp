#ifndef TSFEPETRAVECTORSPACE_HPP
#define TSFEPETRAVECTORSPACE_HPP

#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFHandleable.hpp"
#include "TSFDescribable.hpp"

namespace TSFExtended
{
  using namespace Teuchos;
  /**
   * TSF extension of TSFCore::EpetraVectorSpace, allowing use in handles.
   * This class derives
   * from TSFCore::EpetraVectorSpace, so it can be used seamlessly in any 
   * TSFCore-based code.
   */
  class EpetraVectorSpace : public TSFCore::EpetraVectorSpace,
                            public Handleable<TSFCore::VectorSpace<double> >,
                            public Describable
    {
    public:
      /** */
      EpetraVectorSpace();
      /** */
      EpetraVectorSpace(const RefCountPtr<const Epetra_Map>& map);
      
      /** virtual dtor */
      virtual ~EpetraVectorSpace() {;}

      /** \name Describable interface */
      //@{
      /** Return a short description  */
      string describe() const ;
      //@}

      /** \name Handleable interface */
      virtual RefCountPtr<TSFCore::VectorSpace<double> > getRcp() 
      {return rcp(this);}
    };
  
}

#endif
