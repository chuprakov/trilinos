#ifndef TSFEPETRAVECTOR_HPP
#define TSFEPETRAVECTOR_HPP

#include "TSFPrintable.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFAccessibleVector.hpp"


namespace TSFExtended
{
  using namespace Teuchos;
  /**
   * TSF extension of TSFCore::EpetraVector, implementing the LoadableVector
   * interface allowing an application to access elements. This class derives
   * from TSFCore::EpetraVector, so it can be used seamlessly in any 
   * TSFCore-based code.
   */
  class EpetraVector : public TSFCore::EpetraVector,
                       public VectorExtensions
    {
    public:
      /** Construct with a smart pointer to an Epetra FE vector. */
      EpetraVector(const RefCountPtr<Epetra_FEVector>& vec);

      /** virtual dtor */
      virtual ~EpetraVector() {;}

      /** \name LoadableVector interface */
      //@{
      /** set a single element */
      void setElement(int index, const double& value);

      /** get a single element */
      double getElement(int index) const ;

      /** set a group of elements */
      void setElements(int numElems, const int* globalIndices, 
                       const double* values);


      /** add to a group of elements */
      void addToElements(int numElems, const int* globalIndices, 
                         const double* values);

      /** */
      void finalizeAssembly();
      //@}
      

      /** \name Printable interface */
      //@{
      /** Write to a stream  */
      void print(ostream& os) const 
      {
        epetra_vec()->Print(os);
      }
      //@}
    };
  
}

#endif
