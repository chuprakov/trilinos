// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_InputAdapter.hpp
    \brief Defines the InputAdapter interface.
*/

#ifndef _ZOLTAN2_INPUTADAPTER_HPP_
#define _ZOLTAN2_INPUTADAPTER_HPP_

#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! \brief An enum to identify general types of input adapters.
 *
 *  If you change this, update inputAdapterTypeName().
 */
enum InputAdapterType {
  InvalidAdapterType = 0,    /*!< \brief unused value */
  IdentifierAdapterType,    /*!< \brief plain identifier input, just a list of Ids*/
  VectorAdapterType,    /*!< \brief vector input*/
  CoordinateAdapterType,    /*!< \brief coordinate input */
  GraphAdapterType,    /*!< \brief graph input */
  MeshAdapterType,    /*!< \brief mesh input */
  MatrixAdapterType    /*!< \brief matrix input */
};


/*! \brief InputAdapter defines methods required by all InputAdapters

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    \todo Add add a MeshInput adapter
 */

class InputAdapter {
private:
public:

  /*! \brief Returns the type of adapter.
   */
  virtual enum InputAdapterType inputAdapterType()const = 0;

  /*! \brief Desstructor
   */
  virtual ~InputAdapter() {};

  /*! \brief Returns a descriptive name that identifies the concrete adapter.
   */
  virtual string inputAdapterName() const = 0;

  /*! \brief Returns the number of objects in the input.
   *
   *  Objects may be coordinates, graph vertices, matrix rows, etc.
   *  They are the objects to be partitioned, ordered, or colored.
   */
  virtual size_t getLocalNumberOfObjects() const = 0;

  /*! \brief Returns the number of weights per object.
   *   Number of weights per object should be zero or greater.  If
   *   zero, then it is assumed that all objects are equally weighted.
   */ 
  virtual int getNumberOfWeightsPerObject() const = 0;

  /*! \brief Returns the name of the input adapter
   */
  static string inputAdapterTypeName(InputAdapterType iaType);
};

string InputAdapter::inputAdapterTypeName(InputAdapterType iaType)
{
  string typeName;
  switch (iaType){
    case InvalidAdapterType:
      typeName = string("invalid");
      break;
    case IdentifierAdapterType:
      typeName = string("identifier");
      break;
    case VectorAdapterType:
      typeName = string("vector");
      break;
    case CoordinateAdapterType:
      typeName = string("coordinate");
      break;
    case GraphAdapterType:
      typeName = string("graph");
      break;
    case MeshAdapterType:
      typeName = string("mesh");
      break;
    case MatrixAdapterType:
      typeName = string("matrix");
      break;
    default:
      typeName = string("unknown");
      break;
  }

  return typeName;
}
  
  
}  //namespace Zoltan2
  
#endif
