// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_CONDITIONXMLCONVERTER_HPP
#define TEUCHOS_CONDITIONXMLCONVERTER_HPP

/*! \file Teuchos_ConditionXMLConverter.hpp
 * \brief Converts back and forth between XML
 * and Dependencies.
*/

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Describable.hpp"


namespace Teuchos {


/** \brief An abstract base class for converting Dependencies to
 * and from XML.
 */
class ConditionXMLConverter : public Describable {

public:

  /** \name Converter Functions */
  //@{
  
  /** \brief Converts a given XMLObject to a Condition.
   *
   * @param xmlObj The XMLObject to convert to a Condition.
   * @return The converted Condition.
   */
  RCP<Condition>
  fromXMLtoCondition(const XMLObject& xmlObj) const;

  /** \brief Preforms any and all special xml conversion that is specific to a
   * particular Condition.
   *
   * @param xmlObj The xml to be converted.
   * in which this resulting condition will be inserted.
   * @return The converted Condition.
   */
  virtual RCP<Condition> convertXML(
    const XMLObject& xmlObj) const=0;

  /** \brief Converters a given ParameterEntryValidator to XML.
   *
   * @param condition The Condition to be converted to XML.
   * @return An XML representation of the given Condition.
   */
  XMLObject fromConditiontoXML(const RCP<const Condition> condition);
  
  /** \brief Preforms any and all special condition conversion that is
   * specific to a particlar Condition.
   *
   * @param condition The validator to be converted.
   * @return An XML representation of the given Condition.
   */
  virtual void convertCondition(
    const RCP<const Condition> condition, 
    XMLObject& xmlObj) const = 0;
  
  //@}

  //! \name Attribute/Query Functions
  //@{

  /** \brief gets the value to be used for the type attribute. */
  const std::string& getTypeAttributeValue() const = 0;

  /**
   * \brief Returns the string to be used for the type attribute.
   */
  static const std::string& getTypeAttributeName(){
    static const std::string typeAttributeName = "type";
    return typeAttributeName;
  }
 
  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_CONDITIONXMLCONVERTER_HPP
