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


#ifndef TEUCHOS_CONDITIONXMLCONVERTERDB_HPP
#define TEUCHOS_CONDITIONXMLCONVERTERDB_HPP

/*! \file Teuchos_ConditionXMLConverterDB.hpp
 * \brief A database for ConditionXMLConverters.
*/

#include "Teuchos_ConditionXMLConverter.hpp"


namespace Teuchos {

class Condition;

/** \brief Provides ability to lookup ConditionXMLConverters
 */
class ConditionXMLConverterDB {

public:

  /** \name Modifier Functions */
  //@{
  
  /** \brief Add a converter to the database.
   *
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(Condition& condition,
    RCP<ConditionXMLConverter> converterToAdd);
  
  //@}

  /** \name Converter Functions */
  //@{
  
  /** \brief Get an appropriate ConditionXMLConverter given a 
   *  Condition.
   *
   * \param condition The Condition for which a converter is
   * desired.
   */
  static RCP<const ConditionXMLConverter> 
    getConverter(const Condition& condition);

  /** \brief Get an appropriate ConditionXMLConverter given a XMLObject.
   *
   * @param xmlObject The XMLObject for which a converter is desired.
   */
  static RCP<const ConditionXMLConverter> 
    getConverter(const XMLObject& xmlObject);

  /**
   * \brief Given a condition and ConditiontoIDMap, converts the
   * condition to XML.
   *
   * \return XML representation of the condition.
   */
  static XMLObject convertCondition(
    const RCP<const Condition> condition,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap);

  /**
   * \brief Given an XMLObject and IDtoConditionMap, converts the XMLObject 
   * to a ParameterEntryCondition and inserts the condition into the map.
   *
   * \return A ParameterEntryCondition that was represented by the XML.
   */
  static RCP<Condition> convertXML(
    const XMLObject& xmlObject,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap);
  
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief prints the xml tags associated with all known converters
   *
   * \param out Stream to which tags should be printed.
   */
  static void printKnownConverters(std::ostream& out){
    out << "Known ConditionXMLConverters: " << std::endl;
    for(
      ConverterMap::const_iterator it = getConverterMap().begin();
      it != getConverterMap().end();
      ++it)
    {
      out << "\t" << it->first <<std::endl;
    }
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief convience class. */
  typedef std::map<std::string, RCP<ConditionXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<ConditionXMLConverter> > ConverterPair;

  /** \brief Gets the default converter to be used to convert
   * Conditions.
   */
  static ConverterMap& getConverterMap();
  
  //@}

};


} // end namespace Teuchos


#endif // TEUCHOS_CONDITIONXMLCONVERTERDB_HPP
