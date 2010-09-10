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


#ifndef TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
#define TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP

#include "Teuchos_StandardParameterEntryXMLConverters.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"


/*! \file Teuchos_ParameterEntryXMLCoverterDB.hpp
 * \brief A database for ParameterEntryXMLConverters.
*/


namespace Teuchos {

/** \brief Provides ability to lookup ParameterEntryXMLConverters
 */
class ParameterEntryXMLConverterDB {
public:

  /** \name Modifier Functions */
  //@{
  
  /** \brief Add a converter to the database.
   *
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(RCP<ParameterEntryXMLConverter> converterToAdd){
    getConverterMap().insert(
      ConverterPair(converterToAdd->getTypeAttributeValue(), converterToAdd));
  }
  
  //@}

  /** \name Getter Functions */
  //@{

  /** \brief Get an appropriate ParameterEntryXMLConverter given a ParameterEntry.
   *
   * \param entry The ParameterEntry for which a converter is desired.
   */
  static RCP<const ParameterEntryXMLConverter> 
    getConverter(RCP<const ParameterEntry> entry);

  /** \brief Get an appropriate ParameterEntryXMLConverter given a XMLObject.
   *
   * \param xmlObject The XMLObject for which a converter is desired.
   */
  static RCP<const ParameterEntryXMLConverter> 
    getConverter(const XMLObject& xmlObject);

  /** \brief Gets the default converter to be used on Parameter Entries */
  static RCP<const ParameterEntryXMLConverter> getDefaultConverter();

  //@}

  // 2010/07/30: rabarlt: The above two functions should be moved into
  // Teuchos_ParameterEntryXMLConvergerDB.cpp.  These functions don't need to
  // be inlined and it will be easier to set breakpoints in the debugger if
  // they are in a *.cpp file.

  /** \name Converter Functions */
  //@{
  
  /**
   * \brief Converts the given ParameterEntry to XML.
   */
  static XMLObject convertEntry(
    RCP<const ParameterEntry> entry, 
    const std::string& name,
    const ParameterEntry::ParameterEntryID& id,
    const XMLParameterListWriter::ValidatorIDsMap& validatorIDsMap)
  {
    return getConverter(entry)->fromParameterEntrytoXML(
      entry, name, id, validatorIDsMap);
  }

  /**
   * \brief Converts XML to a ParameterEntry.
   */
  static ParameterEntry convertXML(const XMLObject& xmlObj)
  {
    return getConverter(xmlObj)->fromXMLtoParameterEntry(xmlObj);
  }
  
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief prints the xml tags associated with all known converters
   *
   * \param out Stream to which tags should be printed.
   */
  static void printKnownConverters(std::ostream& out);
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief convience typedef */
  typedef std::map<std::string, RCP<ParameterEntryXMLConverter> > ConverterMap;

  /** \brief convience typedef */
  typedef std::pair<std::string, RCP<ParameterEntryXMLConverter> > ConverterPair;

  /** \brief Gets the map containing all the converters. */
  static ConverterMap& getConverterMap();
  
  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
