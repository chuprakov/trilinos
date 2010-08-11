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


#ifndef TEUCHOS_VALIDATORXMLCONVERTERDB_HPP
#define TEUCHOS_VALIDATORXMLCONVERTERDB_HPP

/*! \file Teuchos_ValidatorXMLConverterDB.hpp
 * \brief A database for ValidatorXMLConverters.
*/

#include "Teuchos_ValidatorXMLConverter.hpp"


namespace Teuchos {

class ParameterEntryValidator;

/** \brief Provides ability to lookup ValidatorXMLConverterDB
 */
class ValidatorXMLConverterDB {
public:

  /** \name Modifier Functions */
  //@{
  
  /** \brief Add a converter to the database.
   *
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(ParameterEntryValidator& validator,
    RCP<ValidatorXMLConverter> converterToAdd);
  
  //@}

  /** \name Converter Functions */
  //@{
  
  /** \brief Get an appropriate ValidatorXMLConverter given a 
   * Validator.
   *
   * \param validator The ParameterEntryValidator for which a converter is
   * desired.
   */
  static RCP<const ValidatorXMLConverter> getConverter(
    const ParameterEntryValidator& validator);

  /** \brief Get an appropriate ValidatorXMLConverter given a XMLObject.
   *
   * @param xmlObject The XMLObject for which a converter is desired.
   */
  static RCP<const ValidatorXMLConverter> 
    getConverter(const XMLObject& xmlObject);

  /**
   * \brief Given a validator converts the
   * validator to XML.
   *
   * \return XML representation of the validator.
   */
  static XMLObject convertValidator(
    RCP<const ParameterEntryValidator> validator); 

  /**
   * \brief Given an XMLObject converts the XMLObject 
   * to a ParameterEntryValidator and inserts the validator into the map.
   *
   * \return A ParameterEntryValidator that was represented by the XML.
   */
  static RCP<ParameterEntryValidator> 
    convertXML(const XMLObject& xmlObject);
  
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
  
  /** \brief convience class. */
  typedef std::map<std::string, RCP<ValidatorXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<ValidatorXMLConverter> > ConverterPair;

  /** \brief Gets the default converter to be used to convert
   * Validators.
   */
  static ConverterMap& getConverterMap();
  
  //@}

};


} // end namespace Teuchos


#endif // TEUCHOS_VALIDATORXMLCONVERTERDB_HPP
