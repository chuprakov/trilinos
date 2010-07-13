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

#ifndef TEUCHOS_VALIDATORXMLCONVERTER_HPP
#define TEUCHOS_VALIDATORXMLCONVERTER_HPP

/*! \file Teuchos_ValidatorXMLConverter.hpp
*/

#include "Teuchos_XMLObject.hpp"

namespace Teuchos {

class ParameterEntryValidator;

/* \class Teuchos::ValidatorXMLConverter
 * \brief An abstract base class for converting ParameterEntryValidators to and from XML.
 */
class ValidatorXMLConverter{
public:
	/* \brief Converts a given XMLObject to a ParameterEntryValidator.
	 *
	 * @param xmlObj The XMLObject to convert to a ParameterEntryValidator.
	 * @return The converted ParameterEntryValidator.
	 */
	virtual RCP<ParameterEntryValidator> fromXMLtoValidator(const XMLObject& xmlObj) const=0;

	/* \brief Converters a given ParameterEntryValidator to XML.
	 *
	 * @param validator The ParameterEntryValidator to be converted to XML.
	 * @return An XML representation of the given ParameterEntryValidator.
	 */
	virtual XMLObject fromValidatortoXML(const RCP<const ParameterEntryValidator> validator) const=0;

	/* \brief Determines whether or not this is the appropriate converter given a ParameterEntryValidator.
	 *
	 * @param validator The validator to test.
	 * @return True if the converter is appropriate for the ParameterEntryValidator, false otherwise.
	 */
	virtual bool isAppropriateConverter(const RCP<const ParameterEntryValidator> validator) const=0;
};

}

#endif // TEUCHOS_VALIDATORXMLCONVERTER_HPP
