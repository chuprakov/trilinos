// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#ifndef TPETRA_SERIALPLATFORM_HPP
#define TPETRA_SERIALPLATFORM_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Object.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include "Tpetra_Platform.hpp"

namespace Tpetra {

	//! Tpetra::SerialPlatform: Serial Implementation of the Platform class.
	template<typename OrdinalType>
	class SerialPlatform : public Teuchos::Object, public virtual Platform<OrdinalType> {
	public:

		//@{ \name Constructor/Destructor Methods

		//! Constructor
		SerialPlatform();

		//! Copy constructor
		SerialPlatform(SerialPlatform<OrdinalType> const& platform);

		//! Destructor
		~SerialPlatform();

		//! Clone constructor
		Teuchos::RCP< Platform<OrdinalType> > clone() const;

		//@}

		//@{ \name Class Creation and Accessor Methods

		//! Comm Instance
		Teuchos::RCP< Teuchos::Comm<OrdinalType> > createComm() const;

		//@}

		//@{ \name I/O Methods

		//! print - implements Teuchos::Object virtual print method.
		void print(ostream& os) const;

		//! printInfo - implements Tpetra::Platform virtual printInfo method.
		void printInfo(ostream& os) const;

		//@}

	}; // SerialPlatform class

  template <typename OrdinalType>
  SerialPlatform<OrdinalType>::SerialPlatform() 
    : Teuchos::Object("Tpetra::SerialPlatform") 
  {}

  template <typename OrdinalType>
  SerialPlatform<OrdinalType>::SerialPlatform(SerialPlatform<OrdinalType> const& platform) 
    : Teuchos::Object(platform.label()) 
  {}

  template <typename OrdinalType>
  SerialPlatform<OrdinalType>::~SerialPlatform() 
  {}

  template <typename OrdinalType>
  Teuchos::RCP< Platform<OrdinalType> > 
  SerialPlatform<OrdinalType>::clone() const 
  {
    Teuchos::RCP< Platform<OrdinalType> > platform;
    platform = Teuchos::rcp(new SerialPlatform<OrdinalType>(*this));
    return platform;
  }

  template <typename OrdinalType>
  Teuchos::RCP< Teuchos::Comm<OrdinalType> > 
  SerialPlatform<OrdinalType>::createComm() const 
  {
    return Teuchos::rcp(new Teuchos::SerialComm<OrdinalType>() );
  }

  template <typename OrdinalType>
  void SerialPlatform<OrdinalType>::print(ostream& os) const 
  {}

  template <typename OrdinalType>
  void SerialPlatform<OrdinalType>::printInfo(ostream& os) const 
  {
    os << *this;
  }


} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
