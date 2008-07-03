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

#ifndef TPETRA_MAP_DECL_HPP
#define TPETRA_MAP_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Object.hpp>
#include "Tpetra_Platform.hpp"

namespace Tpetra {

  template<typename OrdinalType> class MapData;

  //! Tpetra::Map

  template<typename OrdinalType>
  class Map : public Teuchos::Object {

  public:

    //@{ \name Constructor/Destructor Methods

    /*! \brief Map constructor with Tpetra-defined contiguous uniform distribution.
     *   The entries are distributed among nodes so that the subsets of global entries
     *   are non-overlapping and contiguous and as evenly distributed across the nodes as 
     *   possible.
     */
    Map (OrdinalType numGlobalEntries, OrdinalType indexBase, const Platform &platform);

    /*! \brief Map constructor with a user-defined contiguous distribution.
     *  The entries are distributed among the nodes so that the subsets of global entries
     *  are non-overlapping and contiguous 
     *  
     *  If numGlobalEntries == -1, it will be computed via a global communication.
     *  Otherwise, it must be equal to the sum of the local entries across all 
     *  nodes. This will only be verified if Trilinos was compiled with --enable-teuchos-debug.
     *  If this verification fails, a std::invalid_argument exception will be thrown.
     */
    Map (OrdinalType numGlobalEntries, OrdinalType numMyEntries, OrdinalType indexBase, 
         OrdinalType myOffset, const Platform &platform);

    //! Map constructor with user-defined non-contiguous (arbitrary) distribution.
    Map (OrdinalType numGlobalEntries, OrdinalType numMyEntries, 
         const std::vector< OrdinalType > &entryList, OrdinalType indexBase, 
         const Platform &platform);

    //! Map copy constructor.
    Map (const Map<OrdinalType> &Map);

    //! Map destructor. 
    ~Map ();

    //@}


    //@{ \name Map Attribute Methods

    //! Returns the number of entries in this Map.
    OrdinalType getNumGlobalEntries() const;

    //! Returns the number of entries belonging to the calling image.
    OrdinalType getNumMyEntries() const;

    //! Returns the index base for this Map.
    OrdinalType getIndexBase() const;

    //! Returns minimum local index
    OrdinalType getMinLocalIndex() const;

    //! Returns maximum local index
    OrdinalType getMaxLocalIndex() const;

    //! Returns minimum global index owned by this image
    OrdinalType getMinGlobalIndex() const;

    //! Returns maximum global index owned by this image
    OrdinalType getMaxGlobalIndex() const;

    //! Return the local index for a given global index
    OrdinalType getLocalIndex(OrdinalType globalIndex) const;

    //! Return the minimum global index over all images
    OrdinalType getMinAllGlobalIndex() const;

    //! Return the maximum global index over all images
    OrdinalType getMaxAllGlobalIndex() const;

    //! Return the global index for a given local index
    OrdinalType getGlobalIndex(OrdinalType localIndex) const;

    //! Returns true if the local index value passed in is found on the calling image, returns false if it doesn't.
    bool isMyLocalIndex(OrdinalType localIndex) const;

    //! Returns true if the global index value passed in is found the calling image, returns false if it doesn't.
    bool isMyGlobalIndex(OrdinalType globalIndex) const;

    //! Returns true if this Map is distributed contiguously, returns false otherwise.
    bool isContiguous() const;

    //@}

    //@{ \name Boolean Tests

    //! Returns true if \c map is compatible with this Map.
    bool isCompatible (const Map< OrdinalType> &map) const;

    //! Returns true if \c map is identical to this Map.
    bool isSameAs (const Map<OrdinalType> &map) const;

    //! Returns true if \c map is identical to this Map. Implemented in isSameAs().
    bool operator== (const Map< OrdinalType > &map) const;

    //! Returns true if \c map is not identical to this Map. Implemented in isSameAs().
    bool operator!= (const Map< OrdinalType > &map) const;

    //@}

    //@{ Misc. 

    //! Assignment operator
    Map<OrdinalType>& operator = (const Map<OrdinalType> & Source);

    //@}

    //@{ Implements Teuchos::Object 

    //! Prints the Map object to the output stream.
    /*! An << operator is inherited from Teuchos::Object, which uses the print method.*/
    void print(ostream& os) const;

    //@}


  private:

    Teuchos::RCP< Map<OrdinalType> > MapData_;

    // setup the directory
    void directorySetup();

  }; // Map class

} // Tpetra namespace

#endif // TPETRA_MAP_DECL_HPP

