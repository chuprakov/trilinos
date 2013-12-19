// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_BaseAdapter.hpp
    \brief Defines the Adapter interface for accessing user data.
*/

#ifndef _ZOLTAN2_ADAPTER_HPP_
#define _ZOLTAN2_ADAPTER_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

/*! \brief An enum to identify general types of adapters.
 *
 */
enum BaseAdapterType {
  InvalidAdapterType = 0, /*!< \brief unused value */
  IdentifierAdapterType,  /*!< \brief identifier data, just a list of IDs*/
  CoordinateAdapterType,  /*!< \brief coordinate data */
  VectorAdapterType,      /*!< \brief vector data */
  MatrixAdapterType,      /*!< \brief matrix data */
  GraphAdapterType,       /*!< \brief graph data */
  MeshAdapterType         /*!< \brief mesh data */
};


/*! \brief BaseAdapter defines methods required by all Adapters

    Adapters provide access from Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    \todo Add add a MeshAdapter
 */

template <typename User>
  class BaseAdapter {

private:

  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::scalar_t scalar_t;

public:

  /*! \brief Returns the type of adapter.
   */
  virtual enum BaseAdapterType adapterType()const = 0;

  /*! \brief Destructor
   */
  virtual ~BaseAdapter() {};

  /*! \brief Returns the number of objects on this process
   *
   *  Objects may be coordinates, graph vertices, matrix rows, etc.
   *  They are the objects to be partitioned, ordered, or colored.
   */
  virtual size_t getLocalNumIDs() const = 0;

  /*! \brief Provide a pointer to this process' identifiers.

      \param Ids will on return point to the list of the global Ids for 
        this process.
   */
  virtual void getIDsView(const gid_t *&Ids) const = 0;

  /*! \brief Returns the number of weights per object.
   *   Number of weights per object should be zero or greater.  If
   *   zero, then it is assumed that all objects are equally weighted.
   */ 
  virtual int getNumWeightsPerID() const = 0;

  /*! \brief Provide pointer to a weight array with stride.
   *    \param wgt on return a pointer to the weights for this dimension
   *    \param stride on return, the value such that
   *       the \t nth weight should be found at <tt> wgt[n*stride] </tt>.
   *    \param idx  the weight index, zero or greater
   */ 
  virtual void getWeightsView(const scalar_t *&wgt, int &stride,
                              int idx = 0) const = 0;

 /*! \brief Apply a PartitioningSolution to an input.
   *
   *  This is not a required part of the InputAdapter interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *             of global Ids to processes that matches that of the input
   *             data that instantiated this Adapter.
   *  \param out On return this should point to a newly created object
   *             with the specified partitioning.
   *  \param solution  The Solution object created by a Problem should
   *      be supplied as the third argument.  It must have been templated
   *      on user data that has the same global ID distribution as this
   *      user data.
   *  \return   Returns the number of local Ids in the new partition.
   */

  template <typename Adapter>
    size_t applyPartitioningSolution(const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const
  {
//KDDDEC  This function should be void
//KDDDEC  This function should throw an error.
    return 0;
  }

};
  
}  //namespace Zoltan2
  
#endif
