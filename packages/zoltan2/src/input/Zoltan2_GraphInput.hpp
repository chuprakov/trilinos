// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GraphInput.hpp
    \brief Defines the GraphInput interface.
*/


#ifndef _ZOLTAN2_GRAPHINPUT_HPP_
#define _ZOLTAN2_GRAPHINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

/*!  \brief GraphInput defines the interface for graph input adapters.

    InputAdapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t vertex and edge weights 
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c gid_t    application global Ids
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel

    See IdentifierTraits to understand why the user's global ID type (\c gid_t)
    may differ from that used by Zoltan2 (\c gno_t).

    The Kokkos node type can be safely ignored.

    The template parameter \c User is a user-defined data type
    which, through a traits mechanism, provides the actual data types
    with which the Zoltan2 library will be compiled.
    \c User may be the actual class or structure used by application to
    represent a vector, or it may be the helper class BasicUserTypes.
    See InputTraits for more information.
*/

template <typename User>
  class GraphInput : public InputAdapter {
private:

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  enum InputAdapterType inputAdapterType() const {return GraphAdapterType;}

  /*! \brief Destructor
   */
  virtual ~GraphInput() {};

  /*! \brief Returns the number vertices on this process.
   */
  virtual size_t getLocalNumberOfVertices() const = 0;

  /*! \brief Returns the global number vertices.
   *   \todo For GraphInput, should a user have to tell us
   *            global number of vertices and edges?  Zoltan2
   *            can figure that out.
   */
  virtual global_size_t getGlobalNumberOfVertices() const = 0;

  /*! \brief Returns the number edges on this process.
   */
  virtual size_t getLocalNumberOfEdges() const = 0;

  /*! \brief Returns the global number edges.
   */
  virtual global_size_t getGlobalNumberOfEdges() const = 0;

  /*! \brief Returns the dimension (0 or greater) of vertex weights.
   */
  virtual int getVertexWeightDimension() const = 0;

  /*! \brief Returns the dimension (0 or greater) of edge weights.
   */
  virtual int getEdgeWeightDimension() const = 0;

  /*! \brief Returns the dimension of the geometry, if any.
   *
   *  Some algorithms can use geometric vertex coordinate 
   *    information if it is present.
   */
  virtual int getCoordinateDimension() const = 0;

  /*! \brief Sets pointers to this process' graph entries.
      \param vertexIds will on return a pointer to vertex global Ids
      \param offsets is an array of size numVertices + 1.  
         The edge Ids for vertexId[i] begin at edgeIds[offsets[i]].  
          The last element of offsets
          is the size of the edgeIds array.
      \param edgeIds on return will point to the global edge Ids for
         for each vertex.
       \return The number of ids in the vertexIds list.

      Zoltan2 does not copy your data.  The data pointed to by 
      vertexIds, offsets and edgeIds
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getVertexListView(const gid_t *&vertexIds, 
    const lno_t *&offsets, const gid_t *& edgeIds) const = 0; 

  /*! \brief  Provide a pointer to the vertex weights, if any.

      \param weightDim ranges from zero to one less than 
                   getVertexWeightDimension().
      \param weights is the list of weights of the given dimension for
           the vertices returned in getVertexListView().
       \param stride The k'th weight is located at weights[stride*k]
       \return The number of weights listed, which should be the same
               as the number of vertices in getVertexListView().

      Zoltan2 does not copy your data.  The data pointed to by weights
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getVertexWeights(int weightDim,
     const scalar_t *&weights, int &stride) const = 0;

  /*! \brief  Provide a pointer to the edge weights, if any.

      \param weightDim ranges from zero to one less than 
                   getEdgeWeightDimension().
      \param weights is the list of weights of the given dimension for
           the edges returned in getVertexListView().
       \param stride The k'th weight is located at weights[stride*k]
       \return The number of weights listed, which should be the same
               as the number of edges in getVertexListView().

      Zoltan2 does not copy your data.  The data pointed to by weights
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getEdgeWeights(int weightDim,
     const scalar_t *&weights, int &stride) const = 0;

  /*! \brief Provide a pointer to one dimension of vertex coordinates.
      \param coordDim  is a value from 0 to one less than
         getCoordinateDimension() specifying which dimension is
         being provided in the coords list.
      \param coords  points to a list of coordinate values for the dimension.
      \param stride  describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].

       \return The length of the \c coords list.  This may be more than
              getLocalNumberOfVertices() because the \c stride
              may be more than one.

      Zoltan2 does not copy your data.  The data pointed to by coords
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getVertexCoordinates(int coordDim, 
    const scalar_t *&coords, int &stride) const = 0;


  /*! \brief Apply a partitioning problem solution to an input.  
   *
   *  This is not a required part of the GraphInput interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *           of global Ids to processes that matches that of the input
   *           data that instantiated this InputAdapter.
   *  \param out On return this should point to a newly created object 
   *            with the specified partitioning.
   *  \param solution  The Solution object created by a Problem should
   *      be supplied as the third argument.  It must have been templated
   *      on user data that has the same global ID distribution as this
   *      user data.
   *  \return   Returns the number of local Ids in the new partitioning.
   */

  template <typename User2>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<User2> &solution) const
  {
    return 0;
  }
  
};
  
}  //namespace Zoltan2
  
#endif
