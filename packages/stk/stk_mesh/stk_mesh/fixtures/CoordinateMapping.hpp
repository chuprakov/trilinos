/*------------------------------------------------------------------------*/
/*                 Copyright 2014 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_COORDINATE_MAPPING_HPP
#define STK_MESH_FIXTURES_COORDINATE_MAPPING_HPP

#include <math.h>                       // for cos, sin
#include <stddef.h>                     // for size_t, NULL
#include <stk_mesh/base/Field.hpp>      // for Field

namespace stk {
namespace mesh {
namespace fixtures {

/**
 * A mapping for the coordinates as a function of the node indices
 */
class CoordinateMapping
{
public:
  typedef double Scalar;
  CoordinateMapping() {}
  virtual void getNodeCoordinates(Scalar * field, const size_t nx, const size_t ny, const size_t nz) const = 0;
  virtual ~CoordinateMapping() {};
};


/**
 * Standard Cartesian X-Y-Z coordinate mapping
 */
class CartesianCoordinateMapping : public CoordinateMapping
{
public:
  CartesianCoordinateMapping() : CoordinateMapping() {}
  virtual void getNodeCoordinates(Scalar * field, const size_t nx, const size_t ny, const size_t nz) const
  {
    field[0] = nx;
    field[1] = ny;
    field[2] = nz;
  }
};


/**
 * Cylindrical coordinate mapping where the standard X-Y-Z coordinates
 * are mapped onto a wedge slice specified by the inner radius and the
 * maximum angle theta.
 */
class CylindricalCoordinateMapping : public CoordinateMapping
{
public:
  CylindricalCoordinateMapping(Scalar radius, Scalar theta, size_t numTheta)
      : CoordinateMapping(), m_radius(radius), m_theta(theta), m_numTheta(numTheta)
  { }
  virtual void getNodeCoordinates(Scalar * field, const size_t nx, const size_t ny, const size_t nz) const
  {
    Scalar fracTheta = nx/(m_numTheta - 1);

    //we want the angle to go from pi/2 to pi/2 - theta so we do not
    //invert any elements
    Scalar angle = M_PI/2.0 + m_theta*fracTheta;
    field[0] = (m_radius + ny)*std::cos(angle);
    field[1] = (m_radius + ny)*std::sin(angle);
    field[2] = nz;
  }
private:
  Scalar m_radius;
  Scalar m_theta;
  size_t m_numTheta;
};

}}}

#endif
