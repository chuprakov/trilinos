/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFDOUBLE_HPP
#define TSFDOUBLE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFVector.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFVectorType.hpp"
#include "TSFEpetraVector.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFLinearCombination.hpp"

namespace TSFDouble
{
  typedef TSFExtended::VectorSpace<double> VectorSpace;
  typedef TSFExtended::VectorType<double> VectorType;
  typedef TSFExtended::EpetraVectorType EpetraVectorType;
  typedef TSFExtended::Vector<double> Vector;
}

#endif
