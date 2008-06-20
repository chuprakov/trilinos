#ifndef EXAMPLE_FACTORY_TRAITS_HPP
#define EXAMPLE_FACTORY_TRAITS_HPP

// mpl (Meta Programming Library) templates
#include "Sacado_mpl_vector.hpp"

// User Defined FieldEvaluator Types
#include "Evaluator_Constant.hpp"
#include "Evaluator_Density.hpp"
#include "Evaluator_EnergyFlux_Fourier.hpp"
#include "Evaluator_FEInterpolation.hpp"
#include "Evaluator_NonlinearSource.hpp"


#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

/*! \brief Struct to define FieldEvaluators objects for the FieldEvaluatorFactory.
    
    Preconditions:
    - You must provide a Sacado::mpl::vector named FieldEvaluatorTypes that contain all FieldEvaluator objects. 

*/
template<typename Traits>
struct MyFactoryTraits {
  
  static const int id_constant = 0;
  static const int id_density = 1;
  static const int id_fourier = 2;
  static const int id_feinterpolation = 3;
  static const int id_nonlinearsource = 4;

  typedef Sacado::mpl::vector< Constant<_,Traits>,             // 0
 			       Density<_,Traits>,              // 1
 			       Fourier<_,Traits>,              // 2
 			       FEInterpolation<_,Traits>,      // 3
 			       NonlinearSource<_,Traits>       // 4
  > FieldEvaluatorTypes;

};

#endif
