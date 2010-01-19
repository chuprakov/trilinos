
INCLUDE(PackageLibraryMacros)
INCLUDE(PackageAddExecutableAndTest)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  FEApp_AbstractElement.hpp
  FEApp_LinearElement.hpp
  FEApp_ElementFactory.hpp
  FEApp_AbstractQuadrature.hpp
  FEApp_GaussianQuadrature2.hpp
  FEApp_QuadratureFactory.hpp
  FEApp_AbstractSourceFunction.hpp
  FEApp_ConstantSourceFunction.hpp
  FEApp_QuadraticSourceFunction.hpp
  FEApp_CubicSourceFunction.hpp
  FEApp_ExponentialSourceFunction.hpp
  FEApp_MultiVariateExponentialSourceFunction.hpp
  FEApp_SourceFunctionFactory.hpp
  FEApp_SourceFunctionFactoryImpl.hpp
  FEApp_AbstractFunction.hpp
  FEApp_ConstantFunction.hpp
  FEApp_KLExponentialFunction.hpp
  FEApp_KLExponentialFunctionImpl.hpp
  FEApp_FunctionFactory.hpp
  FEApp_FunctionFactoryImpl.hpp
  FEApp_AbstractPDE_NTBase.hpp
  FEApp_AbstractPDE.hpp
  FEApp_HeatNonlinearSourcePDE.hpp
  FEApp_HeatNonlinearSourcePDEImpl.hpp
  FEApp_BrusselatorPDE.hpp
  FEApp_BrusselatorPDEImpl.hpp
  FEApp_AbstractNodeBCStrategy_NTBase.hpp
  FEApp_AbstractNodeBCStrategy.hpp
  FEApp_NodeBC.hpp
  FEApp_BrusselatorNodeBCStrategy.hpp
  FEApp_BrusselatorNodeBCStrategyImpl.hpp
  FEApp_ConstantNodeBCStrategy.hpp
  FEApp_ConstantNodeBCStrategyImpl.hpp
  FEApp_TemplateTypes.hpp
  FEApp_AbstractProblem.hpp
  FEApp_BrusselatorProblem.hpp
  FEApp_HeatNonlinearSourceProblem.hpp
  FEApp_ProblemFactory.hpp
  FEApp_BrusselatorParameters.hpp
  FEApp_Mesh.hpp
  FEApp_AbstractInitPostOp.hpp
  FEApp_InitPostOps.hpp
  FEApp_GlobalFill.hpp
  FEApp_GlobalFillImpl.hpp
  FEApp_AbstractDiscretization.hpp
  FEApp_CZeroDiscretization.hpp
  FEApp_DiscretizationFactory.hpp
  FEApp_AbstractResponseFunction.hpp
  FEApp_BoundaryFlux1DResponseFunction.hpp
  FEApp_SolutionAverageResponseFunction.hpp
  FEApp_SolutionTwoNormResponseFunction.hpp
  FEApp_Application.hpp
  FEApp_ModelEvaluator.hpp
  ) 

APPEND_SET(SOURCES
  FEApp_LinearElement.cpp
  FEApp_ElementFactory.cpp
  FEApp_GaussianQuadrature2.cpp
  FEApp_QuadratureFactory.cpp
  FEApp_SourceFunctionFactory.cpp
  FEApp_FunctionFactory.cpp
  FEApp_HeatNonlinearSourcePDE.cpp
  FEApp_BrusselatorPDE.cpp
  FEApp_BrusselatorNodeBCStrategy.cpp
  FEApp_ConstantNodeBCStrategy.cpp
  FEApp_BrusselatorProblem.cpp
  FEApp_HeatNonlinearSourceProblem.cpp
  FEApp_ProblemFactory.cpp
  FEApp_Mesh.cpp
  FEApp_InitPostOps.cpp
  FEApp_GlobalFill.cpp
  FEApp_CZeroDiscretization.cpp
  FEApp_DiscretizationFactory.cpp
  FEApp_BoundaryFlux1DResponseFunction.cpp
  FEApp_SolutionAverageResponseFunction.cpp
  FEApp_SolutionTwoNormResponseFunction.cpp
  FEApp_Application.cpp
  FEApp_ModelEvaluator.cpp
  )

IF(Sacado_ENABLE_Stokhos)
  APPEND_SET(HEADERS
    FEApp_SGGaussQuadResidualGlobalFill.hpp
    FEApp_SGGaussQuadJacobianGlobalFill.hpp
  )
  APPEND_SET(SOURCES
    FEApp_SGGaussQuadResidualGlobalFill.cpp
    FEApp_SGGaussQuadJacobianGlobalFill.cpp
  )
ENDIF()

IF(Sacado_ENABLE_FEApp)

  PACKAGE_ADD_LIBRARY(
    FEApp
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ENDIF()
