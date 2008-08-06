/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */

/*! \file TMPQualityMetricTest.cpp

Unit testing for the TMPQualityMetric class
\author Jasno Kraftcheck
*/


#include "TMPQualityMetric.hpp"
#include "TargetMetric2D.hpp"
#include "TargetMetric3D.hpp"
#include "UnitWeight.hpp"
#include "IdealTargetCalculator.hpp"
#include "MsqMatrix.hpp"
#include "QualityMetricTester.hpp"
#include "LinearFunctionSet.hpp"
#include "UnitUtil.hpp"
#include "PlanarDomain.hpp"
#include "SamplePoints.hpp"
#include "PatchData.hpp"
#include "Target3DShapeSizeOrient.hpp"
#include "Target2DShapeSizeOrient.hpp"

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <iostream.h>
#else
#include <iostream>
#endif

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;


/** Target metric (templatized by dimension) for use in misc. tests.
 *  'evaluate' method records input values and returns a constant.
 */
template <typename B>
class FauxTarget : public B
{
public:
  int count;
  double value;
  bool rval;
  MsqMatrix<B::MATRIX_DIM,B::MATRIX_DIM> last_A, last_W;

  FauxTarget(double v) : count(0), value(v), rval(true) {}
  
  bool evaluate( const MsqMatrix<B::MATRIX_DIM,B::MATRIX_DIM>& A, 
                 const MsqMatrix<B::MATRIX_DIM,B::MATRIX_DIM>& W, 
                 double& result, MsqError&  )
  {
    last_A = A;
    last_W = W;
    result = value;
    ++count;
    return rval;
  }
};

/** Weight calculator used for testing.  Returns constant weight. */
class ScaleWeight : public WeightCalculator
{
  public:
    ScaleWeight( double s ) : value(s) {}
    double get_weight( PatchData&, size_t, const SamplePoints*, unsigned, MsqError& )
      { return value; }
    double value;
};

/** wrapper class to force numeric approximation of derivatives */
template <class Base>
class NumericalTarget : public Base
{
public:
  
  NumericalTarget( Base* real_metric ) : mMetric(real_metric) {}

  ~NumericalTarget() {}

  bool evaluate( const MsqMatrix<Base::MATRIX_DIM,Base::MATRIX_DIM>& A, 
                 const MsqMatrix<Base::MATRIX_DIM,Base::MATRIX_DIM>& W, 
                 double& result, 
                 MsqError& err )
  { return mMetric->evaluate( A, W, result, err ); }
private:
  Base* mMetric;
};


/** Simple target metric for testing first partial derivatives.  
 *  \f$\mu(A,W) = |A|^2\f$
 *  \f$\frac{\partial\mu}{\partial \A} = 2 A \f$
 */
class TestGradTargetMetric3D : public TargetMetric3D
{
  public:
  
    bool evaluate( const MsqMatrix<3,3>& A, const MsqMatrix<3,3>&, double& result, MsqError& err )
      { result = sqr_Frobenius(A); return true; }
    
    bool evaluate_with_grad( const MsqMatrix<3,3>& A, 
                             const MsqMatrix<3,3>&,
                             double& result,
                             MsqMatrix<3,3>& d,
                             MsqError& err )
    {
      result = sqr_Frobenius(A);
      d = 2*A;
      return true;
    }
};

class TMPQualityMetricTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(TMPQualityMetricTest);
  
  CPPUNIT_TEST (test_negate_flag);
  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_evaluations);
  
  CPPUNIT_TEST (test_evaluate_2D);
  CPPUNIT_TEST (test_evaluate_3D);
  CPPUNIT_TEST (test_evaluate_2D_weight);
  CPPUNIT_TEST (test_evaluate_3D_weight);
  CPPUNIT_TEST (test_evaluate_2D_corner);
  CPPUNIT_TEST (test_evaluate_2D_edge);
  CPPUNIT_TEST (test_evaluate_2D_elem);
  CPPUNIT_TEST (test_evaluate_3D_corner);
  CPPUNIT_TEST (test_evaluate_3D_edge);
  CPPUNIT_TEST (test_evaluate_3D_edge);
  CPPUNIT_TEST (test_evaluate_3D_face);
  
  CPPUNIT_TEST (test_sample_indices);
  CPPUNIT_TEST (test_evaluate_with_indices);
  CPPUNIT_TEST (test_evaluate_fixed_indices);
  
  CPPUNIT_TEST (test_gradient_3D);
  CPPUNIT_TEST (compare_indices_and_gradient);
  CPPUNIT_TEST (test_ideal_element_gradient);
  CPPUNIT_TEST (compare_analytical_and_numerical_gradient);
  CPPUNIT_TEST (test_weighted_gradients);
  CPPUNIT_TEST (test_gradient_with_fixed_vertices);

  CPPUNIT_TEST (compare_indices_and_hessian);
  CPPUNIT_TEST (compare_gradient_and_hessian);
  CPPUNIT_TEST (compare_analytical_and_numerical_hessians);
  CPPUNIT_TEST (test_symmetric_hessian_diagonal);
  CPPUNIT_TEST (test_weighted_hessians);
  CPPUNIT_TEST (test_hessian_with_fixed_vertices);

  CPPUNIT_TEST (compare_indices_and_diagonal);
  CPPUNIT_TEST (compare_gradient_and_diagonal);
  CPPUNIT_TEST (compare_analytical_and_numerical_diagonals);
  CPPUNIT_TEST (test_weighted_diagonals);
  CPPUNIT_TEST (test_diagonal_with_fixed_vertices);

  CPPUNIT_TEST_SUITE_END();
  
  void test_2d_eval_ortho_quad( unsigned dim );
  void test_3d_eval_ortho_hex( unsigned dim );

  QualityMetricTester tester;

  LinearFunctionSet mf;
  SamplePoints corners, center;
  UnitWeight one;
  IdealTargetCalculator ideal;
  ScaleWeight e_weight;

  FauxTarget<TargetMetric2D> faux_2d_pi, faux_2d_zero;
  FauxTarget<TargetMetric3D> faux_3d_zero, faux_3d_two;
  Target3DShapeSizeOrient test_metric_3D;
  Target2DShapeSizeOrient test_metric_2D;
  NumericalTarget<TargetMetric3D> num_metric_3D;
  NumericalTarget<TargetMetric2D> num_metric_2D;
  TMPQualityMetric test_qm, zero_qm, weight_qm, center_qm;
  
public:
  TMPQualityMetricTest() : 
    tester( QualityMetricTester::ALL_FE_EXCEPT_SEPTAHEDRON, &mf ),
    corners( true, false, false, false ),
    center( false, false, false, true ),
    e_weight( 2.7182818284590451 ),
    faux_2d_pi(3.14159), faux_2d_zero(0.0),
    faux_3d_zero(0.0), faux_3d_two(2.0),
    num_metric_3D( &test_metric_3D ),
    num_metric_2D( &test_metric_2D ),
    test_qm( &corners, &ideal, &one, &num_metric_2D, &num_metric_3D ),
    zero_qm( &corners, &ideal, &one, &faux_2d_zero, &faux_3d_zero ),
    weight_qm( &corners, &ideal, &e_weight, &test_metric_2D, &test_metric_3D ),
    center_qm( &center, &ideal, &one, &test_metric_2D, &test_metric_3D )
  {  }
  
  void test_negate_flag()
    { CPPUNIT_ASSERT_EQUAL( 1, zero_qm.get_negate_flag() ); }
  void test_supported_types()     
    { tester.test_supported_element_types( &zero_qm ); }
  void test_get_evaluations();
  void test_get_element_evaluations();
  void test_evaluate_2D();
  void test_evaluate_3D();
  void test_evaluate_2D_weight();
  void test_evaluate_3D_weight();
  void test_evaluate_2D_corner()  { test_2d_eval_ortho_quad(0); }
  void test_evaluate_2D_edge()    { test_2d_eval_ortho_quad(1); }
  void test_evaluate_2D_elem()    { test_2d_eval_ortho_quad(3); }
  void test_evaluate_3D_corner()  { test_3d_eval_ortho_hex(0); }
  void test_evaluate_3D_edge()    { test_3d_eval_ortho_hex(1); }
  void test_evaluate_3D_face()    { test_3d_eval_ortho_hex(2); }
  void test_evaluate_3D_elem()    { test_3d_eval_ortho_hex(3); }
  void test_gradient_3D();

  void test_sample_indices()
    { tester.test_get_sample_indices( &zero_qm, &corners ); }
  void test_evaluate_with_indices()
    { tester.compare_eval_and_eval_with_indices( &zero_qm ); }
  void test_evaluate_fixed_indices() 
    { tester.test_get_indices_fixed( &zero_qm ); }
    
  void compare_indices_and_gradient()
    { tester.compare_eval_with_indices_and_eval_with_gradient( &test_qm ); }
  void test_ideal_element_gradient()
    { tester.test_ideal_element_zero_gradient( &test_qm, true ); }
  void compare_analytical_and_numerical_gradient()
    { tester.compare_analytical_and_numerical_gradients( &test_qm ); }
  void test_weighted_gradients()
    { tester.compare_analytical_and_numerical_gradients( &weight_qm ); }
  void test_gradient_with_fixed_vertices()
    { tester.test_gradient_with_fixed_vertex( &center_qm ); }

  void compare_indices_and_hessian()
    { tester.compare_eval_with_indices_and_eval_with_hessian( &test_qm ); }
  void compare_gradient_and_hessian()
    { tester.compare_eval_with_grad_and_eval_with_hessian( &test_qm );  }
  void compare_analytical_and_numerical_hessians()
    { tester.compare_analytical_and_numerical_hessians( &test_qm ); }
  void test_symmetric_hessian_diagonal()
    { tester.test_symmetric_Hessian_diagonal_blocks( &test_qm ); }
  void test_weighted_hessians()
    { tester.compare_analytical_and_numerical_hessians( &weight_qm ); }
  void test_hessian_with_fixed_vertices()
    { tester.test_hessian_with_fixed_vertex( &center_qm ); }

  void compare_indices_and_diagonal()
    { tester.compare_eval_with_indices_and_eval_with_diagonal( &test_qm ); }
  void compare_gradient_and_diagonal()
    { tester.compare_eval_with_grad_and_eval_with_diagonal( &test_qm ); }
  void compare_analytical_and_numerical_diagonals()
    { tester.compare_analytical_and_numerical_diagonals( &test_qm ); }
  void test_weighted_diagonals()
    { tester.compare_analytical_and_numerical_diagonals( &weight_qm ); }
  void test_diagonal_with_fixed_vertices()
    { tester.test_diagonal_with_fixed_vertex( &center_qm ); }
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TMPQualityMetricTest, "TMPQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TMPQualityMetricTest, "Unit");

void TMPQualityMetricTest::test_get_evaluations()
{
  SamplePoints corners_and_edges( true, true, false, false );
  TMPQualityMetric edge_metric( &corners_and_edges, &ideal, &one, &faux_2d_zero, &faux_3d_zero );
  
  tester.test_get_sample_evaluations( &zero_qm, &corners );
  tester.test_get_sample_evaluations( &edge_metric, &corners_and_edges );
}
  

void TMPQualityMetricTest::test_get_element_evaluations()
{
  SamplePoints edges( false, true, false, false );
  TMPQualityMetric edge_metric( &edges, &ideal, &one, &faux_2d_zero, &faux_3d_zero );
  
  tester.test_get_in_element_evaluations( &zero_qm );
  tester.test_get_in_element_evaluations( &edge_metric );
}

static double col_dot_prod( MsqMatrix<2,2>& m )
  { return m(0,0) * m(0,1) + m(1,0) * m(1,1); }

void TMPQualityMetricTest::test_evaluate_2D()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &one, &faux_2d_pi, &faux_3d_zero );
  
    // test with aligned elements
  faux_2d_pi.count = faux_3d_zero.count = 0;
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_pi.count );
  CPPUNIT_ASSERT_EQUAL( 0, faux_3d_zero.count );
  ASSERT_MATRICES_EQUAL( faux_2d_pi.last_W, faux_2d_pi.last_A, 1e-6 );
  
    // test that columns are orthogonal for ideal quad element
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, col_dot_prod(faux_2d_pi.last_A), 1e-6 );
  
    // test with an element rotated about X-axis
  faux_2d_pi.count = faux_3d_zero.count = 0;
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  // rotate by 90 degrees about X axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    double z = pd.vertex_by_index(i)[2];
    pd.vertex_by_index(i)[2] = pd.vertex_by_index(i)[1];
    pd.vertex_by_index(i)[1] =-z;
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_pi.count );
  CPPUNIT_ASSERT_EQUAL( 0, faux_3d_zero.count );
  ASSERT_MATRICES_EQUAL( faux_2d_pi.last_W, faux_2d_pi.last_A, 1e-6 );
  
    // test with an element rotated about Y-axis
  faux_2d_pi.count = faux_3d_zero.count = 0;
  tester.get_ideal_element( TRIANGLE, true, pd );
  // rotate by -90 degrees about Y axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    double z = pd.vertex_by_index(i)[2];
    pd.vertex_by_index(i)[2] = -pd.vertex_by_index(i)[0];
    pd.vertex_by_index(i)[0] = z;
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_pi.count );
  CPPUNIT_ASSERT_EQUAL( 0, faux_3d_zero.count );
  ASSERT_MATRICES_EQUAL( faux_2d_pi.last_W, faux_2d_pi.last_A, 1e-6 );
}  
 
  
void TMPQualityMetricTest::test_evaluate_3D()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &one, &faux_2d_zero, &faux_3d_two );
  
    // test with aligned elements
  faux_2d_zero.count = faux_3d_two.count = 0;
  tester.get_ideal_element( HEXAHEDRON, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_3d_two.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 0, faux_2d_zero.count );
  CPPUNIT_ASSERT_EQUAL( 1, faux_3d_two.count );
  ASSERT_MATRICES_EQUAL( faux_3d_two.last_W, faux_3d_two.last_A, 1e-6 );
  
    // test that columns are orthogonal for ideal hex element
  MsqMatrix<3,3> A = faux_3d_two.last_A;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(1) % A.column(2), 1e-6 );
  
    // test with rotated element
  faux_2d_zero.count = faux_3d_two.count = 0;
  tester.get_ideal_element( TETRAHEDRON, true, pd );
  // rotate by 90-degrees about X axis
  for (size_t i = 0; i < pd.num_nodes(); ++i) {
    double z = pd.vertex_by_index(i)[2];
    pd.vertex_by_index(i)[2] = pd.vertex_by_index(i)[1];
    pd.vertex_by_index(i)[1] =-z;
  }
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_3d_two.value, value, DBL_EPSILON );
  CPPUNIT_ASSERT_EQUAL( 0, faux_2d_zero.count );
  CPPUNIT_ASSERT_EQUAL( 1, faux_3d_two.count );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( determinant(faux_3d_two.last_W),
                                determinant(faux_3d_two.last_A),
                                1e-6 );
}


void TMPQualityMetricTest::test_evaluate_2D_weight()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &e_weight, &faux_2d_pi, &faux_3d_zero );
  
  tester.get_ideal_element( TRIANGLE, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_2d_pi.value*e_weight.value, value, DBL_EPSILON );
}


void TMPQualityMetricTest::test_evaluate_3D_weight()
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  TMPQualityMetric m( &corners, &ideal, &e_weight, &faux_2d_zero, &faux_3d_two );
  
  tester.get_ideal_element( PRISM, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( faux_3d_two.value*e_weight.value, value, DBL_EPSILON );
}

void TMPQualityMetricTest::test_2d_eval_ortho_quad( unsigned dim )
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  bool locs[4] = {false, false, false, false};
  locs[dim] = true;
  
  SamplePoints pts( locs[0], locs[1], locs[2], locs[3] );
  TMPQualityMetric m( &pts, &ideal, &one, &faux_2d_zero, &faux_3d_zero );
  faux_2d_zero.count = faux_3d_zero.count = 0;
  
  tester.get_ideal_element( QUADRILATERAL, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( 1, faux_2d_zero.count );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, col_dot_prod(faux_2d_zero.last_A), DBL_EPSILON );
}  

void TMPQualityMetricTest::test_3d_eval_ortho_hex( unsigned dim )
{
  MsqPrintError err(cout);
  PatchData pd;
  bool rval;
  double value;
  
  bool locs[4] = {false, false, false, false};
  locs[dim] = true;
  
  SamplePoints pts( locs[0], locs[1], locs[2], locs[3] );
  TMPQualityMetric m( &pts, &ideal, &one, &faux_2d_zero, &faux_3d_zero );
  faux_2d_zero.count = faux_3d_zero.count = 0;
  
  tester.get_ideal_element( HEXAHEDRON, true, pd );
  rval = m.evaluate( pd, 0, value, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  CPPUNIT_ASSERT(rval);
  CPPUNIT_ASSERT_EQUAL( 1, faux_3d_zero.count );
  
    // test that columns are orthogonal for ideal hex element
  MsqMatrix<3,3> A = faux_3d_zero.last_A;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(1), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(0) % A.column(2), 1e-6 );
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, A.column(1) % A.column(2), 1e-6 );
}  



void TMPQualityMetricTest::test_gradient_3D()
{
  MsqPrintError err(msq_stdio::cout);
  
    // check for expected value at center of flattened hex
  
    // construct flattened hex
  const double z = 0.5;
  const double vertices[] = { 0.0, 0.0, 0.0,
                              1.0, 0.0, 0.0,
                              1.0, 1.0, 0.0,
                              0.0, 1.0, 0.0,
                              0.0, 0.0, z,
                              1.0, 0.0, z,
                              1.0, 1.0, z,
                              0.0, 1.0, z };
  size_t conn[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  PatchData pd;
  pd.fill( 8, vertices, 1, HEXAHEDRON, conn, 0, err );
  ASSERT_NO_ERROR(err);
  
    // calculate Jacobian matrix at element center
  const double corner_xi[8][3] = { { -1, -1, -1 },
                                   {  1, -1, -1 },
                                   {  1,  1, -1 },
                                   { -1,  1, -1 },
                                   { -1, -1,  1 },
                                   {  1, -1,  1 },
                                   {  1,  1,  1 },
                                   { -1,  1,  1 } };
  MsqMatrix<8,3> coeff_derivs(&corner_xi[0][0]);
  coeff_derivs *= 0.125;  // derivatives of trilinear map at hex center
  MsqMatrix<8,3> coords( vertices );
  MsqMatrix<3,3> J = transpose(coords) * coeff_derivs;
    // calculate expected metric value
  const double expt_val = sqr_Frobenius( J );
    // calculate derivative for each element vertex
  MsqVector<3> expt_grad[8];
  for (int v = 0; v < 8; ++v)
    expt_grad[v] = 2 * J * transpose( coeff_derivs.row(v) );
    
  
    // construct metric
  LinearFunctionSet lfs;
  pd.set_mapping_functions( &lfs );
  SamplePoints center( false, false, false, true );
  TestGradTargetMetric3D tm;
  IdealTargetCalculator tc;
  UnitWeight wc;
  TMPQualityMetric m( &center, &tc, &wc, 0, &tm );
  
    // evaluate metric
  double act_val;
  msq_std::vector<size_t> indices;
  msq_std::vector<Vector3D> act_grad;
  m.evaluate_with_gradient( pd, 0, act_val, indices, act_grad, err );
  ASSERT_NO_ERROR(err);
  
    // compare values
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expt_val, act_val, 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[0]].data()), act_grad[0], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[1]].data()), act_grad[1], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[2]].data()), act_grad[2], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[3]].data()), act_grad[3], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[4]].data()), act_grad[4], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[5]].data()), act_grad[5], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[6]].data()), act_grad[6], 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[7]].data()), act_grad[7], 1e-10 );

    // check numerical approx of gradient
  m.QualityMetric::evaluate_with_gradient( pd, 0, act_val, indices, act_grad, err );
  ASSERT_NO_ERROR(err);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( expt_val, act_val, 1e-10 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[0]].data()), act_grad[0], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[1]].data()), act_grad[1], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[2]].data()), act_grad[2], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[3]].data()), act_grad[3], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[4]].data()), act_grad[4], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[5]].data()), act_grad[5], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[6]].data()), act_grad[6], 1e-5 );
  CPPUNIT_ASSERT_VECTORS_EQUAL( Vector3D(expt_grad[indices[7]].data()), act_grad[7], 1e-5 );
}
