/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <impl/Kokkos_Preprocessing_macros.hpp>


/********************************************************/

#define numNodesPerElem 8 /* don't change */
#define spatialDim 3

/********************************************************/

//----------------------------------------------------------------------------

template< typename Scalar , unsigned Dim , unsigned N >
struct TensorIntegration ;

template<typename Scalar >
struct TensorIntegration<Scalar,1,1> {
  Scalar pts[1] ;
  Scalar wts[1] ;

  TensorIntegration() { pts[0] = 0 ; wts[0] = 2 ; }
};

template<typename Scalar >
struct TensorIntegration<Scalar,1,2>
{
  Scalar pts[2] ;
  Scalar wts[2] ;

  TensorIntegration()
  {
    const Scalar x2 = 0.577350269 ;
    pts[0] = -x2; wts[0] = 1.0;
    pts[1] =  x2; wts[1] = 1.0;
  }
};

template<typename Scalar >
struct TensorIntegration<Scalar,1,3>
{
  Scalar pts[3] ;
  Scalar wts[3] ;

  TensorIntegration()
  {
    const Scalar x3 = 0.774596669 ;
    const Scalar w1 = 0.555555556 ;
    const Scalar w2 = 0.888888889 ;
    pts[0] =  -x3 ;  wts[0] = w1 ;
    pts[1] =    0 ;  wts[1] = w2 ;
    pts[2] =   x3 ;  wts[2] = w1 ;
  }
};

template< typename Scalar , unsigned Order >
struct TensorIntegration<Scalar,3,Order>
{
  enum { N = Order * Order * Order };

  Scalar pts[N][3] ;
  Scalar wts[N];

  TensorIntegration()
  {
    TensorIntegration<Scalar,1,Order> oneD ;

    unsigned n = 0 ;
    for ( unsigned k = 0 ; k < Order ; ++k ) {
    for ( unsigned j = 0 ; j < Order ; ++j ) {
    for ( unsigned i = 0 ; i < Order ; ++i , ++n ) {
      pts[n][0] = oneD.pts[i] ;
      pts[n][1] = oneD.pts[j] ;
      pts[n][2] = oneD.pts[k] ;
      wts[n] = oneD.wts[i] * oneD.wts[j] * oneD.wts[k] ;
    }}}
  }
};

//----------------------------------------------------------------------------

template< typename Scalar , typename ScalarCoord , class DeviceType >
struct assembleFE;

template<typename Scalar , typename ScalarCoord >
struct assembleFE<Scalar, ScalarCoord, KOKKOS_MACRO_DEVICE> {

  enum { numGaussPointsPerDim = 2 };

  typedef KOKKOS_MACRO_DEVICE                            device_type;
  typedef Scalar                                         scalar_type ;
  typedef device_type::size_type                         index_type ;
  typedef Kokkos::MDArrayView<index_type,device_type>    index_array ;
  typedef Kokkos::MDArrayView<Scalar,      device_type>  scalar_array ;
  typedef Kokkos::MDArrayView<ScalarCoord, device_type>  coord_array ;

  index_array   elem_node_ids ;
  coord_array   node_coords ;
  scalar_array  element_stiffness;
  scalar_array  element_vectors;
  TensorIntegration< Scalar , spatialDim , numGaussPointsPerDim > integration ;
  Scalar coeff_K ;
  Scalar coeff_Q ;

  assembleFE( const index_array  & arg_elem_node_ids ,
              const coord_array  & arg_node_coords ,
              const scalar_array & arg_element_stiffness , 
              const scalar_array & arg_element_vectors ,
              const Scalar       & arg_coeff_K ,
              const Scalar       & arg_coeff_Q )
  : elem_node_ids( arg_elem_node_ids )
  , node_coords(   arg_node_coords )
  , element_stiffness( arg_element_stiffness )
  , element_vectors( arg_element_vectors )
  , integration()
  , coeff_K( arg_coeff_K )
  , coeff_Q( arg_coeff_Q )
  {}


  //  Shape function values and gradients
  //  with respect to master element coordinate system

  KOKKOS_MACRO_DEVICE_FUNCTION
  void shape_fns( const Scalar* x,
                  Scalar * fn_value ,
                  Scalar * fn_grad )const
  {
    const Scalar ONE8TH = 0.125 ;

    const Scalar u = 1.0 - x[0];
    const Scalar v = 1.0 - x[1];
    const Scalar w = 1.0 - x[2];

    const Scalar up1 = 1.0 + x[0];
    const Scalar vp1 = 1.0 + x[1];
    const Scalar wp1 = 1.0 + x[2];

    // Vaues:
    fn_value[0] = ONE8TH *   u *   v *  w ;
    fn_value[1] = ONE8TH * up1 *   v *  w ;
    fn_value[2] = ONE8TH * up1 * vp1 *  w ;
    fn_value[3] = ONE8TH *   u * vp1 *  w ;

    fn_value[4] = ONE8TH *   u *   v *  wp1 ;
    fn_value[5] = ONE8TH * up1 *   v *  wp1 ;
    fn_value[6] = ONE8TH * up1 * vp1 *  wp1 ;
    fn_value[7] = ONE8TH *   u * vp1 *  wp1 ;

    //fn 0 = u * v * w
    fn_grad[ 0] = ONE8TH * -1  *  v  *  w  ;
    fn_grad[ 1] = ONE8TH *  u  * -1  *  w  ;
    fn_grad[ 2] = ONE8TH *  u  *  v  * -1  ;

    //fn 1 = up1 * v * w
    fn_grad[ 3] = ONE8TH *  1  *  v  *  w  ;
    fn_grad[ 4] = ONE8TH * up1 * -1  *  w  ;
    fn_grad[ 5] = ONE8TH * up1 *  v  * -1  ;

    //fn 2 = up1 * vp1 * w
    fn_grad[ 6] = ONE8TH *  1  * vp1 *  w ;
    fn_grad[ 7] = ONE8TH * up1 *  1  *  w ;
    fn_grad[ 8] = ONE8TH * up1 * vp1 * -1 ;

    //fn 3 = u * vp1 * w
    fn_grad[ 9] = ONE8TH * -1 * vp1 *  w ;
    fn_grad[10] = ONE8TH *  u *  1  *  w ;
    fn_grad[11] = ONE8TH *  u * vp1 * -1 ;

    //fn 4 = u * v * wp1
    fn_grad[12] = ONE8TH * -1  *  v  * wp1 ;
    fn_grad[13] = ONE8TH *  u  * -1  * wp1 ;
    fn_grad[14] = ONE8TH *  u  *  v  *  1  ;

    //fn 5 = up1 * v * wp1
    fn_grad[15] = ONE8TH *  1  *  v  * wp1 ;
    fn_grad[16] = ONE8TH * up1 * -1  * wp1 ;
    fn_grad[17] = ONE8TH * up1 *  v  *  1  ;

    //fn 6 = up1 * vp1 * wp1
    fn_grad[18] = ONE8TH *  1  * vp1 * wp1 ;
    fn_grad[19] = ONE8TH * up1 *  1  * wp1 ;
    fn_grad[20] = ONE8TH * up1 * vp1 *  1 ;

    //fn 7 = u * vp1 * wp1
    fn_grad[21] = ONE8TH * -1 * vp1 * wp1 ;
    fn_grad[22] = ONE8TH *  u *  1  * wp1 ;
    fn_grad[23] = ONE8TH *  u * vp1 *  1 ;
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void jacobian( const ScalarCoord * x, 
                 const ScalarCoord * y, 
                 const ScalarCoord * z, 
                 const Scalar * grad_vals, 
                 Scalar * J) const
  {
    J[0] = 0.0;
    J[1] = 0.0;
    J[2] = 0.0;
    J[3] = 0.0;
    J[4] = 0.0;
    J[5] = 0.0;
    J[6] = 0.0;
    J[7] = 0.0;
    J[8] = 0.0;

    int i_X_spatialDim = 0;

    for(int i = 0; i < 8; ++i) {
      J[0] += grad_vals[i_X_spatialDim+0]*x[i];
      J[1] += grad_vals[i_X_spatialDim+0]*y[i];
      J[2] += grad_vals[i_X_spatialDim+0]*z[i];

      J[3] += grad_vals[i_X_spatialDim+1]*x[i];
      J[4] += grad_vals[i_X_spatialDim+1]*y[i];
      J[5] += grad_vals[i_X_spatialDim+1]*z[i];

      J[6] += grad_vals[i_X_spatialDim+2]*x[i];
      J[7] += grad_vals[i_X_spatialDim+2]*y[i];
      J[8] += grad_vals[i_X_spatialDim+2]*z[i];

      i_X_spatialDim += spatialDim;
    }

  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  Scalar inverse_and_determinant3x3( Scalar * const J ) const
  {
    const Scalar J00 = J[0];
    const Scalar J01 = J[1];
    const Scalar J02 = J[2];

    const Scalar J10 = J[3];
    const Scalar J11 = J[4];
    const Scalar J12 = J[5];

    const Scalar J20 = J[6];
    const Scalar J21 = J[7];
    const Scalar J22 = J[8];

    const Scalar term0 = J22*J11 - J21*J12;
    const Scalar term1 = J22*J01 - J21*J02;
    const Scalar term2 = J12*J01 - J11*J02;

    const Scalar detJ = J00*term0 - J10*term1 + J20*term2;

    const Scalar inv_detJ = 1.0/detJ;

    J[0] =  term0*inv_detJ;
    J[1] = -term1*inv_detJ;
    J[2] =  term2*inv_detJ;

    J[3] = -(J22*J10 - J20*J12)*inv_detJ;
    J[4] =  (J22*J00 - J20*J02)*inv_detJ;
    J[5] = -(J12*J00 - J10*J02)*inv_detJ;

    J[6] =  (J21*J10 - J20*J11)*inv_detJ;
    J[7] = -(J21*J00 - J20*J01)*inv_detJ;
    J[8] =  (J11*J00 - J10*J01)*inv_detJ;

    return detJ ;
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void matTransMat3x3_X_3xn(const Scalar * A, int n, const Scalar * B, Scalar * C)const {

    //A is 3x3, B is 3xn. So C is also 3xn.
    //A,B,C are all assumed to be ordered such that columns are contiguous.

    Scalar* Cj = C;
    const Scalar* Bj = B;

    for(int j=0; j<n; ++j) {
      Cj[0] = A[0]*Bj[0] + A[1]*Bj[1] + A[2]*Bj[2];
      Cj[1] = A[3]*Bj[0] + A[4]*Bj[1] + A[5]*Bj[2];
      Cj[2] = A[6]*Bj[0] + A[7]*Bj[1] + A[8]*Bj[2];
      Bj += 3;
      Cj += 3;
    }

  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void contributeDiffusionMatrix(
    const Scalar weight ,
    const Scalar grad_vals[] ,
    const Scalar invJ[] ,
    Scalar elem_stiff[][8] ) const
  {
    Scalar dpsidx[8], dpsidy[8], dpsidz[8];

    for(int i = 0; i < numNodesPerElem ; i++){

      dpsidx[i] = grad_vals[i * 3 + 0] * invJ[0] + 
                  grad_vals[i * 3 + 1] * invJ[1] + 
                  grad_vals[i * 3 + 2] * invJ[2];
      dpsidy[i] = grad_vals[i * 3 + 0] * invJ[3] + 
                  grad_vals[i * 3 + 1] * invJ[4] + 
                  grad_vals[i * 3 + 2] * invJ[5];
      dpsidz[i] = grad_vals[i * 3 + 0] * invJ[6] + 
                  grad_vals[i * 3 + 1] * invJ[7] + 
                  grad_vals[i * 3 + 2] * invJ[8];
    }

    for(int m = 0; m < numNodesPerElem; m++) {
      for(int n = 0; n < numNodesPerElem; n++) {

        elem_stiff[m][n] += weight * 
          ((dpsidx[m] * dpsidx[n]) + 
           (dpsidy[m] * dpsidy[n]) +
           (dpsidz[m] * dpsidz[n]));            
      }
    }
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void contributeSourceVector( const Scalar term ,
                               const Scalar psi[] ,
                               Scalar elem_vec[] ) const
  {
     for(int i=0; i<numNodesPerElem; ++i) {
       elem_vec[i] += psi[i] * term ;
     }
  }

  
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int ielem )const {

    Scalar elem_vec[8] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
    Scalar elem_stiff[8][8] =
      { { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } };

    ScalarCoord x[8], y[8], z[8];

    Scalar psi[numNodesPerElem];
    Scalar grad_vals[numNodesPerElem*spatialDim];

    Scalar J[spatialDim*spatialDim];

    for ( int i = 0 ; i < 8 ; ++i ) {
      const int node_index = elem_node_ids( ielem , i );
      x[i] = node_coords( node_index , 0 );
      y[i] = node_coords( node_index , 1 );
      z[i] = node_coords( node_index , 2 );
    }

    // This loop could be parallelized; however,
    // it would require per-thread temporaries
    // of 'elem_vec' and 'elem_stiff' which would
    // have to be reduced.

    for ( int i = 0 ; i < integration.N ; ++i ) {

      shape_fns(integration.pts[i] , psi, grad_vals);

      jacobian( x, y, z, grad_vals, J );

      // Overwrite J with its inverse
      const Scalar detJ = inverse_and_determinant3x3(J);

      const Scalar k_detJ_w = coeff_K * detJ * integration.wts[i] ;
      const Scalar Q_detJ_w = coeff_Q * detJ * integration.wts[i] ;

      contributeDiffusionMatrix( k_detJ_w , grad_vals , J , elem_stiff );

      contributeSourceVector( Q_detJ_w , psi , elem_vec );
    }

    for(int i=0; i<numNodesPerElem; ++i) {
      element_vectors(ielem, i) = elem_vec[i] ;
    }

    for(int i = 0; i < 8; i++){
      for(int j = 0; j < 8; j++){
        element_stiffness(ielem, i, j) = elem_stiff[i][j] ;
      }
    }
  }


};// assembleFE

