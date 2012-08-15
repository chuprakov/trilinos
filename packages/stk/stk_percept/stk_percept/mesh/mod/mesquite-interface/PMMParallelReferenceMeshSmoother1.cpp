#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/JacobianUtil.hpp>
#include <stk_percept/math/Math.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stdio.h>

#include "mpi.h"
#include <cstdio>

#define DEBUG_PRINT 0
#define PRINT(a) do { if (DEBUG_PRINT && !m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_1(a) do { if (!m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_2(a) do {  std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << " "; } while(0)

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;
    const bool do_tot_test = false;
    bool do_print_elem_val = false;

    double PMMParallelReferenceMeshSmoother1::nodal_gradient(stk::mesh::Entity& node, double alpha, double *coord_current, double *cg_d, bool& valid, double *ng )
    {
      int spatialDim = m_eMesh->get_spatial_dim();
      valid = true;
      double nm=0.0;

      double xc[3]={0,0,0};
      double eps=1.e-6;
      double eps1 = eps*nodal_edge_length_ave(node);

      for (int i=0; i < spatialDim; i++)
        {
          xc[i]=coord_current[i];
          double dt = eps1;
          coord_current[i] += dt;  
          double mp = nodal_metric(node, 0.0, coord_current, cg_d, valid);
          coord_current[i] -= 2.0*dt;
          double mm = nodal_metric(node, 0.0, coord_current, cg_d, valid);
          coord_current[i] = xc[i];
          ng[i] = (mp-mm)/(2.0*eps1);
        }

      return nm;
    }

    double PMMParallelReferenceMeshSmoother1::nodal_metric(stk::mesh::Entity& node, double alpha, double *coord_current, double *cg_d, bool& valid )
    {
      int spatialDim = m_eMesh->get_spatial_dim();
      valid = true;
      double nm=0.0;

      double xc[3]={0,0,0};

      for (int i=0; i < spatialDim; i++)
        {
          xc[i]=coord_current[i];
          double dt = alpha*cg_d[i];
          coord_current[i] += dt;  
        }

      stk::mesh::PairIterRelation node_elems = node.relations(m_eMesh->element_rank());
      for (unsigned i_elem=0; i_elem < node_elems.size(); i_elem++)
        {
          stk::mesh::Entity& element = *node_elems[i_elem].entity();
          bool local_valid = true;
          double val = m_metric->metric(element, local_valid);
          nm += val;
          valid = valid && local_valid;
        }
      for (int i=0; i < spatialDim; i++)
        {
          coord_current[i] = xc[i];
        }

      return nm;
    }

    double PMMParallelReferenceMeshSmoother1::nodal_edge_length_ave(stk::mesh::Entity& node)
    {
      //int spatialDim = m_eMesh->get_spatial_dim();
      double nm=0.0;

      stk::mesh::PairIterRelation node_elems = node.relations(m_eMesh->element_rank());
      for (unsigned i_elem=0; i_elem < node_elems.size(); i_elem++)
        {
          stk::mesh::Entity& element = *node_elems[i_elem].entity();
          double elem_edge_len = m_eMesh->edge_length_ave(element);
          nm += elem_edge_len;
        }
      nm /= double(node_elems.size());
      return nm;
    }


    bool PMMParallelReferenceMeshSmoother1::check_convergence()
    {
      if (m_stage == 0 && (m_dnew == 0.0 || m_total_metric == 0.0))
        {
          return true; // for untangle
        }
      if (m_stage == 0 && m_num_invalid == 0 && m_dmax < gradNorm)
        {
          return true;
        }
      if (m_num_invalid == 0 && (m_dnew < gradNorm*gradNorm*m_d0 && m_dmax < gradNorm))
        {
          return true;
        }      
      return false;
    }

    double PMMParallelReferenceMeshSmoother1::metric(stk::mesh::Entity& element, bool& valid)
    {
      return m_metric->metric(element,valid);
    }

    /// fills cg_g_field with f'(x)
    void PMMParallelReferenceMeshSmoother1::get_gradient( Mesh* mesh, MeshDomain *domain)
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field;
      stk::mesh::FieldBase *cg_g_field    = eMesh->get_field("cg_g");

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = eMesh->get_spatial_dim();
      bool valid=true;

      m_scale = 1.e-10;

      // g=0
      eMesh->nodal_field_set_value(cg_g_field, 0.0);

      // element loop: compute deltas
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (PerceptMesquiteMesh::select_bucket(**k, m_eMesh))
            //if (on_locally_owned_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();
                m_metric->m_topology_data = m_eMesh->get_cell_topology(bucket);

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];

                    const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    double edge_length_ave = m_eMesh->edge_length_ave(element);

                    double metric_0 = metric(element, valid);

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        mesh::Entity & node = * elem_nodes[ inode ].entity();

                        bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                        //VERIFY_OP_ON(isGhostNode, ==, false, "hmmmm");
                        bool node_locally_owned = (eMesh->get_rank() == node.owner_rank());
                        bool fixed = pmm->get_fixed_flag(&node);
                        if (fixed || isGhostNode)
                          continue;

                        double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                        double *cg_g = PerceptMesh::field_data(cg_g_field, node);
                        
                        double eps = 1.e-6;
                        double eps1 = eps*edge_length_ave;

                        double gsav[3]={0,0,0};
                        for (int idim=0; idim < spatialDim; idim++)
                          {
                            double cc = coord_current[idim];
                            coord_current[idim] += eps1;
                            bool pvalid=false, mvalid=false;
                            double mp = metric(element, pvalid);
                            coord_current[idim] -= 2.0*eps1;
                            double mm = metric(element, mvalid);
                            coord_current[idim] = cc;
                            double dd = 0.0;
                            if ((pvalid && mvalid) || m_stage == 0)
                              dd = (mp - mm)/(2*eps1);
                            gsav[idim] = dd;
                            m_scale = std::max(m_scale, std::abs(dd)/edge_length_ave);

                            if (node_locally_owned)
                              {
                                cg_g[idim] += dd;
                              }
                            else
                              {
                                cg_g[idim] = 0.0;
                              }
                          }

                        // FIXME 
                        if (0)
                          {
                            double gsn=0.0;
                            for (int idim=0; idim < spatialDim; idim++)
                              {
                                gsn += gsav[idim]*gsav[idim];
                              }
                            gsn = std::sqrt(gsn);
                            gsn = std::max(gsn,1.e-12);
                            if (gsn > 1.e-4*edge_length_ave)
                              {
                                for (int idim=0; idim < spatialDim; idim++)
                                  {
                                    double dd = gsav[idim]/gsn*edge_length_ave;
                                    coord_current[idim] -= dd*eps;
                                  }
                                double m1=metric(element, valid);
                                for (int idim=0; idim < spatialDim; idim++)
                                  {
                                    double dd = gsav[idim]/gsn*edge_length_ave;
                                    coord_current[idim] += dd*eps;
                                  }
                                if (metric_0 > 1.e-6)
                                  {
                                    if (!(m1 < metric_0*(1.0+eps)))
                                      {
                                        PRINT( "bad grad" << " m1-metric_0 = " << (m1-metric_0) << " gsav= " << gsav[0] << " " << gsav[1] << " " << gsav[2] );
                                      }
                                    VERIFY_OP_ON(m1, <, metric_0*(1.0+eps), "bad gradient");
                                  }
                              }
                          }

                      } // inode...
                  } // i_element
              } // on_locally_owned_part...
          } // buckets
      } // element loop...

      VectorFieldType *cg_g_field_v = static_cast<VectorFieldType *>(cg_g_field);
      stk::mesh::parallel_reduce(*m_eMesh->get_bulk_data(), stk::mesh::sum(*cg_g_field_v));

      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(cg_g_field);

        // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
        stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->shared_aura(), fields); 

        // the shared part (just the shared boundary)
        //stk::mesh::communicate_field_data(*m_eMesh->get_bulk_data()->ghostings()[0], fields);
      }

      {
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMax<1>( & m_scale ) );
        m_scale = (m_scale < 1.0) ? 1.0 : 1.0/m_scale;
        PRINT("tmp srk m_scale= " << m_scale);
      }
        
      //get_scale(mesh, domain);
    }

    /// fills cg_g_field with f'(x)
    void PMMParallelReferenceMeshSmoother1::get_scale( Mesh* mesh, MeshDomain *domain)
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::FieldBase *cg_g_field    = eMesh->get_field("cg_g");

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = eMesh->get_spatial_dim();

      m_scale = 1.e-10;

      // element loop
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (PerceptMesquiteMesh::select_bucket(**k, m_eMesh) && on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];

                    const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    double edge_length_ave = m_eMesh->edge_length_ave(element);

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        mesh::Entity & node = * elem_nodes[ inode ].entity();

                        bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                        VERIFY_OP_ON(isGhostNode, ==, false, "hmmmm");
                        bool fixed = pmm->get_fixed_flag(&node);
                        if (fixed || isGhostNode)
                          continue;

                        double *cg_g = PerceptMesh::field_data(cg_g_field, node);
                        
                        for (int idim=0; idim < spatialDim; idim++)
                          {
                            m_scale = std::max(m_scale, std::abs(cg_g[idim])/edge_length_ave);
                          }
                      }
                  }
              }
          }
      }

      {
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMax<1>( & m_scale ) );
        m_scale = (m_scale < 1.0) ? 1.0 : 1.0/m_scale;
        PRINT("tmp srk m_scale= " << m_scale);
      }

    }

    double PMMParallelReferenceMeshSmoother1::run_one_iteration( Mesh* mesh, MeshDomain *domain,
                                                              MsqError& err )
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      m_pmm  = pmm;
      PerceptMesh *eMesh = pmm->getPerceptMesh();

      stk::mesh::FieldBase *cg_g_field    = eMesh->get_field("cg_g");
      stk::mesh::FieldBase *cg_r_field    = eMesh->get_field("cg_r");
      stk::mesh::FieldBase *cg_d_field    = eMesh->get_field("cg_d");
      stk::mesh::FieldBase *cg_s_field    = eMesh->get_field("cg_s");

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      bool total_valid=true;

      if (m_iter == 0)
        {
          m_dmax = 0.0;

          //PRINT_1("tmp srk get_gradient at m_iter=0");
          get_gradient(mesh, domain);
          /// r = -g
          eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);
          /// s = r  (allows for preconditioning later s = M^-1 r)
          eMesh->copy_field(cg_s_field, cg_r_field);
          /// d = s
          eMesh->copy_field(cg_d_field, cg_s_field);
          /// dnew = r.d
          m_dnew = eMesh->nodal_field_dot(cg_r_field, cg_d_field);
          PRINT("tmp srk m_dnew[0] = " << m_dnew);

          // FIXME
          if (0)
            {
              double coord_mag = eMesh->nodal_field_dot(m_coord_field_current, m_coord_field_current);
              if (!m_eMesh->get_rank()) printf("tmp srk m_dnew[%d] = %30.10g m_dmid= %30.10g coord_mag= %30.10g\n",  m_iter, m_dnew, m_dmid, coord_mag);
            }

          /// d0 = dnew
          m_d0 = m_dnew;
        }

      double metric_check = total_metric(mesh, 0.0, 1.0, total_valid);
      m_total_metric = metric_check;
      if (check_convergence() || metric_check == 0.0)
        {
          PRINT_1( "tmp srk already converged m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " m_d0= " << m_d0 );
          //update_node_positions
          return total_metric(mesh,0.0,1.0, total_valid);
        }

      //m_dd = eMesh->nodal_field_dot(cg_d_field, cg_d_field);
      //PRINT(" tmp srk m_dd= " << m_dd);

      /// line search
      bool restarted=false;
      double alpha = m_scale;
      {
        //double metric_1 = total_metric(mesh, 1.e-6, 1.0, total_valid);

        double metric_0 = total_metric(mesh, 0.0, 1.0, total_valid);
        //double norm_gradient = eMesh->nodal_field_dot(cg_g_field, cg_g_field);
        //PRINT( "tmp srk norm_gradient= " << norm_gradient );
        //PRINT( "tmp srk " << " metric_0= " << metric_0 << " metric(1.e-6) = " << metric_1 << " diff= " << metric_1-metric_0 );
        //metric_1 = total_metric(mesh, -1.e-6, 1.0, total_valid);
        //PRINT( "tmp srk " << " metric_0= " << metric_0 << " metric(-1.e-6)= " << metric_1 << " diff= " << metric_1-metric_0 );
        double metric=0.0;
        //double sigma=0.95;
        double tau = 0.5;
        double c0 = 1.e-4;
        double min_alpha_factor=1.e-12;

        //PRINT_1("tmp srk get_gradient at linesearch");
        //get_gradient(mesh, domain);
        double norm_gradient2 = eMesh->nodal_field_dot(cg_g_field, cg_g_field);

        double armijo_offset_factor = c0*norm_gradient2;
        bool converged = false;
        while (!converged)
          {
            metric = total_metric(mesh, alpha, 1.0, total_valid);

            //converged = (metric > sigma*metric_0) && (alpha > 1.e-16);
            double mfac = alpha*armijo_offset_factor;
            converged = (metric < metric_0 + mfac);
            if (m_untangled) converged = converged && total_valid;
            PRINT(  "tmp srk alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac) 
                    << " m_untangled = " << m_untangled
                    << " total_valid= " << total_valid );
            if (!converged)
              alpha *= tau;
            if (alpha < std::max(min_alpha_factor*m_scale, 1.e-16))
              break;
          }

        //if (metric > sigma*metric_0)
        if (!converged)
          {
            alpha=m_scale;
            /// d = -g
            eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_d_field);
            restarted = true;

            PRINT_1( "can't reduce metric= " << metric << " metric_0 + armijo_offset " << metric_0+alpha*armijo_offset_factor << " norm_gradient = " << std::sqrt(norm_gradient2) );
            double metric_1 = total_metric(mesh, 1.e-6, 1.0, total_valid);
            metric_0 = total_metric(mesh, 0.0, 1.0, total_valid);
            PRINT_1( "tmp srk " << " metric_0= " << metric_0 << " metric(1.e-6) = " << metric_1 << " diff= " << metric_1-metric_0 );
            metric_1 = total_metric(mesh, -1.e-6, 1.0, total_valid);
            PRINT_1( "tmp srk " << " metric_0= " << metric_0 << " metric(-1.e-6)= " << metric_1 << " diff= " << metric_1-metric_0 );

            //throw std::runtime_error("can't reduce metric");

          }

        while (!converged)
          {
            metric = total_metric(mesh, alpha, 1.0, total_valid);

            //converged = (metric > sigma*metric_0) && (alpha > 1.e-16);
            double mfac = alpha*armijo_offset_factor;
            converged = (metric < metric_0 + mfac);
            if (m_untangled) converged = converged && total_valid;
            PRINT(  "tmp srk ### alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac) 
                    << " m_untangled = " << m_untangled
                    << " total_valid= " << total_valid );
            if (!converged)
              alpha *= tau;
            if (alpha < std::max(min_alpha_factor*m_scale, 1.e-16))
              break;
          }

        //if (metric > sigma*metric_0)
        if (!converged)
          {
            PRINT_1( "can't reduce metric= " << metric << " metric_0 + armijo_offset " << metric_0+alpha*armijo_offset_factor << " norm_gradient = " << std::sqrt(norm_gradient2) );
            do_print_elem_val = true;
            double metric_1 = total_metric(mesh, 1.e-6, 1.0, total_valid);
            metric_0 = total_metric(mesh, 0.0, 1.0, total_valid);
            do_print_elem_val = false;
            PRINT_1( "tmp srk " << " metric_0= " << metric_0 << " metric(1.e-6) = " << metric_1 << " diff= " << metric_1-metric_0 );
            metric_1 = total_metric(mesh, -1.e-6, 1.0, total_valid);
            PRINT_1( "tmp srk " << " metric_0= " << metric_0 << " metric(-1.e-6)= " << metric_1 << " diff= " << metric_1-metric_0 );

            throw std::runtime_error("can't reduce metric");
          }

        else
          {
            double a1 = alpha/2.;
            double a2 = alpha;
            double f0 = metric_0, f1 = total_metric(mesh, a1, 1.0, total_valid), f2 = total_metric(mesh, a2, 1.0, total_valid);
            double den = 2.*(a2*(-f0 + f1) + a1*(f0 - f2));
            double num = a2*a2*(f1-f0)+a1*a1*(f0-f2);
            if (std::fabs(den) > 1.e-10)
              {
                double alpha_quadratic = num/den;
                if (alpha_quadratic < 2*alpha)
                  {
                    double fm=total_metric(mesh, alpha_quadratic, 1.0, total_valid);
                    //if (fm < f2 && (!m_untangled || total_valid))
                    if (fm < f2)
                      {
                        alpha = alpha_quadratic;
                        PRINT( "tmp srk alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 );
                      }
                  } 
              }
          }
      }

      /// x = x + alpha*d
      update_node_positions(mesh, alpha);
      //PRINT_1( "tmp srk iter= "<< m_iter << " dmax= " << m_dmax << " alpha= " << alpha);
      bool debug_par = false;

      // FIXME
      double coord_mag = 0.0;
      if (debug_par) coord_mag= eMesh->nodal_field_dot(m_coord_field_current, m_coord_field_current);

      /// f'(x)
      //PRINT_1("tmp srk get_gradient after update_node_positions");
      get_gradient(mesh, domain);
      double d_g = 0.0;
      if (debug_par) d_g = eMesh->nodal_field_dot(cg_g_field, cg_g_field);
      if (debug_par) debug_print(alpha);

      /// r = -g
      eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);
      double d_r = 0.0;
      if (debug_par) d_r = eMesh->nodal_field_dot(cg_r_field, cg_r_field);

      /// dold = dnew
      m_dold = m_dnew;

      /// dmid = r.s
      m_dmid = eMesh->nodal_field_dot(cg_r_field, cg_s_field);

      /// s = r
      eMesh->copy_field(cg_s_field, cg_r_field);

      /// dnew = r.s
      m_dnew = eMesh->nodal_field_dot(cg_r_field, cg_s_field);
      //double f_cur=total_metric(mesh, 0.0, 1.0, total_valid);
      //PRINT("tmp srk m_dnew[n] = " << m_dnew << " f_cur= " << f_cur << " total_valid= " << total_valid);

      double cg_beta = 0.0;
      if (std::fabs(m_dold) < 1.e-12) 
        cg_beta = 0.0;
      else
        cg_beta = (m_dnew - m_dmid) / m_dold;

      PRINT("tmp srk beta = " << cg_beta);


      //int N = num_nodes;
      //if (m_iter == N || cg_beta <= 0.0) 
      if (cg_beta <= 0.0 || restarted) 
        {
          /// d = s
          eMesh->copy_field(cg_d_field, cg_s_field);
        }
      else
        {
          /// d = s + beta * d
          eMesh->nodal_field_axpby(1.0, cg_s_field, cg_beta, cg_d_field);
        }
      //m_dd = eMesh->nodal_field_dot(cg_d_field, cg_d_field);

      double tm = total_metric(mesh,0.0,1.0, total_valid);

      if (debug_par)
        {
          double d_s=eMesh->nodal_field_dot(cg_s_field, cg_s_field);
          if (0 && !m_eMesh->get_rank()) printf("dmax[%3d]=%20.12g tm=%20.12g dnew=%20.12g dmid= %20.12g d_g= %20.12g d_r= %20.12g d_s= %20.12g beta= %20.12g coord= %20.12g\n",  
                                           m_iter, m_dmax, tm, m_dnew, m_dmid, d_g, d_r, d_s, cg_beta, coord_mag);
          if (!m_eMesh->get_rank()) printf("dmax[%3d]=%20.12e tm=%20.12e dnew=%20.12e d_g= %20.12e beta= %20.12e coord= %20.12e scl= %20.12e alp= %20.12e\n",  
                                           m_iter, m_dmax, tm, m_dnew, d_g, cg_beta, coord_mag, m_scale, alpha);
        }

      return tm;
    }
    
    void PMMParallelReferenceMeshSmoother1::debug_print(double alpha)
    {
      if (0)
        {
          //stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
          //stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
          int spatialDim = m_eMesh->get_spatial_dim();
          stk::mesh::FieldBase *cg_d_field    = m_eMesh->get_field("cg_d");
          stk::mesh::FieldBase *cg_g_field    = m_eMesh->get_field("cg_g");
          stk::mesh::FieldBase *cg_r_field    = m_eMesh->get_field("cg_r");
          stk::mesh::FieldBase *cg_s_field    = m_eMesh->get_field("cg_s");

          const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              // update local and globally shared 
              //if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];

                    double *coord_current = PerceptMesh::field_data(m_coord_field_current, node);
                    double *cg_d = PerceptMesh::field_data(cg_d_field, node);
                    double *cg_g = PerceptMesh::field_data(cg_g_field, node);
                    double *cg_r = PerceptMesh::field_data(cg_r_field, node);
                    double *cg_s = PerceptMesh::field_data(cg_s_field, node);

                    bool dopr=false;
                    unsigned nid = node.identifier();
                    if (nid == 71 || nid == 84 || nid == 97)
                      //if (nid == 15)
                      dopr=true;
                    char buf[1024];
                    std::ostringstream ostr;
                    bool owned = (node.owner_rank() == m_eMesh->get_rank());
                    if (dopr)
                      ostr << "P[" << m_eMesh->get_rank() << "] dbp iter= " << m_iter << " nid= " << nid << " o= " << owned;
                    for (int i=0; i < spatialDim; i++)
                      {
                        double dt = alpha*cg_d[i];
                        if (dopr) {
                          sprintf(buf, "dmax i= %d dt= %12g coord= %12g g= %12g r= %12g s= %12g d= %12g", i, dt, coord_current[i], cg_g[i], cg_r[i], cg_s[i], cg_d[i]);
                          ostr << buf;
                        }
                      }
                    if (dopr && owned) std::cout << ostr.str() << std::endl;
                    //if (dopr && m_eMesh->get_rank()==0) std::cout << ostr.str() << std::endl;
                  }
              }
            }
        }

      if (1)
        {
          stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
          stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
          int spatialDim = m_eMesh->get_spatial_dim();
          stk::mesh::FieldBase *cg_d_field    = m_eMesh->get_field("cg_d");
          stk::mesh::FieldBase *cg_g_field    = m_eMesh->get_field("cg_g");
          stk::mesh::FieldBase *cg_r_field    = m_eMesh->get_field("cg_r");
          stk::mesh::FieldBase *cg_s_field    = m_eMesh->get_field("cg_s");

          //if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          unsigned ids[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,71,84,97};
          unsigned idl = sizeof(ids)/sizeof(unsigned);
          idl = 13*13;
          for (unsigned i_node=0; i_node < idl; i_node++)
            {
              MPI_Barrier( MPI_COMM_WORLD );

              //stk::mesh::Entity* node_p = m_eMesh->get_bulk_data()->get_entity(0, ids[i_node]);
              stk::mesh::Entity* node_p = m_eMesh->get_bulk_data()->get_entity(0, i_node+1);
              if (node_p)
                {
                  stk::mesh::Entity& node = *node_p;

                  double *coord_current = PerceptMesh::field_data(m_coord_field_current, node);
                  double *cg_d = PerceptMesh::field_data(cg_d_field, node);
                  double *cg_g = PerceptMesh::field_data(cg_g_field, node);
                  double *cg_r = PerceptMesh::field_data(cg_r_field, node);
                  double *cg_s = PerceptMesh::field_data(cg_s_field, node);

                  bool dopr=true;
                  unsigned nid = node.identifier();
                  char buf[1024];
                  std::ostringstream ostr;
                  bool owned = (node.owner_rank() == m_eMesh->get_rank());
                  if (dopr)
                    ostr << " nid= " << nid << " dmax ";
                  if (0) ostr << " owned= " << owned;
                  //  ostr << "P[" << m_eMesh->get_rank() << "] dbp iter= " << m_iter << " nid= " << nid << " o= " << owned;
                  double ng[3]={0,0,0};
                  bool ng_valid=true;

                  {
                    bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                    bool fixed = m_pmm->get_fixed_flag(&node);
                    if (fixed || isGhostNode)
                      {
                      }
                    else
                      {
                        nodal_gradient(node, 0.0, coord_current, cg_d, ng_valid, ng);
                      }
                  }
                  for (int i=0; i < spatialDim; i++)
                    {
                      double dt = alpha*cg_d[i];
                      if (dopr) {
                        if (0) sprintf(buf, "dt=%12g coord=%12g g=%12g r=%12g s=%12g d=%12g", dt, coord_current[i], cg_g[i], cg_r[i], cg_s[i], cg_d[i]);
                        if (0) sprintf(buf, "coord=%35.18e g=%35.18e ",  coord_current[i], cg_g[i]);
                        sprintf(buf, "g=%20.12e ng=%20.12e diff=%20.12e ",  cg_g[i], ng[i], cg_g[i]-ng[i]);
                        ostr << buf;
                      }
                    }
                  //if (dopr && owned) std::cout << ostr.str() << std::endl;
                  //if (dopr && m_eMesh->get_rank()==0) std::cout << ostr.str() << std::endl;
                  if (dopr && owned) std::cout << ostr.str() << std::endl;
                }
              MPI_Barrier( MPI_COMM_WORLD );

            }
        }        

    }

    void PMMParallelReferenceMeshSmoother1::update_node_positions(Mesh* mesh, double alpha)
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = eMesh->get_spatial_dim();
      stk::mesh::FieldBase *cg_d_field    = eMesh->get_field("cg_d");
      //stk::mesh::FieldBase *cg_g_field    = eMesh->get_field("cg_g");

      m_dmax = 0.0;
      //debug_print(alpha);

      // node loop: update node positions
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            // update local and globally shared 
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    bool fixed = pmm->get_fixed_flag(&node);
                    bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                    if (fixed || isGhostNode)
                      {
                        continue;
                      }

                    double *coord_current = PerceptMesh::field_data(m_coord_field_current, node);
                    double *cg_d = PerceptMesh::field_data(cg_d_field, node);

                    for (int i=0; i < spatialDim; i++)
                      {
                        double dt = alpha*cg_d[i];
                        m_dmax = std::max(std::fabs(dt), m_dmax);
                        coord_current[i] += dt;  
                      }
                  }
              }
          }
      }

      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMax<1>( & m_dmax ) );

      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(m_eMesh->get_coordinates_field());
        //fields.push_back(cg_g_field);

        // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
        stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->shared_aura(), fields); 
        // the shared part (just the shared boundary)
        //stk::mesh::communicate_field_data(*m_eMesh->get_bulk_data()->ghostings()[0], fields);
      }


    }

    double PMMParallelReferenceMeshSmoother1::total_metric(Mesh *mesh, double alpha, double multiplicative_edge_scaling, bool& valid)
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      stk::mesh::FieldBase *coord_field = m_eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field;
      stk::mesh::FieldBase *coord_field_lagged  = m_eMesh->get_field("coordinates_lagged");

      stk::mesh::FieldBase *cg_d_field    = m_eMesh->get_field("cg_d");

      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = m_eMesh->get_spatial_dim();

      double mtot = 0.0;

      // cache coordinates
      m_eMesh->copy_field(coord_field_lagged, coord_field_current);

      // node loop
      {
        const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            // update local and globally shared 
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    bool fixed = pmm->get_fixed_flag(&node);
                    if (fixed)
                      {
                        continue;
                      }

                    double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                    double *cg_d = PerceptMesh::field_data(cg_d_field, node);

                    for (int i=0; i < spatialDim; i++)
                      {
                        double dt = alpha * cg_d[i];
                        coord_current[i] += dt;
                      }
                  }
              }
          }
      }

      valid = true;
      // element loop
      {
        const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (PerceptMesquiteMesh::select_bucket(**k, m_eMesh) && on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                m_metric->m_topology_data = m_eMesh->get_cell_topology(bucket);
            
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];
                    bool local_valid=true;
                    double mm = metric(element, local_valid);
                    valid = valid && local_valid;
                    if (do_print_elem_val) PRINT( "element= " << element.identifier() << " metric= " << mm );
                    if (do_print_elem_val && element.identifier() == 13) { std::cout << element.identifier() << " iter= " << m_iter << " element= " << element.identifier() << " metric= " << mm << std::endl;}
                    mtot += mm;
                  }
              }
          }
      }

      // reset coordinates
      m_eMesh->copy_field(coord_field_current, coord_field_lagged);

      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceSum<1>( & mtot ) );
      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMin<1>( & valid ) );

      return mtot;
      
    }
  }
}


#endif
#endif
