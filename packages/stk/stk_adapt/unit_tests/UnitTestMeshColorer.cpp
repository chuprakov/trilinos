/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/FunctionOperator.hpp>
#include <stk_percept/function/ConstantFunction.hpp>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_adapt/Colorer.hpp>


namespace stk
{
  namespace adapt
  {
    namespace unit_tests
    {
      static int dw_enabled = 1;
      static stk::diag::Writer s_diagWriter(std::cout.rdbuf(), dw_enabled);

      static stk::diag::Writer &
      dw()
      {
        //static stk::diag::Writer s_diagWriter(dwout().rdbuf(), 0);

        s_diagWriter.setPrintMask(percept::LOG_NORM+percept::LOG_ALWAYS);

        return s_diagWriter;
      }
      //static stk::diag::Writer &s_dw_tmp = dw();

#define EXTRA_PRINT 0

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      TEST(mesh_colorer, test1)
      {
        dw().m(percept::LOG_MESH_COLORER) << "TEST::mesh_colorer::test1 " << stk::diag::dendl;

        const size_t numxyz=3;
        const size_t num_x = numxyz;
        const size_t num_y = numxyz;
        const size_t num_z = numxyz;
        std::string config_mesh = 
          Ioss::Utils::to_string(num_x) + "x" + 
          Ioss::Utils::to_string(num_y) + "x" +
          Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
        percept::PerceptMesh eMesh;
        eMesh.newMesh(percept::PerceptMesh::GMeshSpec(config_mesh));
        int vectorDimension = 0;
        FieldBase *element_color_field = eMesh.addField("element_colors", mesh::Element, vectorDimension);
        eMesh.commit();

        Colorer meshColorer;
        unsigned elementType = 0u;
        meshColorer.color(eMesh, &elementType, element_color_field);
        //eMesh.saveAs("./cube_colored.e");

        //std::cout << "Mesh coloring info: " << meshColorer.getElementColors() << std::endl;
        
      }
    }
  }
}
