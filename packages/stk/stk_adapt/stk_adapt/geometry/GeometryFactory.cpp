#include "GeometryFactory.hpp"
#include "PerceptMesh.hpp"
#include "Util.hpp"

using namespace stk;
using namespace mesh;
using namespace percept;

GeometryFactory::GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry)
{
    geomKernel = kernel;
    geomDatabase = geometry;
}

GeometryFactory::~GeometryFactory()
{

}

bool GeometryFactory::read_file(const std::string& filename, PerceptMesh* mesh_data)
{
    std::vector<GeometryHandle> geometry_entities;
    if (!geomKernel->read_file(filename, geometry_entities))
        return false;
    for (size_t i=0; i<geometry_entities.size(); i++)
    {
        std::string str = geomKernel->get_attribute(geometry_entities[i]);
        Part* part = mesh_data->getNonConstPart(str);
        if (part)
            geomDatabase->add_evaluator(new GeometryEvaluator(part));
    }
    return true;
}