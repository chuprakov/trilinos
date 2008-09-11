#ifndef RELAXATION_SMOOTHER_HPP
#define RELAXATION_SMOOTHER_HPP

#include "Mesquite.hpp"
#include "VertexPatches.hpp"
#include "VertexMover.hpp"

namespace Mesquite
{
  /**\brief Base class for LaPlacian and other relaxation smoothers */
  class RelaxationSmoother : public VertexMover
  {
  public:
    /**
     * \param OF  For many relaxation solvers (e.g. Laplacian)
     * this ObjectiveFunction is used only to evaluate user-specified
     * termination criteria that require an objective function.  
     */
    RelaxationSmoother( ObjectiveFunction* OF = NULL ) : VertexMover(OF) {}
    
    virtual ~RelaxationSmoother();

    PatchSet* get_patch_set() { return &patchSet; }
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);

    virtual void optimize_vertex_positions( PatchData &pd,
                                            MsqError &err) = 0;
					 
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);

    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);

    virtual void cleanup();
    
  private:
    VertexPatches patchSet;
  };
  
}

#endif
