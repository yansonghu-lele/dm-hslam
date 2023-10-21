#include "Indirect/Optimizer.h"
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/robust_kernel_impl.h>

#include <g2o/solvers/dense/linear_solver_dense.h>
// #include "g2o/solvers/linear_solver_dense.h"

#include <g2o/solvers/eigen/linear_solver_eigen.h>
// #include <g2o/solvers/linear_solver_eigen.h>

#include <boost/thread.hpp>
#include "Indirect/Frame.h"
#include "util/FrameShell.h"
#include "Indirect/MapPoint.h"
#include "Indirect/Map.h"
#include "FullSystem/FullSystem.h"


#include "Indirect/Sim3_impl.h"
#include "Indirect/Se3.h"
#include <g2o/types/sba/types_six_dof_expmap.h>

#include "g2o/solvers/cholmod/linear_solver_cholmod.h"

// TODO: Replace g2o with GTSAM

namespace dso
{
    using namespace std;
    using namespace OptimizationStructs;

    dso::Sim3 g2oSim3_to_sophusSim3(Sim3Vertex &g2o_sim3)
    {
        Mat44 sim_transform;
        sim_transform.topLeftCorner(3, 3) = g2o_sim3.estimate().rotation().toRotationMatrix();
        sim_transform.topRightCorner(3, 1) = g2o_sim3.estimate().translation();
        sim_transform.topLeftCorner(3, 3) *= g2o_sim3.estimate().scale();
        return Sim3(sim_transform);
    }

    g2o::Sim3 sophusSim3_to_g2oSim3(dso::Sim3 sophus_sim3)
    {
        return g2o::Sim3(sophus_sim3.rotationMatrix(), sophus_sim3.translation(), sophus_sim3.scale());
    }

} // namespace dso