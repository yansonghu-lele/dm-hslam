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

    bool PoseOptimization(std::shared_ptr<Frame> pFrame, CalibHessian *calib, bool updatePose)
    {
        g2o::SparseOptimizer optimizer;
        auto linearSolver = g2o::make_unique<g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>>();
        g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<g2o::BlockSolver_6_3>(std::move(linearSolver)));
        optimizer.setAlgorithm(solver);

        int nInitialCorrespondences = 0;

        // Set Frame vertex
        vertexSE3 *vSE3 = new vertexSE3(); //VertexSE3Expmap Converter::toSE3Quat(pFrame->mTcw));
        vSE3->setEstimate(pFrame->fs->getPoseInverse());
        vSE3->setId(0);
        vSE3->setFixed(false);
        optimizer.addVertex(vSE3);
        // Set MapPoint vertices
        const int N = pFrame->nFeatures;

        vector<edgeSE3XYZPoseOnly *> vpEdgesMono;
        vector<size_t> vnIndexEdgeMono;
        vector<double> initErr;
        double initScale = 1.0;
        initErr.reserve(2 * N);
        vpEdgesMono.reserve(N);
        vnIndexEdgeMono.reserve(N);
        vector<double> vInformation;
        vInformation.reserve(N);
        double maxInfo = 1e7;
        double stdDev = 1e7;

        const float deltaMono = sqrt(1.345); //sqrt(1.345); // sqrt(5.991);
        {
            // boost::unique_lock<boost::mutex> lock(MapPoint::mGlobalMutex); //this would lock ALL map points poses from changing!

            for (int i = 0; i < N; i++)
            {
                std::shared_ptr<MapPoint> pMP = pFrame->tMapPoints[i];
                if (pMP)
                {
                    nInitialCorrespondences++;
                    pFrame->mvbOutlier[i] = false;

                    edgeSE3XYZPoseOnly *e = new edgeSE3XYZPoseOnly();
                    e->setCamera(calib->fxl(), calib->fyl(), calib->cxl(), calib->cyl());
                    e->setXYZ(pMP->getWorldPose().cast<double>());
                    e->setInformation(Eigen::Matrix2d::Identity());
                    e->setMeasurement(Vec2((double)pFrame->mvKeys[i].pt.x, (double)pFrame->mvKeys[i].pt.y));
                    e->setId(i);

                    e->setVertex(0, vSE3);

                    e->computeError();
                    initErr.push_back(e->error()[0]);
                    initErr.push_back(e->error()[1]);

                    double Info = sqrtf(pMP->getidepthHessian());
                    vInformation.push_back(Info);
                    if (Info > maxInfo)
                        maxInfo = Info;

                    g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);

                    rk->setDelta(deltaMono);

                    optimizer.addEdge(e);

                    vpEdgesMono.push_back(e);
                    vnIndexEdgeMono.push_back(i);
                }
            }

            //compute information vector distribution:
            stdDev = getStdDev(vInformation);

            initScale = getStdDev(initErr); //computeScale(initErr);

            for (int i = 0, iend = vpEdgesMono.size(); i < iend; ++i)
            {
                vpEdgesMono[i]->setScale(initScale);

                vpEdgesMono[i]->setInformation(Eigen::Matrix2d::Identity() * (vInformation[i] / (stdDev + 0.00001)));  //* invSigma2); //set this to take into account depth variance!
            }
        }
        if (nInitialCorrespondences < 10 || optimizer.edges().size() < 10)
            return false;

        // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
        // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
        // const float chi2Mono[4] = {5.991, 5.991, 5.991, 5.991}; //1.345
        const float chi2Mono[4] = {1.345, 1.345, 1.345, 1.345}; //5.991
        const int its[4] = {10, 10, 10, 10};

        int nBad = 0;
        for (size_t it = 0; it < 4; it++)
        {
            vSE3->setEstimate(pFrame->fs->getPoseInverse());
            optimizer.initializeOptimization(0);
            optimizer.optimize(its[it]);
            nBad = 0;
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++)
            {
                edgeSE3XYZPoseOnly *e = vpEdgesMono[i];

                const size_t idx = vnIndexEdgeMono[i];

                if (pFrame->mvbOutlier[idx])
                {
                    e->computeError();
                }

                const float chi2 = e->chi2();

                if (chi2 > chi2Mono[it])
                {
                    pFrame->mvbOutlier[idx] = true;
                    e->setLevel(1);
                    nBad++;
                }
                else
                {
                    pFrame->mvbOutlier[idx] = false;
                    e->setLevel(0);
                }

                if (it == 2)
                    e->setRobustKernel(0);
            }

            if (optimizer.edges().size() < 10)
                break;

            if (!updatePose)
                break;
        }

        number_t *hessianData = vSE3->hessianData();
        Vec6 vHessian;
        vHessian<< hessianData[0], hessianData[7], hessianData[14], hessianData[21], hessianData[28], hessianData[35];

        bool isUsable = vHessian.norm() > 1e6;
        // // Recover optimized pose
        if(isUsable && updatePose)
            pFrame->fs->setPose(vSE3->estimate().inverse());
        // bool isUsable = false;
        return isUsable;
    }

    int checkOutliers(std::shared_ptr<Frame> pFrame, CalibHessian* calib)
    {
        g2o::SparseOptimizer optimizer;
        auto linearSolver = g2o::make_unique<g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>>();
        g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<g2o::BlockSolver_6_3>(std::move(linearSolver)));
        optimizer.setAlgorithm(solver);

        vertexSE3 *vSE3 = new vertexSE3(); //VertexSE3Expmap Converter::toSE3Quat(pFrame->mTcw));
        vSE3->setEstimate(pFrame->fs->getPoseInverse());

        vSE3->setId(0);
        vSE3->setFixed(false);
        optimizer.addVertex(vSE3);

        const int N = pFrame->nFeatures;
        vector<edgeSE3XYZPoseOnly *> vpEdgesMono;
        vector<size_t> vnIndexEdgeMono;
        vector<double> initErr;
        double initScale = 1.0;
        initErr.reserve(2 * N);
        vpEdgesMono.reserve(N);
        vnIndexEdgeMono.reserve(N);
        vector<double> vInformation;
        vInformation.reserve(N);
        double stdDev = 1e7;

        for (int i = 0; i < N; i++)
        {
            std::shared_ptr<MapPoint> pMP = pFrame->tMapPoints[i];
            if (pMP)
            {
                pFrame->mvbOutlier[i] = false;

                edgeSE3XYZPoseOnly *e = new edgeSE3XYZPoseOnly();
                e->setCamera(calib->fxl(), calib->fyl(), calib->cxl(), calib->cyl());
                e->setXYZ(pMP->getWorldPose().cast<double>());
                e->setInformation(Eigen::Matrix2d::Identity());
                e->setMeasurement(Vec2((double)pFrame->mvKeys[i].pt.x, (double)pFrame->mvKeys[i].pt.y));
                e->setId(i);

                e->setVertex(0, vSE3);

                e->computeError();
                initErr.push_back(e->error()[0]);
                initErr.push_back(e->error()[1]);

                double Info = pMP->getidepthHessian();
                vInformation.push_back(Info);
                optimizer.addEdge(e);
                vpEdgesMono.push_back(e);
                vnIndexEdgeMono.push_back(i);
            }
        }

        stdDev = getStdDev(vInformation);
        initScale = getStdDev(initErr);
        for (int i = 0, iend = vpEdgesMono.size(); i < iend; ++i)
        {
            // vpEdgesMono[i]->setScale(initScale);
            vpEdgesMono[i]->setInformation(Eigen::Matrix2d::Identity() * (vInformation[i] / (stdDev + 0.00001))); //* invSigma2); //set this to take into account depth variance!
        }

        optimizer.initializeOptimization(0);
        int nBad = 0;
        for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++)
        {
            edgeSE3XYZPoseOnly *e = vpEdgesMono[i];

            const size_t idx = vnIndexEdgeMono[i];

            if (pFrame->mvbOutlier[idx])
            {
                e->computeError();
            }

            const float chi2 = e->chi2();

            if (chi2 > 5.991)
            {
                pFrame->mvbOutlier[idx] = true;
                e->setLevel(1);
                nBad++;
            }
            else
            {
                pFrame->mvbOutlier[idx] = false;
                e->setLevel(0);
            }
        }

        return nBad;
    }

} // namespace dso