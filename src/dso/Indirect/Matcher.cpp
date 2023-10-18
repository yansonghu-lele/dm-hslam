#include "Matcher.h"
#include "Indirect/Frame.h"
#include "util/FrameShell.h"
#include "Indirect/MapPoint.h"
#include "FullSystem/HessianBlocks.h"

#include <cmath>

namespace dso
{
    using namespace std;

    const int Matcher::TH_HIGH = 100;
    const int Matcher::TH_LOW = 50;
    const int Matcher::HISTO_LENGTH = 30;

    int Matcher::SearchLocalMapByProjection(shared_ptr<Frame> F, vector<shared_ptr<MapPoint>> &vpMapPoints, float th, float nnratio) //this is run from tracking thread -> access tMapPoints only!!
    {
        int nmatches = 0;

        const bool bFactor = th != 1.0;

        for (size_t iMP = 0, iend = vpMapPoints.size(); iMP < iend; ++iMP)
        {
            shared_ptr<MapPoint> pMP = vpMapPoints[iMP];
            
            if (!pMP->mbTrackInView)
                continue;

            if (pMP->isBad())
                continue;

            // const int &nPredictedLevel = pMP->mnTrackScaleLevel;

            // The size of the window will depend on the viewing direction
            float r = RadiusByViewingCos(pMP->mTrackViewCos);

            if (bFactor)
                r *= th;
            
            const vector<size_t> vIndices = F->GetFeaturesInArea(pMP->mTrackProjX, pMP->mTrackProjY, r);


            if (vIndices.empty())
                continue;

            const cv::Mat MPdescriptor = pMP->GetDescriptor();

            int bestDist = 256;
            int bestLevel = -1;
            int bestDist2 = 256;
            int bestLevel2 = -1;
            int bestIdx = -1;

            // Get best and second matches with near keypoints
            for (vector<size_t>::const_iterator vit = vIndices.begin(), vend = vIndices.end(); vit != vend; vit++)
            {
                const size_t idx = *vit;

                if (F->tMapPoints[idx])
                    if (F->tMapPoints[idx]->getNObservations() > 0)
                        continue;


                const cv::Mat &d = F->Descriptors.row(idx);

                const int dist = DescriptorDistance(MPdescriptor, d);

                if (dist < bestDist)
                {
                    bestDist2 = bestDist;
                    bestDist = dist;
                    // bestLevel2 = bestLevel;
                    // bestLevel = F.mvKeysUn[idx].octave;
                    bestIdx = idx;
                }
                else if (dist < bestDist2)
                {
                    // bestLevel2 = F.mvKeysUn[idx].octave;
                    bestDist2 = dist;
                }
            }

            // Apply ratio to second match (only if best and second are in the same scale level)
            if (bestDist <= TH_HIGH)
            {
                if (bestLevel == bestLevel2 && bestDist > nnratio * bestDist2)
                    continue;

                F->tMapPoints[bestIdx] = pMP;
                nmatches++;
            }
        }

        return nmatches;
    }

    int Matcher::Fuse(std::shared_ptr<Frame> pKF, const Sim3& Scw, const std::vector<std::shared_ptr<MapPoint>> &vpPoints, float th, std::vector<std::shared_ptr<MapPoint>> &vpReplacePoint) //this is run from mapping thread or loop closure thread -> access mvpMapPoints
    {
        // Get Calibration Parameters for later projection
        const float &fx = pKF->HCalib->fxl();
        const float &fy = pKF->HCalib->fyl();
        const float &cx = pKF->HCalib->cxl();
        const float &cy = pKF->HCalib->cyl();

        // Decompose Scw
        Mat33f sRcw = Scw.rotationMatrix().cast<float>();
        // float scw = sRcw.row(0).dot(sRcw.row(0));
        // if (scw != 1.0)
        //     sRcw = sRcw / scw;
        Vec3f tcw = Scw.translation().cast<float>(); // / Scw.scale();
        Vec3f Ow = -sRcw.transpose() * tcw;

        // cv::Mat sRcw = Scw.rowRange(0, 3).colRange(0, 3);
        // const float scw = sqrt(sRcw.row(0).dot(sRcw.row(0)));
        // cv::Mat Rcw = sRcw / scw;
        // cv::Mat tcw = Scw.rowRange(0, 3).col(3) / scw;
        // cv::Mat Ow = -Rcw.t() * tcw;

        // Set of MapPoints already found in the KeyFrame
        const std::set<std::shared_ptr<MapPoint>> spAlreadyFound = pKF->getMapPointsS();

        int nFused = 0;

        const int nPoints = vpPoints.size();

        // For each candidate MapPoint project and match
        for (int iMP = 0; iMP < nPoints; iMP++)
        {
            std::shared_ptr<MapPoint> pMP = vpPoints[iMP];

            // Discard Bad MapPoints and already found
            if (pMP->isBad() || spAlreadyFound.count(pMP))
                continue;

            // Get 3D Coords.
            Vec3f p3Dw = pMP->getWorldPose();

            // Transform into Camera Coords.
            Vec3f p3Dc = sRcw * p3Dw + tcw;

            // Depth must be positive
            if (p3Dc(2) < 0.0f)
                continue;

            // Project into Image
            const float invz = 1.0 / p3Dc(2);
            const float x = p3Dc(0) * invz;
            const float y = p3Dc(1) * invz;

            const float u = fx * x + cx;
            const float v = fy * y + cy;

            // Point must be inside the image
            if (u < mnMinX || u > mnMaxX || v < mnMinY || v > mnMaxY)
                continue;

            // Depth must be inside the scale pyramid of the image
            // const float maxDistance = pMP->GetMaxDistanceInvariance();
            // const float minDistance = pMP->GetMinDistanceInvariance();
            Vec3f PO = p3Dw - Ow;
            const float dist3D = PO.norm();

            // if (dist3D < minDistance || dist3D > maxDistance)
            //     continue;

            // Viewing angle must be less than 60 deg
            Vec3f Pn = pMP->GetNormal();

            if (PO.dot(Pn) < 0.5 * dist3D)
                continue;

            // Compute predicted scale level
            // const int nPredictedLevel = pMP->PredictScale(dist3D, pKF);

            // Search in a radius
            // const float radius = th * pKF->mvScaleFactors[nPredictedLevel];
            const float radius = th * 1.0f;

            const vector<size_t> vIndices = pKF->GetFeaturesInArea(u, v, radius);

            if (vIndices.empty())
                continue;

            // Match to the most similar keypoint in the radius

            const cv::Mat dMP = pMP->GetDescriptor();

            int bestDist = INT_MAX;
            int bestIdx = -1;
            for (vector<size_t>::const_iterator vit = vIndices.begin(); vit != vIndices.end(); vit++)
            {
                const size_t idx = *vit;
                // const int &kpLevel = pKF->mvKeysUn[idx].octave;

                // if (kpLevel < nPredictedLevel - 1 || kpLevel > nPredictedLevel)
                //     continue;

                const cv::Mat &dKF = pKF->Descriptors.row(idx);

                int dist = DescriptorDistance(dMP, dKF);

                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestIdx = idx;
                }
            }

            // If there is already a MapPoint replace otherwise add new measurement
            if (bestDist <= TH_LOW)
            {
                std::shared_ptr<MapPoint> pMPinKF = pKF->getMapPoint(bestIdx);
                if (pMPinKF)
                {
                    if (!pMPinKF->isBad())
                        vpReplacePoint[iMP] = pMPinKF;
                }
                else
                {
                    pMP->AddObservation(pKF, bestIdx);
                    pKF->addMapPointMatch(pMP, bestIdx);
                }
                nFused++;
            }
        }

        return nFused;
    }

    int Matcher::Fuse(std::shared_ptr<Frame> pKF, const std::vector<std::shared_ptr<MapPoint>> &vpMapPoints, const float th) //this is run from mapping thread or loop closure thread -> access mvpMapPoints
    {
        Mat33f Rcw = pKF->fs->getPoseInverse().rotationMatrix().cast<float>();

        Vec3f tcw = pKF->fs->getPoseInverse().translation().cast<float>();


        const float &fx = pKF->HCalib->fxl();
        const float &fy = pKF->HCalib->fyl();
        const float &cx = pKF->HCalib->cxl();
        const float &cy = pKF->HCalib->cyl();
        // const float &bf = pKF->mbf;

        Vec3f Ow = pKF->fs->getCameraCenter().cast<float>();

        int nFused = 0;

        const int nMPs = vpMapPoints.size();

        for (int i = 0; i < nMPs; i++)
        {
            std::shared_ptr<MapPoint> pMP = vpMapPoints[i];

            if (!pMP)
                continue;

            if (pMP->isBad() || pMP->isInKeyframe(pKF))
                continue;

            Vec3f p3Dw = pMP->getWorldPose();
            Vec3f p3Dc = Rcw * p3Dw + tcw;

            // Depth must be positive
            if (p3Dc(2) < 0.0f)
                continue;

            const float invz = 1 / p3Dc(2);
            const float x = p3Dc(0) * invz;
            const float y = p3Dc(1) * invz;

            const float u = fx * x + cx;
            const float v = fy * y + cy;

            // Point must be inside the image
            if (u < mnMinX || u > mnMaxX || v < mnMinY || v > mnMaxY)
                continue;

            // const float maxDistance = pMP->GetMaxDistanceInvariance();
            // const float minDistance = pMP->GetMinDistanceInvariance();
            Vec3f PO = p3Dw - Ow;
            const float dist3D = PO.norm();

            // Depth must be inside the scale pyramid of the image
            // if (dist3D < minDistance || dist3D > maxDistance)
            //     continue;

            // Viewing angle must be less than 60 deg
            Vec3f Pn = pMP->GetNormal();

            if (PO.dot(Pn) < 0.5 * dist3D)
                continue;

            // int nPredictedLevel = pMP->PredictScale(dist3D, pKF);

            // Search in a radius
            // const float radius = th * pKF->mvScaleFactors[nPredictedLevel];
            const float radius = th * 1.0f;

            const vector<size_t> vIndices = pKF->GetFeaturesInArea(u, v, radius);

            if (vIndices.empty())
                continue;

            // Match to the most similar keypoint in the radius

            const cv::Mat dMP = pMP->GetDescriptor();

            int bestDist = 256;
            int bestIdx = -1;
            for (vector<size_t>::const_iterator vit = vIndices.begin(), vend = vIndices.end(); vit != vend; vit++)
            {
                const size_t idx = *vit;

                const cv::KeyPoint &kp = pKF->mvKeys[idx];

                const int &kpLevel = kp.octave;

                // if (kpLevel < nPredictedLevel - 1 || kpLevel > nPredictedLevel)
                //     continue;

                const float &kpx = kp.pt.x;
                const float &kpy = kp.pt.y;
                const float ex = u - kpx;
                const float ey = v - kpy;
                const float e2 = ex * ex + ey * ey;

                // if (e2 * pKF->mvInvLevelSigma2[kpLevel] > 5.99)
                //     continue;

                if (e2 * 1.0f > 5.99)
                    continue;

                const cv::Mat &dKF = pKF->Descriptors.row(idx);

                const int dist = DescriptorDistance(dMP, dKF);

                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestIdx = idx;
                }
            }

            // If there is already a MapPoint replace otherwise add new measurement
            if (bestDist <= TH_LOW)
            {
                std::shared_ptr<MapPoint> pMPinKF = pKF->getMapPoint(bestIdx);
                if (pMPinKF)
                {
                    if (!pMPinKF->isBad())
                    {

                        // pMP->Replace(pMPinKF); //enable this and comment the if-else to keep the old Mps only
                        if (pMPinKF->getNObservations() > pMP->getNObservations())
                            pMP->Replace(pMPinKF);
                        else
                            pMPinKF->Replace(pMP);
                    }
                }
                else
                {
                    pMP->AddObservation(pKF, bestIdx);
                    pKF->addMapPointMatch(pMP, bestIdx);
                }
                nFused++;
            }
        }

        return nFused;
    }

    int Matcher::SearchByProjectionFrameToFrame(std::shared_ptr<Frame> CurrentFrame, const std::shared_ptr<Frame> LastFrame, const float th, bool mbCheckOrientation) //this is run from tracking thread -> access tMapPoints only!!
    {
        int nmatches = 0;

        // Rotation Histogram (to check rotation consistency)
        vector<int> rotHist[HISTO_LENGTH];
        for (int i = 0; i < HISTO_LENGTH; i++)
            rotHist[i].reserve(500);
        const float factor = 1.0f / HISTO_LENGTH;

        const Mat33f Rcw = CurrentFrame->fs->getPoseInverse().rotationMatrix().cast<float>(); //currFrame not added to map yet (does not have poseOpti!!)  CurrentFrame.mTcw.rowRange(0, 3).colRange(0, 3);
        const Vec3f tcw = CurrentFrame->fs->getPoseInverse().translation().cast<float>();     //CurrentFrame.mTcw.rowRange(0, 3).col(3);

        const Vec3f twc = -Rcw.transpose() * tcw;

        const Mat33f Rlw = LastFrame->fs->getPoseInverse().rotationMatrix().cast<float>(); //mTcw.rowRange(0, 3).colRange(0, 3);
        const Vec3f tlw = LastFrame->fs->getPoseInverse().translation().cast<float>();       //mTcw.rowRange(0, 3).col(3);

   
        for (int i = 0; i < LastFrame->nFeatures; ++i)
        {
            std::shared_ptr<MapPoint> pMP = LastFrame->tMapPoints[i];

            if (pMP)
            {
                if (!LastFrame->mvbOutlier[i])
                {
                    // Project
                    Vec3f x3Dw = pMP->getWorldPose();
                    Vec3f x3Dc = Rcw * x3Dw + tcw;

                    const float xc = x3Dc(0);
                    const float yc = x3Dc(1);
                    const float invzc = 1.0 / x3Dc(2);

                    if (invzc < 0)
                        continue;

                    float u = CurrentFrame->HCalib->fxl() * xc * invzc + CurrentFrame->HCalib->cxl();
                    float v = CurrentFrame->HCalib->fyl() * yc * invzc + CurrentFrame->HCalib->cyl();

                    if (u < mnMinX || u > mnMaxX)
                        continue;
                    if (v < mnMinY || v > mnMaxY)
                        continue;

                    // int nLastOctave = LastFrame.mvKeys[i].octave;

                    // Search in a window. Size depends on scale
                    // float radius = th * CurrentFrame.mvScaleFactors[nLastOctave];
                    float radius = th * 1.0f;


                    vector<size_t> vIndices2;

                    // if (bForward)
                    //     vIndices2 = CurrentFrame.GetFeaturesInArea(u, v, radius, nLastOctave);
                    // else if (bBackward)
                    //     vIndices2 = CurrentFrame.GetFeaturesInArea(u, v, radius, 0, nLastOctave);
                    // else
                    vIndices2 = CurrentFrame->GetFeaturesInArea(u, v, radius);

                    if (vIndices2.empty())
                        continue;

                    const cv::Mat dMP = pMP->GetDescriptor();

                    int bestDist = 256;
                    int bestIdx2 = -1;

                    for (vector<size_t>::const_iterator vit = vIndices2.begin(), vend = vIndices2.end(); vit != vend; vit++)
                    {
                        const size_t i2 = *vit;
                        if (CurrentFrame->tMapPoints[i2])
                            if (CurrentFrame->tMapPoints[i2]->getNObservations() > 0)
                                continue;

                        const cv::Mat &d = CurrentFrame->Descriptors.row(i2);

                        const int dist = DescriptorDistance(dMP, d);

                        if (dist < bestDist)
                        {
                            bestDist = dist;
                            bestIdx2 = i2;
                        }
                    }

                    if (bestDist <= TH_HIGH)
                    {
                        CurrentFrame->tMapPoints[bestIdx2] = pMP;
                        nmatches++;

                        if (mbCheckOrientation)
                        {
                            float rot = LastFrame->mvKeys[i].angle - CurrentFrame->mvKeys[bestIdx2].angle;
                            if (rot < 0.0)
                                rot += 360.0f;
                            int bin = round(rot * factor);
                            if (bin == HISTO_LENGTH)
                                bin = 0;
                            assert(bin >= 0 && bin < HISTO_LENGTH);
                            rotHist[bin].push_back(bestIdx2);
                        }
                    }
                }
            }
        }

        //Apply rotation consistency
        if (mbCheckOrientation)
        {
            int ind1 = -1;
            int ind2 = -1;
            int ind3 = -1;

            ComputeThreeMaxima(rotHist, HISTO_LENGTH, ind1, ind2, ind3);

            for (int i = 0; i < HISTO_LENGTH; i++)
            {
                if (i != ind1 && i != ind2 && i != ind3)
                {
                    for (size_t j = 0, jend = rotHist[i].size(); j < jend; j++)
                    {
                        CurrentFrame->tMapPoints[rotHist[i][j]].reset();
                        nmatches--;
                    }
                }
            }
        }

        return nmatches;
    }

} // namespace dso
