/**
* This file is part of DSO, written by Jakob Engel.
* It has been modified by Lukas von Stumberg for the inclusion in DM-VIO (http://vision.in.tum.de/dm-vio).
*
* Copyright 2022 Lukas von Stumberg <lukas dot stumberg at tum dot de>
* Copyright 2016 Technical University of Munich and Intel.
* Developed by Jakob Engel <engelj at in dot tum dot de>,
* for more information see <http://vision.in.tum.de/dso>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* DSO is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DSO. If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * KFBuffer.cpp
 *
 *  Created on: Jan 7, 2014
 *      Author: engelj
 */



#include "FullSystem/FullSystem.h"

#include "stdio.h"
#include "util/globalFuncs.h"
#include <Eigen/LU>
#include <algorithm>
#include "IOWrapper/ImageDisplay.h"
#include "util/globalCalib.h"
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include "FullSystem/PixelSelector.h"
#include "FullSystem/PixelSelector2.h"
#include "FullSystem/ResidualProjections.h"
#include "FullSystem/ImmaturePoint.h"

#include "FullSystem/CoarseTracker.h"
#include "FullSystem/CoarseInitializer.h"

#include "OptimizationBackend/EnergyFunctional.h"
#include "OptimizationBackend/EnergyFunctionalStructs.h"

#include "IOWrapper/Output3DWrapper.h"
#include "util/ImageAndExposure.h"
#include <cmath>

#include "util/TimeMeasurement.h"
#include "GTSAMIntegration/ExtUtils.h"

#include <iterator>


using dmvio::GravityInitializer;

namespace dso
{

int FrameHessian::instanceCounter=0;
int PointHessian::instanceCounter=0;
unsigned long PointHessian::totalInstantCounter=0;
int CalibHessian::instanceCounter=0;

boost::mutex FrameShell::shellPoseMutex{};

/**
 * @brief Set variables to defaults
 * 
 */
void FullSystem::setDefaults()
{
	int retstat = 0;
	if(globalSettings.setting_logStuff)
	{
		retstat += system("rm -rf logs");
		retstat += system("mkdir logs");

		retstat += system("rm -rf mats");
		retstat += system("mkdir mats");

		calibLog = new std::ofstream();
		calibLog->open("logs/calibLog.txt", std::ios::trunc | std::ios::out);
		calibLog->precision(12);

		numsLog = new std::ofstream();
		numsLog->open("logs/numsLog.txt", std::ios::trunc | std::ios::out);
		numsLog->precision(10);

		coarseTrackingLog = new std::ofstream();
		coarseTrackingLog->open("logs/coarseTrackingLog.txt", std::ios::trunc | std::ios::out);
		coarseTrackingLog->precision(10);

		eigenAllLog = new std::ofstream();
		eigenAllLog->open("logs/eigenAllLog.txt", std::ios::trunc | std::ios::out);
		eigenAllLog->precision(10);

		eigenPLog = new std::ofstream();
		eigenPLog->open("logs/eigenPLog.txt", std::ios::trunc | std::ios::out);
		eigenPLog->precision(10);

		eigenALog = new std::ofstream();
		eigenALog->open("logs/eigenALog.txt", std::ios::trunc | std::ios::out);
		eigenALog->precision(10);

		DiagonalLog = new std::ofstream();
		DiagonalLog->open("logs/diagonal.txt", std::ios::trunc | std::ios::out);
		DiagonalLog->precision(10);

		variancesLog = new std::ofstream();
		variancesLog->open("logs/variancesLog.txt", std::ios::trunc | std::ios::out);
		variancesLog->precision(10);


		nullspacesLog = new std::ofstream();
		nullspacesLog->open("logs/nullspacesLog.txt", std::ios::trunc | std::ios::out);
		nullspacesLog->precision(10);
	}
	else
	{
		nullspacesLog=0;
		variancesLog=0;
		DiagonalLog=0;
		eigenALog=0;
		eigenPLog=0;
		eigenAllLog=0;
		numsLog=0;
		calibLog=0;
	}

	assert(retstat!=293847);

	statistics_lastNumOptIts=0;
	statistics_numDroppedPoints=0;
	statistics_numActivatedPoints=0;
	statistics_numCreatedPoints=0;
	statistics_numForceDroppedResBwd = 0;
	statistics_numForceDroppedResFwd = 0;
	statistics_numMargResFwd = 0;
	statistics_numMargResBwd = 0;
	statistics_lastFineTrackRMSE = 0.0f;

	lastCoarseRMSE.setConstant(100);

	currentMinActDist=2;

	initialized=false;
	isLost=false;
	initFailed=false;
	imuUsedBefore = false;

	needNewKFAfter = -1;
	runMapping=true;
	needToKetchupMapping = false;
	secondKeyframeDone = false;
	lastRefStopID=0;


	minIdJetVisDebug = -1;
	maxIdJetVisDebug = -1;
	minIdJetVisTracker = -1;
	maxIdJetVisTracker = -1;

	firstPose = Sophus::SE3d{};

	framesBetweenKFsRest = 0.0f;
}

/**
 * @brief Init the internal classes
 * 
 */
void FullSystem::setClasses()
{
	coarseDistanceMap = std::make_unique<CoarseDistanceMap> (globalCalib, globalSettings.pyrLevelsUsed);
	coarseTracker = new CoarseTracker(globalCalib, imuIntegration, globalSettings);
	coarseTracker_forNewKF = new CoarseTracker(globalCalib, imuIntegration, globalSettings);
	coarseInitializer = std::make_unique<CoarseInitializer> (globalCalib, globalSettings);
	pixelSelector = std::make_unique<PixelSelector> (globalCalib, globalSettings.setting_minGradHistCut, globalSettings.setting_minGradHistAdd, globalSettings);
}

/**
 * @brief Construct a new FullSystem Object
 * 
 * @param linearizeOperationPassed 
 * @param imuCalibration 
 * @param imuSettings 
 */
FullSystem::FullSystem(bool linearizeOperationPassed, 
						dso::Global_Calib& globalCalib_,
						dmvio::IMUCalibration& imuCalibration,
                       				dmvio::IMUSettings& imuSettings,
					   	GlobalSettings& globalSettings_)
    : globalSettings(globalSettings_), linearizeOperation(linearizeOperationPassed), 
	Hcalib(globalCalib_), globalCalib(globalCalib_), 
	imuIntegration(&Hcalib, imuCalibration, imuSettings, linearizeOperation),
	gravityInit(imuSettings.numMeasurementsGravityInit, imuCalibration),
	shellPoseMutex(FrameShell::shellPoseMutex)
{
	baIntegration = imuIntegration.getBAGTSAMIntegration().get();

	setDefaults();

	selectionMap = new float[globalCalib.wG[0]*globalCalib.hG[0]];

	setClasses();

	ef = new EnergyFunctional(*baIntegration, globalSettings);
	ef->red = &this->treadReduce;

	// MAPPING THREAD IS STARTED HERE
	mappingThread = boost::thread(&FullSystem::mappingLoop, this);
}

/**
 * @brief Destroy the FullSystem Object
 * 
 */
FullSystem::~FullSystem()
{
	blockUntilMappingIsFinished();

	if(globalSettings.setting_logStuff)
	{
		calibLog->close(); delete calibLog;
		numsLog->close(); delete numsLog;
		coarseTrackingLog->close(); delete coarseTrackingLog;
		//errorsLog->close(); delete errorsLog;
		eigenAllLog->close(); delete eigenAllLog;
		eigenPLog->close(); delete eigenPLog;
		eigenALog->close(); delete eigenALog;
		DiagonalLog->close(); delete DiagonalLog;
		variancesLog->close(); delete variancesLog;
		nullspacesLog->close(); delete nullspacesLog;
	}

	delete[] selectionMap;

	for(FrameShell* s : allFrameHistory)
		delete s;
	allFrameHistory.clear();
	for(FrameHessian* fh : unmappedTrackedFrames)
		delete fh;
	unmappedTrackedFrames.clear();

	coarseDistanceMap.reset();
	delete coarseTracker;
	delete coarseTracker_forNewKF;
	coarseInitializer.reset();
	pixelSelector.reset();

	delete ef;
}

void FullSystem::setOriginalCalib(const VecXf &originalCalib, int originalW, int originalH)
{
}

/**
 * @brief Set the gamma
 * 
 * @param BInv 
 */
void FullSystem::setGammaFunction(float* BInv)
{
	if(BInv==0) return;

	// copy BInv.
	memcpy(Hcalib.Binv, BInv, sizeof(float)*256);

	// invert.
	for(int i=1;i<255;i++)
	{
		// find val, such that Binv[val] = i.
		// I dont care about speed for this, so do it the stupid way.

		for(int s=1;s<255;s++)
		{
			if(BInv[s] <= i && BInv[s+1] >= i)
			{
				Hcalib.B[i] = s+(i - BInv[s]) / (BInv[s+1]-BInv[s]);
				break;
			}
		}
	}
	Hcalib.B[0] = 0;
	Hcalib.B[255] = 255;
}

/**
 * @brief Determines best track by tracking possible poses using the coarseTracker class
 * Has three steps
 * 1. Creates a list of likely poses to test
 * 2. Passes possible poses to be tracked by the coarseTracker class till a good track is returned
 * 3. Update variables
 * 
 * If a good track cannot be found, either use the estimated pose from the IMU or assume constant motion
 * If the tracking residual grows too large, this function will kill the program
 * 
 * @param fh 						Current frame
 * @param referenceToFrameHint 		Pose derived from IMU measurements
 * @return std::pair<Vec4, bool> 	Achieved residual, flow vector, and if the tracking is good
 */
std::pair<Vec4, bool> FullSystem::trackNewCoarse(FrameHessian* fh, Sophus::SE3d *referenceToFrameHint)
{
    dmvio::TimeMeasurement timeMeasurement(referenceToFrameHint ? "FullSystem::trackNewCoarse" : "FullSystem::trackNewCoarseNoIMU");
	assert(allFrameHistory.size() > 0);

    for(IOWrap::Output3DWrapper* ow : outputWrapper)
        ow->pushLiveFrame(fh);


	// Frameshell of reference frame
	FrameHessian* lastF = coarseTracker->lastRef;

	AffLight aff_last_2_l = AffLight(0,0);


	// ============== Create a list of possible poses ===================
	// Poses will be relative to the reference frame

    // lastF_2_fh_tries holds list of possible poses
    std::vector<SE3,Eigen::aligned_allocator<SE3>> lastF_2_fh_tries;

	// Have IMU data
	// !inu: Use imu data as initilization motion
    if(referenceToFrameHint)
    {
        // We got a hint (typically from IMU) where our pose is, so we don't need the random initializations below
        lastF_2_fh_tries.push_back(*referenceToFrameHint);
        {
            // Lock on global pose consistency (probably we don't need this for AffineLight, but just to make sure)
            boost::unique_lock<boost::mutex> crlock(shellPoseMutex);

            // Set Affine light to last frame, where tracking was good
            for(int i = allFrameHistory.size() - 2; i >= 0; i--)
            {
                FrameShell* slast = allFrameHistory[i];
                if(slast->trackingWasGood)
                {
                    aff_last_2_l = slast->aff_g2l;
                    break;
                }
                if(slast->trackingRef != lastF->shell)
                {
                    std::cout << "WARNING: No well tracked frame with the same tracking ref available!" << std::endl;
                    aff_last_2_l = lastF->aff_g2l();
                    break;
                }
            }
        }
    }

	// No IMU data
    if(!referenceToFrameHint)
    {
		// Randomly initialize lastF_2_fh_tries with possible poses
        if(allFrameHistory.size() == 2) // Case for intialization
            for(unsigned int i=0;i<lastF_2_fh_tries.size();i++) lastF_2_fh_tries.push_back(SE3());
        else
        {
            FrameShell* slast = allFrameHistory[allFrameHistory.size()-2];
            FrameShell* sprelast = allFrameHistory[allFrameHistory.size()-3];
            SE3 slast_2_sprelast;
            SE3 lastF_2_slast;
            {	
				// lock on global pose consistency!
                boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
                slast_2_sprelast = sprelast->camToWorld.inverse() * slast->camToWorld;	// Constant motion from last frameshell
                lastF_2_slast = slast->camToWorld.inverse() * lastF->shell->camToWorld;	// Changes current motion to reference frame
                aff_last_2_l = slast->aff_g2l;
            }
            SE3 fh_2_slast = slast_2_sprelast;// assumed to be the same as fh_2_slast.

            // get last delta-movement.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast);							// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * fh_2_slast.inverse() * lastF_2_slast);	// assume double motion (frame skipped)
            lastF_2_fh_tries.push_back(SE3::exp(fh_2_slast.log()*0.5).inverse() * lastF_2_slast); 		// assume half motion.
            lastF_2_fh_tries.push_back(lastF_2_slast); 													// assume zero motion.
            lastF_2_fh_tries.push_back(SE3()); 															// assume zero motion from KF.

            // Just try a TON of different initializations (all rotations). In the end,
            // if they don't work they will only be tried on the coarsest level, which is super fast anyway.
            // Tracking failure is really bad, so it needs to be avoided
            for(float rotDelta=0.02; rotDelta < 0.05; rotDelta++)
            {
				// 27 different small rotations in different directions
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,0,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,rotDelta,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,0,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,0,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,-rotDelta,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,0,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,rotDelta,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,rotDelta,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,0,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,rotDelta,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,-rotDelta,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,0,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,-rotDelta,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,rotDelta,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,0,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,-rotDelta,0), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,0,-rotDelta,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,0,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,-rotDelta,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,-rotDelta,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,rotDelta,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,-rotDelta,rotDelta,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,-rotDelta,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,-rotDelta,rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,rotDelta,-rotDelta), Vec3(0,0,0)));
                lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(Sophus::Quaterniond(1,rotDelta,rotDelta,rotDelta), Vec3(0,0,0)));
            }

            if(!slast->poseValid || !sprelast->poseValid || !lastF->shell->poseValid)
            {
                lastF_2_fh_tries.clear();
                lastF_2_fh_tries.push_back(SE3());
            }
        }
    }

	// ============== Test all of the pose guesses ===================
	Vec3 flowVecs = Vec3(100,100,100); // Flow vector of tracked motion
	SE3 lastF_2_fh = SE3();
	AffLight aff_g2l = AffLight(0,0);

	// As long as maxResForImmediateAccept is not reached, it will continue through the options.
	// Keep track of the so-far best achieved residual for each level in achievedRes.
	// If on a coarse level, tracking is WORSE than achievedRes, we will not continue to save time.

	bool trackingGoodRet = false;
	Vec5 achievedRes = Vec5::Constant(NAN);
	bool haveOneGood = false;
	int tryIterations=0;

	for(unsigned int i=0;i<lastF_2_fh_tries.size();i++) // Try tracking for all poses in lastF_2_fh_tries
	{
		AffLight aff_g2l_this = aff_last_2_l;
		SE3 lastF_2_fh_this = lastF_2_fh_tries[i];

		// Tracking is done by the coarseTracker class
		bool trackingIsGood = coarseTracker->trackNewestCoarse(
				fh, lastF_2_fh_this, aff_g2l_this,
				globalSettings.pyrLevelsUsed-1,
				achievedRes);	// in each level this has to be at least as good as the last try.
		
		tryIterations++;


		if(trackingIsGood)
        {
		    trackingGoodRet = true;
        }
		if(!trackingIsGood && imuIntegration.setting_useIMU)
		{
			std::cout << "WARNING: Coarse tracker thinks that tracking was not good!" << std::endl;

			// In IMU mode we can still estimate the pose sufficiently, even if vision is bad.
			trackingIsGood = true;
		}

		if(i != 0)
		{
			if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
				printf("RE-TRACK ATTEMPT %u with initOption %u and start-lvl %u (ab %f %f): %f %f %f %f %f -> %f %f %f %f %f \n",
						i, i, globalSettings.pyrLevelsUsed-1,
						aff_g2l_this.a,
						aff_g2l_this.b,
						achievedRes[0],
						achievedRes[1],
						achievedRes[2],
						achievedRes[3],
						achievedRes[4],
						coarseTracker->lastResidualsStats[0],
						coarseTracker->lastResidualsStats[1],
						coarseTracker->lastResidualsStats[2],
						coarseTracker->lastResidualsStats[3],
						coarseTracker->lastResidualsStats[4]);
		}

		// ============== Update variables if there is a good track ===================
		// Track for given motion is sucessful is:
		// 1. Tracking does not return invalid values
		// 2. Residual is finite
		// 3. Residual is smaller than last achieved residual
		if(trackingIsGood && std::isfinite((float)coarseTracker->lastResidualsStats[0]) && !(coarseTracker->lastResidualsStats[0] >=  achievedRes[0]))
		{
			flowVecs = coarseTracker->lastFlowIndicators;
			aff_g2l = aff_g2l_this;
			lastF_2_fh = lastF_2_fh_this;
			haveOneGood = true;
		}

		// Set achieved residuals to new values
		if(haveOneGood)
		{
			for(int j=0;j<5;j++)
			{
				if(!std::isfinite((float)achievedRes[j]) || achievedRes[j] > coarseTracker->lastResidualsStats[j])
					achievedRes[j] = coarseTracker->lastResidualsStats[j];
			}
		}

		// Stop testing other tracks if result is good enough
        if(haveOneGood &&  achievedRes[0] < lastCoarseRMSE[0]*globalSettings.setting_reTrackThreshold)
            break;
	}

	// Case if no good track is found
	// Use IMU motion if possible else assume constant motion 
	if(!haveOneGood)
	{
        printf("BIG ERROR! tracking failed entirely. Take predicted pose and hope we may somehow recover.\n");
		flowVecs = Vec3(0,0,0);
		aff_g2l = aff_last_2_l;
		lastF_2_fh = lastF_2_fh_tries[0]; // imu!: IMU motion if IMU is active, else this is constant motion
        std::cout << "Predicted pose:\n" << lastF_2_fh.matrix() << std::endl;

		// Degenerate case where tracking completely fails
		if(lastF_2_fh.translation().norm() > 100000 || lastF_2_fh.matrix().hasNaN())
        {
            std::cout << "TRACKING FAILED ENTIRELY, NO HOPE TO RECOVER" << std::endl;
		    std::cerr << "TRACKING FAILED ENTIRELY, NO HOPE TO RECOVER" << std::endl;
		    exit(1);
        }
	}

	lastCoarseRMSE = achievedRes;

	// no lock required, as fh is not used anywhere yet.
	fh->shell->camToTrackingRef = lastF_2_fh.inverse();
	fh->shell->trackingRef = lastF->shell;
	fh->shell->aff_g2l = aff_g2l;
	fh->shell->camToWorld = fh->shell->trackingRef->camToWorld * fh->shell->camToTrackingRef;
	fh->shell->trackingWasGood = trackingGoodRet;


	if(coarseTracker->firstCoarseRMSE < 0)
		coarseTracker->firstCoarseRMSE = achievedRes[0];

    if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        printf("Coarse Tracker tracked ab = %f %f (exp %f). Res %f!\n", aff_g2l.a, aff_g2l.b, fh->ab_exposure, achievedRes[0]);


	if(globalSettings.setting_logStuff)
	{
		(*coarseTrackingLog) << std::setprecision(16)
						<< fh->shell->id << " "
						<< fh->shell->timestamp << " "
						<< fh->ab_exposure << " "
						<< fh->shell->camToWorld.log().transpose() << " "
						<< aff_g2l.a << " "
						<< aff_g2l.b << " "
						<< achievedRes[0] << " "
						<< tryIterations << "\n";
	}

	return std::make_pair(Vec4(achievedRes[0], flowVecs[0], flowVecs[1], flowVecs[2]), trackingGoodRet);
}


/**
 * @brief Traces all the immature points
 * 
 * Because immature points are only made for new frames,
 * the number of immature points per frame will decrease over time
 * 
 * @param fh Current frame
 */
void FullSystem::traceNewCoarse(FrameHessian* fh)
{
    dmvio::TimeMeasurement timeMeasurement("traceNewCoarse");
	boost::unique_lock<boost::mutex> lock(mapMutex);

	// Note the number of immature points whose trace ended with the respective conditions for debugging purposes
	int trace_total=0, trace_good=0, trace_oob=0, trace_out=0, trace_skip=0, trace_badcondition=0, trace_uninitialized=0;

	Mat33f K = Mat33f::Identity();
	// Get K, the camera calibration matrix
	K(0,0) = Hcalib.fxl();
	K(1,1) = Hcalib.fyl();
	K(0,2) = Hcalib.cxl();
	K(1,2) = Hcalib.cyl();

	for(FrameHessian* host : frameHessians)	// Go through all active frames
	{
		// Get SE(3) matrix from old to new position
		SE3 hostToNew = fh->PRE_worldToCam * host->PRE_camToWorld;
		// Seperate into rotation and translation parts, multiply by K at this stage for calculational efficency
		Mat33f KRKi = K * hostToNew.rotationMatrix().cast<float>() * K.inverse();
		Vec3f Kt = K * hostToNew.translation().cast<float>();
		// Get matrix that transforms points from one frame to another
		Vec2f aff = AffLight::fromToVecExposure(host->ab_exposure, fh->ab_exposure, host->aff_g2l(), fh->aff_g2l()).cast<float>();

		for(ImmaturePoint* ph : host->immaturePoints) // For all immature points in active host frame
		{
			ph->traceOn(fh, KRKi, Kt, aff, &Hcalib, !setting_debugout_runquiet ); // trace the immature point

			if(ph->lastTraceStatus==ImmaturePointStatus::IPS_GOOD) trace_good++;
			if(ph->lastTraceStatus==ImmaturePointStatus::IPS_BADCONDITION) trace_badcondition++;
			if(ph->lastTraceStatus==ImmaturePointStatus::IPS_OOB) trace_oob++;
			if(ph->lastTraceStatus==ImmaturePointStatus::IPS_OUTLIER ||
				ph->lastTraceStatus==ImmaturePointStatus::IPS_OUTLIER_OUT) trace_out++;
			if(ph->lastTraceStatus==ImmaturePointStatus::IPS_SKIPPED) trace_skip++;
			if(ph->lastTraceStatus==ImmaturePointStatus::IPS_UNINITIALIZED) trace_uninitialized++;
			trace_total++;
		}
	}

/*
	printf("ADD: TRACE: %'d points. %'d (%.0f%%) good. %'d (%.0f%%) skip. %'d (%.0f%%) badcond. %'d (%.0f%%) oob. %'d (%.0f%%) out. %'d (%.0f%%) uninit.\n",
			trace_total,
			trace_good, 100*trace_good/(float)trace_total,
			trace_skip, 100*trace_skip/(float)trace_total,
			trace_badcondition, 100*trace_badcondition/(float)trace_total,
			trace_oob, 100*trace_oob/(float)trace_total,
			trace_out, 100*trace_out/(float)trace_total,
			trace_uninitialized, 100*trace_uninitialized/(float)trace_total);
*/
}

/**
 * @brief Helper function for activatePointsMT
 * 
 * Does the type conversion from ImmaturePoints to PointHessians
 * 
 * @param optimized 
 * @param toOptimize 
 * @param min 
 * @param max 
 * @param stats 
 * @param tid 
 */
void FullSystem::activatePointsMT_Reductor(
		std::vector<PointHessian*>* optimized,
		std::vector<ImmaturePoint*>* toOptimize,
		int min, int max, Vec10* stats, int tid)
{
	ImmaturePointTemporaryResidual* tr = new ImmaturePointTemporaryResidual[frameHessians.size()];
	for(int k=min;k<max;k++)
	{
		(*optimized)[k] = optimizeImmaturePoint((*toOptimize)[k],1,tr);
	}
	delete[] tr;
}

/**
 * @brief Converts immature points to active points
 * 
 */
void FullSystem::activatePointsMT()
{
    dmvio::TimeMeasurement timeMeasurement("activatePointsMT");

	// ============== Update variables for desired point density ===================
	// Change the currentMinActDist in order to achieve the desired point number
    if(ef->nPoints < globalSettings.setting_desiredPointDensity*0.66)
		currentMinActDist -= 0.8;
	if(ef->nPoints < globalSettings.setting_desiredPointDensity*0.8)
		currentMinActDist -= 0.5;
	else if(ef->nPoints < globalSettings.setting_desiredPointDensity*0.9)
		currentMinActDist -= 0.2;
	else if(ef->nPoints < globalSettings.setting_desiredPointDensity)
		currentMinActDist -= 0.1;

	if(ef->nPoints > globalSettings.setting_desiredPointDensity*1.5)
		currentMinActDist += 0.8;
	if(ef->nPoints > globalSettings.setting_desiredPointDensity*1.3)
		currentMinActDist += 0.5;
	if(ef->nPoints > globalSettings.setting_desiredPointDensity*1.15)
		currentMinActDist += 0.2;
	if(ef->nPoints > globalSettings.setting_desiredPointDensity)
		currentMinActDist += 0.1;

	// Max and min currentMinActDist values
	if(currentMinActDist < 0) currentMinActDist = 0;
	if(currentMinActDist > 4) currentMinActDist = 4;

    if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        printf("SPARSITY:  MinActDist %f (need %d points, have %d points)!\n",
                currentMinActDist, (int)(globalSettings.setting_desiredPointDensity), ef->nPoints);

	// ============== Loop through all of the immature points in active frames ===================
	FrameHessian* newestHs = frameHessians.back();

	// Make dist map.
	coarseDistanceMap->makeK(&Hcalib);
	coarseDistanceMap->makeDistanceMap(frameHessians, newestHs);
	//coarseTracker->debugPlotDistMap("distMap");

	unsigned int max_points = globalSettings.setting_desiredImmatureDensity;
	std::vector<ImmaturePoint*> toOptimize; toOptimize.reserve(max_points);

	for(FrameHessian* host : frameHessians)	// go through all active frames
	{
		if(host == newestHs) continue; // exclude newest frame

		SE3 fhToNew = newestHs->PRE_worldToCam * host->PRE_camToWorld; // Transformation matrix from host to newest frame
		// Seperate transformation matrix to rotation and translation parts
		Mat33f KRKi = (coarseDistanceMap->K[1] * fhToNew.rotationMatrix().cast<float>() * coarseDistanceMap->Ki[0]);
		Vec3f Kt = (coarseDistanceMap->K[1] * fhToNew.translation().cast<float>());

		for(unsigned int i=0;i<host->immaturePoints.size();i+=1) // for every immature point in the frame
		{
			ImmaturePoint* ph = host->immaturePoints[i];
			ph->idxInImmaturePoints = i;

			// ============== Delete invalid immature points ===================
			// Delete points that have never been traced successfully, or that are outlier on the last trace.
			if(!std::isfinite(ph->idepth_max) || ph->lastTraceStatus == IPS_OUTLIER ||
			ph->lastTraceStatus == IPS_OUTLIER_OUT)
			{
				// immature_invalid_deleted++;
				delete ph;
				host->immaturePoints[i]=0;
				continue;
			}

			// Activate only if this is true.
			// Don't activate if immature point outlier or unintialized
			// Activate if immature point was
			// traced along a good path, has good energy quality, and realistic depth
			bool canActivate = (ph->lastTraceStatus == IPS_GOOD
					|| ph->lastTraceStatus == IPS_SKIPPED
					|| ph->lastTraceStatus == IPS_BADCONDITION
					|| ph->lastTraceStatus == IPS_OOB )
							&& ph->lastTracePixelInterval < 8
							&& ph->quality > globalSettings.setting_minTraceQuality
							&& (ph->idepth_max+ph->idepth_min) > 0;

			// if cannot activate the point, skip it. Maybe also delete it.
			if(!canActivate)
			{
				// if point will be out afterwards due 
				// to being oob or in a deleted frame, delete it instead
				if(ph->host->flaggedForMarginalization || ph->lastTraceStatus == IPS_OOB)
				{
					// immature_notReady_deleted++;
					delete ph;
					host->immaturePoints[i]=0;
				}
				// immature_notReady_skipped++;
				continue;
			}

			// ============== Add immature point to activation list if it meets conditions ===================
			// Determine if the point is far away enough from other points
			// Points are only activated if it has good spacing compared to other active points
			// once projected on the newest keyframe
			Vec3f ptp = KRKi * Vec3f(ph->u, ph->v, 1) + Kt*(0.5f*(ph->idepth_max+ph->idepth_min));
			int u = ptp[0] / ptp[2] + 0.5f;
			int v = ptp[1] / ptp[2] + 0.5f;

			if((u > 0 && v > 0 && u < globalCalib.wG[1] && v < globalCalib.hG[1])) // delete point if it is out of bounds
			{
				// Find distance to closest point
				float dist = coarseDistanceMap->fwdWarpedIDDistFinal[u+globalCalib.wG[1]*v] + (ptp[0]-floorf((float)(ptp[0])));

				if(dist>=currentMinActDist * ph->my_type) // check if point meets density requirements
				{
					coarseDistanceMap->addIntoDistFinal(u,v);
					toOptimize.push_back(ph); // add point to list of points to be optimized
				}
			}
			else
			{
				delete ph;
				host->immaturePoints[i]=0;
			}
		}
	}


//	printf("ACTIVATE: %d. (del %d, notReady %d, marg %d, good %d, marg-skip %d)\n",
//			(int)toOptimize.size(), immature_deleted, immature_notReady, immature_needMarg, immature_want, immature_margskip);

	std::vector<PointHessian*> optimized; optimized.resize(toOptimize.size());

	// ============== Activate points in activation list ===================
	// Activate points by optimizing them and converting them into PointHessians
	if(!globalSettings.settings_no_multiThreading)
		treadReduce.reduce(boost::bind(&FullSystem::activatePointsMT_Reductor, this, &optimized, &toOptimize, _1, _2, _3, _4), 0, toOptimize.size(), 50);
	else
		activatePointsMT_Reductor(&optimized, &toOptimize, 0, toOptimize.size(), 0, 0);

	// Check if all the points are valid after optimization
	for(unsigned k=0;k<toOptimize.size();k++)
	{
		PointHessian* newpoint = optimized[k];
		ImmaturePoint* ph = toOptimize[k];

		if(newpoint != 0 && newpoint != (PointHessian*)((long)(-1)))
		{
			// Add new point into the optimization and active point list
			// Remove point from immature point list
			newpoint->host->immaturePoints[ph->idxInImmaturePoints]=0;
			// Add active point
			newpoint->host->pointHessians.push_back(newpoint);
			ef->insertPoint(newpoint);
			for(PointFrameResidual* r : newpoint->residuals)
				ef->insertResidual(r);
			assert(newpoint->efPoint != 0);
			delete ph;
		}
		else if(newpoint == (PointHessian*)((long)(-1)) || ph->lastTraceStatus==IPS_OOB)
		{
			ph->host->immaturePoints[ph->idxInImmaturePoints]=0;
            delete ph;
		}
		else
		{
			assert(newpoint == 0 || newpoint == (PointHessian*)((long)(-1)));
		}
	}

	// Removes immature points that have been deleted from the frames
	for(FrameHessian* host : frameHessians)
	{
		for(int i=0;i<(int)host->immaturePoints.size();i++)
		{
			if(host->immaturePoints[i]==0)
			{
				host->immaturePoints[i] = host->immaturePoints.back();
				host->immaturePoints.pop_back();
				i--;
			}
		}
	}
}

void FullSystem::activatePointsOldFirst()
{
	assert(false);
}


/**
 * @brief Flags bad points for removal
 * 
 */
void FullSystem::flagPointsForRemoval()
{
	assert(EFIndicesValid);

	std::vector<FrameHessian*> fhsToKeepPoints;
	std::vector<FrameHessian*> fhsToMargPoints;

	//if(globalSettings.setting_margPointVisWindow>0)
	{
		for(int i=((int)frameHessians.size())-1;i>=0 && i >= ((int)frameHessians.size());i--)
			if(!frameHessians[i]->flaggedForMarginalization) fhsToKeepPoints.push_back(frameHessians[i]);

		for(int i=0; i< (int)frameHessians.size();i++)
			if(frameHessians[i]->flaggedForMarginalization) fhsToMargPoints.push_back(frameHessians[i]);
	}

	//ef->setAdjointsF();
	//ef->setDeltaF(&Hcalib);
	int flag_oob=0, flag_in=0, flag_inin=0, flag_nores=0;

	for(FrameHessian* host : frameHessians)	// go through all active frames
	{
		for(unsigned int i=0;i<host->pointHessians.size();i++)
		{
			PointHessian* ph = host->pointHessians[i];
			if(ph==0) continue;

			// Remove points that have invalid depths or no more residuals
			if(ph->idepth_scaled < globalSettings.setting_minIdepth || ph->residuals.size()==0)
			{
				host->pointHessiansOut.push_back(ph);
				ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
				host->pointHessians[i]=0;
				flag_nores++;
			}
			// Remove points that are out of bounds or are in frames that are to be marginalized
			else if(ph->isOOB(fhsToKeepPoints, fhsToMargPoints) || host->flaggedForMarginalization)
			{
				flag_oob++;

				if(ph->isInlierNew()) // Case for points with low residual values
				{
					flag_in++;
					int ngoodRes=0;

					// Linearize all
					for(PointFrameResidual* r : ph->residuals)
					{
						r->resetOOB();
						r->linearize(&Hcalib);
						r->efResidual->isLinearized = false;
						r->applyRes(true);
						if(r->efResidual->isActive())
						{
							r->efResidual->fixLinearizationF(ef);
							ngoodRes++;
						}
					}

                    if(ph->idepth_hessian > globalSettings.setting_minIdepthH_marg)
					{
						flag_inin++;
						ph->efPoint->stateFlag = EFPointStatus::PS_MARGINALIZE;
						host->pointHessiansMarginalized.push_back(ph);
					}
					else
					{
						ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
						host->pointHessiansOut.push_back(ph);
					}
				}
				else
				{
					host->pointHessiansOut.push_back(ph);
					ph->efPoint->stateFlag = EFPointStatus::PS_DROP;

					//printf("drop point in frame %d (%d goodRes, %d activeRes)\n", ph->host->idx, ph->numGoodResiduals, (int)ph->residuals.size());
				}
				host->pointHessians[i]=0;
			}
		}

		// Set the pointHessians list to exclude the removed points
		for(int i=0;i<(int)host->pointHessians.size();i++)
		{
			if(host->pointHessians[i]==0)
			{
				host->pointHessians[i] = host->pointHessians.back();
				host->pointHessians.pop_back();
				i--;
			}
		}
	}
}

/**
 * @brief Executes odometry when new frame is added
 * 
 * This function is the starting point for most of the other functions
 * 
 * The function is passed the IMU-data from the previous frame until the current frame.
 * 
 * Steps:
 * If no frames given, set up first frame
 * After given first frame, it does init process till init finished
 * After init finished, does standard tracking
 * Will create keyframes and update map as required
 * 
 * @param image 
 * @param id 
 * @param imuData 
 * @param gtData 
 */
void FullSystem::addActiveFrame(ImageAndExposure* image, int id, dmvio::IMUData* imuData, dmvio::GTData* gtData)
{
    // Measure Time of the time measurement.
    dmvio::TimeMeasurement timeMeasurementMeasurement("timeMeasurement");
    dmvio::TimeMeasurement timeMeasurementZero("zero");
    timeMeasurementZero.end();
    timeMeasurementMeasurement.end();

    dmvio::TimeMeasurement timeMeasurement("addActiveFrame");
	boost::unique_lock<boost::mutex> lock(trackMutex);


	dmvio::TimeMeasurement measureInit("initObjectsAndMakeImage");
	// =========================== Add new frame into allFrameHistory =========================
	// Initialize variables
	FrameHessian* fh = new FrameHessian(globalCalib, globalSettings);
	FrameShell* shell = new FrameShell();
	shell->camToWorld = SE3();
	shell->aff_g2l = AffLight(0,0);
    shell->marginalizedAt = shell->id = allFrameHistory.size();
    shell->timestamp = image->timestamp;
    shell->incoming_id = id;
	fh->shell = shell;
	allFrameHistory.push_back(shell);


    // =========================== Place image and image settings into frame =========================
	fh->ab_exposure = image->exposure_time;
	fh->makeImages(image->image, &Hcalib); // Image derivative and gradient is also calculated
	if(image->useColour){
		fh->makeColourImages(image->r_image, image->g_image, image->b_image);
	}

    measureInit.end();

	// =========================== Process Image =========================
	if(!initialized) // Need initialization case
	{
		// =========================== Init using image if needed =========================
		if(coarseInitializer->frameID<0) // first frame set. fh is kept by coarseInitializer.
		{
            // Only in this case no IMU-data is accumulated for the BA as this is the first frame.
		    dmvio::TimeMeasurement initMeasure("InitializerFirstFrame");

			coarseInitializer->setFirst(&Hcalib, fh);
		// imu!: Start calculating gravity
            if(imuIntegration.setting_useIMU)
            {
                gravityInit.addMeasure(*imuData, Sophus::SE3d());
            }
            for(IOWrap::Output3DWrapper* ow : outputWrapper)
                ow->publishSystemStatus(dmvio::VISUAL_INIT);
        }
		else // runs till initialization is complete
        {
            dmvio::TimeMeasurement initMeasure("InitializerOtherFrames");

			// Run the initialization through the coarseInitializer class
			bool initDone = coarseInitializer->trackFrame(fh, outputWrapper);
			
			// imu!: Add initial data to imu optimization
			if(imuIntegration.setting_useIMU)
			{
                imuIntegration.addIMUDataToBA(*imuData);
				Sophus::SE3d imuToWorld = gravityInit.addMeasure(*imuData, Sophus::SE3d());
				if(initDone)
				{
					firstPose = imuToWorld * imuIntegration.TS_cam_imu.inverse();
				}
			}

            if (initDone)
            {
				// Use the values from the coarseInitializer to finalize the initialization
                initializeFromInitializer(fh);

		// imu!: Set GT data
                if(imuIntegration.setting_useIMU && linearizeOperation)
                {
                    imuIntegration.setGTData(gtData, fh->shell->id);
                }

                lock.unlock();
                initMeasure.end();
                for(IOWrap::Output3DWrapper* ow : outputWrapper)
                    ow->publishSystemStatus(dmvio::VISUAL_ONLY);

                deliverTrackedFrame(fh, true);
            } 
			else  // if still initializing
            {
                // Maybe change first frame.
                double timeBetweenFrames = fh->shell->timestamp - coarseInitializer->firstFrame->shell->timestamp;
                if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage){
					std::cout << "InitTimeBetweenFrames: " << timeBetweenFrames << std::endl;
				}

		// imu!: Time between frames cannot be too high or else the imu data will be inaccurate
                if(timeBetweenFrames > imuIntegration.getImuSettings().maxTimeBetweenInitFrames)
                {
                    // Do full reset so that the next frame becomes the first initializer frame.
                    setting_fullResetRequested = true;
                }
				else
                {
                    fh->shell->poseValid = false;
                    delete fh;
                }
            }
        }

		return;
	}
	else // do standard front-end operation.
	{
        dmvio::TimeMeasurement coarseTrackingTime("fullCoarseTracking");
		int lastFrameId = -1;

		// =========================== Swap tracking reference =========================
		bool trackingRefChanged = false;
		// If there is a newer keyframe
		// Set that as the reference frame for tracking
		if(coarseTracker_forNewKF->refFrameID > coarseTracker->refFrameID)
		{
            dmvio::TimeMeasurement referenceSwapTime("swapTrackingRef");

			boost::unique_lock<boost::mutex> crlock(coarseTrackerSwapMutex);
			CoarseTracker* tmp = coarseTracker; coarseTracker=coarseTracker_forNewKF; coarseTracker_forNewKF=tmp;

			// imu!: Set up the imu optimzation
			if(imuIntegration.setting_useIMU)
			{
			    // BA for new keyframe has finished and we have a new tracking reference.
                if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
                {
                    std::cout << "New ref frame id: " << coarseTracker->refFrameID << " prepared keyframe id: "
                              << imuIntegration.getPreparedKeyframe() << std::endl;
                }

                lastFrameId = coarseTracker->refFrameID;

				assert(coarseTracker->refFrameID == imuIntegration.getPreparedKeyframe());
				SE3 lastRefToNewRef = imuIntegration.initCoarseGraph();

				trackingRefChanged = true;
			}
		}

        SE3 *referenceToFramePassed = 0;

	// imu!: Add imu data
        if(imuIntegration.setting_useIMU)
        {
			SE3 referenceToFrame = imuIntegration.addIMUData(*imuData, fh->shell->id,
                                                                fh->shell->timestamp, trackingRefChanged, lastFrameId);
           
		    // If initialized we use the prediction from IMU data as initialization for the coarse tracking.
            referenceToFramePassed = &referenceToFrame;
			if(!imuIntegration.isCoarseInitialized())
            {
			    referenceToFramePassed = nullptr;
            }
            imuIntegration.addIMUDataToBA(*imuData);
        }

		// =========================== Do coarse tracking =========================
		// Coarse Tracking is only done on only the new frame and reference frame
		// The reference frame should be the latest keyframe
        std::pair<Vec4, bool> pair = trackNewCoarse(fh, referenceToFramePassed);

        dso::Vec4 tres = std::move(pair.first); // Achieved residual and flow vectors
        bool forceNoKF = !pair.second; // If coarse tracking was bad don't make KF.
        bool forceKF = false;

		// Case if tracking failed
		if(!std::isfinite((double)tres[0]) || !std::isfinite((double)tres[1]) || !std::isfinite((double)tres[2]) || !std::isfinite((double)tres[3]))
        {
            if(imuIntegration.setting_useIMU)
            {
                // If completely Nan, don't force noKF!
                forceNoKF = false;
                forceKF = true; // actually we force a KF in that situation as there are no points to track.
            }else
            {
                printf("Initial Tracking failed: LOST!\n");
                isLost=true;
                return;
            }
        }

		// =========================== Make keyframe or non-keyframe =========================
        double timeSinceLastKeyframe = fh->shell->timestamp - allKeyFramesHistory.back()->timestamp;
		bool needToMakeKF = false;

		// Decide if keyframe or non-keyframe needs to be made
		// If globalSettings.setting_keyframesPerSecond is set, the keyframe is made at a specific frequency
		// Otherwise, the keyframe is made depending on specific conditions

		if(globalSettings.setting_keyframesPerSecond > 0) // fixed keyframe rate
		{
			// Make keyframe at a set frequency
			needToMakeKF = allFrameHistory.size()== 1 ||
					(fh->shell->timestamp - allKeyFramesHistory.back()->timestamp) > 0.95f/globalSettings.setting_keyframesPerSecond;
		}
		else // makes keyframe under specific conditions
		{
			Vec2 refToFh=AffLight::fromToVecExposure(coarseTracker->lastRef->ab_exposure, fh->ab_exposure,
					coarseTracker->lastRef_aff_g2l, fh->shell->aff_g2l);

			// Make keyframe if:
			// 1. If the field of view changes too much. FOV change is measured by the mean optical flow
			// 2. If motion causes occlusions and disocclusions. This is measured by mean translation flow
			// 3. If camera exposure is changed significantly.
			// 4. Residual is too high.
			// 5. If max time between keyframes is surpassed.
			// 6. The IMU system needs a keyframe.
			needToMakeKF = allFrameHistory.size()== 1 ||
					globalSettings.setting_kfGlobalWeight*globalSettings.setting_maxShiftWeightT *  sqrtf((double)tres[1]) / (globalCalib.wG[0]+globalCalib.hG[0]) + 			// Translation flow
					globalSettings.setting_kfGlobalWeight*globalSettings.setting_maxShiftWeightR *  sqrtf((double)tres[2]) / (globalCalib.wG[0]+globalCalib.hG[0]) + 			// Rotation flow
					globalSettings.setting_kfGlobalWeight*globalSettings.setting_maxShiftWeightRT * sqrtf((double)tres[3]) / (globalCalib.wG[0]+globalCalib.hG[0]) + 			// Overal flow
					globalSettings.setting_kfGlobalWeight*globalSettings.setting_maxAffineWeight * fabs(logf((float)refToFh[0])) > 1 || 				// Exposure changes
					2*coarseTracker->firstCoarseRMSE < tres[0] ||														// Residual value
                    (globalSettings.setting_maxTimeBetweenKeyframes > 0 && timeSinceLastKeyframe > globalSettings.setting_maxTimeBetweenKeyframes) || // Time between keyframes
                    forceKF;																							// IMU system needs keyframe

			if(needToMakeKF && !setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
            {
                std::cout << "Time since last keyframe: " << timeSinceLastKeyframe << std::endl;
            }
		}

		// Do not make keyframe if translation is too low
		double transNorm = fh->shell->camToTrackingRef.translation().norm() * imuIntegration.getCoarseScale();
		if(imuIntegration.isCoarseInitialized() && transNorm < globalSettings.setting_forceNoKFTranslationThresh)
        {
		    forceNoKF = true;
        }
		if(forceNoKF)
        {
		    std::cout << "Forcing NO KF!" << std::endl;
		    needToMakeKF = false;
        }

        if(needToMakeKF)
        {
            int prevKFId = fh->shell->trackingRef->id;
            // In non real-time mode this will always be accurate, 
			// but in real time mode the printout in makeKeyframe is correct 
			// (because some of these KFs do not end up getting created).
            int framesBetweenKFs = fh->shell->id - prevKFId - 1;

            // Enforce globalSettings.setting_minFramesBetweenKeyframes.
            if(framesBetweenKFs < (int) globalSettings.setting_minFramesBetweenKeyframes) // if integer value is smaller we just skip.
            {
                if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
					std::cout << "Skipping KF because of minFramesBetweenKeyframes." << std::endl;
                needToMakeKF = false;
            }
			else if(framesBetweenKFs < globalSettings.setting_minFramesBetweenKeyframes) // Enforce it for non-integer values.
            {
                double fractionalPart = globalSettings.setting_minFramesBetweenKeyframes - (int) globalSettings.setting_minFramesBetweenKeyframes;
                framesBetweenKFsRest += fractionalPart;
                if(framesBetweenKFsRest >= 1.0)
                {
                    if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
						std::cout << "Skipping KF because of minFramesBetweenKeyframes." << std::endl;
                    needToMakeKF = false;
                    framesBetweenKFsRest--;
                }
            }
        }

	// imu!: Finish handling the imu data
        if(imuIntegration.setting_useIMU)
        {
            imuIntegration.finishCoarseTracking(*(fh->shell), needToMakeKF);
        }

        if(needToMakeKF && imuIntegration.setting_useIMU && linearizeOperation)
        {
            imuIntegration.setGTData(gtData, fh->shell->id);
        }

        dmvio::TimeMeasurement timeLastStuff("afterCoarseTracking");

        for(IOWrap::Output3DWrapper* ow : outputWrapper)
            ow->publishCamPose(fh->shell, &Hcalib);

        lock.unlock();
        timeLastStuff.end();
        coarseTrackingTime.end();
		deliverTrackedFrame(fh, needToMakeKF);
		return;
	}
}


/**
 * @brief Creates a keyframe or non-keyframe
 * 
 * Also parses outputted frame for GUI and debugging
 * Helper function for addActiveFrame
 * 
 * @param fh 
 * @param needKF 
 */
void FullSystem::deliverTrackedFrame(FrameHessian* fh, bool needKF)
{
    dmvio::TimeMeasurement timeMeasurement("deliverTrackedFrame");
	// There seems to be exactly one instance where needKF is false but the mapper 
	// creates a keyframe nevertheless: if it is the second tracked frame (so it will become the third keyframe in total)
	// There are also some cases where needKF is true but the mapper does not create a keyframe.

	bool alreadyPreparedKF = imuIntegration.setting_useIMU && imuIntegration.getPreparedKeyframe() != -1 && !linearizeOperation;

    if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
    {
        std::cout << "Frame history size: " << allFrameHistory.size() << std::endl;
    }

    if((needKF || (!secondKeyframeDone && !linearizeOperation)) && imuIntegration.setting_useIMU && !alreadyPreparedKF)
    {
        // imu!: PrepareKeyframe tells the IMU-Integration that this frame will become a keyframe. -> don' marginalize it during addIMUData.
        // Also resets the IMU preintegration for the BA.
        if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        {
            std::cout << "Preparing keyframe: " << fh->shell->id << std::endl;
        }
        imuIntegration.prepareKeyframe(fh->shell->id);

		if(!needKF)
		{
			secondKeyframeDone = true;
		}
    }
	else
	{
        if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        {
            std::cout << "Creating a non-keyframe: " << fh->shell->id << std::endl;
        }
    }

	if(linearizeOperation) // play as fast as possible
	{
		
#ifdef GRAPHICAL_DEBUG
		if(globalSettings.setting_goStepByStep && lastRefStopID != coarseTracker->refFrameID)
		{
			MinimalImageF3 img(globalCalib.wG[0], globalCalib.hG[0], fh->dI);
			if (!globalSettings.setting_disableAllDisplay){
				IOWrap::displayImage("frameToTrack", &img);
				while(true)
				{
					char k=IOWrap::waitKey(0);
					if(k==' ') break;
					globalSettings.handleKey( k );
				}
			}
			lastRefStopID = coarseTracker->refFrameID;
		}
		else{
			if(!globalSettings.setting_disableAllDisplay){
				globalSettings.handleKey( IOWrap::waitKey(1) );
			}
		}
#endif
		// Saves the frame as a keyframe or treats it as a non-keyframe
		if(needKF)
		{
            if(imuIntegration.setting_useIMU)
            {
		// imu!: Create keyframe in imu framework
                imuIntegration.keyframeCreated(fh->shell->id);
            }
            makeKeyFrame(fh);
		}
		else makeNonKeyFrame(fh);
	}
	else // wait for mapper to synchronize
	{
		boost::unique_lock<boost::mutex> lock(trackMapSyncMutex);
		unmappedTrackedFrames.push_back(fh);

		// imu!: If the prepared KF is still in the queue right now the current frame will become a KF instead.
		if(alreadyPreparedKF && !imuIntegration.isPreparedKFCreated())
		{
			imuIntegration.prepareKeyframe(fh->shell->id);
			needKF = true;
		}

		if(imuIntegration.setting_useIMU)
        {
            if(needKF) needNewKFAfter=imuIntegration.getPreparedKeyframe();
        }else
        {
             if(needKF) needNewKFAfter=fh->shell->trackingRef->id;
        }
		trackedFrameSignal.notify_all();

		while(coarseTracker_forNewKF->refFrameID == -1 && coarseTracker->refFrameID == -1 )
		{
			mappedFrameSignal.wait(lock);
		}

		lock.unlock();
	}
}

/**
 * @brief Map using the set of active frames
 * 
 * Mapping is done when a keyframe is created
 * Runs on its own thread
 * 
 */
void FullSystem::mappingLoop()
{
	boost::unique_lock<boost::mutex> lock(trackMapSyncMutex);

	while(runMapping)
	{
		while(unmappedTrackedFrames.size()==0) // do not run if there are no frames
		{
			trackedFrameSignal.wait(lock);
			if(!runMapping) return;
		}

		FrameHessian* fh = unmappedTrackedFrames.front();
		unmappedTrackedFrames.pop_front();

        if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        {
            std::cout << "Current mapping id: " << fh->shell->id << " create KF after: " << needNewKFAfter << std::endl;
        }

        // Guaranteed to make a KF for the very first two tracked frames.
		if(allKeyFramesHistory.size() <= 2)
		{
            if(imuIntegration.setting_useIMU)
            {
                imuIntegration.keyframeCreated(fh->shell->id);
            }
            lock.unlock();
			makeKeyFrame(fh);
			lock.lock();
			mappedFrameSignal.notify_all();
			continue;
		}

		if(unmappedTrackedFrames.size() > 3)
			needToKetchupMapping=true;

		// if there are other frames to track, do that first.
		if(unmappedTrackedFrames.size() > 0)
		{
			if(imuIntegration.setting_useIMU && needNewKFAfter == fh->shell->id)
			{
                if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
                {
                    std::cout << "WARNING: Prepared keyframe got skipped!" << std::endl;
                }
                imuIntegration.skipPreparedKeyframe();
				assert(false);
			}

			lock.unlock();
			makeNonKeyFrame(fh);
			lock.lock();

			if(needToKetchupMapping && unmappedTrackedFrames.size() > 0)
			{
				FrameHessian* fh_catch = unmappedTrackedFrames.front();
				unmappedTrackedFrames.pop_front();
				{
					boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
					assert(fh_catch->shell->trackingRef != 0);
					fh_catch->shell->camToWorld = fh_catch->shell->trackingRef->camToWorld * fh_catch->shell->camToTrackingRef;
					fh_catch->setEvalPT_scaled(fh_catch->shell->camToWorld.inverse(),fh_catch->shell->aff_g2l);
				}
				delete fh_catch;
			}

		}
		else // Do mapping
		{
		    bool createKF = imuIntegration.setting_useIMU ? needNewKFAfter==fh->shell->id : needNewKFAfter >= frameHessians.back()->shell->id;
			if(globalSettings.setting_realTimeMaxKF || createKF)
			{
                if(imuIntegration.setting_useIMU)
                {
			// imu: Make keyframe for imu framework
                    imuIntegration.keyframeCreated(fh->shell->id);
                }
                lock.unlock();
				makeKeyFrame(fh);
				needToKetchupMapping=false;
				lock.lock();
			}
			else
			{
				lock.unlock();
				makeNonKeyFrame(fh);
				lock.lock();
			}
		}
		mappedFrameSignal.notify_all();
	}
	printf("MAPPING FINISHED!\n");
}

void FullSystem::blockUntilMappingIsFinished()
{
	boost::unique_lock<boost::mutex> lock(trackMapSyncMutex);
	runMapping = false;
	trackedFrameSignal.notify_all();
	lock.unlock();

	mappingThread.join();
}

/**
 * @brief Creates a non-keyframe
 * 
 * Only sets the new poses and traces immature keypoints
 * 
 * @param fh 
 */
void FullSystem::makeNonKeyFrame( FrameHessian* fh)
{
    dmvio::TimeMeasurement timeMeasurement("makeNonKeyframe");
	// needs to be set by mapping thread. no lock required since we are in mapping thread.
	{
		boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
		assert(fh->shell->trackingRef != 0);
		fh->shell->camToWorld = fh->shell->trackingRef->camToWorld * fh->shell->camToTrackingRef;
		fh->setEvalPT_scaled(fh->shell->camToWorld.inverse(),fh->shell->aff_g2l);
	}

	traceNewCoarse(fh); // trace the immature points
	delete fh;
}

/**
 * @brief Create a keyframe and optimize all of the active frames
 * 
 * Steps:
 * 1. Sets new pose
 * 2. Traces the current immature points
 * 3. Flags frames for marginalization
 * 4. Creates new keyframe
 * 5. Adds current active points into the new keyframe
 * 6. Converts immature points in new keyframe into active points
 * 7. Do full map optimization
 * 8. Remove outlier points
 * 9. Set current keyframe as tracking frame
 * 10. Remove other bad points
 * 11. Make new immature points for keyframe
 * 12. Marginalize flagged frames
 * 
 * @param fh 
 */
void FullSystem::makeKeyFrame(FrameHessian* fh)
{
    dmvio::TimeMeasurement timeMeasurement("makeKeyframe");

	// needs to be set by mapping thread
	{
		boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
		assert(fh->shell->trackingRef != 0);
		fh->shell->camToWorld = fh->shell->trackingRef->camToWorld * fh->shell->camToTrackingRef;
		fh->setEvalPT_scaled(fh->shell->camToWorld.inverse(),fh->shell->aff_g2l);
		int prevKFId = fh->shell->trackingRef->id;
		int framesBetweenKFs = fh->shell->id - prevKFId - 1;
        if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        {
            std::cout << "Frames between KFs: " << framesBetweenKFs << std::endl;
        }
    }

	traceNewCoarse(fh); // trace the immature points

	boost::unique_lock<boost::mutex> lock(mapMutex);


	// =========================== Flag Frames to be Marginalized. =========================
	flagFramesForMarginalization(fh);


	// =========================== add New Frame to Hessian Struct. =========================
    dmvio::TimeMeasurement timeMeasurementAddFrame("newFrameAndNewResidualsForOldPoints");
	fh->idx = frameHessians.size();
	frameHessians.push_back(fh);
	fh->frameID = allKeyFramesHistory.size();
    fh->shell->keyframeId = fh->frameID;
	allKeyFramesHistory.push_back(fh->shell);
	ef->insertFrame(fh, &Hcalib);

	setPrecalcValues();


	// =========================== add new residuals for old points =========================
	int numFwdResAdde=0;
	for(FrameHessian* fh1 : frameHessians) // go through all active frames
	{
		if(fh1 == fh) continue;
		for(PointHessian* ph : fh1->pointHessians)
		{
			PointFrameResidual* r = new PointFrameResidual(std::addressof(globalCalib), ph, fh1, fh, std::addressof(globalSettings));
			r->setState(ResState::IN);
			ph->residuals.push_back(r);
			ef->insertResidual(r);
			ph->lastResiduals[1] = ph->lastResiduals[0];
			ph->lastResiduals[0] = std::pair<PointFrameResidual*, ResState>(r, ResState::IN);
			numFwdResAdde+=1;
		}
	}

    timeMeasurementAddFrame.end();


	// =========================== Activate Points (& flag for marginalization). =========================
	activatePointsMT();
	ef->makeIDX();

    if(imuIntegration.setting_useGTSAMIntegration)
    {
        // Adds new keyframe to the BA graph, together with matching factors (e.g. IMUFactors).
        baIntegration->addKeyframeToBA(fh->shell->id, fh->shell->camToWorld, ef->frames);
    }


	// =========================== OPTIMIZE ALL!!!!! =========================
	fh->frameEnergyTH = frameHessians.back()->frameEnergyTH;
	float rmse = optimize(globalSettings.setting_maxOptIterations);


	// =========================== Figure Out if INITIALIZATION FAILED =========================
	if(allKeyFramesHistory.size() <= 4)
	{
		if(allKeyFramesHistory.size()==2 && rmse > 20*globalSettings.benchmark_initializerSlackFactor)
		{
			printf("I THINK INITIALIZATINO FAILED! Resetting.\n");
			initFailed=true;
		}
		if(allKeyFramesHistory.size()==3 && rmse > 13*globalSettings.benchmark_initializerSlackFactor)
		{
			printf("I THINK INITIALIZATINO FAILED! Resetting.\n");
			initFailed=true;
		}
		if(allKeyFramesHistory.size()==4 && rmse > 9*globalSettings.benchmark_initializerSlackFactor)
		{
			printf("I THINK INITIALIZATINO FAILED! Resetting.\n");
			initFailed=true;
		}
	}


	// =========================== REMOVE OUTLIER =========================
	removeOutliers();

	// imu!: Handle imu post optimization
	if(imuIntegration.setting_useIMU)
    {
	    imuIntegration.postOptimization(fh->shell->id);
    }

    bool imuReady = false;
	{
        dmvio::TimeMeasurement timeMeasurement2("makeKeyframeChangeTrackingRef");
		boost::unique_lock<boost::mutex> crlock(coarseTrackerSwapMutex);

        if(imuIntegration.setting_useIMU)
        {
            imuReady = imuIntegration.finishKeyframeOptimization(fh->shell->id);
        }

        coarseTracker_forNewKF->makeK(&Hcalib);
		// Set the reference frame for coarse tracking
		coarseTracker_forNewKF->setCoarseTrackingRef(frameHessians);

        coarseTracker_forNewKF->debugPlotIDepthMap(&minIdJetVisTracker, &maxIdJetVisTracker, outputWrapper);
        coarseTracker_forNewKF->debugPlotIDepthMapFloat(outputWrapper);
	}

#ifdef GRAPHICAL_DEBUG
	debugPlot("post Optimize");
#endif

    for(auto* ow : outputWrapper)
    {
        if(imuReady && !imuUsedBefore)
        {
            // Update state if this is the first time after IMU init.
            // VIO is now initialized the next published scale will be useful.
            ow->publishSystemStatus(dmvio::VISUAL_INERTIAL);
        }
        ow->publishTransformDSOToIMU(imuIntegration.getTransformDSOToIMU());
    }
    imuUsedBefore = imuReady;


    // =========================== (Activate-)Marginalize Points =========================
    dmvio::TimeMeasurement timeMeasurementMarginalizePoints("marginalizeAndRemovePoints");
	flagPointsForRemoval();
	ef->dropPointsF();
	getNullspaces(
			ef->lastNullspaces_pose,
			ef->lastNullspaces_scale,
			ef->lastNullspaces_affA,
			ef->lastNullspaces_affB);
	ef->marginalizePointsF();
	timeMeasurementMarginalizePoints.end();


	// =========================== add new Immature points & new residuals =========================
	makeNewTraces(fh, 0);

    dmvio::TimeMeasurement timeMeasurementPublish("publishInMakeKeyframe");
    for(IOWrap::Output3DWrapper* ow : outputWrapper)
    {
        ow->publishGraph(ef->connectivityMap);
        ow->publishKeyframes(frameHessians, false, &Hcalib);
    }
    timeMeasurementPublish.end();


    // =========================== Marginalize Frames =========================
    dmvio::TimeMeasurement timeMeasurementMargFrames("marginalizeFrames");
	for(unsigned int i=0;i<frameHessians.size();i++)
		if(frameHessians[i]->flaggedForMarginalization)
        {
		    marginalizeFrame(frameHessians[i]);
		    i=0;
            if(imuIntegration.setting_useGTSAMIntegration)
            {
                baIntegration->updateBAOrdering(ef->frames);
            }
        }
	timeMeasurementMargFrames.end();

	printLogLine();
	printEigenValLine();

    if(imuIntegration.setting_useGTSAMIntegration)
    {
        baIntegration->updateBAValues(ef->frames);
    }

    if(imuIntegration.setting_useIMU)
    {
        imuIntegration.finishKeyframeOperations(fh->shell->id);
    }
}

/**
 * @brief Initializes the system if coarseInitializer finds good initial systems
 * 
 * @param newFrame 
 */
void FullSystem::initializeFromInitializer(FrameHessian* newFrame)
{
	boost::unique_lock<boost::mutex> lock(mapMutex);

    // Add firstframe.
	FrameHessian* firstFrame = coarseInitializer->firstFrame;
	firstFrame->idx = frameHessians.size();
	frameHessians.push_back(firstFrame);
	firstFrame->frameID = allKeyFramesHistory.size();
	allKeyFramesHistory.push_back(firstFrame->shell);
	ef->insertFrame(firstFrame, &Hcalib);
	setPrecalcValues();

	// imu!: Start BA with first frame
	baIntegration->addFirstBAFrame(firstFrame->shell->id);

	firstFrame->pointHessians.reserve(globalCalib.wG[0]*globalCalib.hG[0]);
	firstFrame->pointHessiansMarginalized.reserve(globalCalib.wG[0]*globalCalib.hG[0]);
	firstFrame->pointHessiansOut.reserve(globalCalib.wG[0]*globalCalib.hG[0]);


	float sumID=1e-5, numID=1e-5;
    double sumFirst = 0.0;
    double sumSecond = 0.0;
    int num = 0;
	for(int i=0;i<coarseInitializer->numPoints[0];i++)
	{
		sumID += coarseInitializer->points[0][i].iR;
		numID++;
	}
    sumFirst /= num;
    sumSecond /= num;

    float rescaleFactor = 1;
	rescaleFactor = 1 / (sumID / numID);
    SE3 firstToNew = coarseInitializer->thisToNext;
    if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
		std::cout << "Scaling with rescaleFactor: " << rescaleFactor << std::endl;
    firstToNew.translation() /= rescaleFactor;

	// randomly sub-select the points
	float keepPercentage = globalSettings.setting_desiredPointDensity / coarseInitializer->numPoints[0];

    if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        printf("Initialization: keep %.1f%% (need %d, have %d)!\n", 100*keepPercentage,
                (int)(globalSettings.setting_desiredPointDensity), coarseInitializer->numPoints[0] );

	for(int i=0;i<coarseInitializer->numPoints[0];i++)
	{
		if(rand()/(float)RAND_MAX > keepPercentage) continue;

		Pnt* point = coarseInitializer->points[0]+i;
		ImmaturePoint* pt = new ImmaturePoint(point->u+0.5f,point->v+0.5f, globalCalib, firstFrame,point->my_type, &Hcalib, globalSettings);

		if(!std::isfinite(pt->energyTH)) { delete pt; continue; }

		pt->idepth_max=pt->idepth_min=1;
		PointHessian* ph = new PointHessian(pt, &Hcalib, globalSettings);
		delete pt;
		if(!std::isfinite(ph->energyTH)) {delete ph; continue;}

        ph->setIdepthScaled(point->iR * rescaleFactor);
		ph->setIdepthZero(ph->idepth);
		ph->hasDepthPrior=true;
		ph->setPointStatus(PointHessian::ACTIVE);

		firstFrame->pointHessians.push_back(ph);
		ef->insertPoint(ph);
	}


	// really no lock required, as we are initializing.
	{
		boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
        firstFrame->shell->camToWorld = firstPose;
		firstFrame->shell->aff_g2l = AffLight(0,0);
		firstFrame->setEvalPT_scaled(firstFrame->shell->camToWorld.inverse(),firstFrame->shell->aff_g2l);
		firstFrame->shell->trackingRef=0;
		firstFrame->shell->camToTrackingRef = SE3();
		firstFrame->shell->keyframeId = 0;

		newFrame->shell->camToWorld = firstPose * firstToNew.inverse();
		newFrame->shell->aff_g2l = AffLight(0,0);
        newFrame->setEvalPT_scaled(newFrame->shell->camToWorld.inverse(),newFrame->shell->aff_g2l);
		newFrame->shell->trackingRef = firstFrame->shell;
		newFrame->shell->camToTrackingRef = firstToNew.inverse();

    }
    imuIntegration.finishCoarseTracking(*(newFrame->shell), true);

    initialized=true;
	if (!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
		printf("INITIALIZE FROM INITIALIZER (%d pts)!\n", (int)firstFrame->pointHessians.size());
}

/**
 * @brief Create new immature points
 * 
 * Only made for the newest keyframe
 * 
 * @param newFrame 
 * @param gtDepth 
 */
void FullSystem::makeNewTraces(FrameHessian* newFrame, float* gtDepth)
{
    dmvio::TimeMeasurement timeMeasurement("makeNewTraces");

	int numPointsTotal = pixelSelector->makeMaps(newFrame, selectionMap,globalSettings.setting_desiredImmatureDensity);

	newFrame->pointHessians.reserve(numPointsTotal*1.2f);
	//fh->pointHessiansInactive.reserve(numPointsTotal*1.2f);
	newFrame->pointHessiansMarginalized.reserve(numPointsTotal*1.2f);
	newFrame->pointHessiansOut.reserve(numPointsTotal*1.2f);

	// Makes immature points at the selectionMap points
	for(int y=PATTERNPADDING+1;y<globalCalib.hG[0]-PATTERNPADDING-2;y++)
	for(int x=PATTERNPADDING+1;x<globalCalib.wG[0]-PATTERNPADDING-2;x++)
	{
		int i = x+y*globalCalib.wG[0];
		if(selectionMap[i]==0) continue;

		ImmaturePoint* impt = new ImmaturePoint(x, y, globalCalib, newFrame, selectionMap[i], &Hcalib, globalSettings);
		if(!std::isfinite(impt->energyTH)) delete impt;
		else newFrame->immaturePoints.push_back(impt);

	}
	//printf("MADE %d IMMATURE POINTS!\n", (int)newFrame->immaturePoints.size());
}


/**
 * @brief Sets pre-calculation values for the active frames
 * 
 */
void FullSystem::setPrecalcValues()
{
	// For all active frames to all other active frames
	for(FrameHessian* fh : frameHessians)
	{
		fh->targetPrecalc.resize(frameHessians.size());
		for(unsigned int i=0;i<frameHessians.size();i++)
			fh->targetPrecalc[i].set(fh, frameHessians[i], &Hcalib);
	}
	ef->setDeltaF(&Hcalib);
}

}
