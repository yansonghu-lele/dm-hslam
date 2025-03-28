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



#include <util/TimeMeasurement.h>
#include "FullSystem/FullSystem.h"

#include "stdio.h"
#include <algorithm>

#include "util/globalFuncs.h"
#include "util/globalCalib.h"

#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

#include "FullSystem/ResidualProjections.h"
#include "FullSystem/ImmaturePoint.h"

#include "OptimizationBackend/EnergyFunctional.h"
#include "OptimizationBackend/EnergyFunctionalStructs.h"

#include "IOWrapper/Output3DWrapper.h"
#include "IOWrapper/ImageDisplay.h"

#include "FullSystem/CoarseTracker.h"



namespace dso
{


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
 * @brief Flags frames for marginalization
 * 
 * Conditions for marginalization:
 * 1. Always keep latest keyframe
 * 2. Marginalize if less than 5% of points are visable
 * 3. If the frame window is full, marginalize the frame with a highest distance score
 * 
 * @param newFH 
 */
void FullSystem::flagFramesForMarginalization(FrameHessian* newFH)
{
    dmvio::TimeMeasurement timeMeasurement("flagFramesForMarginalization");

	if(globalSettings.setting_minFrameAge > globalSettings.setting_maxFrames)
	{
		for(int i=globalSettings.setting_maxFrames;i<(int)frameHessians.size();i++)
		{
			FrameHessian* fh = frameHessians[i-globalSettings.setting_maxFrames];
			fh->flaggedForMarginalization = true;
		}
		return;
	}


	int flagged = 0;
	// marginalize all frames that have not enough points.
	for(int i=0;i<(int)frameHessians.size();i++)
	{
		FrameHessian* fh = frameHessians[i];
		int in = fh->pointHessians.size() + fh->immaturePoints.size();
		int out = fh->pointHessiansMarginalized.size() + fh->pointHessiansOut.size();


		Vec2 refToFh=AffLight::fromToVecExposure(frameHessians.back()->ab_exposure, fh->ab_exposure,
				frameHessians.back()->aff_g2l(), fh->aff_g2l());


		if( (in < globalSettings.setting_minPointsRemaining *(in+out) || fabs(logf((float)refToFh[0])) > globalSettings.setting_maxLogAffFacInWindow)
				&& ((int)frameHessians.size())-flagged > globalSettings.setting_minFrames)
		{
//			printf("MARGINALIZE frame %d, as only %'d/%'d points remaining (%'d %'d %'d %'d). VisInLast %'d / %'d. traces %d, activated %d!\n",
//					fh->frameID, in, in+out,
//					(int)fh->pointHessians.size(), (int)fh->immaturePoints.size(),
//					(int)fh->pointHessiansMarginalized.size(), (int)fh->pointHessiansOut.size(),
//					visInLast, outInLast,
//					fh->statistics_tracesCreatedForThisFrame, fh->statistics_pointsActivatedForThisFrame);
			fh->flaggedForMarginalization = true;
			flagged++;
		}
		else
		{
//			printf("May Keep frame %d, as %'d/%'d points remaining (%'d %'d %'d %'d). VisInLast %'d / %'d. traces %d, activated %d!\n",
//					fh->frameID, in, in+out,
//					(int)fh->pointHessians.size(), (int)fh->immaturePoints.size(),
//					(int)fh->pointHessiansMarginalized.size(), (int)fh->pointHessiansOut.size(),
//					visInLast, outInLast,
//					fh->statistics_tracesCreatedForThisFrame, fh->statistics_pointsActivatedForThisFrame);
		}
	}

	// Marginalize one if active window size is too large
	// Chose frame to be marginalized using the distance score
	// The distance score tries to keep keyframes well distributed in 3D space
	if((int)frameHessians.size()-flagged >= globalSettings.setting_maxFrames)
	{
		double smallestScore = 1;
		FrameHessian* toMarginalize=0;
		FrameHessian* latest = frameHessians.back();


		for(FrameHessian* fh : frameHessians)
		{
			if(fh->frameID > latest->frameID-globalSettings.setting_minFrameAge || fh->frameID == 0) continue;
			//if(fh==frameHessians.front() == 0) continue;

			// Calculate distance score
			double distScore = 0;
			for(FrameFramePrecalc &ffh : fh->targetPrecalc)
			{
				if(ffh.target->frameID > latest->frameID-globalSettings.setting_minFrameAge+1 || ffh.target == ffh.host) continue;
				distScore += 1/(1e-5+ffh.distanceLL);
			}
			distScore *= -sqrtf(fh->targetPrecalc.back().distanceLL);


			if(distScore < smallestScore)
			{
				smallestScore = distScore;
				toMarginalize = fh;
			}
		}

//		printf("MARGINALIZE frame %d, as it is the closest (score %.2f)!\n",
//				toMarginalize->frameID, smallestScore);
		toMarginalize->flaggedForMarginalization = true;
		flagged++;
	}

//	printf("FRAMES LEFT: ");
//	for(FrameHessian* fh : frameHessians)
//		printf("%d ", fh->frameID);
//	printf("\n");
}

/**
 * @brief Create Point Cloud output for Point
 * 
 * @param p 
 * @param camToWorld 
 * @return PC_output 
 */
PC_output FullSystem::createPCOutput(Point* p, SE3 camToWorld)
{
	PC_output tmp_point_output;

	Eigen::Vector3d worldPoint = p->getWorldPosition(1 / Hcalib.fxl(), 1 / Hcalib.fyl(), 
								-Hcalib.cxl() / Hcalib.fxl(), -Hcalib.cyl() / Hcalib.fyl(), 
								camToWorld);
	tmp_point_output.x = worldPoint[0];
	tmp_point_output.y = worldPoint[1];
	tmp_point_output.z = worldPoint[2];

	Eigen::Vector3f color = p->getColourRGBfloat();
	tmp_point_output.r = color[0];
	tmp_point_output.g = color[1];
	tmp_point_output.b = color[2];

	return tmp_point_output;
}

/**
 * @brief Marginalizes frame
 * 
 * Marginalizes the points in the frame and the frame itself
 * 
 * @param frame 
 */
void FullSystem::marginalizeFrame(FrameHessian* frame)
{
    dmvio::TimeMeasurement timeMeasurement("marginalizeFrame");
	// marginalize or remove all this frames points.

	assert((int)frame->pointHessians.size()==0);

	// This actually does the marginalization math!!!
	ef->marginalizeFrame(frame->efFrame, imuIntegration->setting_useGTSAMIntegration);

	dmvio::TimeMeasurement timeMeasurementEnd("marginalizeFrameOverhead");

	// Marginalize all the points represented in the frame
	// by dropping all observations of existing points in that frame
	for(FrameHessian* fh : frameHessians)
	{
		if(fh==frame) continue;

		for(PointHessian* ph : fh->pointHessians)
		{
			for(unsigned int i=0;i<ph->residuals.size();i++)
			{
				PointFrameResidual* r = ph->residuals[i];
				if(r->target == frame)
				{
					if(ph->lastResiduals[0].first == r)
						ph->lastResiduals[0].first=0;
					else if(ph->lastResiduals[1].first == r)
						ph->lastResiduals[1].first=0;


					if(r->host->frameID < r->target->frameID)
						statistics_numForceDroppedResFwd++;
					else
						statistics_numForceDroppedResBwd++;

					ef->dropResidual(r->efResidual);
					deleteOut<PointFrameResidual>(ph->residuals,i);
					break;
				}
			}
		}
	}

    // ef->marginalizeFrame does not delete efFrame anymore, because it's field frameID was needed in dropResidual.
    delete frame->efFrame;
	frame->efFrame = nullptr;

	{
        std::vector<FrameHessian*> v;
        v.push_back(frame);
        for(IOWrap::Output3DWrapper* ow : outputWrapper){
            ow->publishKeyframes(v, true, &Hcalib);
		}
	} 

	frame->shell->marginalizedAt = frameHessians.back()->shell->id;
	frame->shell->movedByOpt = frame->w2c_leftEps().norm();

	auto frameID = frame->frameID;

	// Only marginalized points are included in the final output point cloud
	if (globalSettings.setting_outputPC){
		for(PointHessian* p : frame->pointHessiansMarginalized){
			if (allMargPointsHistory.find(p->point_id) == allMargPointsHistory.end())
					allMargPointsHistory[p->point_id] = createPCOutput(p, firstPose.inverse() * frame->shell->camToWorld);
		}
	}


	deleteOutOrder<FrameHessian>(frameHessians, frame);
	for(unsigned int i=0;i<frameHessians.size();i++)
		frameHessians[i]->idx = i;


	// Culling the connectivity map is incompatible with viewing full constraints in GUI
	// Culled connectivity map information is unneeded for operation and culling results in performance improvement
	if (globalSettings.setting_disableAllDisplay){
		int numDel = 0;
		for(auto it = ef->connectivityMap.begin(); it != ef->connectivityMap.end();)
		{
			int host = (int)(it->first >> 32);
			int target = (int)(it->first & (uint64_t)0xFFFFFFFF);
			if(host == frameID || target == frameID)
			{
				numDel++;
				it = ef->connectivityMap.erase(it);
			}else
			{
				it++;
			}
		}
	}

    setPrecalcValues();
	ef->setAdjointsF(&Hcalib);
}

void FullSystem::cleanFrame(FrameHessian* frame)
{
	// Only marginalized points are included in the final output point cloud
	if (globalSettings.setting_outputPC){
		for(PointHessian* p : frame->pointHessiansMarginalized){
			if (allMargPointsHistory.find(p->point_id) == allMargPointsHistory.end())
					allMargPointsHistory[p->point_id] = createPCOutput(p,firstPose.inverse() * frame->shell->camToWorld);
		}
	}
}

}
