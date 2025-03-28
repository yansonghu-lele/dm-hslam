/**
* This file is part of DSO.
*
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
#include <algorithm>

#include "util/globalFuncs.h"
#include "util/globalCalib.h"
#include "util/TimeMeasurement.h"

#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

#include "FullSystem/CoarseTracker.h"
#include "FullSystem/ImmaturePoint.h"
#include "IOWrapper/ImageDisplay.h"

#include "math.h"



namespace dso
{

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
 * @brief Do optimization calculations for immature points that are being activated
 * 
 * For all of the active frames
 * - Starts by calculating an optimized depth value
 * - Checks the energy from the optimization to determine if the point should count as visible in a frame
 * - Creates a PointHessian struct for the activated point and frame with the immature point values and new depth
 * 
 * @param point 			List of immature points that are to be optimized
 * @param minObs 
 * @param residuals 		Residual calculations
 * @return PointHessian* 	Output array of activated points
 */
PointHessian* FullSystem::optimizeImmaturePoint(
		ImmaturePoint* point, int minObs,
		ImmaturePointTemporaryResidual* residuals)
{
	int nres = 0;
	for(FrameHessian* fh : frameHessians) // for all active frames
	{
		// Initialize variables for all connected frames
		if(fh != point->host)
		{
			residuals[nres].state_NewEnergy = residuals[nres].state_energy = 0;
			residuals[nres].state_NewState = ResState::OUTLIER;
			residuals[nres].state_state = ResState::IN;
			residuals[nres].target = fh;
			nres++;
		}
	}
	assert(nres == ((int)frameHessians.size())-1);

	bool print = false;

	float lastEnergy = 0;
	float lastHdd=0;
	float lastbd=0;
	float currentIdepth=(point->idepth_max+point->idepth_min)*0.5f; // Initial depth from immature point

	// Initial calculations
	for(int i=0;i<nres;i++)
	{
		// Calculate residual of immature point for all connected frames
		lastEnergy += point->linearizeResidual(&Hcalib, 1000, residuals+i, lastHdd, lastbd, currentIdepth);
		residuals[i].state_state = residuals[i].state_NewState;
		residuals[i].state_energy = residuals[i].state_NewEnergy;
	}

	if(!std::isfinite(lastEnergy) || lastHdd < globalSettings.setting_minIdepthH_act)
	{
		if(print)
			printf("OptPoint: Not well-constrained (%d res, H=%.1f). E=%f. SKIP!\n",
				nres, lastHdd, lastEnergy);
		return 0;
	}

	if(print) printf("Activate point. %d residuals. H=%f. Initial Energy: %f. Initial Id=%f\n" ,
			nres, lastHdd, lastEnergy, currentIdepth);

	// Optimize new values (depth) for activated point
	float lambda = 0.1;
	for(int iteration=0;iteration<globalSettings.setting_GNItsOnPointActivation;iteration++)
	{
		float H = lastHdd;
		H *= 1+lambda;
		float step = (1.0/H) * lastbd;
		float newIdepth = currentIdepth - step;

		float newHdd=0; float newbd=0; float newEnergy=0;
		for(int i=0;i<nres;i++)
			newEnergy += point->linearizeResidual(&Hcalib, 1, residuals+i, newHdd, newbd, newIdepth);

		if(!std::isfinite(lastEnergy) || newHdd < globalSettings.setting_minIdepthH_act)
		{
			if(print) printf("OptPoint: Not well-constrained (%d res, H=%.1f). E=%f. SKIP!\n",
					nres,
					newHdd,
					lastEnergy);
			return 0;
		}

		if(print) printf("%s %d (L %.2f) %s: %f -> %f (idepth %f)!\n",
				(true || newEnergy < lastEnergy) ? "ACCEPT" : "REJECT",
				iteration,
				log10(lambda),
				"",
				lastEnergy, newEnergy, newIdepth);

		if(newEnergy < lastEnergy) // use new values and increase step
		{
			currentIdepth = newIdepth;
			lastHdd = newHdd;
			lastbd = newbd;
			lastEnergy = newEnergy;
			for(int i=0;i<nres;i++)
			{
				residuals[i].state_state = residuals[i].state_NewState;
				residuals[i].state_energy = residuals[i].state_NewEnergy;
			}

			lambda *= 0.5;
		}
		else // decrease step
		{
			lambda *= 5;
		}

		if(fabsf(step) < 0.0001*currentIdepth)
			break;
	}

	if(!std::isfinite(currentIdepth))
	{
		printf("MAJOR ERROR! point idepth is nan after initialization (%f).\n", currentIdepth);
		return (PointHessian*)((long)(-1));		// yeah I'm like 99% sure this is OK on 32bit systems.
	}


	int numGoodRes=0;
	for(int i=0;i<nres;i++)
		if(residuals[i].state_state == ResState::IN) numGoodRes++;

	if(numGoodRes < minObs)
	{
		if(print) printf("OptPoint: OUTLIER!\n");
		return (PointHessian*)((long)(-1));		// 99% sure this is OK on 32bit systems.
	}


	// Set new PointHessian and PointFrameResidual structs for activated points
	PointHessian* p = new PointHessian(point, globalSettings);
	if(!std::isfinite(p->energyTH)) {delete p; return (PointHessian*)((long)(-1));}

	p->lastResiduals[0].first = 0;
	p->lastResiduals[0].second = ResState::OOB;
	p->lastResiduals[1].first = 0;
	p->lastResiduals[1].second = ResState::OOB;
	p->setIdepthZero(currentIdepth);
	p->setIdepth(currentIdepth);
	p->setPointStatus(PointHessian::ACTIVE);

	// Do all of the required optimization calculations for the new points
	for(int i=0;i<nres;i++)
		if(residuals[i].state_state == ResState::IN)
		{
			PointFrameResidual* r = new PointFrameResidual(std::addressof(globalCalib), p, p->host, residuals[i].target, std::addressof(globalSettings));
			r->state_NewEnergy = r->state_energy = 0;
			r->state_NewState = ResState::OUTLIER;
			r->setState(ResState::IN);
			p->residuals.push_back(r);

			if(r->target == frameHessians.back())
			{
				p->lastResiduals[0].first = r;
				p->lastResiduals[0].second = ResState::IN;
			}
			else if(r->target == (frameHessians.size()<2 ? 0 : frameHessians[frameHessians.size()-2]))
			{
				p->lastResiduals[1].first = r;
				p->lastResiduals[1].second = ResState::IN;
			}
		}

	if(print) printf("Point activated!\n");

	statistics_numActivatedPoints++;
	return p;
}
}
