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



#include "FullSystem/FullSystem.h"

#include "stdio.h"
#include "util/globalFuncs.h"
#include <Eigen/LU>
#include <algorithm>
#include "IOWrapper/ImageDisplay.h"
#include "util/globalCalib.h"
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include "FullSystem/ResidualProjections.h"

#include "OptimizationBackend/EnergyFunctional.h"
#include "OptimizationBackend/EnergyFunctionalStructs.h"
#include "util/TimeMeasurement.h"

#include <cmath>

#include <algorithm>



namespace dso
{

/**
 * @brief Actually linearizes all activeResiduals
 * 
 * Assists linearizeAll
 * 
 * @param fixLinearization 
 * @param toRemove 
 * @param min 
 * @param max 
 * @param stats 
 * @param tid 
 */
void FullSystem::linearizeAll_Reductor(bool fixLinearization, std::vector<PointFrameResidual*>* toRemove, int min, int max, Vec10* stats, int tid)
{
	for(int k=min;k<max;k++)
	{
		PointFrameResidual* r = activeResiduals[k];

		// The linearization is done by the residual class here!!!
		(*stats)[0] += r->linearize(&Hcalib);

		if(fixLinearization)
		{
			r->applyRes(true);

			if(r->efResidual->isActive())
			{
				PointHessian* p = r->point;
				Vec3f ptp_inf = r->host->targetPrecalc[r->target->idx].PRE_KRKiTll * Vec3f(p->u,p->v, 1);	// projected point assuming infinite depth.
				Vec3f ptp = ptp_inf + r->host->targetPrecalc[r->target->idx].PRE_KtTll*p->idepth_scaled;	// projected point with real depth.
				float relBS = 0.01*((ptp_inf.head<2>() / ptp_inf[2])-(ptp.head<2>() / ptp[2])).norm();		// 0.01 = one pixel.

				if(relBS > p->maxRelBaseline)
					p->maxRelBaseline = relBS;

				p->numGoodResiduals++;
			}
			else
			{
				toRemove[tid].push_back(activeResiduals[k]);
			}
		}
	}
}

/**
 * @brief Makes all of the active residuals set their new states
 * 
 * @param copyJacobians 
 * @param min 
 * @param max 
 * @param stats 
 * @param tid 
 */
void FullSystem::applyRes_Reductor(bool copyJacobians, int min, int max, Vec10* stats, int tid)
{
	for(int k=min;k<max;k++)
		activeResiduals[k]->applyRes(true);
}

/**
 * @brief Sets new frame energy threshold
 * 
 */
void FullSystem::setNewFrameEnergyTH()
{
	// collect all residuals and make decision on TH.
	allResVec.clear();
	allResVec.reserve(activeResiduals.size()*2);
	FrameHessian* newFrame = frameHessians.back();

	for(PointFrameResidual* r : activeResiduals)
		if(r->state_NewEnergyWithOutlier >= 0 && r->target == newFrame)
		{
			allResVec.push_back(r->state_NewEnergyWithOutlier);

		}

	if(allResVec.size()==0) // should never happen, but lets make sure
	{
		newFrame->frameEnergyTH = 12*12*PATTERNNUM;
		return;
	}


	int nthIdx = globalSettings.setting_frameEnergyTHN*allResVec.size();

	assert(nthIdx < (int)allResVec.size());
	assert(globalSettings.setting_frameEnergyTHN < 1);

	std::nth_element(allResVec.begin(), allResVec.begin()+nthIdx, allResVec.end());
	float nthElement = sqrtf(allResVec[nthIdx]);


    newFrame->frameEnergyTH = nthElement*globalSettings.setting_frameEnergyTHFacMedian;
	newFrame->frameEnergyTH = 26.0f*globalSettings.setting_frameEnergyTHConstWeight + newFrame->frameEnergyTH*(1-globalSettings.setting_frameEnergyTHConstWeight);
	newFrame->frameEnergyTH = newFrame->frameEnergyTH*newFrame->frameEnergyTH;
	newFrame->frameEnergyTH *= globalSettings.setting_overallEnergyTHWeight*globalSettings.setting_overallEnergyTHWeight;

	// imu!: Used to enforce a maximum energy threshold.
	if(imuIntegration.setting_useIMU)
    {
	    imuIntegration.newFrameEnergyTH(newFrame->frameEnergyTH);
    }

//
//	int good=0,bad=0;
//	for(float f : allResVec) if(f<newFrame->frameEnergyTH) good++; else bad++;
//	printf("EnergyTH: mean %f, median %f, result %f (in %d, out %d)! \n",
//			meanElement, nthElement, sqrtf(newFrame->frameEnergyTH),
//			good, bad);

}

/**
 * @brief Linearizes all activeResiduals and other steps
 * 
 * Assists optimize
 * 
 * @param fixLinearization 
 * @return Vec3 
 */
Vec3 FullSystem::linearizeAll(bool fixLinearization)
{
	double lastEnergyP = 0;
	double lastEnergyR = 0;
	double num = 0;


	std::vector<PointFrameResidual*> toRemove[NUM_THREADS];
	for(int i=0;i<NUM_THREADS;i++) toRemove[i].clear();

	// Linearize all of the residuals
	if(!globalSettings.settings_no_multiThreading)
	{
		treadReduce.reduce(boost::bind(&FullSystem::linearizeAll_Reductor, this, fixLinearization, toRemove, _1, _2, _3, _4), 0, activeResiduals.size(), 0);
		lastEnergyP = treadReduce.stats[0];
	}
	else
	{
		Vec10 stats;
		linearizeAll_Reductor(fixLinearization, toRemove, 0,activeResiduals.size(),&stats,0);
		lastEnergyP = stats[0];
	}


	// Set new frame energy threshold
	setNewFrameEnergyTH();


	if(fixLinearization)
	{
		for(PointFrameResidual* r : activeResiduals)
		{
			PointHessian* ph = r->point;
			if(ph->lastResiduals[0].first == r)
				ph->lastResiduals[0].second = r->state_state;
			else if(ph->lastResiduals[1].first == r)
				ph->lastResiduals[1].second = r->state_state;
		}

		// Remove residuals
		int nResRemoved=0;
		for(int i=0;i<NUM_THREADS;i++)
		{
			for(PointFrameResidual* r : toRemove[i])
			{
				PointHessian* ph = r->point;

				if(ph->lastResiduals[0].first == r)
					ph->lastResiduals[0].first=0;
				else if(ph->lastResiduals[1].first == r)
					ph->lastResiduals[1].first=0;

				for(unsigned int k=0; k<ph->residuals.size();k++)
					if(ph->residuals[k] == r)
					{
						ef->dropResidual(r->efResidual);
						deleteOut<PointFrameResidual>(ph->residuals,k);
						nResRemoved++;
						break;
					}
			}
		}
		//printf("FINAL LINEARIZATION: removed %d / %d residuals!\n", nResRemoved, (int)activeResiduals.size());
	}

	return Vec3(lastEnergyP, lastEnergyR, num);
}


/**
 * @brief Does step from backup
 * 
 * Applies step to linearization point
 * 
 * @param stepfacC 
 * @param stepfacT 
 * @param stepfacR 
 * @param stepfacA 
 * @param stepfacD 
 * @return true 
 * @return false 
 */
bool FullSystem::doStepFromBackup(float stepfacC,float stepfacT,float stepfacR,float stepfacA,float stepfacD)
{
//	float meanStepC=0,meanStepP=0,meanStepD=0;
//	meanStepC += Hcalib.step.norm();

	Vec10 pstepfac;
	pstepfac.segment<3>(0).setConstant(stepfacT);
	pstepfac.segment<3>(3).setConstant(stepfacR);
	pstepfac.segment<4>(6).setConstant(stepfacA);


	float sumA=0, sumB=0, sumT=0, sumR=0, numID=0;

	float sumNID=0;

	if(globalSettings.setting_solverMode & SOLVER_MOMENTUM)
	{
		Hcalib.setValue(Hcalib.value_backup + Hcalib.step);
		for(FrameHessian* fh : frameHessians)
		{
			Vec10 step = fh->step;
			step.head<6>() += 0.5f*(fh->step_backup.head<6>());

			fh->setState(fh->state_backup + step);
			sumA += step[6]*step[6];
			sumB += step[7]*step[7];
			sumT += step.segment<3>(0).squaredNorm();
			sumR += step.segment<3>(3).squaredNorm();

			for(PointHessian* ph : fh->pointHessians)
			{
				float step_ph = ph->step+0.5f*(ph->step_backup);
				ph->setIdepth(ph->idepth_backup + step_ph);
				sumNID += fabsf(ph->idepth_backup);
				numID++;

                ph->setIdepthZero(ph->idepth_backup + step_ph);
			}
		}
	}
	else
	{
		Hcalib.setValue(Hcalib.value_backup + stepfacC*Hcalib.step);
		for(FrameHessian* fh : frameHessians)
		{
			fh->setState(fh->state_backup + pstepfac.cwiseProduct(fh->step));
			sumA += fh->step[6]*fh->step[6];
			sumB += fh->step[7]*fh->step[7];
			sumT += fh->step.segment<3>(0).squaredNorm();
			sumR += fh->step.segment<3>(3).squaredNorm();

			for(PointHessian* ph : fh->pointHessians)
			{
				ph->setIdepth(ph->idepth_backup + stepfacD*ph->step);
				sumNID += fabsf(ph->idepth_backup);
				numID++;

                ph->setIdepthZero(ph->idepth_backup + stepfacD*ph->step);
			}
		}
	}

	sumA /= frameHessians.size();
	sumB /= frameHessians.size();
	sumR /= frameHessians.size();
	sumT /= frameHessians.size();
	sumNID /= numID;

    if(!setting_debugout_runquiet)
        printf("STEPS: A %.1f; B %.1f; R %.1f; T %.1f. \t",
                sqrtf(sumA) / (0.0005*globalSettings.setting_thOptIterations),
                sqrtf(sumB) / (0.00005*globalSettings.setting_thOptIterations),
                sqrtf(sumR) / (0.00005*globalSettings.setting_thOptIterations),
                sqrtf(sumT)*sumNID / (0.00005*globalSettings.setting_thOptIterations));

	EFDeltaValid=false;
	setPrecalcValues();

	return sqrtf(sumA) < 0.0005*globalSettings.setting_thOptIterations &&
			sqrtf(sumB) < 0.00005*globalSettings.setting_thOptIterations &&
			sqrtf(sumR) < 0.00005*globalSettings.setting_thOptIterations &&
			sqrtf(sumT)*sumNID < 0.00005*globalSettings.setting_thOptIterations;
//
//	printf("mean steps: %f %f %f!\n",
//			meanStepC, meanStepP, meanStepD);
}

// sets linearization point.
/**
 * @brief Set up backup state if the optimization goes too far
 * 
 * @param backupLastStep 
 */
void FullSystem::backupState(bool backupLastStep)
{
	if(globalSettings.setting_solverMode & SOLVER_MOMENTUM)
	{
		if(backupLastStep)
		{
			Hcalib.step_backup = Hcalib.step;
			Hcalib.value_backup = Hcalib.value;
			for(FrameHessian* fh : frameHessians)
			{
				fh->step_backup = fh->step;
				fh->state_backup = fh->get_state();
				for(PointHessian* ph : fh->pointHessians)
				{
					ph->idepth_backup = ph->idepth;
					ph->step_backup = ph->step;
				}
			}
		}
		else
		{
			Hcalib.step_backup.setZero();
			Hcalib.value_backup = Hcalib.value;
			for(FrameHessian* fh : frameHessians)
			{
				fh->step_backup.setZero();
				fh->state_backup = fh->get_state();
				for(PointHessian* ph : fh->pointHessians)
				{
					ph->idepth_backup = ph->idepth;
					ph->step_backup=0;
				}
			}
		}
	}
	else
	{
		Hcalib.value_backup = Hcalib.value;
		for(FrameHessian* fh : frameHessians)
		{
			fh->state_backup = fh->get_state();
			for(PointHessian* ph : fh->pointHessians)
				ph->idepth_backup = ph->idepth;
		}
	}
}

// sets linearization point.
/**
 * @brief Loads backup values
 * 
 */
void FullSystem::loadSateBackup()
{
	Hcalib.setValue(Hcalib.value_backup);
	for(FrameHessian* fh : frameHessians)
	{
		fh->setState(fh->state_backup);
		for(PointHessian* ph : fh->pointHessians)
		{
			ph->setIdepth(ph->idepth_backup);

            ph->setIdepthZero(ph->idepth_backup);
		}

	}

	EFDeltaValid=false;
	setPrecalcValues();
}

/**
 * @brief Wrapper for calcMEnergyF
 * 
 * @param useNewValues 
 * @return double 
 */
double FullSystem::calcMEnergy(bool useNewValues)
{
	if(globalSettings.setting_forceAceptStep) return 0;

	// calculate (x-x0)^T * [2b + H * (x-x0)] for everything saved in L.
	//ef->makeIDX();
	//ef->setDeltaF(&Hcalib);

	return ef->calcMEnergyF(useNewValues, imuIntegration.setting_useGTSAMIntegration);

}


void FullSystem::printOptRes(const Vec3 &res, double resL, double resM, double resPrior, double LExact, float a, float b)
{
	printf("A(%f)=(AV %.3f). Num: A(%'d) + M(%'d); ab %f %f!\n",
			res[0],
			sqrtf((float)(res[0] / (PATTERNNUM*ef->resInA))),
			ef->resInA,
			ef->resInM,
			a,
			b
	);

}

/**
 * @brief Optimizes the window
 * 
 * Main function for full optimization
 * 
 * @param mnumOptIts 
 * @return float 
 */
float FullSystem::optimize(int mnumOptIts)
{
    dmvio::TimeMeasurement timeMeasurement("FullSystemOptimize");

	if(frameHessians.size() < 2) return 0;
	if(frameHessians.size() < 3) mnumOptIts = 20;
	if(frameHessians.size() < 4) mnumOptIts = 15;

	// get statistics and active residuals.

	activeResiduals.clear();
	int numPoints = 0;
	int numLRes = 0;
	// for all points in active frames
	for(FrameHessian* fh : frameHessians)
		for(PointHessian* ph : fh->pointHessians)
		{
			// for all frames the point is visible in
			for(PointFrameResidual* r : ph->residuals)
			{
				if(!r->efResidual->isLinearized)
				{
					activeResiduals.push_back(r);
					r->resetOOB();
				}
				else
					numLRes++;
			}
			numPoints++;
		}

    if(!setting_debugout_runquiet)
        printf("OPTIMIZE %d pts, %d active res, %d lin res!\n",ef->nPoints,(int)activeResiduals.size(), numLRes);


	// Initial calculation
	// Do optimization process
	Vec3 lastEnergy = linearizeAll(false);
	double lastEnergyL = calcLEnergy();
	double lastEnergyM = calcMEnergy(false);

	// Set new states
	if(!globalSettings.settings_no_multiThreading)
		treadReduce.reduce(boost::bind(&FullSystem::applyRes_Reductor, this, true, _1, _2, _3, _4), 0, activeResiduals.size(), 50);
	else
		applyRes_Reductor(true,0,activeResiduals.size(),0,0);

    if(!setting_debugout_runquiet)
    {
        printf("Initial Error \t");
        printOptRes(lastEnergy, lastEnergyL, lastEnergyM, 0, 0, frameHessians.back()->aff_g2l().a, frameHessians.back()->aff_g2l().b);
    }

	debugPlotTracking();

    double dynamicGTSAMWeight = 1.0;
    double minLambda = 1e-5;

//	double lambda = 1e-1;
    double lambda = minLambda;
	float stepsize=1;
	VecX previousX = VecX::Constant(CPARS+ 8*frameHessians.size(), NAN);

	// Iterative optimization process
	int numIterations = 0;
	for(int iteration=0;iteration<mnumOptIts;iteration++)
	{
	    dmvio::TimeMeasurement timeMeasurement_("baIteration");
		// Set up backup
		backupState(iteration!=0);

		// Solve!
	// imu!: Update the dynamic weight
        if(imuIntegration.getImuSettings().updateDynamicWeightDuringOptimization || iteration==0)
        {
            // Update dynamic weight before solving (where the active DSO factor will be scaled accordingly).
            dynamicGTSAMWeight = baIntegration->updateDynamicWeight(lastEnergy[0], sqrtf((float) (lastEnergy[0] /
                                                                                                  (PATTERNNUM *
                                                                                                   ef->resInA))),
                                                                    frameHessians.back()->shell->trackingWasGood);
            if(!setting_debugout_runquiet)
            {
                std::cout << "Dynamic weight: " << dynamicGTSAMWeight << std::endl;
            }
        }

		// Solves the Hessian
        solveSystem(iteration, lambda);

		// Increment optimization
		double incDirChange = (1e-20 + previousX.dot(ef->lastX)) / (1e-20 + previousX.norm() * ef->lastX.norm());
		previousX = ef->lastX;

		if(std::isfinite(incDirChange) && (globalSettings.setting_solverMode & SOLVER_STEPMOMENTUM))
		{
			float newStepsize = exp(incDirChange*1.4);
			if(incDirChange<0 && stepsize>1) stepsize=1;

			stepsize = sqrtf(sqrtf(newStepsize*stepsize*stepsize*stepsize));
			if(stepsize > 2) stepsize=2;
			if(stepsize <0.25) stepsize=0.25;
		}

		// Check if the optimization is good enough to stop
		bool canbreak = doStepFromBackup(stepsize,stepsize,stepsize,stepsize,stepsize);

		canbreak = canbreak && baIntegration->canBreak();


		// Eval new energy!
		// Do optimization process
		Vec3 newEnergy = linearizeAll(false);
		double newEnergyL = calcLEnergy();
		double newEnergyM = calcMEnergy(true);

	// imu!: Update dynamic weight before deciding whether to accept the step
		if(imuIntegration.getImuSettings().updateDynamicWeightDuringOptimization)
        {
            dynamicGTSAMWeight = baIntegration->updateDynamicWeight(lastEnergy[0], sqrtf((float)(lastEnergy[0] / (PATTERNNUM*ef->resInA))), frameHessians.back()->shell->trackingWasGood);
        }

        if(!setting_debugout_runquiet)
        {
            printf("%s %d (L %.2f, dir %.2f, ss %.1f): \t",
				(newEnergy[0] +  newEnergy[1] +  newEnergyL + newEnergyM / dynamicGTSAMWeight <
						lastEnergy[0] + lastEnergy[1] + lastEnergyL + lastEnergyM / dynamicGTSAMWeight) ? "ACCEPT" : "REJECT",
				iteration,
				log10(lambda),
				incDirChange,
				stepsize);
            printOptRes(newEnergy, newEnergyL, newEnergyM , 0, 0, frameHessians.back()->aff_g2l().a, frameHessians.back()->aff_g2l().b);
        }

		if(globalSettings.setting_forceAceptStep || (newEnergy[0] +  newEnergy[1] +  newEnergyL + newEnergyM / dynamicGTSAMWeight <
				lastEnergy[0] + lastEnergy[1] + lastEnergyL + lastEnergyM / dynamicGTSAMWeight))
		{
			if(!globalSettings.settings_no_multiThreading)
				treadReduce.reduce(boost::bind(&FullSystem::applyRes_Reductor, this, true, _1, _2, _3, _4), 0, activeResiduals.size(), 50);
			else
				applyRes_Reductor(true,0,activeResiduals.size(),0,0);

			lastEnergy = newEnergy;
			lastEnergyL = newEnergyL;
			lastEnergyM = newEnergyM;

			lambda *= 0.25;
            lambda = std::max(lambda, minLambda);

			// imu!: Accept the update to the bundle adjustment
			if(imuIntegration.setting_useGTSAMIntegration)
			{
				baIntegration->acceptBAUpdate(lastEnergy[0]);
			}
		}
		else // Undo step
		{
			loadSateBackup();

			// Do optimization process
			lastEnergy = linearizeAll(false);
			lastEnergyL = calcLEnergy();
			lastEnergyM = calcMEnergy(false);

			lambda *= 1e2;
		}
		numIterations++;

		if(canbreak && iteration >= globalSettings.setting_minOptIterations) break;
	}

    if(!setting_debugout_runquiet)
    {
        std::cout << "Num BA Iterations done: " << numIterations << "\n";
    }


    // Update IMU again
    baIntegration->updateDynamicWeight(lastEnergy[0], sqrtf((float)(lastEnergy[0] / (PATTERNNUM*ef->resInA))), frameHessians.back()->shell->trackingWasGood);

    Vec10 newStateZero = Vec10::Zero();
	newStateZero.segment<2>(6) = frameHessians.back()->get_state().segment<2>(6);

	// Set new poses
	frameHessians.back()->setEvalPT(frameHessians.back()->PRE_worldToCam,
			newStateZero);

	EFDeltaValid=false;
	EFAdjointsValid=false;
	ef->setAdjointsF(&Hcalib);

	setPrecalcValues();

	// Do a fixed linearization
	lastEnergy = linearizeAll(true);


	if(!std::isfinite((double)lastEnergy[0]) || !std::isfinite((double)lastEnergy[1]) || !std::isfinite((double)lastEnergy[2]))
	{
		std::cout << "Tracking lost after bundle adjustment!" << std::endl;
		isLost = true;
	}


	statistics_lastFineTrackRMSE = sqrtf((float)(lastEnergy[0] / (PATTERNNUM*ef->resInA)));

	if(calibLog != 0)
	{
		(*calibLog) << Hcalib.value_scaled.transpose() <<
				" " << frameHessians.back()->get_state_scaled().transpose() <<
				" " << sqrtf((float)(lastEnergy[0] / (PATTERNNUM*ef->resInA))) <<
				" " << ef->resInM << "\n";
		calibLog->flush();
	}

	{
		boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
		// Set new poses
		for(FrameHessian* fh : frameHessians)
		{
			fh->shell->camToWorld = fh->PRE_camToWorld;
			fh->shell->aff_g2l = fh->aff_g2l();
		}
	}

	// imu!: Post optimization for BA
    baIntegration->postOptimization(ef->frames);

	debugPlotTracking();

	return statistics_lastFineTrackRMSE;
}


/**
 * @brief Wrapper for solveSystemF
 * 
 * @param iteration 
 * @param lambda 
 */
void FullSystem::solveSystem(int iteration, double lambda)
{
	ef->lastNullspaces_forLogging = getNullspaces(
			ef->lastNullspaces_pose,
			ef->lastNullspaces_scale,
			ef->lastNullspaces_affA,
			ef->lastNullspaces_affB);

	ef->solveSystemF(iteration, lambda,&Hcalib,imuIntegration.setting_useGTSAMIntegration);
}


/**
 * @brief Wrapper for calcLEnergyF_MT
 * 
 * @return double 
 */
double FullSystem::calcLEnergy()
{
	if(globalSettings.setting_forceAceptStep) return 0;

	double Ef = ef->calcLEnergyF_MT();
	return Ef;
}


/**
 * @brief Removes outliers
 * 
 */
void FullSystem::removeOutliers()
{
    dmvio::TimeMeasurement timeMeasurement("removeOutliers");
	int numPointsDropped=0;
	for(FrameHessian* fh : frameHessians)
	{
		for(unsigned int i=0;i<fh->pointHessians.size();i++)
		{
			PointHessian* ph = fh->pointHessians[i];
			if(ph==0) continue;

			// Remove point
			if(ph->residuals.size() == 0)
			{
				fh->pointHessiansOut.push_back(ph);
				ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
				fh->pointHessians[i] = fh->pointHessians.back();
				fh->pointHessians.pop_back();
				i--;
				numPointsDropped++;
			}
		}
	}
	ef->dropPointsF();
}


/**
 * @brief Get nullspaces
 * 
 * @param nullspaces_pose 
 * @param nullspaces_scale 
 * @param nullspaces_affA 
 * @param nullspaces_affB 
 * @return std::vector<VecX> 
 */
std::vector<VecX> FullSystem::getNullspaces(
		std::vector<VecX> &nullspaces_pose,
		std::vector<VecX> &nullspaces_scale,
		std::vector<VecX> &nullspaces_affA,
		std::vector<VecX> &nullspaces_affB)
{
	nullspaces_pose.clear();
	nullspaces_scale.clear();
	nullspaces_affA.clear();
	nullspaces_affB.clear();


	int n=CPARS+frameHessians.size()*8;
	std::vector<VecX> nullspaces_x0_pre;
	for(int i=0;i<6;i++)
	{
		VecX nullspace_x0(n);
		nullspace_x0.setZero();
		for(FrameHessian* fh : frameHessians)
		{
			nullspace_x0.segment<6>(CPARS+fh->idx*8) = fh->nullspaces_pose.col(i);
			nullspace_x0.segment<3>(CPARS+fh->idx*8) *= SCALE_XI_TRANS_INVERSE;
			nullspace_x0.segment<3>(CPARS+fh->idx*8+3) *= SCALE_XI_ROT_INVERSE;
		}
		nullspaces_x0_pre.push_back(nullspace_x0);
		nullspaces_pose.push_back(nullspace_x0);
	}
	for(int i=0;i<2;i++)
	{
		VecX nullspace_x0(n);
		nullspace_x0.setZero();
		for(FrameHessian* fh : frameHessians)
		{
			nullspace_x0.segment<2>(CPARS+fh->idx*8+6) = fh->nullspaces_affine.col(i).head<2>();
			nullspace_x0[CPARS+fh->idx*8+6] *= SCALE_A_INVERSE;
			nullspace_x0[CPARS+fh->idx*8+7] *= SCALE_B_INVERSE;
		}
		nullspaces_x0_pre.push_back(nullspace_x0);
		if(i==0) nullspaces_affA.push_back(nullspace_x0);
		if(i==1) nullspaces_affB.push_back(nullspace_x0);
	}

	VecX nullspace_x0(n);
	nullspace_x0.setZero();
	for(FrameHessian* fh : frameHessians)
	{
		nullspace_x0.segment<6>(CPARS+fh->idx*8) = fh->nullspaces_scale;
		nullspace_x0.segment<3>(CPARS+fh->idx*8) *= SCALE_XI_TRANS_INVERSE;
		nullspace_x0.segment<3>(CPARS+fh->idx*8+3) *= SCALE_XI_ROT_INVERSE;
	}
	nullspaces_x0_pre.push_back(nullspace_x0);
	nullspaces_scale.push_back(nullspace_x0);

	return nullspaces_x0_pre;
}


dmvio::IMUIntegration &FullSystem::getImuIntegration()
{
    return imuIntegration;
}

}
