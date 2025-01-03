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



#include "FullSystem/CoarseInitializer.h"
#include "FullSystem/FullSystem.h"
#include "FullSystem/HessianBlocks.h"
#include "FullSystem/Residuals.h"
#include "FullSystem/PixelSelector.h"
#include "FullSystem/PixelSelector2.h"
#include "util/nanoflann.h"

#include <opencv2/highgui.hpp>

#if !defined(__SSE3__) && !defined(__SSE2__) && !defined(__SSE1__)
#include "SSE2NEON.h"
#endif



namespace dso
{

/**
 * @brief Construct a new Coarse Initializer:: Coarse Initializer object
 * 
 * @param globalCalib_ 
 * @param globalSettings_ 
 */
CoarseInitializer::CoarseInitializer(Global_Calib& globalCalib_, GlobalSettings& globalSettings_)
		: thisToNext_aff(0, 0), thisToNext(SE3()), globalCalib(globalCalib_), globalSettings(globalSettings_)
{
	for(int lvl=0; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		points[lvl] = 0;
		numPoints[lvl] = 0;
	}

	wG0 = globalCalib.wG[0];
	hG0 = globalCalib.hG[0];

	pixelSelectorSparsity = globalSettings.setting_sparsityFactor;

	JbBuffer = new Vec10f[wG0*hG0];
	JbBuffer_new = new Vec10f[wG0*hG0];


	frameID=-1;
	fixAffine=true;
	printDebug=!setting_debugout_runquiet;

	wM.diagonal()[0] = wM.diagonal()[1] = wM.diagonal()[2] = SCALE_XI_ROT;
	wM.diagonal()[3] = wM.diagonal()[4] = wM.diagonal()[5] = SCALE_XI_TRANS;
	wM.diagonal()[6] = SCALE_A;
	wM.diagonal()[7] = SCALE_B;
}

/**
 * @brief Destroy the Coarse Initializer
 * 
 */
CoarseInitializer::~CoarseInitializer()
{
	for(int lvl=0; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		if(points[lvl] != 0) delete[] points[lvl];
	}

	delete[] JbBuffer;
	delete[] JbBuffer_new;
}


/**
 * @brief Track frame without initialization
 * 
 * The main purpose of this function is to calculate initial point depths
 * setFirst() should be called first to get the initial pixel list
 * 
 * This function works by using optimizations and heuristics till a good set of depths are found
 * 
 * This function should find a good frame and then finished the init after five frames
 * If this function cannot find a good frame, it will trigger a reset
 * 
 * @param newFrameHessian 
 * @param wraps 
 * @return true 
 * @return false 
 */
bool CoarseInitializer::trackFrame(FrameHessian *newFrameHessian, std::vector<IOWrap::Output3DWrapper*> &wraps)
{
	// Attach frame
	newFrame = newFrameHessian;

	for(IOWrap::Output3DWrapper* ow : wraps)
		ow->pushLiveFrame(newFrameHessian);

	int maxIterations[] = {5,5,10,30,50};


	alphaK = 2.5*2.5;
	alphaW = 150*150;
	regWeight = 0.8;
	couplingWeight = 1;

	if(!snapped)
	{
		// Set all the points to default values
		thisToNext.translation().setZero();
		for(int lvl=0;lvl<globalSettings.pyrLevelsUsed;lvl++)
		{
			int npts = numPoints[lvl];
			Pnt* ptsl = points[lvl];
			for(int i=0;i<npts;i++)
			{
				ptsl[i].iR = 1;
				ptsl[i].idepth_new = 1;
				ptsl[i].lastHessian = 0;
			}
		}
	}


	SE3 refToNew_current = thisToNext;
	AffLight refToNew_aff_current = thisToNext_aff;
	if(firstFrame->ab_exposure>0 && newFrame->ab_exposure>0)
		refToNew_aff_current = AffLight(logf(newFrame->ab_exposure /  firstFrame->ab_exposure),0); // coarse approximation


	// Start calculating values
	Vec3f latestRes = Vec3f::Zero();
	for(int lvl=globalSettings.pyrLevelsUsed-1; lvl>=0; lvl--)
	{
		// Applys the valid values calculated from the level below (lower to higher) to the current level
		if(lvl<globalSettings.pyrLevelsUsed-1)
			propagateDown(lvl+1);

		Mat88f H,Hsc; Vec8f b,bsc;
		// Prepare points for energy calculations
		resetPoints(lvl);

		// Calculate initial energy
		// refToNew_current for the first frame will have no movement
		// refToNew_aff_current for the first frame will be the default or a coarse approximation
		Vec3f resOld = calcResAndGS(lvl, H, b, Hsc, bsc, refToNew_current, refToNew_aff_current, false);
		// Set point values to new values
		applyStep(lvl);

		float lambda = 0.1;
		float eps = 1e-4;
		int fails=0;

		if(printDebug && !globalSettings.no_CoarseInit_debugMessage)
		{
			printf("lvl %d, it %d (l=%f) %s: %.3f+%.5f -> %.3f+%.5f (%.3f->%.3f) (|inc| = %f)! \t",
					lvl, 0, lambda,
					"INITIA",
					sqrtf((float)(resOld[0] / resOld[2])),
					sqrtf((float)(resOld[1] / resOld[2])),
					sqrtf((float)(resOld[0] / resOld[2])),
					sqrtf((float)(resOld[1] / resOld[2])),
					(resOld[0]+resOld[1]) / resOld[2],
					(resOld[0]+resOld[1]) / resOld[2],
					0.0f);
			std::cout << refToNew_current.log().transpose() << " AFF " << refToNew_aff_current.vec().transpose() <<"\n";
		}

		// Iterate till good values are found or too many iterations
		int iteration=0;
		while(true)
		{
			Mat88f Hl = H;
			for(int i=0;i<8;i++) Hl(i,i) *= (1+lambda);
			Hl -= Hsc*(1/(1+lambda));
			Vec8f bl = b - bsc*(1/(1+lambda));

			Hl = wM * Hl * wM * (0.01f/(w[lvl]*h[lvl]));
			bl = wM * bl * (0.01f/(w[lvl]*h[lvl]));


			Vec8f inc;
			SE3 refToNew_new;

			// Calculate the increment from the given H and b matrices
			if (fixAffine)
			{
				// Note as we set the weights of rotation and translation to 1 the wM is just the identity in this case.
				inc.head<6>() = -(wM.toDenseMatrix().topLeftCorner<6, 6>() *
								  (Hl.topLeftCorner<6, 6>().ldlt().solve(bl.head<6>())));
				inc.tail<2>().setZero();
			} else
				inc = -(wM * (Hl.ldlt().solve(bl)));    // = -H^-1 * b.
			double incNorm = inc.norm();

			// Calculate new values
			// First size values of inc are the translation and rotation (SE3 matrix)
			refToNew_new = SE3::exp(inc.head<6>().cast<double>()) * refToNew_current;

			// Calculate new photometric values
			// Last two values of the inc are the photometric values
			AffLight refToNew_aff_new = refToNew_aff_current;
			refToNew_aff_new.a += inc[6];
			refToNew_aff_new.b += inc[7];

			// Update point depths
			doStep(lvl, lambda, inc);

			// Repeat optimization with new values
			Mat88f H_new, Hsc_new; Vec8f b_new, bsc_new;
			Vec3f resNew = calcResAndGS(lvl, H_new, b_new, Hsc_new, bsc_new, refToNew_new, refToNew_aff_new, false);
			Vec3f regEnergy = calcEC(lvl);

			// Decide to accept or reject
			float eTotalNew = (resNew[0]+resNew[1]+regEnergy[1]);
			float eTotalOld = (resOld[0]+resOld[1]+regEnergy[0]);

			bool accept = eTotalOld > eTotalNew;

			if(printDebug && !globalSettings.no_CoarseInit_debugMessage)
			{
				printf("lvl %d, it %d (l=%f) %s: %.5f + %.5f + %.5f -> %.5f + %.5f + %.5f (%.2f->%.2f) (|inc| = %f)! \t",
						lvl, iteration, lambda,
						(accept ? "ACCEPT" : "REJECT"),
						sqrtf((float)(resOld[0] / resOld[2])),
						sqrtf((float)(regEnergy[0] / regEnergy[2])),
						sqrtf((float)(resOld[1] / resOld[2])),
						sqrtf((float)(resNew[0] / resNew[2])),
						sqrtf((float)(regEnergy[1] / regEnergy[2])),
						sqrtf((float)(resNew[1] / resNew[2])),
						eTotalOld / resNew[2],
						eTotalNew / resNew[2],
						incNorm);
				std::cout << refToNew_new.log().transpose() << " AFF " << refToNew_aff_new.vec().transpose() <<"\n";
			}

			if(accept)
			{
				if(resNew[1] == alphaK*numPoints[lvl]) // Stop init if values good enough
					snapped = true;

				// Update H, b, and residual
				H = H_new;
				b = b_new;
				Hsc = Hsc_new;
				bsc = bsc_new;
				resOld = resNew;
				// Update transform
				refToNew_aff_current = refToNew_aff_new;
				refToNew_current = refToNew_new;
				// Update points
				applyStep(lvl);
				optReg(lvl);

				lambda *= 0.5;
				fails=0;
				if(lambda < 0.0001) lambda = 0.0001;
			}
			else // reduce step size and try again
			{
				fails++;
				lambda *= 4;

				if(lambda > 10000) lambda = 10000;
			}


			bool quitOpt = false;
			if(!(incNorm > eps) || iteration >= maxIterations[lvl] || fails >= 2)
			{
				quitOpt = true;
			}

			if(quitOpt) break;
			iteration++;
		}
		latestRes = resOld;
	}

	// Update transformation with optimized values
	thisToNext = refToNew_current;
	thisToNext_aff = refToNew_aff_current;

	for(int i=0;i<globalSettings.pyrLevelsUsed-1;i++)
		propagateUp(i);

	frameID++;
	if(!snapped) snappedAt=0;

	if(snapped && snappedAt==0)
		snappedAt = frameID;

	debugPlot(0,wraps);

	return snapped && frameID > snappedAt+5;
}

/**
 * @brief Creates depth image
 * 
 * @param lvl 
 * @param wraps 
 */
void CoarseInitializer::debugPlot(int lvl, std::vector<IOWrap::Output3DWrapper*> &wraps)
{
	bool needCall = false;
	for(IOWrap::Output3DWrapper* ow : wraps)
		needCall = needCall || ow->needPushDepthImage();
	if(!needCall) return;


	int wl = w[lvl], hl = h[lvl];
	Eigen::Vector3f* colorRef = firstFrame->dIp[lvl];

	MinimalImageB3 iRImg(wl,hl);

	for(int i=0;i<wl*hl;i++)
		iRImg.at(i) = Vec3b(colorRef[i][0],colorRef[i][0],colorRef[i][0]);

	int npts = numPoints[lvl];

	float nid = 0, sid=0;
	for(int i=0;i<npts;i++)
	{
		Pnt* point = points[lvl]+i;
		if(point->isGood)
		{
			nid++;
			sid += point->iR;
		}
	}
	float fac = nid / sid;


	for(int i=0;i<npts;i++)
	{
		Pnt* point = points[lvl]+i;

		if(!point->isGood)
			iRImg.setPixel9(point->u+0.5f,point->v+0.5f,Vec3b(0,0,0));
		else
			iRImg.setPixel9(point->u+0.5f,point->v+0.5f,makeRainbow3B(point->iR*fac));
	}

	for(IOWrap::Output3DWrapper* ow : wraps)
		ow->pushDepthImage(&iRImg);
}

/**
 * @brief Calculates residual, Hessian and Hessian-block needed for re-substituting depth.
 * 
 * @param lvl 
 * @param H_out 
 * @param b_out 
 * @param H_out_sc 
 * @param b_out_sc 
 * @param refToNew 
 * @param refToNew_aff 
 * @param plot 
 * @return Vec3f 		(Energy Residual, alphaEnergy, number of points)
 */
Vec3f CoarseInitializer::calcResAndGS(
		int lvl, Mat88f &H_out, Vec8f &b_out,
		Mat88f &H_out_sc, Vec8f &b_out_sc,
		const SE3 &refToNew, AffLight refToNew_aff,
		bool plot)
{
	int wl = w[lvl], hl = h[lvl];
	// Get pixel intensities
	Eigen::Vector3f* colorRef = firstFrame->dIp[lvl];
	Eigen::Vector3f* colorNew = newFrame->dIp[lvl];

	// Set transformation matrix
	Mat33f RKi = (refToNew.rotationMatrix() * Ki[lvl]).cast<float>();
	Vec3f t = refToNew.translation().cast<float>();
	Eigen::Vector2f r2new_aff = Eigen::Vector2f(exp(refToNew_aff.a), refToNew_aff.b);

	// Set calibration matrix
	float fxl = fx[lvl];
	float fyl = fy[lvl];
	float cxl = cx[lvl];
	float cyl = cy[lvl];

	// Intialize accumulators
	for(auto&& acc9 : acc9s)
	{
		acc9.initialize();
	}
	for(auto&& E : accE)
	{
		E.initialize();
	}

	int npts = numPoints[lvl];
	Pnt* ptsl = points[lvl];

	// lambda function so it can be parallelized
	// This part takes most of the time for this method --> parallelize this only
	auto processPointsForReduce = [&](int min=0, int max=1, double* stats=0, int tid=0)
	{
		auto& acc9 = acc9s[tid];
		auto& E = accE[tid];

		// For every active point
		for(int i = min; i < max; i++)
		{
			Pnt* point = ptsl + i;

			point->maxstep = 1e10;

			// use previous valid value of energy of point if bad point
			if(!point->isGood)
			{
				E.updateSingle((float) (point->energy[0]));
				point->energy_new = point->energy;
				point->isGood_new = false;
				continue;
			}

			VecNRf dp0;
			VecNRf dp1;
			VecNRf dp2;
			VecNRf dp3;
			VecNRf dp4;
			VecNRf dp5;
			VecNRf dp6;
			VecNRf dp7;
			VecNRf dd;
			VecNRf r;
			JbBuffer_new[i].setZero();

			// sum over all residuals.
			bool isGood = true;
			float energy = 0;
			// Account for the pixel pattern used
			for(int idx = 0; idx < PATTERNNUM; idx++)
			{
				int dx = PATTERNP[idx][0];
				int dy = PATTERNP[idx][1];

				// Convert to Cartesian coordinates and apply transform
				Vec3f pt = RKi * Vec3f(point->u + dx, point->v + dy, 1) + t * point->idepth_new;
				// Calculate new position and depth
				// Transform back to projective Image/Sensor coordinates
				float u_im_j = pt[0] / pt[2];
				float v_im_j = pt[1] / pt[2];
				// Transform back to projective Pixel coordinates
				float Ku_pc_j = fxl * u_im_j + cxl;
				float Kv_pc_j = fyl * v_im_j + cyl;
				float new_idepth = point->idepth_new / pt[2];

				if(!(Ku_pc_j > 1 && Kv_pc_j > 1 && Ku_pc_j < wl - 2 && Kv_pc_j < hl - 2 && new_idepth > 0)) // out of bounds
				{
					isGood = false;
					break;
				}

				// Get intensitives of the pixel values in the new and first frames
				Vec3f hitColor = getInterpolatedElement33(colorNew, Ku_pc_j, Kv_pc_j, wl, hl);
				float rlR = getInterpolatedElement31(colorRef, point->u + dx, point->v + dy, wl, hl);

				if(!std::isfinite(rlR) || !std::isfinite((float) hitColor[0])) // infinite residual
				{
					isGood = false;
					break;
				}


				// Calculate residual
				// Note that exposure is not included because first frame
				// has no exposure ratio to compare to
				float residual = hitColor[0] - r2new_aff[0] * rlR - r2new_aff[1];

				// Huber loss
				float hw = fabs(residual) < globalSettings.setting_huberTH ? 1 : globalSettings.setting_huberTH / fabs(residual);
				energy += hw * residual * residual * (2 - hw);


				// Calculate Jacobian
				// Pixel depth partial derivative
				// du/dd
				float dxdd = (t[0] - t[2] * u_im_j) / pt[2];
				// dv/dd
				float dydd = (t[1] - t[2] * v_im_j) / pt[2];

				if(hw < 1) hw = sqrtf(hw); // huber weight
				// Pixel Intensity derivatives
				float dxInterp = hw * hitColor[1] * fxl;
				float dyInterp = hw * hitColor[2] * fyl;
				
				// d_I[u,v]/d_tx
				dp0[idx] = new_idepth * dxInterp;						// inverse_depth * dx
				// d_I[u,v]/d_ty
				dp1[idx] = new_idepth * dyInterp;						// inverse_depth * dy
				// d_I[u,v]/d_tz
				dp2[idx] = -new_idepth * (u_im_j * dxInterp + v_im_j * dyInterp);	// inverse_depth* (u*dx + v*dy)
				dp3[idx] = -u_im_j * v_im_j * dxInterp - (1 + v_im_j * v_im_j) * dyInterp;	// -dx * (u*v) - dy * (1 + v^2)
				dp4[idx] = (1 + u_im_j * u_im_j) * dxInterp + u_im_j * v_im_j * dyInterp;	// dy * (u*v) + dx * (1 + u^2)
				dp5[idx] = -v_im_j * dxInterp + u_im_j * dyInterp;				// dy * u - dx * v
				//dE/da
				dp6[idx] = -hw * r2new_aff[0] * rlR;					// a * (b0 - I)
				//dE/db
				dp7[idx] = -hw * 1;							// -1
				// dE/dd
				dd[idx] = dxInterp * dxdd + dyInterp * dydd;
				r[idx] = hw * residual;									// residual energy		

				float maxstep = 1.0f / Vec2f(dxdd * fxl, dydd * fyl).norm();
				if(maxstep < point->maxstep) point->maxstep = maxstep;

				// Calculate point depths
				// immediately compute dp*dd' and dd*dd' in JbBuffer1.
				JbBuffer_new[i][0] += dp0[idx] * dd[idx];
				JbBuffer_new[i][1] += dp1[idx] * dd[idx];
				JbBuffer_new[i][2] += dp2[idx] * dd[idx];
				JbBuffer_new[i][3] += dp3[idx] * dd[idx];
				JbBuffer_new[i][4] += dp4[idx] * dd[idx];
				JbBuffer_new[i][5] += dp5[idx] * dd[idx];
				JbBuffer_new[i][6] += dp6[idx] * dd[idx];
				JbBuffer_new[i][7] += dp7[idx] * dd[idx];
				JbBuffer_new[i][8] += r[idx] * dd[idx];
				JbBuffer_new[i][9] += dd[idx] * dd[idx];
			}

			 // use previous valid value of energy of point if outlier
			if(!isGood || energy > point->outlierTH * 20)
			{
				E.updateSingle((float) (point->energy[0]));
				point->isGood_new = false;
				point->energy_new = point->energy;
				continue;
			}


			// Add into energy.
			E.updateSingle(energy);
			point->isGood_new = true;
			point->energy_new[0] = energy;

			// Update Hessian matrix.
			for(int j = 0; j + 3 < PATTERNNUM; j += 4){
				acc9.updateSSE(
						_mm_load_ps(((float*) (&dp0)) + j),
						_mm_load_ps(((float*) (&dp1)) + j),
						_mm_load_ps(((float*) (&dp2)) + j),
						_mm_load_ps(((float*) (&dp3)) + j),
						_mm_load_ps(((float*) (&dp4)) + j),
						_mm_load_ps(((float*) (&dp5)) + j),
						_mm_load_ps(((float*) (&dp6)) + j),
						_mm_load_ps(((float*) (&dp7)) + j),
						_mm_load_ps(((float*) (&r)) + j));
			}

			for(int j = ((PATTERNNUM >> 2) << 2); j < PATTERNNUM; j++){
				acc9.updateSingle(
						(float) dp0[j], (float) dp1[j], (float) dp2[j], (float) dp3[j],
						(float) dp4[j], (float) dp5[j], (float) dp6[j], (float) dp7[j],
						(float) r[j]);
			}

		}
	};

	// Calculate Hessian (acc9) and Energy (E)
	reduce.reduce(processPointsForReduce, 0, npts, 50);

	for(auto&& acc9 : acc9s)
	{
		acc9.finish();
	}
	for(auto&& E : accE)
	{
		E.finish();
	}


	// Calculate alpha energy (EAlpha), and decide if we cap it.
	Accumulator11 EAlpha;
	EAlpha.initialize();
	for(int i=0;i<npts;i++)
	{
		Pnt* point = ptsl+i;
		if(!point->isGood_new)
		{
			// This should actually be EAlpha, but it seems like fixing this might change the optimal values of some
			// parameters, so it's kept like it is (see https://github.com/JakobEngel/dso/issues/52)
			// At the moment, this code will not change the value of E.A (because E.finish() is not called again after
			// this. It will however change E.num.
			accE[0].updateSingle((float)(point->energy[1]));
		}
		else
		{
			point->energy_new[1] = (point->idepth_new-1)*(point->idepth_new-1);
			accE[0].updateSingle((float)(point->energy_new[1]));
		}
	}
	EAlpha.finish();
	float alphaEnergy = alphaW*(EAlpha.A + refToNew.translation().squaredNorm() * npts);


	// compute alpha opt
	float alphaOpt;
	if(alphaEnergy > alphaK*npts)
	{
		alphaOpt = 0;
		alphaEnergy = alphaK*npts;
	}
	else
	{
		alphaOpt = alphaW;
	}

	
	// Calculate the HessianScaled
	acc9SC.initialize();
	for(int i=0;i<npts;i++)
	{
		Pnt* point = ptsl+i;
		if(!point->isGood_new)
			continue;

		// Set the value of (dE/dd)^2
		point->lastHessian_new = JbBuffer_new[i][9];
		JbBuffer_new[i][9] += alphaOpt;

		// Set the value of R*(dE/dd)
		JbBuffer_new[i][8] += alphaOpt*(point->idepth_new - 1);

		if(alphaOpt==0)
		{
			JbBuffer_new[i][8] += couplingWeight*(point->idepth_new - point->iR);
			JbBuffer_new[i][9] += couplingWeight;
		}

		// Set (dE/dd)^2 to 1/(1+(dE/dd)^2) for computational effciency
		JbBuffer_new[i][9] = 1/(1+JbBuffer_new[i][9]);
		acc9SC.updateSingleWeighted(
				(float)JbBuffer_new[i][0],(float)JbBuffer_new[i][1],(float)JbBuffer_new[i][2],(float)JbBuffer_new[i][3],
				(float)JbBuffer_new[i][4],(float)JbBuffer_new[i][5],(float)JbBuffer_new[i][6],(float)JbBuffer_new[i][7],
				(float)JbBuffer_new[i][8],(float)JbBuffer_new[i][9]);
	}
	acc9SC.finish();

	// Extract the H and b matrices
	H_out.setZero();
	b_out.setZero();
	// This needs to sum up the acc9s from all the workers!
	for(auto&& acc9 : acc9s)
	{
		H_out += acc9.H.topLeftCorner<8,8>();// / acc9.num;
		b_out += acc9.H.topRightCorner<8,1>();// / acc9.num;
	}
	H_out_sc = acc9SC.H.topLeftCorner<8,8>();// / acc9.num;
	b_out_sc = acc9SC.H.topRightCorner<8,1>();// / acc9.num;


	H_out(0,0) += alphaOpt*npts;
	H_out(1,1) += alphaOpt*npts;
	H_out(2,2) += alphaOpt*npts;

	Vec3f tlog = refToNew.log().head<3>().cast<float>();
	b_out[0] += tlog[0]*alphaOpt*npts;
	b_out[1] += tlog[1]*alphaOpt*npts;
	b_out[2] += tlog[2]*alphaOpt*npts;


	// Add zero prior to translation.
	// setting_weightZeroPriorDSOInitY is the squared weight of the prior residual.
	H_out(1, 1) += globalSettings.setting_weightZeroPriorDSOInitY;
	b_out(1) += globalSettings.setting_weightZeroPriorDSOInitY * refToNew.translation().y();

	H_out(0, 0) += globalSettings.setting_weightZeroPriorDSOInitX;
	b_out(0) += globalSettings.setting_weightZeroPriorDSOInitX * refToNew.translation().x();


	double A = 0;
	int num = 0;
	for(auto&& E : accE)
	{
		A += E.A;
		num += E.num;
	}

	return Vec3f(A, alphaEnergy, num);
}

float CoarseInitializer::rescale()
{
	float factor = 20*thisToNext.translation().norm();

	return factor;
}


Vec3f CoarseInitializer::calcEC(int lvl)
{
	if(!snapped) return Vec3f(0,0,numPoints[lvl]);
	AccumulatorX<2> E;
	E.initialize();
	int npts = numPoints[lvl];
	for(int i=0;i<npts;i++)
	{
		Pnt* point = points[lvl]+i;
		if(!point->isGood_new) continue;
		float rOld = (point->idepth-point->iR);
		float rNew = (point->idepth_new-point->iR);
		E.updateNoWeight(Vec2f(rOld*rOld,rNew*rNew));
	}
	E.finish();

	return Vec3f(couplingWeight*E.A1m[0], couplingWeight*E.A1m[1], E.num);
}

/**
 * @brief Smooths out the iR values between point neighbours
 * 
 * Smooths with closest N points at same level
 * 
 * @param lvl 
 */
void CoarseInitializer::optReg(int lvl)
{
	int npts = numPoints[lvl];
	Pnt* ptsl = points[lvl];

	if(!snapped)
	{
		return;
	}

	for(size_t i=0;i<npts;i++)
	{
		Pnt* point = ptsl+i;
		if(!point->isGood) continue;

		// Get all of the iR values for the good neighbours
		float idnn[10];
		unsigned int nnn = 0;
		for(unsigned int j=0;j<10;j++) // for all point neighbours
		{
			if(point->neighbours[j] == -1) continue;

			Pnt* other = ptsl+point->neighbours[j];

			if(!other->isGood) continue;

			idnn[nnn] = other->iR;
			nnn++;
		}

		if(nnn > 2)
		{
			// Sort to get median value of iR for the good neighbours
			// idnn[nnn/2] will be the median
			std::nth_element(idnn,idnn+nnn/2,idnn+nnn);

			// Set new iR to account for current inverse depth and median iR value of neighbours
			point->iR = (1-regWeight)*point->idepth + regWeight*idnn[nnn/2];
		}
	}
}


/**
 * @brief Apply information from children point to parent points for given level
 * 
 * @param srcLvl 
 */
void CoarseInitializer::propagateUp(int srcLvl)
{
	assert(srcLvl+1<globalSettings.pyrLevelsUsed);
	// set idepth of target

	int nptss= numPoints[srcLvl];
	int nptst= numPoints[srcLvl+1];
	Pnt* ptss = points[srcLvl];
	Pnt* ptst = points[srcLvl+1];

	// set to zero.
	for(int i=0;i<nptst;i++) // for all parent points
	{
		Pnt* parent = ptst+i;
		parent->iR=0;
		parent->iRSumNum=0;
	}

	// Set parent point parameters with child values
	for(int i=0;i<nptss;i++) // for all points
	{
		Pnt* point = ptss+i;
		if(!point->isGood) continue;

		Pnt* parent = ptst + point->parent;
		parent->iR += point->iR * point->lastHessian;
		parent->iRSumNum += point->lastHessian;
	}

	for(int i=0;i<nptst;i++) // for all parent points
	{
		Pnt* parent = ptst+i;
		if(parent->iRSumNum > 0)
		{
			parent->idepth = parent->iR = (parent->iR / parent->iRSumNum);
			parent->isGood = true;
		}
	}

	optReg(srcLvl+1);
}

/**
 * @brief Apply information from parent (lower) point to children (higher) points for given level
 * 
 * This function sets a better starting point for the optimization than the default value (1)
 * by using the valid values calculated in a pyramid level above
 * 
 * @param srcLvl 
 */
void CoarseInitializer::propagateDown(int srcLvl)
{
	assert(srcLvl>0);
	// set idepth of target

	unsigned int nptst = numPoints[srcLvl-1];
	Pnt* ptst = points[srcLvl-1];
	Pnt* ptss = points[srcLvl];

	for(size_t i=0;i<nptst;i++)
	{
		Pnt* point = ptst+i;
		Pnt* parent = ptss+(point->parent);

		if(!parent->isGood || parent->lastHessian < 0.1) continue;

		if(!point->isGood) // Use the value of the parent if depth is not good
		{
			// Set values between point and parent to be equal
			point->iR = point->idepth = point->idepth_new = parent->iR;

			point->isGood=true;
			point->lastHessian=0;
		}
		else // Take weighted average is both depth values are good
		{
			// Calculate new iR weighted by the lastHessians
			float newiR = (point->iR * point->lastHessian * 2 + parent->iR * parent->lastHessian) 
						/ (point->lastHessian * 2 + parent->lastHessian);
			// Set new value
			point->iR = point->idepth = point->idepth_new = newiR;
		}
	}

	// Smooth out depth values on child layer
	optReg(srcLvl-1);
}


void CoarseInitializer::makeGradients(Eigen::Vector3f** data)
{
	for(int lvl=1; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		int lvlm1 = lvl-1;
		int wl = w[lvl], hl = h[lvl], wlm1 = w[lvlm1];

		Eigen::Vector3f* dINew_l = data[lvl];
		Eigen::Vector3f* dINew_lm = data[lvlm1];

		for(int y=0;y<hl;y++)
			for(int x=0;x<wl;x++)
				dINew_l[x + y*wl][0] = 0.25f * (dINew_lm[2*x   + 2*y*wlm1][0] +
													dINew_lm[2*x+1 + 2*y*wlm1][0] +
													dINew_lm[2*x   + 2*y*wlm1+wlm1][0] +
													dINew_lm[2*x+1 + 2*y*wlm1+wlm1][0]);

		for(int idx=wl;idx < wl*(hl-1);idx++)
		{
			dINew_l[idx][1] = 0.5f*(dINew_l[idx+1][0] - dINew_l[idx-1][0]);
			dINew_l[idx][2] = 0.5f*(dINew_l[idx+wl][0] - dINew_l[idx-wl][0]);
		}
	}
}


/**
 * @brief Sets the first frame
 * 
 * Starting point of the CoarseInitializer
 * Chooses the pixel values used of initial tracking
 * 
 * @param HCalib 
 * @param newFrameHessian 
 */
void CoarseInitializer::setFirst(CalibHessian* HCalib, FrameHessian* newFrameHessian)
{
	// Set calibration matrix
	makeK(HCalib);
	// Attrach frame
	firstFrame = newFrameHessian;

	// setting_minGradHistCut and setting_minGradHistAdd can change
	PixelSelector sel(globalCalib, globalSettings.setting_minGradHistCut, globalSettings.setting_minGradHistAdd, globalSettings);

	// Stores positiions where points were selected to be
	float* statusMap = new float[w[0]*h[0]];
	bool* statusMapB = new bool[w[0]*h[0]];
	float densities[] = {0.03,0.05,0.15,0.5,1};

	// Select pixels on all levels using the pixel selector
	for(unsigned int lvl=0; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		// Select positions for points in first frame
		sel.currentPotential = 3;
		int npts;
		if(lvl == 0)
			npts = sel.makeMaps(firstFrame, statusMap, densities[lvl]*w[0]*h[0], 1, 2);
		else
			npts = makePixelStatus(firstFrame->dIp[lvl], statusMapB, w[lvl], h[lvl], densities[lvl]*w[0]*h[0], pixelSelectorSparsity);


		if(points[lvl] != 0) delete[] points[lvl];
		points[lvl] = new Pnt[npts];

		// Create the points using the postions given by the pixel selector
		// set idepth map to initially 1 everywhere
		int wl = w[lvl], hl = h[lvl];
		Pnt* pl = points[lvl];
		int nl = 0;

		// For all pixels excluding the border
		for(size_t y=PATTERNPADDING+1;y<hl-PATTERNPADDING-2;y++)
		for(size_t x=PATTERNPADDING+1;x<wl-PATTERNPADDING-2;x++)
		{
			//if(x==2) printf("y=%d!\n",y);
			if((lvl!=0 && statusMapB[x+y*wl]) || (lvl==0 && statusMap[x+y*wl] != 0))
			{
				pl[nl].u = x+0.1;
				pl[nl].v = y+0.1;
				pl[nl].idepth = 1;
				pl[nl].iR = 1;
				pl[nl].isGood=true;
				pl[nl].energy.setZero();
				pl[nl].lastHessian=0;
				pl[nl].lastHessian_new=0;
				pl[nl].my_type= (lvl!=0) ? 1 : statusMap[x+y*wl];

/* 				
				// Experimental adjustable outlier threshold code
				Eigen::Vector3f* cpt = firstFrame->dIp[lvl] + x + y*w[lvl];
				float sumGrad2=0; // sum of gradients

				for(int idx=0;idx<PATTERNNUM;idx++)
				{
					int dx = PATTERNP[idx][0];
					int dy = PATTERNP[idx][1];
					float absgrad = cpt[dx + dy*w[lvl]].tail<2>().squaredNorm();
					sumGrad2 += absgrad;
				}

				float gth = setting_outlierTH * (sqrtf(sumGrad2)+setting_outlierTHSumComponent);
				pl[nl].outlierTH = PATTERNNUM*gth*gth;
*/

				pl[nl].outlierTH = PATTERNNUM*globalSettings.setting_outlierTH;

				nl++;
				assert(nl <= npts);
			}
		}
		numPoints[lvl]=nl;
	}
	delete[] statusMap;
	delete[] statusMapB;

	// For each point, find the ten nearest neighbours and their distance to the point
	// Also find the closest point one pyramid level above and it's distance
	makeNN();

	// Set flags and initial movement
	thisToNext = SE3();
	snapped = false;
	frameID = snappedAt = 0;

	for(unsigned int i=0;i<globalSettings.pyrLevelsUsed;i++)
		dGrads[i].setZero();
}

/**
 * @brief Set the values of bad points and resets the energy of all points to zero
 * 
 * @param lvl 
 */
void CoarseInitializer::resetPoints(int lvl)
{
	Pnt* pts = points[lvl];
	int npts = numPoints[lvl];

	for(size_t i=0;i<npts;i++)
	{
		pts[i].energy.setZero();
		pts[i].idepth_new = pts[i].idepth;

		if(lvl==globalSettings.pyrLevelsUsed-1 && !pts[i].isGood)  // Set the values of bad point using it's neighbours
		{
			float snd=0, sn=0;
			for(unsigned int n = 0; n<10; n++)
			{
				if(pts[i].neighbours[n] == -1 || !pts[pts[i].neighbours[n]].isGood) continue;

				snd += pts[pts[i].neighbours[n]].iR;
				sn += 1;
			}

			if(sn > 0) // The new iR of the point will be the average of the neighbouring values
			{
				pts[i].isGood=true;
				pts[i].iR = pts[i].idepth = pts[i].idepth_new = snd/sn;
			}
		}
	}
}

/**
 * @brief Update point depths
 * 
 * Apply optimized point depth increments
 * 
 * @param lvl 
 * @param lambda 
 * @param inc 
 */
void CoarseInitializer::doStep(int lvl, float lambda, Vec8f inc)
{
	const float maxPixelStep = 0.25;
	const float idMaxStep = 1e10;
	Pnt* pts = points[lvl];
	int npts = numPoints[lvl];

	for(size_t i=0;i<npts;i++)
	{
		if(!pts[i].isGood) continue;
		// Calculate increment of depth for point

		// b = (J_T_d * r_d + H_dξ * δd)
		float b = JbBuffer[i][8] + JbBuffer[i].head<8>().dot(inc);
		// δd = H_-1_dd * b
		// Note that H_-1_dd = 1/(1+JbBuffer_new[i][9])
		// JbBuffer_new[i][9] was inverted at a prior step
		float step = - b * JbBuffer[i][9] / (1+lambda);


		float maxstep = maxPixelStep*pts[i].maxstep;
		if(maxstep > idMaxStep) maxstep=idMaxStep;

		if(step >  maxstep) step = maxstep;
		if(step < -maxstep) step = -maxstep;

		float newIdepth = pts[i].idepth + step;
		if(newIdepth < 1e-3 ) newIdepth = 1e-3;
		if(newIdepth > 50) newIdepth = 50;

		pts[i].idepth_new = newIdepth;
	}
}

/**
 * @brief Set the point parameters to the new values
 * 
 * @param lvl 
 */
void CoarseInitializer::applyStep(int lvl)
{
	Pnt* pts = points[lvl];
	int npts = numPoints[lvl];

	for(size_t i=0;i<npts;i++)
	{
		if(!pts[i].isGood)
		{
			pts[i].idepth = pts[i].idepth_new = pts[i].iR;
			continue;
		}

		pts[i].energy = pts[i].energy_new;
		pts[i].isGood = pts[i].isGood_new;
		pts[i].idepth = pts[i].idepth_new;
		pts[i].lastHessian = pts[i].lastHessian_new;
	}
	std::swap<Vec10f*>(JbBuffer, JbBuffer_new);
}

/**
 * @brief Set width, height, and camera parameters for all pyramid levels
 * 
 * Helper function
 * 
 * @param HCalib 
 */
void CoarseInitializer::makeK(CalibHessian* HCalib)
{
	w[0] = wG0;
	h[0] = hG0;

	// Set calibration matrix
	fx[0] = HCalib->fxl();
	fy[0] = HCalib->fyl();
	cx[0] = HCalib->cxl();
	cy[0] = HCalib->cyl();

	for (unsigned int level = 1; level < globalSettings.pyrLevelsUsed; ++ level)
	{
		w[level] = w[0] >> level;
		h[level] = h[0] >> level;
		fx[level] = fx[level-1] * 0.5;
		fy[level] = fy[level-1] * 0.5;
		cx[level] = (cx[0] + 0.5) / ((int)1<<level) - 0.5;
		cy[level] = (cy[0] + 0.5) / ((int)1<<level) - 0.5;
	}

	// Set inverse calibration matrix
	for (unsigned int level = 0; level < globalSettings.pyrLevelsUsed; ++ level)
	{
		K[level] << fx[level], 0.0, cx[level], 
					0.0, fy[level], cy[level], 
					0.0, 0.0, 1.0;
		Ki[level] = K[level].inverse();
		fxi[level] = Ki[level](0,0);
		fyi[level] = Ki[level](1,1);
		cxi[level] = Ki[level](0,2);
		cyi[level] = Ki[level](1,2);
	}
}


/**
 * @brief For performing fast approximate nearest neighbor searches
 * 
 * Finds the values for:
 * - parent (idx (x+y*w) of closest point one pyramid level above)
 * - parent distance
 * - neighbours (idx (x+y*w) of up to N nearest points in pixel space)
 * - neighbours distances
 * 
 * Knowing the NN of each point helps with calculating point depths
 */
void CoarseInitializer::makeNN()
{
	const float NNDistFactor=0.05;

	typedef nanoflann::KDTreeSingleIndexAdaptor
			<nanoflann::L2_Simple_Adaptor<float, FLANNPointcloud>,FLANNPointcloud,2>
			KDTree;

	// Build indices
	FLANNPointcloud pcs[PYR_LEVELS];
	KDTree* indexes[PYR_LEVELS];

	for(unsigned int i=0;i<globalSettings.pyrLevelsUsed;i++)
	{
		// Build 2D KD Tree to help with NN searches
		pcs[i] = FLANNPointcloud(numPoints[i], points[i]);
		indexes[i] = new KDTree(2, pcs[i], nanoflann::KDTreeSingleIndexAdaptorParams(5) );
		indexes[i]->buildIndex();
	}

	// Number of nearest points to find
	const unsigned int nn = 10;

	// Find NN & parents
	for(unsigned int lvl=0;lvl<globalSettings.pyrLevelsUsed;lvl++)
	{
		Pnt* pts = points[lvl];
		int npts = numPoints[lvl];

		int ret_index[nn];
		float ret_dist[nn];
		nanoflann::KNNResultSet<float, int, int> resultSet(nn);
		nanoflann::KNNResultSet<float, int, int> resultSet1(1);

		for(size_t i=0;i<npts;i++)
		{
			resultSet.init(ret_index, ret_dist);
			Vec2f pt = Vec2f(pts[i].u,pts[i].v);

			// Find neighbours and neighbour distances
			indexes[lvl]->findNeighbors(resultSet, (float*)&pt, nanoflann::SearchParams());

			int myidx=0;
			float sumDF = 0;
			// For up to N points
			for(int k=0;k<nn;k++)
			{
				// Store the N nearest points and their distances
				pts[i].neighbours[myidx]=ret_index[k];

				// Distances are stored in exponential form
				float df = expf(-ret_dist[k]*NNDistFactor);
				// Total distance
				sumDF += df;

				pts[i].neighboursDist[myidx]=df;
				assert(ret_index[k]>=0 && ret_index[k] < npts);

				myidx++;
			}

			// Normalize distances
			for(int k=0;k<nn;k++)
				pts[i].neighboursDist[k] *= 10/sumDF;

			// Find parent point and parent distance
			if(lvl < globalSettings.pyrLevelsUsed-1 )
			{
				resultSet1.init(ret_index, ret_dist);
				pt = pt*0.5f-Vec2f(0.25f,0.25f);
				indexes[lvl+1]->findNeighbors(resultSet1, (float*)&pt, nanoflann::SearchParams());

				pts[i].parent = ret_index[0];
				// Distances are stored in exponential form
				pts[i].parentDist = expf(-ret_dist[0]*NNDistFactor);

				assert(ret_index[0]>=0 && ret_index[0] < numPoints[lvl+1]);
			}
			else
			{
				pts[i].parent = -1;
				pts[i].parentDist = -1;
			}
		}
	}

	// Cleanup
	for(unsigned int i=0;i<globalSettings.pyrLevelsUsed;i++)
		delete indexes[i];
}
}

