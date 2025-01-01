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



#include "FullSystem/CoarseTracker.h"
#include "FullSystem/FullSystem.h"
#include "FullSystem/HessianBlocks.h"
#include "FullSystem/Residuals.h"

#include "OptimizationBackend/EnergyFunctionalStructs.h"
#include "IOWrapper/ImageRW.h"
#include "util/TimeMeasurement.h"

#include <algorithm>

#if !defined(__SSE3__) && !defined(__SSE2__) && !defined(__SSE1__)
#include "SSE2NEON.h"
#endif



namespace dso
{

/**
 * @brief Allocates memory aligned to a specific number of bits
 * 
 * Allows for the usage of instructions that need aligned memory
 * All allocated pointers are stored in a vector
 * 
 * @tparam b 
 * @tparam T 
 * @param size 
 * @param rawPtrVec 
 * @return T* 
 */
template<int b, typename T>
T* allocAligned(int size, std::vector<T*> &rawPtrVec)
{
	const int padT = 1 + ((1 << b)/sizeof(T));
	T* ptr = new T[size + padT];
	rawPtrVec.push_back(ptr);
	T* alignedPtr = (T*)((((uintptr_t)(ptr+padT)) >> b) << b);
	return alignedPtr;
}

/**
 * @brief Construct a new Coarse Tracker:: Coarse Tracker object
 * 
 * @param globalCalib_ 
 * @param imuIntegration_ 
 * @param globalSettings_ 
 */
CoarseTracker::CoarseTracker(Global_Calib& globalCalib_, dmvio::IMUIntegration &imuIntegration_, GlobalSettings& globalSettings_) 
: lastRef_aff_g2l(0, 0), imuIntegration(imuIntegration_), globalCalib(globalCalib_), globalSettings(globalSettings_)
{
	// Set width and height
	wG0 = globalCalib.wG[0];
	hG0 = globalCalib.hG[0];

	// Make coarse tracking templates.
	// All arrays are aligned to allow for the use of high performance instructions
	// All pointers are stored in ptrToDelete for the destructor
	for(int lvl=0; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		int wl = wG0>>lvl;
		int hl = hG0>>lvl;

		idepth[lvl] = allocAligned<4,float>(wl*hl, ptrToDelete);
		weightSums[lvl] = allocAligned<4,float>(wl*hl, ptrToDelete);
		weightSums_bak[lvl] = allocAligned<4,float>(wl*hl, ptrToDelete);

		pc_u[lvl] = allocAligned<4,float>(wl*hl, ptrToDelete);
		pc_v[lvl] = allocAligned<4,float>(wl*hl, ptrToDelete);
		pc_idepth[lvl] = allocAligned<4,float>(wl*hl, ptrToDelete);
		pc_color[lvl] = allocAligned<4,float>(wl*hl, ptrToDelete);
	}

	unsigned int wh = wG0*hG0;

	// Warped buffers
	buf_warped_idepth = allocAligned<4,float>(wh, ptrToDelete);
	buf_warped_u = allocAligned<4,float>(wh, ptrToDelete);
	buf_warped_v = allocAligned<4,float>(wh, ptrToDelete);
	buf_warped_dx = allocAligned<4,float>(wh, ptrToDelete);
	buf_warped_dy = allocAligned<4,float>(wh, ptrToDelete);
	buf_warped_residual = allocAligned<4,float>(wh, ptrToDelete);
	buf_warped_weight = allocAligned<4,float>(wh, ptrToDelete);
	buf_warped_refColor = allocAligned<4,float>(wh, ptrToDelete);


	newFrame = 0;
	lastRef = 0;
	debugPrint = true;
	w[0]=h[0]=0;
	refFrameID=-1;
}

/**
 * @brief Destroy the Coarse Tracker:: Coarse Tracker object
 * 
 */
CoarseTracker::~CoarseTracker()
{
	for(float* ptr : ptrToDelete)
		delete[] ptr;
	ptrToDelete.clear();
}

/**
 * @brief Set width, height, and camera parameters for all pyramid levels
 * 
 * @param HCalib 
 */
void CoarseTracker::makeK(CalibHessian* HCalib)
{
	w[0] = wG0;
	h[0] = hG0;

	fx[0] = HCalib->fxl();
	fy[0] = HCalib->fyl();
	cx[0] = HCalib->cxl();
	cy[0] = HCalib->cyl();

	for (int level = 1; level < globalSettings.pyrLevelsUsed; ++ level)
	{
		w[level] = w[0] >> level;
		h[level] = h[0] >> level;

		fx[level] = fx[level-1] * 0.5;
		fy[level] = fy[level-1] * 0.5;
		cx[level] = (cx[0] + 0.5) / ((int)1<<level) - 0.5;
		cy[level] = (cy[0] + 0.5) / ((int)1<<level) - 0.5;
	}

	for (int level = 0; level < globalSettings.pyrLevelsUsed; ++ level)
	{
		K[level]  << fx[level], 0.0, cx[level], 
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
 * @brief Create a depth and weight map for the points in frames
 * 
 * @param frameHessians 
 */
void CoarseTracker::makeCoarseDepthL0(std::vector<FrameHessian*> frameHessians)
{
	// Make coarse tracking templates for latstRef.
	memset(idepth[0], 0, sizeof(float)*w[0]*h[0]);
	memset(weightSums[0], 0, sizeof(float)*w[0]*h[0]);

	// Set the depth and weight for the points involved in coarse tracking
	for(FrameHessian* fh : frameHessians) // for all active frames
	{
		for(PointHessian* ph : fh->pointHessians) // for all points in frame
		{
			if(ph->lastResiduals[0].first != 0 && ph->lastResiduals[0].second == ResState::IN)
			{
				PointFrameResidual* r = ph->lastResiduals[0].first;
				assert(r->efResidual->isActive() && r->target == lastRef);
				int u = r->centerProjectedTo[0] + 0.5f;
				int v = r->centerProjectedTo[1] + 0.5f;
				float new_idepth = r->centerProjectedTo[2];
				float weight = sqrtf(1e-3 / (ph->efPoint->HdiF+1e-12));

				idepth[0][u+w[0]*v] += new_idepth * weight;
				weightSums[0][u+w[0]*v] += weight;
			}
		}
	}

	// Set values on every pytamid level and filtering

	// Do a 2 by 2 box kernal convolution on the depth and weight values
	for(int lvl=1; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		int lvlm1 = lvl-1;
		int wl = w[lvl], hl = h[lvl], wlm1 = w[lvlm1];

		// Current level is set
		float* idepth_l = idepth[lvl];
		float* weightSums_l = weightSums[lvl];

		// Values are obtained from the last level
		float* idepth_lm = idepth[lvlm1];
		float* weightSums_lm = weightSums[lvlm1];

		for(int y=0;y<hl;y++)
			for(int x=0;x<wl;x++)
			{
				int bidx = 2*x + 2*y*wlm1;
				idepth_l[x + y*wl] = 		idepth_lm[bidx] +
											idepth_lm[bidx+1] +
											idepth_lm[bidx+wlm1] +
											idepth_lm[bidx+wlm1+1];

				weightSums_l[x + y*wl] = 	weightSums_lm[bidx] +
											weightSums_lm[bidx+1] +
											weightSums_lm[bidx+wlm1] +
											weightSums_lm[bidx+wlm1+1];
			}
	}


	// Interpolate values to areas without depth
	// Only done at higher resolutions
	for(int lvl=0; lvl<2; lvl++)
	{
		int numIts = 1;

		for(int it=0;it<numIts;it++)
		{
			int wh = w[lvl]*h[lvl]-w[lvl];
			int wl = w[lvl];

			float* weightSumsl = weightSums[lvl];
			float* weightSumsl_bak = weightSums_bak[lvl];
			memcpy(weightSumsl_bak, weightSumsl, w[lvl]*h[lvl]*sizeof(float));
			float* idepthl = idepth[lvl];	// don't need to make a temp copy of depth, since the code only
											// reads values with weightSumsl>0, and write ones with weightSumsl<=0.
			
			for(int i=w[lvl]+1;i<wh-1;i++) // frame area excluding border
			{
				if(weightSumsl_bak[i] <= 0)
				{
					float sum=0, num=0, numn=0;
					// Averaging matrix is a X shape
					if(weightSumsl_bak[i+1+wl] > 0) { sum += idepthl[i+1+wl]; num+=weightSumsl_bak[i+1+wl]; numn++;}
					if(weightSumsl_bak[i-1-wl] > 0) { sum += idepthl[i-1-wl]; num+=weightSumsl_bak[i-1-wl]; numn++;}
					if(weightSumsl_bak[i+wl-1] > 0) { sum += idepthl[i+wl-1]; num+=weightSumsl_bak[i+wl-1]; numn++;}
					if(weightSumsl_bak[i-wl+1] > 0) { sum += idepthl[i-wl+1]; num+=weightSumsl_bak[i-wl+1]; numn++;}
					if(numn>0) {idepthl[i] = sum/numn; weightSumsl[i] = num/numn;}
				}
			}
		}
	}


	// Interpolate values to areas without depth part 2
	// Only done at higher resolutions
	for(int lvl=2; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		int wh = w[lvl]*h[lvl]-w[lvl];
		int wl = w[lvl];

		float* weightSumsl = weightSums[lvl];
		float* weightSumsl_bak = weightSums_bak[lvl];
		memcpy(weightSumsl_bak, weightSumsl, w[lvl]*h[lvl]*sizeof(float));
		float* idepthl = idepth[lvl];	// don't need to make a temp copy of depth, since the code only
										// reads values with weightSumsl>0, and write ones with weightSumsl<=0.
		
		for(int i=w[lvl]+1;i<wh-1;i++) // frame area excluding border
		{
			if(weightSumsl_bak[i] <= 0)
			{
				float sum=0, num=0, numn=0;
				// Averaging matrix is a + shape
				if(weightSumsl_bak[i+1] > 0) { sum += idepthl[i+1]; num+=weightSumsl_bak[i+1]; numn++;}
				if(weightSumsl_bak[i-1] > 0) { sum += idepthl[i-1]; num+=weightSumsl_bak[i-1]; numn++;}
				if(weightSumsl_bak[i+wl] > 0) { sum += idepthl[i+wl]; num+=weightSumsl_bak[i+wl]; numn++;}
				if(weightSumsl_bak[i-wl] > 0) { sum += idepthl[i-wl]; num+=weightSumsl_bak[i-wl]; numn++;}
				if(numn>0) {idepthl[i] = sum/numn; weightSumsl[i] = num/numn;}
			}
		}
	}


	// Sets the semi-dense depth and pixel map for the reference frame
	// normalize idepths and weights
	for(int lvl=0; lvl<globalSettings.pyrLevelsUsed; lvl++)
	{
		float* weightSumsl = weightSums[lvl];
		float* idepthl = idepth[lvl];
		Eigen::Vector3f* dIRefl = lastRef->dIp[lvl];

		int wl = w[lvl], hl = h[lvl];

		int lpc_n=0;
		float* lpc_u = pc_u[lvl];
		float* lpc_v = pc_v[lvl];
		float* lpc_idepth = pc_idepth[lvl];
		float* lpc_color = pc_color[lvl];

		// For all of the frame area excluding the border
		for(int y=2;y<hl-2;y++)
			for(int x=2;x<wl-2;x++)
			{
				int i = x+y*wl;

				// Add the active points to the point list
				if(weightSumsl[i] > 0)
				{
					idepthl[i] /= weightSumsl[i];
					lpc_u[lpc_n] = x;
					lpc_v[lpc_n] = y;
					lpc_idepth[lpc_n] = idepthl[i];
					lpc_color[lpc_n] = dIRefl[i][0];

					if(!std::isfinite(lpc_color[lpc_n]) || !(idepthl[i]>0))
					{
						idepthl[i] = -1;
						continue;	// just skip if something is wrong.
					}
					lpc_n++;
				}
				else
					idepthl[i] = -1;

				weightSumsl[i] = 1;
			}
		pc_n[lvl] = lpc_n;
	}
}


/**
 * @brief Effciently accumulates the values for all of the valid points
 * 
 * Calculates H = J^T*W*J and b = -J^T*W*r
 * where J is the jacobian, W is the weights, and r is the stack residual vector
 * 
 * The varaibles being optimized are [w1, w2, w3, d1, d2, d3, a, b]
 * Point depths and calibration matrix variables are considered fixed
 * 
 * @param lvl 
 * @param H_out 
 * @param b_out 
 * @param refToNew 
 * @param aff_g2l 
 */
void CoarseTracker::calcGSSSE(int lvl, Mat88 &H_out, Vec8 &b_out, const SE3 &refToNew, AffLight aff_g2l)
{
	acc9.initialize(); // Initialize accumulater

	// Set values
	// Camera focus
	__m128 fxl = _mm_set1_ps(fx[lvl]);
	__m128 fyl = _mm_set1_ps(fy[lvl]);
	// Photometric values
	__m128 b0 = _mm_set1_ps(lastRef_aff_g2l.b);
	__m128 a = _mm_set1_ps((float)(AffLight::fromToVecExposure(lastRef->ab_exposure, newFrame->ab_exposure, lastRef_aff_g2l, aff_g2l)[0]));

	// Set 1, -1, and 0
	__m128 one = _mm_set1_ps(1);
	__m128 minusOne = _mm_set1_ps(-1);
	__m128 zero = _mm_set1_ps(0);

	// Make sure point buffer is the right size
	int n = buf_warped_n;
	assert(n%4==0);

	for(int i=0;i<n;i+=4) // for all valid points
	{
		// Pixel intensity derivative
		__m128 dx = _mm_mul_ps(_mm_load_ps(buf_warped_dx+i), fxl);
		__m128 dy = _mm_mul_ps(_mm_load_ps(buf_warped_dy+i), fyl);
		// New pixel position
		__m128 u = _mm_load_ps(buf_warped_u+i);
		__m128 v = _mm_load_ps(buf_warped_v+i);
		__m128 id = _mm_load_ps(buf_warped_idepth+i);

		// Accumulate matrix
		// Sum of all of the values multiplied with each other in every combination
		acc9.updateSSE_eighted(
				// inverse_depth * dx
				_mm_mul_ps(id,dx),
				// inverse_depth * dy
				_mm_mul_ps(id,dy),
				// inverse_depth* (u*dx + v*dy)
				_mm_sub_ps(zero, _mm_mul_ps(id,_mm_add_ps(_mm_mul_ps(u,dx), _mm_mul_ps(v,dy)))),
				// -dx * (u*v) - dy * (1 + v^2)
				_mm_sub_ps(zero, _mm_add_ps(
						_mm_mul_ps(_mm_mul_ps(u,v),dx),
						_mm_mul_ps(dy,_mm_add_ps(one, _mm_mul_ps(v,v))))),
				// dy * (u*v) + dx * (1 + u^2)
				_mm_add_ps(
						_mm_mul_ps(_mm_mul_ps(u,v),dy),
						_mm_mul_ps(dx,_mm_add_ps(one, _mm_mul_ps(u,u)))),
				// dy * u - dx * v
				_mm_sub_ps(_mm_mul_ps(u,dy), _mm_mul_ps(v,dx)),
				// a * (b0 - I) = (dE/da)
				_mm_mul_ps(a,_mm_sub_ps(b0, _mm_load_ps(buf_warped_refColor+i))),
				// -1  = (dE/db)
				minusOne,
				// Residual Energy before huber loss
				_mm_load_ps(buf_warped_residual+i),
				// Huber loss is implemented by multiplying by the huber weight
				_mm_load_ps(buf_warped_weight+i));
	}

	acc9.finish(); // Set the H and b matrix from the accumulated values

	// Extract H and b from the Hessian Structure matrix
	H_out = acc9.H.topLeftCorner<8,8>().cast<double>() * (1.0f/n);
	b_out = acc9.H.topRightCorner<8,1>().cast<double>() * (1.0f/n);

	// Scale H and b
	H_out.block<8,3>(0,0) *= SCALE_XI_ROT;
	H_out.block<8,3>(0,3) *= SCALE_XI_TRANS;
	H_out.block<8,1>(0,6) *= SCALE_A;
	H_out.block<8,1>(0,7) *= SCALE_B;
	H_out.block<3,8>(0,0) *= SCALE_XI_ROT;
	H_out.block<3,8>(3,0) *= SCALE_XI_TRANS;
	H_out.block<1,8>(6,0) *= SCALE_A;
	H_out.block<1,8>(7,0) *= SCALE_B;

	b_out.segment<3>(0) *= SCALE_XI_ROT;
	b_out.segment<3>(3) *= SCALE_XI_TRANS;
	b_out.segment<1>(6) *= SCALE_A;
	b_out.segment<1>(7) *= SCALE_B;
}

/**
 * @brief Calculates the coarse tracking residual
 * 
 * Only two-frame residual is done for coarse tracking
 * Note that the pattern is not used for coarse tracking
 * 
 * @param lvl 
 * @param refToNew 
 * @param aff_g2l 
 * @param cutoffTH 
 * @return Vec6 		Energy, number of points, Estimate of translation flow, 0 , Estimate of transformation flow, Percentage of high energy points
 */
Vec6 CoarseTracker::calcRes(int lvl, const SE3 &refToNew, AffLight aff_g2l, float cutoffTH)
{
	// Initialize variables
	float E = 0;
	int numTermsInE = 0;
	int numTermsInWarped = 0;
	int numSaturated=0;

	int wl = w[lvl];
	int hl = h[lvl];
	float fxl = fx[lvl];
	float fyl = fy[lvl];
	float cxl = cx[lvl];
	float cyl = cy[lvl];
	
	// Target frame is the newly inserted frame
	Eigen::Vector3f* dINewl = newFrame->dIp[lvl];


	Mat33f RKi = (refToNew.rotationMatrix().cast<float>() * Ki[lvl]);
	Vec3f t = (refToNew.translation()).cast<float>();
	Vec2f affLL = AffLight::fromToVecExposure(lastRef->ab_exposure, newFrame->ab_exposure, lastRef_aff_g2l, aff_g2l).cast<float>();


	float sumSquaredShiftT=0;
	float sumSquaredShiftRT=0;
	float sumSquaredShiftNum=0;

	// energy for r=setting_coarseCutoffTH
	float maxEnergy = 2*globalSettings.setting_huberTH*cutoffTH-globalSettings.setting_huberTH*globalSettings.setting_huberTH;


	MinimalImageB3* resImage = 0;

#ifdef GRAPHICAL_DEBUG
	if(globalSettings.setting_render_displayCoarseTrackingFull)
	{
		resImage = new MinimalImageB3(wl,hl);
		resImage->setConst(Vec3b(255,255,255));
	}
#endif

	int nl = pc_n[lvl];						// Number of points
	float* lpc_u = pc_u[lvl];				// u coordinates of points
	float* lpc_v = pc_v[lvl];				// v coordinates of points
	float* lpc_idepth = pc_idepth[lvl];		// depth of points
	float* lpc_color = pc_color[lvl];		// pixel intensity of the points


	// Start accumulating residuals
	for(int i=0;i<nl;i++) // for every point in the semi-dense depth and point map
	{
		// Get depth and position of point
		// Every active point from the last keyframe is projected onto the new frame
		float id = lpc_idepth[i];
		float u_pc_i = lpc_u[i];
		float v_pc_i = lpc_v[i];

		// Transform point and reproject
		// Transform to Cartesian coordinates and apply matrixes
		Vec3f pt = RKi * Vec3f(u_pc_i, v_pc_i, 1) + t*id;
		// Transform to Projective Image/Sensor coordinates
		float u_im_j = pt[0] / pt[2];
		float v_im_j = pt[1] / pt[2];
		// Transform to Projective Pixel coordinates
		float Ku_pc_j = fxl * u_im_j + cxl;
		float Kv_pc_j = fyl * v_im_j + cyl;
		float new_idepth = id/pt[2];

		// Cacluate overall flow indicators by sampling some points
		// Flow estimates are for keyframe selection
		if(lvl==0 && i%32==0)
		{
			// translation only (positive)
			Vec3f ptT = Ki[lvl] * Vec3f(u_pc_i, v_pc_i, 1) + t*id;
			float uT = ptT[0] / ptT[2];
			float vT = ptT[1] / ptT[2];
			float KuT = fxl * uT + cxl;
			float KvT = fyl * vT + cyl;

			// translation only (negative)
			Vec3f ptT2 = Ki[lvl] * Vec3f(u_pc_i, v_pc_i, 1) - t*id;
			float uT2 = ptT2[0] / ptT2[2];
			float vT2 = ptT2[1] / ptT2[2];
			float KuT2 = fxl * uT2 + cxl;
			float KvT2 = fyl * vT2 + cyl;

			// translation and rotation (negative)
			Vec3f pt3 = RKi * Vec3f(u_pc_i, v_pc_i, 1) - t*id;
			float u3 = pt3[0] / pt3[2];
			float v3 = pt3[1] / pt3[2];
			float Ku3 = fxl * u3 + cxl;
			float Kv3 = fyl * v3 + cyl;

			// translation and rotation (positive)
			// Is default transformation

			sumSquaredShiftT += (KuT-u_pc_i)*(KuT-u_pc_i) + (KvT-v_pc_i)*(KvT-v_pc_i);
			sumSquaredShiftT += (KuT2-u_pc_i)*(KuT2-u_pc_i) + (KvT2-v_pc_i)*(KvT2-v_pc_i);
			sumSquaredShiftRT += (Ku_pc_j-u_pc_i)*(Ku_pc_j-u_pc_i) + (Kv_pc_j-v_pc_i)*(Kv_pc_j-v_pc_i);
			sumSquaredShiftRT += (Ku3-u_pc_i)*(Ku3-u_pc_i) + (Kv3-v_pc_i)*(Kv3-v_pc_i);
			sumSquaredShiftNum+=2;
		}

		// Ignore if out of bounds
		if(!(Ku_pc_j > 2 && Kv_pc_j > 2 && Ku_pc_j < wl-3 && Kv_pc_j < hl-3 && new_idepth > 0)) continue;


		// Cacluate residual
		// Orignal pixel intensity
		float refColor = lpc_color[i];
		// Pixel intensity of transformed location
		Vec3f hitColor = getInterpolatedElement33(dINewl, Ku_pc_j, Kv_pc_j, wl, hl);
		if(!std::isfinite((float)hitColor[0])) continue;
		
		float residual = hitColor[0] - (float)(affLL[0] * refColor + affLL[1]);
		float hw = fabs(residual) < globalSettings.setting_huberTH ? 1 : globalSettings.setting_huberTH / fabs(residual);

		// Accumulate residuals
		if(fabs(residual) > cutoffTH) // energy is too high
		{
#ifdef GRAPHICAL_DEBUG
			if(globalSettings.setting_render_displayCoarseTrackingFull) resImage->setPixel4(lpc_u[i], lpc_v[i], Vec3b(0,0,255));
#endif
			// Add energy cutoff and note high energy point
			E += maxEnergy;
			numTermsInE++;
			numSaturated++;
		}
		else // energy is within range
		{
#ifdef GRAPHICAL_DEBUG
			if(globalSettings.setting_render_displayCoarseTrackingFull) resImage->setPixel4(lpc_u[i], lpc_v[i], Vec3b(residual+128,residual+128,residual+128));
#endif

			E += hw*residual*residual*(2-hw); // Accumulate huber norm
			numTermsInE++;

			// Record values of good point
			buf_warped_idepth[numTermsInWarped] = new_idepth;		// inverse depth
			buf_warped_u[numTermsInWarped] = u_im_j;						// new horizontal position
			buf_warped_v[numTermsInWarped] = v_im_j;						// new vertical position
			buf_warped_dx[numTermsInWarped] = hitColor[1];			// pixel intensity x derivative
			buf_warped_dy[numTermsInWarped] = hitColor[2];			// pixel intensity y derivative
			buf_warped_residual[numTermsInWarped] = residual;		// non huber energy residual
			buf_warped_weight[numTermsInWarped] = hw;				// huber weight
			buf_warped_refColor[numTermsInWarped] = lpc_color[i];	// pixel intensity in reference image
			numTermsInWarped++;
		}
	}

	// Align buffer
	while(numTermsInWarped%4!=0)
	{
		buf_warped_idepth[numTermsInWarped] = 0;
		buf_warped_u[numTermsInWarped] = 0;
		buf_warped_v[numTermsInWarped] = 0;
		buf_warped_dx[numTermsInWarped] = 0;
		buf_warped_dy[numTermsInWarped] = 0;
		buf_warped_residual[numTermsInWarped] = 0;
		buf_warped_weight[numTermsInWarped] = 0;
		buf_warped_refColor[numTermsInWarped] = 0;
		numTermsInWarped++;
	}
	buf_warped_n = numTermsInWarped;

#ifdef GRAPHICAL_DEBUG
	if(globalSettings.setting_render_displayCoarseTrackingFull&& !globalSettings.setting_disableAllDisplay)
	{
		IOWrap::displayImage("RES", resImage, 8);
		IOWrap::waitKey(0);
		delete resImage;
	}
#endif


	Vec6 rs;
	rs[0] = E;												// Residual energy of transform estimation
	rs[1] = numTermsInE;									// Number of points
	rs[2] = sumSquaredShiftT/(sumSquaredShiftNum+0.1);		// Estimate of translation flow
	rs[3] = 0;
	rs[4] = sumSquaredShiftRT/(sumSquaredShiftNum+0.1);		// Estimate of transformation flow
	rs[5] = numSaturated / (float)numTermsInE;				// Percentage of high energy points

	return rs;
}


/**
 * @brief Setup
 * 
 * Sets the last reference frame as the newest keyframe
 * 
 * @param frameHessians 
 */
void CoarseTracker::setCoarseTrackingRef(
		std::vector<FrameHessian*> frameHessians)
{
	assert(frameHessians.size()>0);
	lastRef = frameHessians.back();
	makeCoarseDepthL0(frameHessians);


	refFrameID = lastRef->shell->id;
	lastRef_aff_g2l = lastRef->aff_g2l();

	firstCoarseRMSE=-1;
}

/**
 * @brief Coarse frame tracking function
 * 
 * Only two-frame tracking is done for coarse tracking
 * 
 * @param newFrameHessian 	New frame structure
 * @param lastToNew_out 	Motion prediction to initalize
 * @param aff_g2l_out 		Lighting prediction to intialize
 * @param coarsestLvl 
 * @param minResForAbort 	Best previously achieved residual for current track
 * @param wrap 
 * 
 * @return true 
 * @return false 
 */
bool CoarseTracker::trackNewestCoarse(
		FrameHessian* newFrameHessian,
		SE3 &lastToNew_out, AffLight &aff_g2l_out,
		int coarsestLvl,
		Vec5 minResForAbort,
		IOWrap::Output3DWrapper* wrap)
{
	debugPrint = !setting_debugout_runquiet;

	assert(coarsestLvl < 5 && coarsestLvl < globalSettings.pyrLevelsUsed);

	lastResidualsStats.setConstant(NAN);
	lastFlowIndicators.setConstant(1000);

	newFrame = newFrameHessian;
	int maxIterations[] = {10,20,50,50,50};
	float lambdaExtrapolationLimit = 0.001;

	SE3 refToNew_current = lastToNew_out;
	AffLight aff_g2l_current = aff_g2l_out;

	bool haveRepeated = false;


	Mat88 H; Vec8 b;
	int lastLvl = -1;
	// Do all the pyramid levels from the coarsest to the normal image
	for(int lvl=coarsestLvl; lvl>=0; lvl--)
	{
		// Do initial calculation

		float levelCutoffRepeat=1;
		// Calulate residual of transformation estimate
		Vec6 resOld = calcRes(lvl, refToNew_current, aff_g2l_current, globalSettings.setting_coarseCutoffTH*levelCutoffRepeat);

		while(resOld[5] > 0.6 && (levelCutoffRepeat < 50 || resOld[5] > 0.99) ) // If too many high energy points
		{
			// Try again with lower high energy point threshold
			levelCutoffRepeat*=2;
			resOld = calcRes(lvl, refToNew_current, aff_g2l_current, globalSettings.setting_coarseCutoffTH*levelCutoffRepeat);

			if(!setting_debugout_runquiet && !globalSettings.no_CoarseTracker_debugMessage)
				printf("INCREASING cutoff to %f (ratio is %f)!\n", globalSettings.setting_coarseCutoffTH*levelCutoffRepeat, resOld[5]);
		}

		// Calculate H and b for Gauss Newton Optimization
		// Coarse tracker only optimizes the pose and photometric variables
		// point depths and the calibration matrix are not optimized
		calcGSSSE(lvl, H, b, refToNew_current, aff_g2l_current);

		float lambda = 0.01;

		if(debugPrint && !globalSettings.no_CoarseTracker_debugMessage)
		{
			Vec2f relAff = AffLight::fromToVecExposure(lastRef->ab_exposure, newFrame->ab_exposure, lastRef_aff_g2l, aff_g2l_current).cast<float>();
			printf("lvl%d, it %d (l=%f / %f) %s: %.3f->%.3f (%d -> %d) (|inc| = %f)! \t",
					lvl, -1, lambda, 1.0f,
					"INITIA",
					0.0f,
					resOld[0] / resOld[1],
					 0,(int)resOld[1],
					0.0f);
			std::cout << refToNew_current.log().transpose() << " AFF " << aff_g2l_current.vec().transpose() <<" (rel " << relAff.transpose() << ")\n";
		}

		
		// Do iterations
		for(int iteration=0; iteration < maxIterations[lvl]; iteration++)
		{
			dmvio::TimeMeasurement timeMeasurement("coarseTrackingIteration");
			Mat88 Hl = H;

			// Multiply diagonal by current lambda value
			for(int i=0;i<8;i++) Hl(i,i) *= (1+lambda);

			float extrapFac = 1;
			if(lambda < lambdaExtrapolationLimit) extrapFac = sqrt(sqrt(lambdaExtrapolationLimit / lambda));

			SE3 refToNew_new;
			AffLight aff_g2l_new = aff_g2l_current;
			double incNorm;
			if(imuIntegration.setting_useIMU && imuIntegration.isCoarseInitialized())
			{
				// imu!: The idea of the integration of the IMU (and GTSAM) into the coarse tracking is to replace the line
				// Vec8 inc = Hl.ldlt().solve(-b);
				// with a call to computeCoarseUpdate, which will add GTSAM factors before calculating the update.

				double incA, incB;
				// Note that we pass H instead of Hl as the lambda multiplication is done inside...
				refToNew_new = imuIntegration.computeCoarseUpdate(H, b, extrapFac, lambda, incA, incB, incNorm);

				SE3 oldVal = refToNew_current;
				SE3 newVal = refToNew_new;
				dso::Vec6 increment = (newVal * oldVal.inverse()).log();

				dso::Vec8 totalIncrement;
				totalIncrement.segment(0, 6) = increment;

				totalIncrement(6) = incA;
				totalIncrement(7) = incB;

				incA *= SCALE_A;
				incB *= SCALE_B;

				aff_g2l_new.a += incA;
				aff_g2l_new.b += incB;
			}else
			{
				// Calculate increment
				// inc = [w1, w2, w3, d1, d2, d3, a, b]
				Vec8 inc = Hl.ldlt().solve(-b);

				if(globalSettings.setting_affineOptModeA < 0 && globalSettings.setting_affineOptModeB < 0)	// fix a, b
				{
					inc.head<6>() = Hl.topLeftCorner<6,6>().ldlt().solve(-b.head<6>());
					inc.tail<2>().setZero();
				}
				if(!(globalSettings.setting_affineOptModeA < 0) && globalSettings.setting_affineOptModeB < 0)	// fix b
				{
					inc.head<7>() = Hl.topLeftCorner<7,7>().ldlt().solve(-b.head<7>());
					inc.tail<1>().setZero();
				}
				if(globalSettings.setting_affineOptModeA < 0 && !(globalSettings.setting_affineOptModeB < 0))	// fix a
				{
					Mat88 HlStitch = Hl;
					Vec8 bStitch = b;
					HlStitch.col(6) = HlStitch.col(7);
					HlStitch.row(6) = HlStitch.row(7);
					bStitch[6] = bStitch[7];
					Vec7 incStitch = HlStitch.topLeftCorner<7,7>().ldlt().solve(-bStitch.head<7>());
					inc.setZero();
					inc.head<6>() = incStitch.head<6>();
					inc[6] = 0;
					inc[7] = incStitch[6];
				}

				// Scale increment
				inc *= extrapFac;
				Vec8 incScaled = inc;
				incScaled.segment<3>(0) *= SCALE_XI_ROT;
				incScaled.segment<3>(3) *= SCALE_XI_TRANS;
				incScaled.segment<1>(6) *= SCALE_A;
				incScaled.segment<1>(7) *= SCALE_B;

				if(!std::isfinite(incScaled.sum())) incScaled.setZero();

				// Calculate new pose
				// exp: first three: translational part, last three: rotational part.
				// Note: gtsam::Pose3 contains first rotational and then translational part!
				refToNew_new = SE3::exp((Vec6) (incScaled.head<6>())) * refToNew_current;
				aff_g2l_new = aff_g2l_current;
				aff_g2l_new.a += incScaled[6];
				aff_g2l_new.b += incScaled[7];

				incNorm = inc.norm();
			}

			// Cacluate residual for new pose
			Vec6 resNew = calcRes(lvl, refToNew_new, aff_g2l_new, globalSettings.setting_coarseCutoffTH*levelCutoffRepeat);

			// Accept if residual energy per point lowers
			bool accept = (resNew[0] / resNew[1]) < (resOld[0] / resOld[1]);

			if(debugPrint && !globalSettings.no_CoarseTracker_debugMessage)
			{
				Vec2f relAff = AffLight::fromToVecExposure(lastRef->ab_exposure, newFrame->ab_exposure, lastRef_aff_g2l, aff_g2l_new).cast<float>();
				printf("lvl %d, it %d (l=%f / %f) %s: %.3f->%.3f (%d -> %d) (|inc| = %f)! \t",
						lvl, iteration, lambda,
						extrapFac,
						(accept ? "ACCEPT" : "REJECT"),
						resOld[0] / resOld[1],
						resNew[0] / resNew[1],
						(int)resOld[1], (int)resNew[1],
						incNorm);
				std::cout << refToNew_new.log().transpose() << " AFF " << aff_g2l_new.vec().transpose() <<" (rel " << relAff.transpose() << ")\n";
			}

			if(accept) // if energy lowers, do another iteration with lower lambda
			{
				calcGSSSE(lvl, H, b, refToNew_new, aff_g2l_new);

				resOld = resNew;
				aff_g2l_current = aff_g2l_new;
				refToNew_current = refToNew_new;
				
				// imu!: Calculate the update with IMU factors
				if(imuIntegration.setting_useIMU)
					imuIntegration.acceptCoarseUpdate();

				lambda *= 0.5;
			}
			else // if energy increases or stays the same, increase lambda
			{
				lambda *= 4;
				if(lambda < lambdaExtrapolationLimit) lambda = lambdaExtrapolationLimit;
			}

			lastLvl = lvl;

			if(!(incNorm > 1e-3)) // break if increment is too small
			{
				if(debugPrint && !globalSettings.no_CoarseTracker_debugMessage)
					printf("inc too small, break!\n");
				break;
			}
		}

		// Set last residual for that level, as well as flow indicators.
		lastResidualsStats[lvl] = sqrtf((float)(resOld[0] / resOld[1]));
		lastFlowIndicators = resOld.segment<3>(2);

		// Track is bad if the residual is invalid or more than many times the previously achieved residual
		if(std::isnan(lastResidualsStats[lvl])) return false;
		if(lastResidualsStats[lvl] > 1.5*minResForAbort[lvl]) return false;


		if(levelCutoffRepeat > 1 && !haveRepeated)
		{
			lvl++;
			haveRepeated=true;
			if(!setting_debugout_runquiet && !globalSettings.no_CoarseTracker_debugMessage) printf("REPEAT LEVEL!\n");
		}
	}

	// Set new pose and photometric values
	lastToNew_out = refToNew_current;
	aff_g2l_out = aff_g2l_current;

	// Set tracking as bad if any of the final values are too large or small
	bool trackingGood = true;

	if((globalSettings.setting_affineOptModeA != 0 && (fabsf(aff_g2l_out.a) > 1.2))
	|| (globalSettings.setting_affineOptModeB != 0 && (fabsf(aff_g2l_out.b) > 200)))
		trackingGood = false;

	Vec2f relAff = AffLight::fromToVecExposure(lastRef->ab_exposure, newFrame->ab_exposure, lastRef_aff_g2l, aff_g2l_out).cast<float>();

	if((globalSettings.setting_affineOptModeA == 0 && (fabsf(logf((float)relAff[0])) > 1.5))
	|| (globalSettings.setting_affineOptModeB == 0 && (fabsf((float)relAff[1]) > 200)))
		trackingGood = false;

	if(globalSettings.setting_affineOptModeA < 0) aff_g2l_out.a=0;
	if(globalSettings.setting_affineOptModeB < 0) aff_g2l_out.b=0;

	// imu!: Add visual to coarse graph of last level is zero
	if(lastLvl == 0)
	{
		if(imuIntegration.setting_useIMU)
			imuIntegration.addVisualToCoarseGraph(H, b, trackingGood);
	}

	return trackingGood;
}


/**
 * @brief Creates the depth map for the GUI
 * 
 * @param minID_pt 
 * @param maxID_pt 
 * @param wraps 
 */
void CoarseTracker::debugPlotIDepthMap(float* minID_pt, float* maxID_pt, std::vector<IOWrap::Output3DWrapper*> &wraps) const
{
	dmvio::TimeMeasurement timeMeasurement("debugPlotIDepthMap");
	if(wraps.empty() && !globalSettings.setting_debugSaveImages)
	{
		return;
	}
	if(w[1] == 0) return;

	{
		int lvl = 0;

		std::vector<float> allID;
		for(int i=0;i<h[lvl]*w[lvl];i++)
		{
			if(idepth[lvl][i] > 0)
				allID.push_back(idepth[lvl][i]);
		}
		std::sort(allID.begin(), allID.end());
		int n = allID.size()-1;
		if(n <= 0)
		{
			return;
		}

		float minID_new = allID[(int)(n*0.05)];
		float maxID_new = allID[(int)(n*0.95)];

		float minID, maxID;
		minID = minID_new;
		maxID = maxID_new;
		if(minID_pt!=0 && maxID_pt!=0)
		{
			if(*minID_pt < 0 || *maxID_pt < 0)
			{
				*maxID_pt = maxID;
				*minID_pt = minID;
			}
			else
			{
				// slowly adapt: change by maximum 10% of old span.
				float maxChange = 0.3*(*maxID_pt - *minID_pt);

				if(minID < *minID_pt - maxChange)
					minID = *minID_pt - maxChange;
				if(minID > *minID_pt + maxChange)
					minID = *minID_pt + maxChange;


				if(maxID < *maxID_pt - maxChange)
					maxID = *maxID_pt - maxChange;
				if(maxID > *maxID_pt + maxChange)
					maxID = *maxID_pt + maxChange;

				*maxID_pt = maxID;
				*minID_pt = minID;
			}
		}

		MinimalImageB3 mf(w[lvl], h[lvl]);
		mf.setBlack();
		for(int i=0;i<h[lvl]*w[lvl];i++)
		{
			int c = lastRef->dIp[lvl][i][0]*0.9f;
			if(c>255) c=255;
			mf.at(i) = Vec3b(c,c,c);
		}
		int wl = w[lvl];
		for(int y=3;y<h[lvl]-3;y++)
			for(int x=3;x<wl-3;x++)
			{
				int idx=x+y*wl;
				float sid=0, nid=0;
				float* bp = idepth[lvl]+idx;

				if(bp[0] > 0) {sid+=bp[0]; nid++;}
				if(bp[1] > 0) {sid+=bp[1]; nid++;}
				if(bp[-1] > 0) {sid+=bp[-1]; nid++;}
				if(bp[wl] > 0) {sid+=bp[wl]; nid++;}
				if(bp[-wl] > 0) {sid+=bp[-wl]; nid++;}

				if(bp[0] > 0 || nid >= 3)
				{
					float id = ((sid / nid)-minID) / ((maxID-minID));
					mf.setPixelCirc(x,y,makeJet3B(id));
					//mf.at(idx) = makeJet3B(id);
				}
			}


		for(IOWrap::Output3DWrapper* ow : wraps)
			ow->pushDepthImage(&mf);

		if(globalSettings.setting_debugSaveImages)
		{
			char buf[1024];
			snprintf(buf, 1024, "images_out/predicted_%05d_%05d.png", lastRef->shell->id, refFrameID);
			IOWrap::writeImage(buf,&mf);
		}
	}
}

/**
 * @brief Creates the depth map for the GUI
 * 
 * @param wraps 
 */
void CoarseTracker::debugPlotIDepthMapFloat(std::vector<IOWrap::Output3DWrapper*> &wraps)
{
	dmvio::TimeMeasurement timeMeasurement("debugPlotIDepthMapFloat");
	if(w[1] == 0) return;
	int lvl = 0;
	MinimalImageF mim(w[lvl], h[lvl], idepth[lvl]);
	for(IOWrap::Output3DWrapper* ow : wraps)
		ow->pushDepthImageFloat(&mim, lastRef);
}



/**
 * @brief Construct a new Coarse Distance Map:: Coarse Distance Map object
 * 
 * The coarse distance map shows how far each pixel in a frame is from a valid point
 * The coarse distance map is used to maintain immature point density
 * 
 * @param globalCalib_ 
 * @param pyrLevelsUsed_ 
 * 
 */
CoarseDistanceMap::CoarseDistanceMap(Global_Calib& globalCalib_, int pyrLevelsUsed_):
globalCalib(globalCalib_), pyrLevelsUsed(pyrLevelsUsed_)
{
	wG0 = globalCalib_.wG[0];
	hG0 = globalCalib_.hG[0];
	unsigned int wh = wG0*hG0;

	fwdWarpedIDDistFinal = new float[wh/4];

	bfsList1 = new Eigen::Vector2i[wh/4];
	bfsList2 = new Eigen::Vector2i[wh/4];

	int fac = 1 << (pyrLevelsUsed-1);


	coarseProjectionGrid = new PointFrameResidual*[2048*(wh/(fac*fac))];
	coarseProjectionGridNum = new int[wh/(fac*fac)];

	w[0]=h[0]=0;
}

/**
 * @brief Destroy the Coarse Distance Map
 * 
 */
CoarseDistanceMap::~CoarseDistanceMap()
{
	delete[] fwdWarpedIDDistFinal;
	delete[] bfsList1;
	delete[] bfsList2;
	delete[] coarseProjectionGrid;
	delete[] coarseProjectionGridNum;
}

/**
 * @brief Makes the distance map
 * 
 * @param frameHessians 
 * @param frame 
 */
void CoarseDistanceMap::makeDistanceMap(
		std::vector<FrameHessian*> frameHessians,
		FrameHessian* frame)
{
	int w1 = w[1];
	int h1 = h[1];
	int wh1 = w1*h1;
	for(int i=0;i<wh1;i++)
		fwdWarpedIDDistFinal[i] = 1024;


	// make coarse tracking templates for latstRef.
	int numItems = 0;

	for(FrameHessian* fh : frameHessians) // for all active frames
	{
		if(frame == fh) continue;

		// Get transformation matrix
		SE3 fhToNew = frame->PRE_worldToCam * fh->PRE_camToWorld;
		Mat33f KRKi = (K[1] * fhToNew.rotationMatrix().cast<float>() * Ki[0]);
		Vec3f Kt = (K[1] * fhToNew.translation().cast<float>());

		for(PointHessian* ph : fh->pointHessians) // for all points in frame
		{
			assert(ph->status == PointHessian::ACTIVE);
			// Get transformed position of points
			Vec3f ptp = KRKi * Vec3f(ph->u, ph->v, 1) + Kt*ph->idepth_scaled;
			int u = ptp[0] / ptp[2] + 0.5f;
			int v = ptp[1] / ptp[2] + 0.5f;

			if(!(u > 0 && v > 0 && u < w[1] && v < h[1])) continue;
			fwdWarpedIDDistFinal[u+w1*v]=0;
			bfsList1[numItems] = Eigen::Vector2i(u,v);
			numItems++;
		}
	}

	growDistBFS(numItems);
}

void CoarseDistanceMap::makeInlierVotes(std::vector<FrameHessian*> frameHessians)
{
}

/**
 * @brief Grows distance map
 * 
 * @param bfsNum 
 */
void CoarseDistanceMap::growDistBFS(int bfsNum)
{
	assert(w[0] != 0);
	int w1 = w[1], h1 = h[1];

	for(int k=1;k<64;k++) // Only update pixels a certain number of pixels away from the point
	{
		int bfsNum2 = bfsNum;
		std::swap<Eigen::Vector2i*>(bfsList1,bfsList2);
		bfsNum=0;

		if(k%2==0) // Update only top, left, right and bottom pixels
		{
			for(int i=0;i<bfsNum2;i++)
			{
				int x = bfsList2[i][0];
				int y = bfsList2[i][1];

				if(x==0 || y== 0 || x==w1-1 || y==h1-1) continue;
				int idx = x + y * w1;

				if(fwdWarpedIDDistFinal[idx+1] > k)
				{
					fwdWarpedIDDistFinal[idx+1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x+1,y); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx-1] > k)
				{
					fwdWarpedIDDistFinal[idx-1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x-1,y); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx+w1] > k)
				{
					fwdWarpedIDDistFinal[idx+w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x,y+1); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx-w1] > k)
				{
					fwdWarpedIDDistFinal[idx-w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x,y-1); bfsNum++;
				}
			}
		}
		else // Update all surrounding pixels
		{
			for(int i=0;i<bfsNum2;i++)
			{
				int x = bfsList2[i][0];
				int y = bfsList2[i][1];

				if(x==0 || y== 0 || x==w1-1 || y==h1-1) continue;
				int idx = x + y * w1;

				if(fwdWarpedIDDistFinal[idx+1] > k)
				{
					fwdWarpedIDDistFinal[idx+1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x+1,y); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx-1] > k)
				{
					fwdWarpedIDDistFinal[idx-1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x-1,y); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx+w1] > k)
				{
					fwdWarpedIDDistFinal[idx+w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x,y+1); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx-w1] > k)
				{
					fwdWarpedIDDistFinal[idx-w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x,y-1); bfsNum++;
				}

				if(fwdWarpedIDDistFinal[idx+1+w1] > k)
				{
					fwdWarpedIDDistFinal[idx+1+w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x+1,y+1); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx-1+w1] > k)
				{
					fwdWarpedIDDistFinal[idx-1+w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x-1,y+1); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx-1-w1] > k)
				{
					fwdWarpedIDDistFinal[idx-1-w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x-1,y-1); bfsNum++;
				}
				if(fwdWarpedIDDistFinal[idx+1-w1] > k)
				{
					fwdWarpedIDDistFinal[idx+1-w1] = k;
					bfsList1[bfsNum] = Eigen::Vector2i(x+1,y-1); bfsNum++;
				}
			}
		}
	}
}

/**
 * @brief Add new distance to map
 * 
 * @param u 
 * @param v 
 */
void CoarseDistanceMap::addIntoDistFinal(int u, int v)
{
	if(w[0] == 0) return;
	bfsList1[0] = Eigen::Vector2i(u,v);
	fwdWarpedIDDistFinal[u+w[1]*v] = 0;
	growDistBFS(1);
}

/**
 * @brief Set width, height, and camera parameters for all pyramid levels
 * 
 * @param HCalib 
 */
void CoarseDistanceMap::makeK(CalibHessian* HCalib)
{
	w[0] = wG0;
	h[0] = hG0;

	fx[0] = HCalib->fxl();
	fy[0] = HCalib->fyl();
	cx[0] = HCalib->cxl();
	cy[0] = HCalib->cyl();

	for (int level = 1; level < pyrLevelsUsed; ++ level)
	{
		w[level] = w[0] >> level;
		h[level] = h[0] >> level;
		fx[level] = fx[level-1] * 0.5;
		fy[level] = fy[level-1] * 0.5;
		cx[level] = (cx[0] + 0.5) / ((int)1<<level) - 0.5;
		cy[level] = (cy[0] + 0.5) / ((int)1<<level) - 0.5;
	}

	for (int level = 0; level < pyrLevelsUsed; ++ level)
	{
		K[level]  << fx[level], 0.0, cx[level], 
					0.0, fy[level], cy[level], 
					0.0, 0.0, 1.0;
		Ki[level] = K[level].inverse();
		fxi[level] = Ki[level](0,0);
		fyi[level] = Ki[level](1,1);
		cxi[level] = Ki[level](0,2);
		cyi[level] = Ki[level](1,2);
	}
}
}
