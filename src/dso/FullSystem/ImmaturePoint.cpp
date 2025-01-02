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



#include "FullSystem/ImmaturePoint.h"
#include "util/FrameShell.h"
#include "FullSystem/ResidualProjections.h"



namespace dso
{

/**
 * @brief Construct a new Immature Point:: Immature Point object
 * 
 * @param u_ 		x position of point
 * @param v_ 		y position of point
 * @param globalCalib_ 
 * @param host_ 	Frame containing point
 * @param type 		Pyramid level point is detected in
 * @param HCalib 
 * @param globalSettings_ 
 */
ImmaturePoint::ImmaturePoint(int u_, int v_, Global_Calib& globalCalib_, FrameHessian* host_, float type, CalibHessian* HCalib, GlobalSettings& globalSettings_)
: globalSettings(globalSettings_), globalCalib(globalCalib_), host(host_), my_type(type), idepth_min(0), idepth_max(NAN), lastTraceStatus(IPS_UNINITIALIZED), Point()
{
	u = u_;
	v = v_;
	hostFrameID = host_->frameID;

	wG0 = globalCalib.wG[0];
	hG0 = globalCalib.hG[0];

	gradH.setZero();

	colourValid = false;
	if(host->colourValid) colourValid = true;

	// Add all pixels in the pattern to the point
	for(int idx=0;idx<PATTERNNUM;idx++)
	{
		int dx = PATTERNP[idx][0];
		int dy = PATTERNP[idx][1];

		// ptc is (pixel intensity, dx, dy)
        Vec3f ptc = getInterpolatedElement33BiLin(host->dI, u+dx, v+dy, wG0, hG0);
		color[idx] = ptc[0]; // Set pixel internsity

		if(colourValid){
			int u_int = static_cast<int>(u);
			int v_int = static_cast<int>(v);
			Vec3f ptc_3 = host->dI_c[u_int+dx+(v_int+dy)*wG0];
			colour3[idx] = ptc_3;
		}

		if(!std::isfinite(color[idx])) {energyTH=NAN; return;}

		// 2 by 2 outer product matrix
		// [I_dx^2 I_dx*I_dxy]
		// [I_dx*I_dxy I_dy^2]
		gradH += ptc.tail<2>() * ptc.tail<2>().transpose(); 
		// Weights for gradient-dependent weighting
		weights[idx] = sqrtf(globalSettings.setting_outlierTHSumComponent / (globalSettings.setting_outlierTHSumComponent + ptc.tail<2>().squaredNorm()));
	}

	energyTH = PATTERNNUM*globalSettings.setting_outlierTH;
	energyTH *= globalSettings.setting_overallEnergyTHWeight*globalSettings.setting_overallEnergyTHWeight;

	idepth_GT=0;
	quality=10000;
}

/**
 * @brief Destroy the Immature object
 * 
 */
ImmaturePoint::~ImmaturePoint()
{
}

/**
 * @brief Traces the points as the frames progress
 * 
 * Only the depth of the point is updated
 * 
 * Steps:
 * 1. Project (u,v) point assuming min and max depth
 * 2. Search along line between the max and min depth points uniformly for lowest energy
 * 3. Do Gauss Newton optimization with found point to get a more precise point
 * 
 * @param frame 				// Frame that point in contained in
 * @param hostToFrame_KRKi		// Rotation matrix from frame movement
 * @param hostToFrame_Kt 		// Translation matrix from frame movement
 * @param hostToFrame_affine	// Constants for photogrammetric conversion
 * @param HCalib 
 * @param debugPrint 
 * @return ImmaturePointStatus: 
 * OOB -> point is optimized and marginalized.
 * UPDATED -> point has been updated.
 * SKIP -> point has not been updated.
 */
ImmaturePointStatus ImmaturePoint::traceOn(FrameHessian* frame,const Mat33f &hostToFrame_KRKi, const Vec3f &hostToFrame_Kt, const Vec2f& hostToFrame_affine, CalibHessian* HCalib, bool debugPrint)
{
	// Do not trace immature points that are OOB
	if(lastTraceStatus == ImmaturePointStatus::IPS_OOB ||
		lastTraceStatus == ImmaturePointStatus::IPS_OUTLIER_OUT) return lastTraceStatus;

	float maxPixSearch = (wG0+hG0)*globalSettings.setting_maxPixSearch;

	if(debugPrint && !globalSettings.no_Immature_debugMessage)
		printf("trace pt (%.1f %.1f) from frame %d to %d. Range %f -> %f. t %f %f %f!\n",
				u,v,
				host->shell->id, frame->shell->id,
				idepth_min, idepth_max,
				hostToFrame_Kt[0],hostToFrame_Kt[1],hostToFrame_Kt[2]);


	// ============== Check if points are OOB. Project assumming min and max depth and then check bounds ===================
	// Find the maximum size of the pattern
    Mat22f Rplane = hostToFrame_KRKi.topLeftCorner<2,2>();
    int maxRotPatX = 0;
    int maxRotPatY = 0;
    Vec2f rotatetPattern[MAX_RES_PER_POINT];
    for(int idx=0;idx<PATTERNNUM;idx++)
    {
		// Calculate the new u,v positions of the pattern
		// Calculations done on viewing plane, so depth is assumed to be zero
        rotatetPattern[idx] = Rplane * Vec2f(PATTERNP[idx][0], PATTERNP[idx][1]);
        int absX = (int) abs(rotatetPattern[idx][0]);
        int absY = (int) abs(rotatetPattern[idx][1]);
        maxRotPatX = std::max(absX, maxRotPatX);
        maxRotPatY = std::max(absY, maxRotPatY);
    }
    int realBoundU = maxRotPatX + 2;
    int realBoundV = maxRotPatY + 2;
    int boundU = 4;
    int boundV = 4;
    boundU = std::max(boundU, realBoundU);
    boundV = std::max(boundV, realBoundV);
	
	// Unproject (u,v) (depth does not matter for rotation), rotate, and the re-project
	Vec3f pr = hostToFrame_KRKi * Vec3f(u,v,1); 
	
	// Do min depth first
	// Translation is dependent on depth, which is unknown, so depth must be estimated
	Vec3f ptpMin = pr + hostToFrame_Kt*idepth_min;	// Translate rotated (u,v) point assuming depth is idepth_min

	float uMin = ptpMin[0] / ptpMin[2];
	float vMin = ptpMin[1] / ptpMin[2];

	if(!(uMin > boundU && vMin > boundV && uMin < wG0-boundU-1 && vMin < hG0-boundV-1)) // Pattern is OOB
	{
		if(debugPrint && !globalSettings.no_Immature_debugMessage)
			printf("OOB uMin %f %f - %f %f %f (id %f-%f)!\n", u,v,uMin, vMin,  ptpMin[2], idepth_min, idepth_max);
		lastTraceUV = Vec2f(-1,-1);
		lastTracePixelInterval=0;
		return lastTraceStatus = ImmaturePointStatus::IPS_OOB;
	}

	// Now do max depth
	float dist;
	float uMax;
	float vMax;
	Vec3f ptpMax;
	if(std::isfinite(idepth_max)) // If last max depth was valid
	{
		ptpMax = pr + hostToFrame_Kt*idepth_max; // Repeat, but assuming the depth is idepth_max
		uMax = ptpMax[0] / ptpMax[2];
		vMax = ptpMax[1] / ptpMax[2];

		if(!(uMax > boundU && vMax > boundV && uMax < wG0-boundU-1 && vMax < hG0-boundV-1)) // Pattern is OOB
		{
			if(debugPrint && !globalSettings.no_Immature_debugMessage)
				printf("OOB uMax  %f %f - %f %f!\n",u,v, uMax, vMax);
			lastTraceUV = Vec2f(-1,-1);
			lastTracePixelInterval=0;
			return lastTraceStatus = ImmaturePointStatus::IPS_OOB;
		}

		// ============== check their distance. everything below 2px is OK (-> skip). ===================
		dist = (uMin-uMax)*(uMin-uMax) + (vMin-vMax)*(vMin-vMax);
		dist = sqrtf(dist);
		if(dist < globalSettings.setting_trace_slackInterval) // Tracing is not needed if the distance between max and min is low
		{
			if(debugPrint && !globalSettings.no_Immature_debugMessage)
				printf("TOO CERTAIN ALREADY (dist %f)!\n", dist);

			lastTraceUV = Vec2f(uMax+uMin, vMax+vMin)*0.5; // Take average
			lastTracePixelInterval=dist;
			return lastTraceStatus = ImmaturePointStatus::IPS_SKIPPED;
		}
		assert(dist>0);
	}
	else // Intial search will start off with max depth at infinite
	{
		dist = maxPixSearch; // Set length of search line to maxPixSearch

		// Project to arbitrary depth to get direction
		ptpMax = pr + hostToFrame_Kt*0.01;

		uMax = ptpMax[0] / ptpMax[2];
		vMax = ptpMax[1] / ptpMax[2];

		float dx = uMax-uMin;
		float dy = vMax-vMin;
		float d = 1.0f / sqrtf(dx*dx+dy*dy);

		// Set to [setting_maxPixSearch].
		uMax = uMin + dist*dx*d;
		vMax = vMin + dist*dy*d;

		if(!(uMax > boundU && vMax > boundV && uMax < wG0-boundU-1 && vMax < hG0-boundV-1)) // Pattern may still be OOB
		{
			if(debugPrint && !globalSettings.no_Immature_debugMessage)
				printf("OOB uMax-coarse %f %f %f!\n", uMax, vMax,  ptpMax[2]);
			lastTraceUV = Vec2f(-1,-1);
			lastTracePixelInterval=0;
			return lastTraceStatus = ImmaturePointStatus::IPS_OOB;
		}
		assert(dist>0);
	}

	// ptpMin[2] is approximately 1+t_z/d
	// Do not accept transformations where the z change is large compared to the depth
	if(!(idepth_min<0 || (ptpMin[2]>0.75 && ptpMin[2]<1.5))) 
	{
		if(debugPrint && !globalSettings.no_Immature_debugMessage)
			printf("OOB SCALE %f %f %f!\n", uMax, vMax,  ptpMin[2]);
		lastTraceUV = Vec2f(-1,-1);
		lastTracePixelInterval=0;
		return lastTraceStatus = ImmaturePointStatus::IPS_OOB;
	}


	// ============== compute error-bounds on result in pixel ===================
	float dx = globalSettings.setting_trace_stepsize*(uMax-uMin);
	float dy = globalSettings.setting_trace_stepsize*(vMax-vMin);

	// I_dx^2*dx^2 + 2*dx*dy*I_dx*I_dy + I_dy^2*dy^2
	float a_err = (Vec2f(dx,dy).transpose() * gradH * Vec2f(dx,dy));
	// I_dx^2*dy^2 + I_dx*I_dy(dx^2+dy^2) + I_dy^2*dx^2
	float b_err = (Vec2f(dy,-dx).transpose() * gradH * Vec2f(dy,-dx));
	// errorInPixel is a function that grows extremely fast if the dx/dy isn't
	// aligned with I_dx/Id_y or it's opposite direction
	// It measures the alignment with the pixel derivative
	// Searching along the pixel edge should be difficult, 
	// so the point is considered IPS_BADCONDITION if the alignment too off
	float errorInPixel = 0.2f + 0.2f * (a_err+b_err) / a_err;

	if(errorInPixel*globalSettings.setting_trace_minImprovementFactor > dist && std::isfinite(idepth_max))
	{
		if(debugPrint && !globalSettings.no_Immature_debugMessage)
			printf("NO SIGNIFICANT IMPROVMENT (%f)!\n", errorInPixel);
		lastTraceUV = Vec2f(uMax+uMin, vMax+vMin)*0.5;
		lastTracePixelInterval=dist;
		return lastTraceStatus = ImmaturePointStatus::IPS_BADCONDITION;
	}

	if(errorInPixel > 10) errorInPixel = 10;


	// ============== do the discrete search ===================
	// Find a good starting point for the Gauss Newton optimization
	dx /= dist;
	dy /= dist;

	if(debugPrint && !globalSettings.no_Immature_debugMessage)
		printf("trace pt (%.1f %.1f) from frame %d to %d. Range %f (%.1f %.1f) -> %f (%.1f %.1f)! ErrorInPixel %.1f!\n",
				u,v,
				host->shell->id, frame->shell->id,
				idepth_min, uMin, vMin,
				idepth_max, uMax, vMax,
				errorInPixel
				);

	if(dist > maxPixSearch)
	{
		uMax = uMin + maxPixSearch*dx;
		vMax = vMin + maxPixSearch*dy;
		dist = maxPixSearch;
	}

	int numSteps = 1.9999f + dist / globalSettings.setting_trace_stepsize;

	// Add some randomness
	float randShift = uMin*1000-floorf(uMin*1000);
	float ptx = uMin-randShift*dx;
	float pty = vMin-randShift*dy;

    if(!std::isfinite(dx) || !std::isfinite(dy)) // OOD if dx or dy is infinite
	{
		printf("CAUGHT INF / NAN dxdy!\n");

		lastTraceUV = Vec2f(-1,-1);
		lastTracePixelInterval=0;
		return lastTraceStatus = ImmaturePointStatus::IPS_OOB;
	}

	// Start linearly searching starting from the uMin, vMin values
	// Get point between the min and max depth points with the lowest energy
	float errors[100];
	float bestU=0, bestV=0, bestEnergy=1e10;
	int bestIdx=-1;
	if(numSteps >= 100) numSteps = 99;
	for(int i=0;i<numSteps;i++)
	{
		float energy=0;

		// Calculate residual of pattern at each searched point
		for(int idx=0;idx<PATTERNNUM;idx++)
		{
			float hitColor = getInterpolatedElement31(frame->dI,
										(float)(ptx+rotatetPattern[idx][0]),
										(float)(pty+rotatetPattern[idx][1]),
										wG0, hG0);

			if(!std::isfinite(hitColor)) {energy+=1e5; continue;}
			float residual = hitColor - (float)(hostToFrame_affine[0] * color[idx] + hostToFrame_affine[1]);
			float hw = fabs(residual) < globalSettings.setting_huberTH ? 1 : globalSettings.setting_huberTH / fabs(residual);
			energy += hw*residual*residual*(2-hw);
		}

		if(debugPrint && !globalSettings.no_Immature_debugMessage)
			printf("step %.1f %.1f (id %f): energy = %f!\n",
					ptx, pty, 0.0f, energy);

		errors[i] = energy;
		if(energy < bestEnergy)
		{
			bestU = ptx; bestV = pty; bestEnergy = energy; bestIdx = i;
		}

		ptx+=dx;
		pty+=dy;
	}

	// Find best score outside the min trace test radius.
	float secondBest=1e10;
	for(int i=0;i<numSteps;i++)
	{
		if((i < bestIdx-globalSettings.setting_minTraceTestRadius 
			|| i > bestIdx+globalSettings.setting_minTraceTestRadius) 
			&& errors[i] < secondBest)
			secondBest = errors[i];
	}
	float newQuality = secondBest / bestEnergy;
	if(newQuality < quality || numSteps > 10) quality = newQuality;


	// ============== do Gauss Newton optimization ===================
	// Use starting point from linear search
	// 2D search is done instead of 3D projected search
	float uBak=bestU, vBak=bestV, gnstepsize=1, stepBack=0;
	if(globalSettings.setting_trace_GNIterations>0) bestEnergy = 1e5;
	int gnStepsGood=0, gnStepsBad=0;

	for(int it=0;it<globalSettings.setting_trace_GNIterations;it++)
	{
		float H=1, b=0, energy=0;
		for(int idx=0;idx<PATTERNNUM;idx++)
		{
            float posU = (float)(bestU + rotatetPattern[idx][0]);
            float posV = (float)(bestV + rotatetPattern[idx][1]);

            if(posU < 0 || posV < 0 || posU >= wG0 - 1 || posV >= hG0 - 1) // Check if pattern is OOB
            {
                if(debugPrint && !globalSettings.no_Immature_debugMessage)
					printf("OOB uMax  %f %f - %f %f!\n", posU, posV, uMax, vMax);
                lastTraceUV = Vec2f(-1,-1);
                lastTracePixelInterval=0;

                return lastTraceStatus = ImmaturePointStatus::IPS_OOB;
            }

			Vec3f hitColor = getInterpolatedElement33(frame->dI, posU, posV, wG0, hG0);

			if(!std::isfinite((float)hitColor[0])) {energy+=1e5; continue;}
			// r = (I_j[u_j,v_j] - b_j) - a_j/a_i * (I_i[u_i,v_i] - b_i)
			// u_j = u_i + dx*d, v_j = v_i + dy*d
			float residual = hitColor[0] - (hostToFrame_affine[0] * color[idx] + hostToFrame_affine[1]);
			// dr/dd = dI_j/dx*dx + dI_j/dy*dy
			float dResdDist = dx*hitColor[1] + dy*hitColor[2]; // Calculate derivative for GN
			float hw = fabs(residual) < globalSettings.setting_huberTH ? 1 : globalSettings.setting_huberTH / fabs(residual);

			// Calculate GN variables
			H += hw*dResdDist*dResdDist; 
			b += hw*residual*dResdDist;

			// Calculate energy
			energy += weights[idx]*weights[idx]*hw*residual*residual*(2-hw);
		}

		if(energy > bestEnergy)
		{
			gnStepsBad++;

			// Do a smaller step from old point.
			stepBack*=0.5;
			bestU = uBak + stepBack*dx;
			bestV = vBak + stepBack*dy;
			if(debugPrint && !globalSettings.no_Immature_debugMessage)
				printf("GN BACK %d: E %f, H %f, b %f. id-step %f. UV %f %f -> %f %f.\n",
						it, energy, H, b, stepBack,
						uBak, vBak, bestU, bestV);
		}
		else
		{
			// Step towards to optimal point
			gnStepsGood++;

			// Solve GN
			float step = -gnstepsize*b/H;
			if(step < -0.5) step = -0.5;
			else if(step > 0.5) step = 0.5;

			if(!std::isfinite(step)) step=0;

			uBak=bestU;
			vBak=bestV;
			stepBack=step;

			bestU += step*dx;
			bestV += step*dy;
			bestEnergy = energy;

			if(debugPrint && !globalSettings.no_Immature_debugMessage)
				printf("GN step %d: E %f, H %f, b %f. id-step %f. UV %f %f -> %f %f.\n",
						it, energy, H, b, step,
						uBak, vBak, bestU, bestV);
		}

		if(fabsf(stepBack) < globalSettings.setting_trace_GNThreshold) break;
	}


	// ============== detect energy-based outlier. ===================
	if(!(bestEnergy < energyTH*globalSettings.setting_trace_extraSlackOnTH))
	{
		if(debugPrint && !globalSettings.no_Immature_debugMessage)
			printf("OUTLIER!\n");

		lastTracePixelInterval=0;
		lastTraceUV = Vec2f(-1,-1);
		if(lastTraceStatus == ImmaturePointStatus::IPS_OUTLIER)
			return lastTraceStatus = ImmaturePointStatus::IPS_OUTLIER_OUT;
		else
			return lastTraceStatus = ImmaturePointStatus::IPS_OUTLIER;
	}


	// ============== set new interval for idepth_min and idepth_max===================
	// depth can be calculated from either the x or y direction
	// u_j = (r_x+t_x*p_-1) / (r_z+t_z*p_-1)
	// p_-1 = (u_j*r_z - r_x) / (t_x - u_j*t_z)

	// v_j = (r_y+t_y*p_-1) / (r_z+t_z*p_-1)
	// p_-1 = (v_j*r_z - r_y) / (t_y - v_j*t_z)

	float dx_dy_ratio = (dx*dx/(dy*dy+dx*dx));
	idepth_min = (dx_dy_ratio)*
		(pr[2]*(bestU-errorInPixel*dx) - pr[0]) / 
		(hostToFrame_Kt[0] - hostToFrame_Kt[2]*(bestU-errorInPixel*dx)) +
		(1-dx_dy_ratio)*
		(pr[2]*(bestV-errorInPixel*dy) - pr[1]) / 
		(hostToFrame_Kt[1] - hostToFrame_Kt[2]*(bestV-errorInPixel*dy));

	idepth_max = (dx_dy_ratio)*
		(pr[2]*(bestU+errorInPixel*dx) - pr[0]) / 
		(hostToFrame_Kt[0] - hostToFrame_Kt[2]*(bestU+errorInPixel*dx)) +
		(1-dx_dy_ratio)*
		(pr[2]*(bestV+errorInPixel*dy) - pr[1]) / 
		(hostToFrame_Kt[1] - hostToFrame_Kt[2]*(bestV+errorInPixel*dy));

	if(idepth_min > idepth_max) std::swap<float>(idepth_min, idepth_max);


	if(!std::isfinite(idepth_min) || !std::isfinite(idepth_max) || (idepth_max<0))
	{
		if(!setting_debugout_runquiet && !globalSettings.no_Immature_debugMessage)
			printf("CAUGHT INF / NAN minmax depth!\n");

		lastTracePixelInterval=0;
		lastTraceUV = Vec2f(-1,-1);
		return lastTraceStatus = ImmaturePointStatus::IPS_OUTLIER;
	}

	lastTracePixelInterval=2*errorInPixel;
	lastTraceUV = Vec2f(bestU, bestV);

	return lastTraceStatus = ImmaturePointStatus::IPS_GOOD;
}


/**
 * @brief Calculate some residual values for the immature points
 * 
 * Will only calculate the depth Hessian
 * Calculates the energy residual and values that are helpful for point optimizations
 * 
 * @param HCalib 			Calibration matrices
 * @param outlierTHSlack 	Outlier threshold
 * @param tmpRes 			Struct storing residual information
 * @param Hdd 				
 * @param bd 				
 * @param idepth 			Input inverse depth
 * @return double 			EnergyLeft
 */
double ImmaturePoint::linearizeResidual(
		CalibHessian *  HCalib, const float outlierTHSlack,
		ImmaturePointTemporaryResidual* tmpRes,
		float &Hdd, float &bd,
		float idepth)
{
	if(tmpRes->state_state == ResState::OOB)
		{ tmpRes->state_NewState = ResState::OOB; return tmpRes->state_energy; }

	FrameFramePrecalc* precalc = &(host->targetPrecalc[tmpRes->target->idx]);

	float energyLeft=0;
	// Get image
	const Eigen::Vector3f* dIl = tmpRes->target->dI;

	// Set transforms
	const Mat33f &PRE_RTll = precalc->PRE_RTll;
	const Vec3f &PRE_tTll = precalc->PRE_tTll;
	//const float * const Il = tmpRes->target->I;

	// Get photometric values
	Vec2f affLL = precalc->PRE_aff_mode;

	for(int idx=0;idx<PATTERNNUM;idx++)
	{
		int dx = PATTERNP[idx][0];
		int dy = PATTERNP[idx][1];

		float drescale, u, v, new_idepth;
		float Ku, Kv;
		Vec3f KliP;

		// Get intensity of point at transformed location
		if(!projectPoint(this->u,this->v, idepth, dx, dy, HCalib,
				PRE_RTll, PRE_tTll, drescale, u, v, Ku, Kv, KliP, new_idepth, wG0, hG0))
			{tmpRes->state_NewState = ResState::OOB; return tmpRes->state_energy;}
		Vec3f hitColor = (getInterpolatedElement33(dIl, Ku, Kv, wG0, hG0));

		if(!std::isfinite((float)hitColor[0])) {tmpRes->state_NewState = ResState::OOB; return tmpRes->state_energy;}
		
		// Calculate residual
		float residual = hitColor[0] - (affLL[0] * color[idx] + affLL[1]);
		// Huber loss
		float hw = fabsf(residual) < globalSettings.setting_huberTH ? 1 : globalSettings.setting_huberTH / fabsf(residual);
		energyLeft += weights[idx]*weights[idx]*hw *residual*residual*(2-hw);

		// Depth derivative
		// r = (I_j[u_j,v_j] - b_j) - a_j/a_i * (I_i[u_i,v_i] - b_i)
		// ptp = RK[u,v,1]^T+tp
		// u_j = f_x*ptp_x/ptp_z+c_x, v_j = f_y*ptp_y/ptp_z+c_y

		// dr/dp = dI_j/du_j*du_j/dp + dI_j/dv_j*dv_j/dp
		// du_j/dp = f_x/ptp_z * (t_x+t_z*(ptp_x/ptp_z))
		// dv_j/dp = f_y/ptp_z * (t_y+t_z*(ptp_y/ptp_z))
		float dxInterp = hitColor[1]*HCalib->fxl();
		float dyInterp = hitColor[2]*HCalib->fyl();
		float d_idepth = derive_idepth(PRE_tTll, u, v, dx, dy, dxInterp, dyInterp, drescale);

		// Huber weight
		hw *= weights[idx]*weights[idx];

		Hdd += (hw*d_idepth)*d_idepth;
		bd += (hw*residual)*d_idepth;
	}


	if(energyLeft > energyTH*outlierTHSlack)
	{
		energyLeft = energyTH*outlierTHSlack;
		tmpRes->state_NewState = ResState::OUTLIER;
	}
	else
	{
		tmpRes->state_NewState = ResState::IN;
	}

	tmpRes->state_NewEnergy = energyLeft;
	return energyLeft;
}

}
