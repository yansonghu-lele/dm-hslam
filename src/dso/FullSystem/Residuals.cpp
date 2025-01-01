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

#include "FullSystem/HessianBlocks.h"



namespace dso
{
int PointFrameResidual::instanceCounter = 0;
long runningResID=0;


PointFrameResidual::PointFrameResidual(){assert(false); instanceCounter++; globalSettings = nullptr;}

PointFrameResidual::~PointFrameResidual(){assert(efResidual==0); instanceCounter--; delete J;}

/**
 * @brief Construct a new Point Frame Residual:: Point Frame Residual object
 * 
 * @param globalCalib_ 
 * @param point_ 
 * @param host_ 
 * @param target_ 
 * @param globalSettings_ 
 */
PointFrameResidual::PointFrameResidual(Global_Calib* globalCalib_, PointHessian* point_, FrameHessian* host_, FrameHessian* target_,  GlobalSettings* globalSettings_) :
	point(point_),
	host(host_),
	target(target_),
	globalCalib(globalCalib_),
	globalSettings(globalSettings_)
{
	wG0 = globalCalib->wG[0];
	hG0 = globalCalib->hG[0];
	efResidual=0;
	instanceCounter++;
	resetOOB();
	J = new RawResidualJacobian();
	assert(((long)J)%16==0);
}

/**
 * @brief Linearize point
 * 
 * Calculates the residual and Jacobians of the point
 * Note that the linearization is only done at one pyramid level
 * 
 * @param HCalib 
 * @return double 
 */
double PointFrameResidual::linearize(CalibHessian* HCalib)
{
	state_NewEnergyWithOutlier=-1;

	if(state_state == ResState::OOB)
		{ state_NewState = ResState::OOB; return state_energy; }


	FrameFramePrecalc* precalc = &(host->targetPrecalc[target->idx]);
	float energyLeft=0;
	// Get image
	const Eigen::Vector3f* dIl = target->dI;
	//const float* const Il = target->I;
	// Get K matrix and transforms from the pre-calculation
	// K matrix and Transformation matrix are multiplied here for effciency
	const Mat33f &PRE_KRKiTll = precalc->PRE_KRKiTll; 	// K * rotationMatrix * K^-1
	const Vec3f &PRE_KtTll = precalc->PRE_KtTll; 		// K * translationMatrix
	const Mat33f &PRE_RTll_0 = precalc->PRE_RTll_0; 	// rotationMatrix_0
	const Vec3f &PRE_tTll_0 = precalc->PRE_tTll_0; 		// translationMatrix_0
	// Get info about point
	const float * const color = point->color;
	const float * const weights = point->weights;

	// Photogrammetric values
	Vec2f affLL = precalc->PRE_aff_mode;
	float b0 = precalc->PRE_b0_mode;

	Vec6f d_xi_x, d_xi_y;
	Vec4f d_C_x, d_C_y;
	float d_d_x, d_d_y;
	
	// Calculate the geometric residual for the point
	{
		float drescale, u, v, new_idepth;
		float Ku, Kv;
		Vec3f KliP;

		// Project point from 3D to 2D
		// We only project the middle point because all the pixels
		// in a pattern share the same geometric residual

		// The calibration and transformation matrices are applied by the projection
		if(!projectPoint(point->u, point->v, point->idepth_zero_scaled, 0, 0, HCalib,
				PRE_RTll_0, PRE_tTll_0, drescale, u, v, Ku, Kv, KliP, new_idepth, wG0, hG0))
			{ state_NewState = ResState::OOB; return state_energy; }

		// Get the new coordinates
		centerProjectedTo = Vec3f(Ku, Kv, new_idepth);

		// Calculate the values of d_point'/d_J_geometric
		// These are evaluated using First Estimate Jacobians
		// therefore they are evaluated at x=0 instead of the point position

		// Partial differential of d_Ku/d_idepth
		d_d_x = drescale * (PRE_tTll_0[0]-PRE_tTll_0[2]*u)*SCALE_IDEPTH*HCalib->fxl();
		// Partial differential of d_Kv/d_idepth
		d_d_y = drescale * (PRE_tTll_0[1]-PRE_tTll_0[2]*v)*SCALE_IDEPTH*HCalib->fyl();


		// Partial differential of calib matrix
		d_C_x[2] = drescale*(PRE_RTll_0(2,0)*u-PRE_RTll_0(0,0)); //d_Ku/d_c_x
		d_C_x[3] = HCalib->fxl() * drescale*(PRE_RTll_0(2,1)*u-PRE_RTll_0(0,1)) * HCalib->fyli(); //d_Ku/d_c_y
		d_C_x[0] = KliP[0]*d_C_x[2]; //d_Ku/d_f_x
		d_C_x[1] = KliP[1]*d_C_x[3]; //d_Ku/d_f_y

		d_C_y[2] = HCalib->fyl() * drescale*(PRE_RTll_0(2,0)*v-PRE_RTll_0(1,0)) * HCalib->fxli(); //d_Kv/d_c_x
		d_C_y[3] = drescale*(PRE_RTll_0(2,1)*v-PRE_RTll_0(1,1)); //d_Kv/d_c_y
		d_C_y[0] = KliP[0]*d_C_y[2]; //d_Kv/d_f_x
		d_C_y[1] = KliP[1]*d_C_y[3]; //d_Kv/d_f_y

		d_C_x[0] = (d_C_x[0]+u)*SCALE_F;
		d_C_x[1] *= SCALE_F;
		d_C_x[2] = (d_C_x[2]+1)*SCALE_C;
		d_C_x[3] *= SCALE_C;

		d_C_y[0] *= SCALE_F;
		d_C_y[1] = (d_C_y[1]+v)*SCALE_F;
		d_C_y[2] *= SCALE_C;
		d_C_y[3] = (d_C_y[3]+1)*SCALE_C;


		// Partial differential of transformations
		// Translation
		d_xi_x[0] = new_idepth*HCalib->fxl(); //d_Ku/d_d_x
		d_xi_x[1] = 0; //d_Ku/d_d_y
		d_xi_x[2] = -new_idepth*u*HCalib->fxl(); //d_Ku/d_d_z
		// Rotation (three lie algebra variables)
		d_xi_x[3] = -u*v*HCalib->fxl();
		d_xi_x[4] = (1+u*u)*HCalib->fxl();
		d_xi_x[5] = -v*HCalib->fxl();

		// Translation
		d_xi_y[0] = 0; //d_Ku/d_d_x
		d_xi_y[1] = new_idepth*HCalib->fyl(); //d_Ku/d_d_y
		d_xi_y[2] = -new_idepth*v*HCalib->fyl(); //d_Ku/d_d_z
		// Rotation (three lie algebra variables)
		d_xi_y[3] = -(1+v*v)*HCalib->fyl();
		d_xi_y[4] = u*v*HCalib->fyl();
		d_xi_y[5] = u*HCalib->fyl();
	}

	{
		// Set the values of the Jacobian matrix
		J->Jpdxi[0] = d_xi_x;
		J->Jpdxi[1] = d_xi_y;

		J->Jpdc[0] = d_C_x;
		J->Jpdc[1] = d_C_y;

		J->Jpdd[0] = d_d_x;
		J->Jpdd[1] = d_d_y;
	}


	// Calculate the projective and image residual for the pattern
	float JIdxJIdx_00=0, JIdxJIdx_11=0, JIdxJIdx_10=0;
	float JabJIdx_00=0, JabJIdx_01=0, JabJIdx_10=0, JabJIdx_11=0;
	float JabJab_00=0, JabJab_01=0, JabJab_11=0;

	float wJI2_sum = 0;

	for(int idx=0;idx<PATTERNNUM;idx++)
	{
		float Ku, Kv;
		if(!projectPoint(point->u+PATTERNP[idx][0], point->v+PATTERNP[idx][1], point->idepth_scaled, PRE_KRKiTll, PRE_KtTll, Ku, Kv, wG0, hG0))
			{ state_NewState = ResState::OOB; return state_energy; }

		projectedTo[idx][0] = Ku;
		projectedTo[idx][1] = Kv;
        Vec3f hitColor = (getInterpolatedElement33(dIl, Ku, Kv, wG0, hG0));

		// Calculate residual
        float residual = hitColor[0] - (float)(affLL[0] * color[idx] + affLL[1]);


		float drdA = (color[idx]-b0);
		if(!std::isfinite((float)hitColor[0]))
			{ state_NewState = ResState::OOB; return state_energy; }

		// Adjust residual
		// Apply gradient-dependent weighting
		float w = sqrtf(globalSettings->setting_outlierTHSumComponent / (globalSettings->setting_outlierTHSumComponent + hitColor.tail<2>().squaredNorm()));
        w = 0.5f*(w + weights[idx]);
		// Huber weight
		float hw = fabsf(residual) < globalSettings->setting_huberTH ? 1 : globalSettings->setting_huberTH / fabsf(residual);
		// Use Huber Norm
		// hw*(2-hw)*residual*residual results in Huber loss for the residual
		energyLeft += w*w*hw*residual*residual*(2-hw);

		{
			// Calculate Jacobian of image variables
			if(hw < 1) hw = sqrtf(hw);
			hw = hw*w;

			// Pixel derivatives
			hitColor[1]*=hw;
			hitColor[2]*=hw;

			J->resF[idx] = residual*hw;

			// Calculate Jacobian of pixel intensities
			J->JIdx[0][idx] = hitColor[1];
			J->JIdx[1][idx] = hitColor[2];

			// Calculate Jacobian of photometric variables
			J->JabF[0][idx] = drdA*hw;
			J->JabF[1][idx] = hw;

			// Pre-calculate values for effciency
			JIdxJIdx_00+=hitColor[1]*hitColor[1];
			JIdxJIdx_11+=hitColor[2]*hitColor[2];
			JIdxJIdx_10+=hitColor[1]*hitColor[2];

			JabJIdx_00+= drdA*hw * hitColor[1];
			JabJIdx_01+= drdA*hw * hitColor[2];
			JabJIdx_10+= hw * hitColor[1];
			JabJIdx_11+= hw * hitColor[2];

			JabJab_00+= drdA*drdA*hw*hw;
			JabJab_01+= drdA*hw*hw;
			JabJab_11+= hw*hw;

			wJI2_sum += hw*hw*(hitColor[1]*hitColor[1]+hitColor[2]*hitColor[2]);


			if(globalSettings->setting_affineOptModeA < 0) J->JabF[0][idx]=0;
			if(globalSettings->setting_affineOptModeB < 0) J->JabF[1][idx]=0;
		}
	}

	// Set pre-calculated values
	J->JIdx2(0,0) = JIdxJIdx_00;
	J->JIdx2(0,1) = JIdxJIdx_10;
	J->JIdx2(1,0) = JIdxJIdx_10;
	J->JIdx2(1,1) = JIdxJIdx_11;
	J->JabJIdx(0,0) = JabJIdx_00;
	J->JabJIdx(0,1) = JabJIdx_01;
	J->JabJIdx(1,0) = JabJIdx_10;
	J->JabJIdx(1,1) = JabJIdx_11;
	J->Jab2(0,0) = JabJab_00;
	J->Jab2(0,1) = JabJab_01;
	J->Jab2(1,0) = JabJab_01;
	J->Jab2(1,1) = JabJab_11;


	state_NewEnergyWithOutlier = energyLeft;
	// Check if point is outlier
	if(energyLeft > std::max<float>(host->frameEnergyTH, target->frameEnergyTH) || wJI2_sum < 2)
	{
		energyLeft = std::max<float>(host->frameEnergyTH, target->frameEnergyTH);
		state_NewState = ResState::OUTLIER;
	}
	else
	{
		state_NewState = ResState::IN;
	}

	state_NewEnergy = energyLeft;
	return energyLeft;
}

#ifdef GRAPHICAL_DEBUG
void PointFrameResidual::debugPlot()
{
	if(state_state==ResState::OOB) return;
	Vec3b cT = Vec3b(0,0,0);

	if(freeDebugParam5==0)
	{
		float rT = 20*sqrt(state_energy/9);
		if(rT<0) rT=0; if(rT>255)rT=255;
		cT = Vec3b(0,255-rT,rT);
	}
	else
	{
		if(state_state == ResState::IN) cT = Vec3b(255,0,0);
		else if(state_state == ResState::OOB) cT = Vec3b(255,255,0);
		else if(state_state == ResState::OUTLIER) cT = Vec3b(0,0,255);
		else cT = Vec3b(255,255,255);
	}

	for(int i=0;i<PATTERNNUM;i++)
	{
		if((projectedTo[i][0] > PIXEL_BORDER && projectedTo[i][1] > PIXEL_BORDER && projectedTo[i][0] < wG0-PIXEL_BORDER-1 && projectedTo[i][1] < hG0-PIXEL_BORDER-1 ))
			target->debugImage->setPixel1((float)projectedTo[i][0], (float)projectedTo[i][1],cT);
	}
}
#endif

/**
 * @brief Set the state and energy to the new values
 * 
 * @param copyJacobians 
 */
void PointFrameResidual::applyRes(bool copyJacobians)
{
	if(copyJacobians)
	{
		if(state_state == ResState::OOB)
		{
			assert(!efResidual->isActiveAndIsGoodNEW);
			return;	// can never go back from OOB
		}
		if(state_NewState == ResState::IN)// && )
		{
			efResidual->isActiveAndIsGoodNEW=true;
			efResidual->takeDataF();
		}
		else
		{
			efResidual->isActiveAndIsGoodNEW=false;
		}
	}

	setState(state_NewState);
	state_energy = state_NewEnergy;
}
}
