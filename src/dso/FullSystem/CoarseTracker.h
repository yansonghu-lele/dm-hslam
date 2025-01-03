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



#pragma once

#include "util/NumType.h"
#include "vector"
#include <math.h>
#include "util/settings.h"
#include "OptimizationBackend/MatrixAccumulators.h"
#include "IOWrapper/Output3DWrapper.h"

#include "IMU/IMUIntegration.hpp"



namespace dso
{

struct CalibHessian;
struct FrameHessian;
struct PointFrameResidual;

class CoarseTracker {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	CoarseTracker(Global_Calib& globalCalib_, std::shared_ptr<dmvio::IMUIntegration> imuIntegration_, GlobalSettings& globalSettings_);
	~CoarseTracker();

	bool trackNewestCoarse(
			FrameHessian* newFrameHessian,
			SE3 &lastToNew_out, AffLight &aff_g2l_out,
			int coarsestLvl, Vec5 minResForAbort,
			IOWrap::Output3DWrapper* wrap=0);

	void setCoarseTrackingRef(
			std::vector<FrameHessian*> frameHessians);

	void makeK(
			CalibHessian* HCalib);

	bool debugPrint;

	int wG0;
	int hG0;
	Mat33f K[PYR_LEVELS];
	Mat33f Ki[PYR_LEVELS];
	float fx[PYR_LEVELS];
	float fy[PYR_LEVELS];
	float fxi[PYR_LEVELS];
	float fyi[PYR_LEVELS];
	float cx[PYR_LEVELS];
	float cy[PYR_LEVELS];
	float cxi[PYR_LEVELS];
	float cyi[PYR_LEVELS];
	int w[PYR_LEVELS];
	int h[PYR_LEVELS];

    void debugPlotIDepthMap(float* minID, float* maxID, std::vector<IOWrap::Output3DWrapper*> &wraps) const;
    void debugPlotIDepthMapFloat(std::vector<IOWrap::Output3DWrapper*> &wraps);

	FrameHessian* lastRef;
	AffLight lastRef_aff_g2l;
	FrameHessian* newFrame;
	int refFrameID;

	// act as pure ouptut
	Vec5 lastResidualsStats;
	Vec3 lastFlowIndicators;
	double firstCoarseRMSE;


private:
	GlobalSettings& globalSettings;
	Global_Calib& globalCalib;

	void makeCoarseDepthL0(std::vector<FrameHessian*> frameHessians);
	float* idepth[PYR_LEVELS];
	float* weightSums[PYR_LEVELS];
	float* weightSums_bak[PYR_LEVELS];


	Vec6 calcRes(int lvl, const SE3 &refToNew, AffLight aff_g2l, float cutoffTH);
	void calcGSSSE(int lvl, Mat88 &H_out, Vec8 &b_out, const SE3 &refToNew, AffLight aff_g2l);

	// pc buffers
	float* pc_u[PYR_LEVELS];
	float* pc_v[PYR_LEVELS];
	float* pc_idepth[PYR_LEVELS];
	float* pc_color[PYR_LEVELS];
	int pc_n[PYR_LEVELS];

	// warped buffers
	float* buf_warped_idepth;
	float* buf_warped_u;
	float* buf_warped_v;
	float* buf_warped_dx;
	float* buf_warped_dy;
	float* buf_warped_residual;
	float* buf_warped_weight;
	float* buf_warped_refColor;
	int buf_warped_n;

    std::vector<float*> ptrToDelete;
	Accumulator9 acc9;

    std::shared_ptr<dmvio::IMUIntegration> imuIntegration;
};


class CoarseDistanceMap {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	CoarseDistanceMap(Global_Calib& globalCalib_, int pyrLevelsUsed_);
	~CoarseDistanceMap();

	void makeDistanceMap(
			std::vector<FrameHessian*> frameHessians,
			FrameHessian* frame);

	void makeInlierVotes(
			std::vector<FrameHessian*> frameHessians);

	void makeK( CalibHessian* HCalib);


	float* fwdWarpedIDDistFinal;

	int wG0;
	int hG0;
	Mat33f K[PYR_LEVELS];
	Mat33f Ki[PYR_LEVELS];
	float fx[PYR_LEVELS];
	float fy[PYR_LEVELS];
	float fxi[PYR_LEVELS];
	float fyi[PYR_LEVELS];
	float cx[PYR_LEVELS];
	float cy[PYR_LEVELS];
	float cxi[PYR_LEVELS];
	float cyi[PYR_LEVELS];
	int w[PYR_LEVELS];
	int h[PYR_LEVELS];

	void addIntoDistFinal(int u, int v);


private:
	Global_Calib& globalCalib;

	int pyrLevelsUsed;

	PointFrameResidual** coarseProjectionGrid;
	int* coarseProjectionGridNum;
	Eigen::Vector2i* bfsList1;
	Eigen::Vector2i* bfsList2;

	void growDistBFS(int bfsNum);
};

}

