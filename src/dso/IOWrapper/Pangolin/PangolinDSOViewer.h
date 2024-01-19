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
#include <pangolin/pangolin.h>
#include "boost/thread.hpp"
#include "util/MinimalImage.h"
#include "IOWrapper/Output3DWrapper.h"
#include <map>
#include <deque>
#include "util/SettingsUtil.h"
#include "FollowCamMode.h"


namespace dmvio
{
class TransformDSOToIMU;
}

namespace dso
{

class FrameHessian;
class CalibHessian;
class FrameShell;


namespace IOWrap
{

class KeyFrameDisplay;

struct GraphConnection
{
	KeyFrameDisplay* from;
	KeyFrameDisplay* to;
	int fwdMarg, bwdMarg, fwdAct, bwdAct;
};

class InternalImageB3: public MinimalImage<Vec3b> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	pangolin::GlTexture FeatureFrameTexture;
	bool IsTextureGood;
	bool HaveNewImage;

	InternalImageB3(int w_, int h_) : MinimalImage(w_, h_)
	{
		IsTextureGood = false;
		HaveNewImage = false;
	}

	InternalImageB3(int w_, int h_, Vec3b* data_) : MinimalImage(w_, h_, data_)
	{
		IsTextureGood = false;
		HaveNewImage = false;
	}

	~InternalImageB3()
	{
	}
};

class PangolinDSOViewer : public Output3DWrapper
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    PangolinDSOViewer(int w, int h, bool startRunThread=true, std::shared_ptr<dmvio::SettingsUtil> settingsUtil =
            nullptr, std::shared_ptr<double> normalizeCamSize = nullptr);
	virtual ~PangolinDSOViewer();

	void run();
	void close();

	void addImageToDisplay(std::string name, MinimalImageB3* image);
	void clearAllImagesToDisplay();


	// ==================== Output3DWrapper Functionality ======================
    virtual void publishTransformDSOToIMU(const dmvio::TransformDSOToIMU& transformDSOToIMU) override;
    virtual void publishGraph(const std::map<uint64_t, Eigen::Vector2i, std::less<uint64_t>, Eigen::aligned_allocator<std::pair<const uint64_t, Eigen::Vector2i>>> &connectivity) override;
    virtual void publishKeyframes( std::vector<FrameHessian*> &frames, bool final, CalibHessian* HCalib) override;
    virtual void publishCamPose(FrameShell* frame, CalibHessian* HCalib) override;
    virtual void publishSystemStatus(dmvio::SystemStatus systemStatus) override;

    void addGTCamPose(const Sophus::SE3d& gtPose);

    virtual void pushLiveFrame(FrameHessian* image) override;
    virtual void pushDepthImage(MinimalImageB3* image) override;
    virtual bool needPushDepthImage() override;

    bool shouldQuit();

    virtual void join() override;

    virtual void reset() override;
private:

	bool needReset;
	void reset_internal();
	void drawConstraints();

	boost::thread runThread;
	bool running;
	bool shouldQuitVar{false};
	int w,h;

	std::vector<int> active_frame_IDs;

	// images rendering
	boost::mutex openImagesMutex;
	std::unique_ptr<InternalImageB3> internalVideoImg;
	std::unique_ptr<InternalImageB3> internalKFImg;
	std::unique_ptr<InternalImageB3> internalResImg;

	pangolin::View *d_kfDepth;
	pangolin::View *d_video;
	pangolin::View *d_residual;

    CalibHessian *HCalib;

	// 3D model rendering
	boost::mutex model3DMutex;
	KeyFrameDisplay* currentCam, *currentGTCam;
	std::vector<KeyFrameDisplay*> keyframes;
	std::vector<Vec3f,Eigen::aligned_allocator<Vec3f>> allFramePoses;
	std::map<int, KeyFrameDisplay*> keyframesByKFID;
	std::vector<GraphConnection,Eigen::aligned_allocator<GraphConnection>> connections;



	// render settings
	bool setting_render_display3D;
	bool setting_render_displayDepth;
	bool setting_render_displayResidual;
	bool setting_render_displayVideo;

	bool settings_showKFCameras;
	bool settings_showCurrentCamera;
	bool settings_showTrajectory;
	bool settings_showFullTrajectory;
	bool settings_showActiveConstraints;
	bool settings_showAllConstraints;

	float settings_scaledVarTH;
	float settings_absVarTH;
	int settings_pointCloudMode;
	bool settings_showDrawPC;
	bool settings_showOnlyActive;
	float settings_minRelBS;
	int settings_sparsity;


	// timings
	struct timeval last_track;
	struct timeval last_map;


	std::deque<float> lastNTrackingMs;
	std::deque<float> lastNMappingMs;

	// GT cam poses
	SE3 gtCamPoseMetric; // in metric frame, but will be displayed in DSO frame using.
	SE3 firstCamPoseDSO; // DSO pose corresponding to the first groundtruth pose.
	SE3 firstGTCamPoseMetric;
	bool gtCamPoseSet = false;
	std::unique_ptr<dmvio::TransformDSOToIMU> transformDSOToIMU;
    dmvio::SystemStatus systemStatus;
	void updateDisplayedCamPose();
    std::shared_ptr<double> normalizeCamSize;

    FollowCamMode followCam;

	std::shared_ptr<dmvio::SettingsUtil> settingsUtil;
};



}



}
