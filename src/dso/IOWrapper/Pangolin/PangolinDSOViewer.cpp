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



#include "PangolinDSOViewer.h"
#include "KeyFrameDisplay.h"

#include "util/settings.h"
#include "util/globalCalib.h"
#include "FullSystem/HessianBlocks.h"
#include "FullSystem/FullSystem.h"
#include "FullSystem/ImmaturePoint.h"
#include "OverwritePangolin.h"

namespace dso
{
namespace IOWrap
{



PangolinDSOViewer::PangolinDSOViewer(int w_, int h_, GlobalSettings& globalSettings_, bool startRunThread, std::shared_ptr<dmvio::SettingsUtil>
        settingsUtilPassed, std::shared_ptr<double> normalizeCamSize)
        : globalSettings(globalSettings_), HCalib(0), settingsUtil(std::move(settingsUtilPassed)), normalizeCamSize(normalizeCamSize)
{
	this->w = w_;
	this->h = h_;
	running=true;


	{
		boost::unique_lock<boost::mutex> lk(openImagesMutex);
		internalVideoImg = std::unique_ptr<InternalImageB3>(new InternalImageB3(w,h));
		internalKFImg = std::unique_ptr<InternalImageB3>(new InternalImageB3(w,h));
		internalResImg = std::unique_ptr<InternalImageB3>(new InternalImageB3(w,h));

		internalVideoImg->HaveNewImage = true;
		internalKFImg->HaveNewImage = true;
		internalResImg->HaveNewImage = true;

		internalVideoImg->setBlack();
		internalKFImg->setBlack();
		internalResImg->setBlack();
	}


	{
		currentCam = new KeyFrameDisplay(w, h);
        currentGTCam = new KeyFrameDisplay(w, h);
	}

	needReset = false;
	active_frame_IDs.resize(globalSettings.setting_maxFrames+1, -1);


    if(startRunThread)
        runThread = boost::thread(&PangolinDSOViewer::run, this);

}


PangolinDSOViewer::~PangolinDSOViewer()
{
	close();
	if(runThread.joinable())
        runThread.join();
}


void PangolinDSOViewer::run()
{
	printf("START PANGOLIN!\n");

	pangolin::CreateWindowAndBind("Main",2*w,2*h);
	const int UI_WIDTH = 180;

	glEnable(GL_DEPTH_TEST);

	// 3D visualization
	pangolin::OpenGlRenderState Visualization3D_camera(
		pangolin::ProjectionMatrix(w,h,400,400,w/2,h/2,0.1,1000),
		pangolin::ModelViewLookAt(-0,-5,-10, 0,0,0, pangolin::AxisNegY)
		);

	pangolin::View& Visualization3D_display = pangolin::CreateDisplay()
		.SetBounds(0.0, 1.0, pangolin::Attach::Pix(UI_WIDTH), 1.0, -w/(float)h)
		.SetHandler(new pangolin::Handler3D(Visualization3D_camera));


	// 3 images
	d_kfDepth = &pangolin::Display("imgKFDepth").SetAspect(w/(float)h);
	d_video = &pangolin::Display("imgVideo").SetAspect(w/(float)h);
	d_residual = &pangolin::Display("imgResidual").SetAspect(w/(float)h);

	internalVideoImg->FeatureFrameTexture.Reinitialise(w,h,GL_RGB,false,0,GL_RGB,GL_UNSIGNED_BYTE);
	internalKFImg->FeatureFrameTexture.Reinitialise(w,h,GL_RGB,false,0,GL_RGB,GL_UNSIGNED_BYTE);
	internalResImg->FeatureFrameTexture.Reinitialise(w,h,GL_RGB,false,0,GL_RGB,GL_UNSIGNED_BYTE);

    pangolin::CreateDisplay()
		  .SetBounds(0.7, 1.0, pangolin::Attach::Pix(UI_WIDTH), 1.0)
		  .SetLayout(pangolin::LayoutEqual)
		  .SetLock(pangolin::LockLeft, pangolin::LockTop)
		  .AddDisplay(*d_kfDepth)
		  .AddDisplay(*d_video)
		  .AddDisplay(*d_residual)
		  .SetHandler(new pangolin::HandlerResize());

	// parameter reconfigure gui
	pangolin::CreatePanel("ui").SetBounds(0.0, 1.0, 0.0, pangolin::Attach::Pix(UI_WIDTH));

	pangolin::Var<int> settings_pointCloudMode("ui.PC_mode",2,0,3,false);
	pangolin::Var<bool> settings_showDrawPC("ui.PC_OnOff",true,true);
	pangolin::Var<bool> settings_showOnlyActive("ui.PC_OnlyActive",true,true);

	pangolin::Var<bool> settings_showKFCameras("ui.KFCam",false,true);
	pangolin::Var<bool> settings_showCurrentCamera("ui.CurrCam",true,true);
	pangolin::Var<bool> settings_showTrajectory("ui.KFTrajectory",false,true);
	pangolin::Var<bool> settings_showFullTrajectory("ui.FullTrajectory",true,true);
	pangolin::Var<bool> settings_showActiveConstraints("ui.ActiveConst",true,true);
	pangolin::Var<bool> settings_showAllConstraints("ui.AllConst",false,true);


	pangolin::Var<bool> settings_show3D("ui.show3D",true,true);
	pangolin::Var<bool> settings_showLiveDepth("ui.showDepth",true,true);
	pangolin::Var<bool> settings_showLiveVideo("ui.showVideo",true,true);
    pangolin::Var<bool> settings_showLiveResidual("ui.showResidual",false,true);

#ifdef GRAPHICAL_DEBUG
	pangolin::Var<bool> settings_showFramesWindow("ui.BD-FramesWindow",false,true);
	pangolin::Var<bool> settings_showFullTracking("ui.BD-FullTracking",false,true);
	pangolin::Var<bool> settings_showCoarseTracking("ui.BD-CoarseTracking",false,true);
	pangolin::Var<bool> settings_showImmatureTracking("ui.BD-ImmatureTracking",false,true);
	pangolin::Var<bool> settings_showFramebyFrame("ui.BD-FramebyFrame",false,true);
#endif

	pangolin::Var<int> settings_sparsity("ui.sparsity",1,1,20,false);
	pangolin::Var<double> settings_scaledVarTH("ui.relVarTH",0.001,1e-10,1e10, true);
	pangolin::Var<double> settings_absVarTH("ui.absVarTH",0.001,1e-10,1e10, true);
	pangolin::Var<double> settings_minRelBS("ui.minRelativeBS",0.1,0,1, false);


	pangolin::Var<bool> settings_resetButton("ui.Reset",false,false);
#ifdef GRAPHICAL_DEBUG
	pangolin::Var<bool> settings_pauseButton("ui.Pause",false,false);
#endif


	pangolin::Var<int> settings_nPts("ui.activePoints",globalSettings.setting_desiredPointDensity, 50,5000, false);
	pangolin::Var<int> settings_nCandidates("ui.pointCandidates",globalSettings.setting_desiredImmatureDensity, 50,5000, false);
	pangolin::Var<int> settings_nMaxFrames("ui.maxFrames",globalSettings.setting_maxFrames, 4,10, false);
	pangolin::Var<double> settings_kfFrequency("ui.kfFrequency",globalSettings.setting_kfGlobalWeight,0.1,3, false);
	pangolin::Var<double> settings_gradHistAdd("ui.minGradAdd",globalSettings.setting_minGradHistAdd,0,1, false);
	pangolin::Var<double> settings_gradHistCut("ui.minGradCut",globalSettings.setting_minGradHistCut,0,1, false);

	pangolin::Var<double> settings_trackFps("ui.Track fps",0,0,0,false);
	pangolin::Var<double> settings_mapFps("ui.KF fps",0,0,0,false);

    pangolin::Var<double> settings_Scale("ui.Scale", 1, 0, 0, false);
    pangolin::Var<std::string> setting_SystemStatus("ui.Status", "");



    if(settingsUtil)
    {
        settingsUtil->createPangolinSettings();
    }

    followCam.createPangolinSettings();


    // Default hooks for exiting (Esc) and fullscreen (tab).
	while(running)
	{
		// Clear entire screen
		glClearColor(0.36, 0.36, 0.55, 0.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if(setting_render_display3D)
		{
            double sizeFactor = 1.0;

            // I need model3DMutex lock before I get currentCam->camToWorld.
            boost::unique_lock<boost::mutex> lk3d(model3DMutex);

            auto updatedVisualization3DCam = followCam.updateVisualizationCam(currentCam->camToWorld,
                                                                         Visualization3D_camera);

			// Activate efficiently by object
            Visualization3D_display.Activate(*updatedVisualization3DCam);

            // Normalize cam size with scale.
            if(transformDSOToIMU && normalizeCamSize && *normalizeCamSize > 0)
            {
                sizeFactor = *normalizeCamSize / transformDSOToIMU->getScale();
            }

			int refreshed=0;

			if(!settings_showOnlyActive){
				for (std::vector<KeyFrameDisplay*>::reverse_iterator i = keyframes.rbegin(); 
					i != keyframes.rend(); ++i ) {
						KeyFrameDisplay* fh = *i;

						if(this->settings_showKFCameras) {
							float blue[3] = {0.32,0.7,0.8};
							fh->drawCam(1,blue,0.1 * sizeFactor);
						}

						refreshed += (int)(fh->refreshPC(refreshed < 10, this->settings_scaledVarTH, this->settings_absVarTH,
								this->settings_pointCloudMode, this->settings_minRelBS, this->settings_sparsity));
						
						if (settings_showDrawPC) fh->drawPC(1);
				} 
			} else {
				for(int id : active_frame_IDs)
				{
					if (id != -1){
						if(keyframesByKFID.find(id) == keyframesByKFID.end()) continue;

						KeyFrameDisplay* fh = keyframesByKFID[id];

						if(this->settings_showKFCameras) {
							float blue[3] = {0.32,0.7,0.8};
							fh->drawCam(1,blue,0.1 * sizeFactor);
						}

						refreshed += (int)(fh->refreshPC(refreshed < 10, this->settings_scaledVarTH, this->settings_absVarTH,
								this->settings_pointCloudMode, this->settings_minRelBS, this->settings_sparsity));
						
						if (settings_showDrawPC) fh->drawPC(1);
					}
				}
			}

			if(this->settings_showCurrentCamera) currentCam->drawCam(2,0,0.2 * sizeFactor);

            if(gtCamPoseSet)
            {
				float yellow[3] = {0.94,0.82,0.0};
                currentGTCam->drawCam(2, yellow, 0.2 * sizeFactor);
            }

            drawConstraints();
			lk3d.unlock();
		}



		openImagesMutex.lock();

		if(internalVideoImg->HaveNewImage) 	internalVideoImg->FeatureFrameTexture.Upload(internalVideoImg->data,GL_BGR,GL_UNSIGNED_BYTE);
		if(internalKFImg->HaveNewImage) 		internalKFImg->FeatureFrameTexture.Upload(internalKFImg->data,GL_BGR,GL_UNSIGNED_BYTE);
		if(internalResImg->HaveNewImage) 		internalResImg->FeatureFrameTexture.Upload(internalResImg->data,GL_BGR,GL_UNSIGNED_BYTE);
		
		internalVideoImg->HaveNewImage = false;
		internalKFImg->HaveNewImage = false;
		internalResImg->HaveNewImage = false;

		openImagesMutex.unlock();




		// update fps counters
		{
			openImagesMutex.lock();
			float sd = std::accumulate(lastNMappingMs.begin(), lastNMappingMs.end(), 0);
			settings_mapFps=lastNMappingMs.size()*1000.0f / sd;
			openImagesMutex.unlock();
		}
		{
			model3DMutex.lock();
			float sd = std::accumulate(lastNTrackingMs.begin(), lastNTrackingMs.end(), 0);
			settings_trackFps = lastNTrackingMs.size()*1000.0f / sd;
			model3DMutex.unlock();
		}

        // Update scale and status text.
        {
            boost::unique_lock<boost::mutex> lk(model3DMutex);
            if(transformDSOToIMU)
            {
                settings_Scale = transformDSOToIMU->getScale();
            }
            switch(systemStatus)
            {
                case dmvio::VISUAL_INIT:
                    setting_SystemStatus = "Visual-init";
                    break;
                case dmvio::VISUAL_ONLY:
                    setting_SystemStatus = "Visual-only";
                    break;
                case dmvio::VISUAL_INERTIAL:
                    setting_SystemStatus = "VIO";
                    break;
            }
        }

		if(setting_render_displayVideo)
		{
			d_video->Activate();
			glColor4f(1.0f,1.0f,1.0f,1.0f);
			internalVideoImg->FeatureFrameTexture.RenderToViewportFlipY();
		}

		if(setting_render_displayDepth)
		{
			d_kfDepth->Activate();
			glColor4f(1.0f,1.0f,1.0f,1.0f);
			internalKFImg->FeatureFrameTexture.RenderToViewportFlipY();
		}

		if(setting_render_displayResidual)
		{
			d_residual->Activate();
			glColor4f(1.0f,1.0f,1.0f,1.0f);
			internalResImg->FeatureFrameTexture.RenderToViewportFlipY();
		}


	    // update parameters
	    this->settings_pointCloudMode = settings_pointCloudMode.Get();
		this->settings_showDrawPC = settings_showDrawPC.Get();
		this->settings_showOnlyActive = settings_showOnlyActive.Get();

	    this->settings_showActiveConstraints = settings_showActiveConstraints.Get();
	    this->settings_showAllConstraints = settings_showAllConstraints.Get();
	    this->settings_showCurrentCamera = settings_showCurrentCamera.Get();
	    this->settings_showKFCameras = settings_showKFCameras.Get();
	    this->settings_showTrajectory = settings_showTrajectory.Get();
	    this->settings_showFullTrajectory = settings_showFullTrajectory.Get();

		setting_render_display3D = settings_show3D.Get();
		setting_render_displayDepth = settings_showLiveDepth.Get();
		setting_render_displayVideo =  settings_showLiveVideo.Get();
		setting_render_displayResidual = settings_showLiveResidual.Get();

#ifdef GRAPHICAL_DEBUG
		globalSettings.setting_render_renderWindowFrames = settings_showFramesWindow.Get();
		globalSettings.setting_render_plotTrackingFull = settings_showFullTracking.Get();
		globalSettings.setting_render_displayCoarseTrackingFull = settings_showCoarseTracking.Get();
		globalSettings.setting_render_displayImmatureTracking = settings_showImmatureTracking.Get();
		globalSettings.setting_goStepByStep = settings_showFramebyFrame.Get();
#endif

	    this->settings_absVarTH = settings_absVarTH.Get();
	    this->settings_scaledVarTH = settings_scaledVarTH.Get();
	    this->settings_minRelBS = settings_minRelBS.Get();
	    this->settings_sparsity = settings_sparsity.Get();

	    globalSettings.setting_desiredPointDensity = settings_nPts.Get();
	    globalSettings.setting_desiredImmatureDensity = settings_nCandidates.Get();
	    globalSettings.setting_maxFrames = settings_nMaxFrames.Get();
	    globalSettings.setting_kfGlobalWeight = settings_kfFrequency.Get();
	    globalSettings.setting_minGradHistAdd = settings_gradHistAdd.Get();
		globalSettings.setting_minGradHistCut = settings_gradHistCut.Get();

        if(settingsUtil)
        {
            settingsUtil->updatePangolinSettings();
        }


        if(settings_resetButton.Get())
	    {
	    	printf("RESET!\n");
	    	settings_resetButton.Reset();
	    	setting_fullResetRequested = true;
	    }

#ifdef GRAPHICAL_DEBUG
		if(settings_pauseButton.Get())
	    {
	    	printf("PAUSE OR UNPAUSE!\n");
	    	settings_pauseButton.Reset();
	    	globalSettings.global_Pause = !globalSettings.global_Pause;
	    }
#endif


		// Swap frames and Process Events
		pangolin::FinishFrame();


	    if(needReset) reset_internal();

	    if(pangolin::ShouldQuit())
        {
	        shouldQuitVar = true;
        }
	}


	printf("QUIT Pangolin thread!\n");
	printf("So Long, and Thanks for All the Fish!\n");
}


void PangolinDSOViewer::close()
{
	running = false;
}

void PangolinDSOViewer::join()
{
    close();
    if(runThread.joinable())
        runThread.join();
	printf("JOINED Pangolin thread!\n");
}

void PangolinDSOViewer::reset()
{
	needReset = true;
}

void PangolinDSOViewer::reset_internal()
{
	model3DMutex.lock();
	for(size_t i=0; i<keyframes.size();i++) delete keyframes[i];
	keyframes.clear();
	allFramePoses.clear();
	keyframesByKFID.clear();
	connections.clear();
	model3DMutex.unlock();


	openImagesMutex.lock();
	internalVideoImg->setBlack();
	internalKFImg->setBlack();
	internalResImg->setBlack();

	internalVideoImg->HaveNewImage = true;
	internalKFImg->HaveNewImage = true;
	internalResImg->HaveNewImage = true;

	openImagesMutex.unlock();

	needReset = false;
}


void PangolinDSOViewer::drawConstraints()
{
	if(settings_showAllConstraints)
	{
		// draw constraints
		glLineWidth(1);
		glColor3f(0.625,0.85,0.88);
		glBegin(GL_LINES);

		for(unsigned int i=0;i<connections.size();i++)
		{
			if(connections[i].to == 0 || connections[i].from==0) continue;
			int nAct = connections[i].bwdAct + connections[i].fwdAct;
			int nMarg = connections[i].bwdMarg + connections[i].fwdMarg;
			if(nAct==0 && nMarg>0  )
			{
				Sophus::Vector3f t = connections[i].from->camToWorld.translation().cast<float>();
				glVertex3f((GLfloat) t[0],(GLfloat) t[1], (GLfloat) t[2]);
				t = connections[i].to->camToWorld.translation().cast<float>();
				glVertex3f((GLfloat) t[0],(GLfloat) t[1], (GLfloat) t[2]);
			}
		}
		glEnd();
	}

	if(settings_showActiveConstraints)
	{
		glLineWidth(3);
		glColor3f(0.32,0.7,0.8);
		glBegin(GL_LINES);

		for(unsigned int i=0;i<connections.size();i++)
		{
			if(connections[i].to == 0 || connections[i].from==0) continue;
			int nAct = connections[i].bwdAct + connections[i].fwdAct;

			if(nAct>0)
			{
				Sophus::Vector3f t = connections[i].from->camToWorld.translation().cast<float>();
				glVertex3f((GLfloat) t[0],(GLfloat) t[1], (GLfloat) t[2]);
				t = connections[i].to->camToWorld.translation().cast<float>();
				glVertex3f((GLfloat) t[0],(GLfloat) t[1], (GLfloat) t[2]);
			}
		}
		glEnd();
	}

	if(settings_showTrajectory)
	{
        float colorYellow[3] = {0.94,0.82,0.0};
		glColor3f(colorYellow[0],colorYellow[1],colorYellow[2]);
		glLineWidth(3);

		glBegin(GL_LINE_STRIP);
		for(unsigned int i=0;i<keyframes.size();i++)
		{
			glVertex3f((float)keyframes[i]->camToWorld.translation()[0],
					(float)keyframes[i]->camToWorld.translation()[1],
					(float)keyframes[i]->camToWorld.translation()[2]);
		}
		glEnd();
	}

	if(settings_showFullTrajectory)
	{
        float colorRed[3] = {0.875,0.1,0.3};
		glColor3f(colorRed[0],colorRed[1],colorRed[2]);
		glLineWidth(3);

		glBegin(GL_LINE_STRIP);
		for(unsigned int i=0;i<allFramePoses.size();i++)
		{
			glVertex3f((float)allFramePoses[i][0],
					(float)allFramePoses[i][1],
					(float)allFramePoses[i][2]);
		}
		glEnd();
	}
}



void PangolinDSOViewer::publishGraph(const std::map<uint64_t, Eigen::Vector2i, std::less<uint64_t>, Eigen::aligned_allocator<std::pair<const uint64_t, Eigen::Vector2i>>> &connectivity)
{
    if(!setting_render_display3D) return;
    if(globalSettings.setting_disableAllDisplay) return;

	model3DMutex.lock();
    connections.resize(connectivity.size());
	int runningID=0;
	int totalActFwd=0, totalActBwd=0, totalMargFwd=0, totalMargBwd=0;
    for(std::pair<uint64_t,Eigen::Vector2i> p : connectivity)
	{
		int host = (int)(p.first >> 32);
        int target = (int)(p.first & (uint64_t)0xFFFFFFFF);

		assert(host >= 0 && target >= 0);
		if(host == target)
		{
			assert(p.second[0] == 0 && p.second[1] == 0);
			continue;
		}

		if(host > target) continue;

		connections[runningID].from = keyframesByKFID.count(host) == 0 ? 0 : keyframesByKFID[host];
		connections[runningID].to = keyframesByKFID.count(target) == 0 ? 0 : keyframesByKFID[target];
		connections[runningID].fwdAct = p.second[0];
		connections[runningID].fwdMarg = p.second[1];
		totalActFwd += p.second[0];
		totalMargFwd += p.second[1];

        uint64_t inverseKey = (((uint64_t)target) << 32) + ((uint64_t)host);
		Eigen::Vector2i st = connectivity.at(inverseKey);
		connections[runningID].bwdAct = st[0];
		connections[runningID].bwdMarg = st[1];

		totalActBwd += st[0];
		totalMargBwd += st[1];

		runningID++;
	}


	model3DMutex.unlock();
}
void PangolinDSOViewer::publishKeyframes(
		std::vector<FrameHessian*> &frames,
		bool final,
		CalibHessian* HCalib)
{
	if(!setting_render_display3D) return;
    if(globalSettings.setting_disableAllDisplay) return;

	boost::unique_lock<boost::mutex> lk(model3DMutex);

	if (!final) std::fill(active_frame_IDs.begin(), active_frame_IDs.end(), -1);
	int active_frame_count = 0;

	for(FrameHessian* fh : frames)
	{
		if(keyframesByKFID.find(fh->frameID) == keyframesByKFID.end())
		{
			KeyFrameDisplay* kfd = new KeyFrameDisplay(w,h);
			keyframesByKFID[fh->frameID] = kfd;
			keyframes.push_back(kfd);
		}
		keyframesByKFID[fh->frameID]->setFromKF(fh, HCalib);

		if (!final){
			active_frame_IDs[active_frame_count] = fh->frameID;
			active_frame_count++;
		}
    }
}

void PangolinDSOViewer::publishCamPose(FrameShell* frame,
		CalibHessian* HCalib)
{
    if(!setting_render_display3D) return;
    if(globalSettings.setting_disableAllDisplay) return;

	boost::unique_lock<boost::mutex> lk(model3DMutex);
	struct timeval time_now;
	gettimeofday(&time_now, NULL);
	lastNTrackingMs.push_back(((time_now.tv_sec-last_track.tv_sec)*1000.0f + (time_now.tv_usec-last_track.tv_usec)/1000.0f));
	if(lastNTrackingMs.size() > 10) lastNTrackingMs.pop_front();
	last_track = time_now;

	if(!setting_render_display3D) return;

	this->HCalib = HCalib;

	currentCam->setFromF(frame, HCalib);
	allFramePoses.push_back(frame->camToWorld.translation().cast<float>());
}


void PangolinDSOViewer::pushLiveFrame(FrameHessian* image)
{
	if(!setting_render_displayVideo) return;
    if(globalSettings.setting_disableAllDisplay) return;

	boost::unique_lock<boost::mutex> lk(openImagesMutex);

	if (!image->colourValid){
		for(int i=0;i<w*h;i++)
			internalVideoImg->data[i][0] =
			internalVideoImg->data[i][1] =
			internalVideoImg->data[i][2] =
				image->dI[i][0]*0.8 > 255.0f ? 255.0 : image->dI[i][0]*0.8;
	} else {
		for(int i=0;i<w*h;i++){
				internalVideoImg->data[i][0] = image->dI_c[i][0]*0.8 > 255.0f ? 255.0 : image->dI_c[i][0]*0.8;
				internalVideoImg->data[i][1] = image->dI_c[i][1]*0.8 > 255.0f ? 255.0 : image->dI_c[i][1]*0.8;
				internalVideoImg->data[i][2] = image->dI_c[i][2]*0.8 > 255.0f ? 255.0 : image->dI_c[i][2]*0.8;
			}
	}

	internalVideoImg->HaveNewImage = true;
}

bool PangolinDSOViewer::needPushDepthImage()
{
    return setting_render_displayDepth;
}
void PangolinDSOViewer::pushDepthImage(MinimalImageB3* image)
{

    if(!setting_render_displayDepth) return;
    if(globalSettings.setting_disableAllDisplay) return;

	boost::unique_lock<boost::mutex> lk(openImagesMutex);

	struct timeval time_now;
	gettimeofday(&time_now, NULL);
	lastNMappingMs.push_back(((time_now.tv_sec-last_map.tv_sec)*1000.0f + (time_now.tv_usec-last_map.tv_usec)/1000.0f));
	if(lastNMappingMs.size() > 10) lastNMappingMs.pop_front();
	last_map = time_now;

	memcpy(internalKFImg->data, image->data, w*h*3);
	internalKFImg->HaveNewImage = true;
}

void PangolinDSOViewer::publishTransformDSOToIMU(const dmvio::TransformDSOToIMU& transformDSOToIMUPassed)
{
    if(!setting_render_display3D) return;
    if(globalSettings.setting_disableAllDisplay) return;

    boost::unique_lock<boost::mutex> lk(model3DMutex);
    transformDSOToIMU = std::make_unique<dmvio::TransformDSOToIMU>(transformDSOToIMUPassed, std::make_shared<bool>(false),
            std::make_shared<bool>(false), std::make_shared<bool>(false));
}

void PangolinDSOViewer::publishSystemStatus(dmvio::SystemStatus systemStatus)
{
    boost::unique_lock<boost::mutex> lk(model3DMutex);
    this->systemStatus = systemStatus;
}

void PangolinDSOViewer::addGTCamPose(const Sophus::SE3d& gtPose)
{
	boost::unique_lock<boost::mutex> lk(model3DMutex);

	if(!setting_render_display3D || !HCalib) return;

	if (!setting_debugout_runquiet && !globalSettings.no_GT_debugMessage) std::cout << "GTPose: " << gtPose.translation().transpose() << std::endl;

	if(!gtCamPoseSet)
    {
	    firstGTCamPoseMetric = gtPose;
	    firstCamPoseDSO = currentCam->camToWorld; // Needed to make sure the first pose is the same for both.
    }

	gtCamPoseMetric = gtPose;
    gtCamPoseSet = true;
	updateDisplayedCamPose();

}

// Caller should aquire model lock for us.
void PangolinDSOViewer::updateDisplayedCamPose()
{
    if(!gtCamPoseSet || !transformDSOToIMU) return;
    if(!setting_render_display3D || !HCalib) return;

    // The visualizer shows cam to world in dso scale. The groundtruth pose is imu to world in metric scale.
    // This transforms to worldToCam
    SE3 worldToCam(transformDSOToIMU->transformPoseInverse(gtCamPoseMetric.matrix()));

    SE3 firstGTWorldToCam(transformDSOToIMU->transformPoseInverse(firstGTCamPoseMetric.matrix()));

    // We want the first pose to stay the same:
    // firstPose = offset * gtFirstPose;
    // --> offset = firstPose * gtFirstPose^(-1)
    SE3 offset = firstCamPoseDSO * firstGTWorldToCam;
    SE3 gtPoseTransformed =  offset * worldToCam.inverse();

    currentGTCam->setFromPose(gtPoseTransformed, HCalib);

}

bool PangolinDSOViewer::shouldQuit()
{
    return shouldQuitVar;
}

}
}
