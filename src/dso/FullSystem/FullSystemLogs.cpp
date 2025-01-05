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

#include "FullSystem/FullSystem.h"

#include "stdio.h"
#include "util/globalFuncs.h"
#include <Eigen/LU>
#include <algorithm>
#include "IOWrapper/ImageDisplay.h"
#include "IOWrapper/ImageRW.h"
#include "util/globalCalib.h"
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <algorithm>

#include "FullSystem/ImmaturePoint.h"
#include "util/TimeMeasurement.h"


namespace dso
{
#ifdef GRAPHICAL_DEBUG
	void FullSystem::debugPlotTracking()
	{
		// Shows the points of each active frame in each active frame
		if(globalSettings.setting_disableAllDisplay) return;
		if(!globalSettings.setting_render_plotTrackingFull) return;
		int wh = globalCalib.hG[0]*globalCalib.wG[0];

		int idx=0;
		for(FrameHessian* f : frameHessians)
		{
			std::vector<MinimalImageB3* > images;

			// make images for all frames. will be deleted by the FrameHessian's destructor.
			for(FrameHessian* f2 : frameHessians)
				if(f2->debugImage == 0) f2->debugImage = new MinimalImageB3( globalCalib.wG[0],  globalCalib.hG[0]);

			for(FrameHessian* f2 : frameHessians)
			{
				MinimalImageB3* debugImage=f2->debugImage;
				images.push_back(debugImage);

				Eigen::Vector3f* fd = f2->dI;

				Vec2 affL = AffLight::fromToVecExposure(f2->ab_exposure, f->ab_exposure, f2->aff_g2l(), f->aff_g2l());

				for(int i=0;i<wh;i++)
				{
					// BRIGHTNESS TRANSFER
					float colL = affL[0] * fd[i][0] + affL[1];
					if(colL<0) colL=0; if(colL>255) colL =255;
					debugImage->at(i) = Vec3b(colL, colL, colL);
				}
			}


			for(PointHessian* ph : f->pointHessians)
			{
				assert(ph->status == PointHessian::ACTIVE);
				if(ph->status == PointHessian::ACTIVE || ph->status == PointHessian::MARGINALIZED)
				{
					for(PointFrameResidual* r : ph->residuals)
						r->debugPlot();
					f->debugImage->setPixel9(ph->u+0.5, ph->v+0.5, makeRainbow3B(ph->idepth_scaled));
				}
			}

			char buf[100];
			snprintf(buf, 100, "IMG %d", idx);
			if(!globalSettings.setting_disableAllDisplay) IOWrap::displayImageStitch(buf, images, globalSettings.setting_maxFrames);
			idx++;
		}
		if(!globalSettings.setting_disableAllDisplay) IOWrap::waitKey(0);
	}
#endif

#ifdef GRAPHICAL_DEBUG
	void FullSystem::debugPlot(std::string name)
	{
        dmvio::TimeMeasurement timeMeasurement("debugPlot");
		if(globalSettings.setting_disableAllDisplay) return;
		if(!globalSettings.setting_render_renderWindowFrames) return;
		std::vector<MinimalImageB3* > images;

		float minID=0, maxID=0;
		if((int)(freeDebugParam5+0.5f) == 7 || (globalSettings.setting_debugSaveImages&&false))
		{
			std::vector<float> allID;
			for(unsigned int f=0;f<frameHessians.size();f++)
			{
				for(PointHessian* ph : frameHessians[f]->pointHessians)
					if(ph!=0) allID.push_back(ph->idepth_scaled);

				for(PointHessian* ph : frameHessians[f]->pointHessiansMarginalized)
					if(ph!=0) allID.push_back(ph->idepth_scaled);

				for(PointHessian* ph : frameHessians[f]->pointHessiansOut)
					if(ph!=0) allID.push_back(ph->idepth_scaled);
			}
			std::sort(allID.begin(), allID.end());
			int n = allID.size()-1;
			minID = allID[(int)(n*0.05)];
			maxID = allID[(int)(n*0.95)];


			// slowly adapt: change by maximum 10% of old span.
			float maxChange = 0.1*(maxIdJetVisDebug - minIdJetVisDebug);
			if(maxIdJetVisDebug < 0  || minIdJetVisDebug < 0 ) maxChange = 1e5;


			if(minID < minIdJetVisDebug - maxChange)
				minID = minIdJetVisDebug - maxChange;
			if(minID > minIdJetVisDebug + maxChange)
				minID = minIdJetVisDebug + maxChange;


			if(maxID < maxIdJetVisDebug - maxChange)
				maxID = maxIdJetVisDebug - maxChange;
			if(maxID > maxIdJetVisDebug + maxChange)
				maxID = maxIdJetVisDebug + maxChange;

			maxIdJetVisDebug = maxID;
			minIdJetVisDebug = minID;
		}


		int wh = globalCalib.hG[0]*globalCalib.wG[0];
		for(unsigned int f=0;f<frameHessians.size();f++)
		{
			MinimalImageB3* img = new MinimalImageB3(globalCalib.wG[0],globalCalib.hG[0]);
			images.push_back(img);
			//float* fd = frameHessians[f]->I;
			Eigen::Vector3f* fd = frameHessians[f]->dI;


			for(int i=0;i<wh;i++)
			{
				int c = fd[i][0]*0.9f;
				if(c>255) c=255;
				img->at(i) = Vec3b(c,c,c);
			}

			if((int)(freeDebugParam5+0.5f) == 0)
			{
				// Mode 0: points are coloured by depth, marginalized points are coloured by depth, outlier points are white
				for(PointHessian* ph : frameHessians[f]->pointHessians)
				{
					if(ph==0) continue;

					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, makeRainbow3B(ph->idepth_scaled));
				}
				for(PointHessian* ph : frameHessians[f]->pointHessiansMarginalized)
				{
					if(ph==0) continue;
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, makeRainbow3B(ph->idepth_scaled));
				}
				for(PointHessian* ph : frameHessians[f]->pointHessiansOut)
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(255,255,255));
			}
			else if((int)(freeDebugParam5+0.5f) == 1)
			{
				// Mode 1: points are coloured by depth, marginalized points are black, outlier points are white
				for(PointHessian* ph : frameHessians[f]->pointHessians)
				{
					if(ph==0) continue;
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, makeRainbow3B(ph->idepth_scaled));
				}

				for(PointHessian* ph : frameHessians[f]->pointHessiansMarginalized)
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,0,0));

				for(PointHessian* ph : frameHessians[f]->pointHessiansOut)
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(255,255,255));
			}
			else if((int)(freeDebugParam5+0.5f) == 2)
			{
				// Mode 2: No points
			}
			else if((int)(freeDebugParam5+0.5f) == 3)
			{
				// Mode 3: Immature points are coloured by depth
				for(ImmaturePoint* ph : frameHessians[f]->immaturePoints)
				{
					if(ph==0) continue;
					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_GOOD ||
							ph->lastTraceStatus==ImmaturePointStatus::IPS_SKIPPED ||
							ph->lastTraceStatus==ImmaturePointStatus::IPS_BADCONDITION)
					{
						if(!std::isfinite(ph->idepth_max))
							img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,0,0));
						else
						{
							img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, makeRainbow3B((ph->idepth_min + ph->idepth_max)*0.5f));
						}
					}
				}
			}
			else if((int)(freeDebugParam5+0.5f) == 4)
			{
				// Mode 4: Immature points are coloured by type
				for(ImmaturePoint* ph : frameHessians[f]->immaturePoints)
				{
					if(ph==0) continue;

					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_GOOD)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,255,0)); // Green
					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_OOB)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(255,0,0)); // Red (will not be seen)
					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_OUTLIER ||
						ph->lastTraceStatus==ImmaturePointStatus::IPS_OUTLIER_OUT)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,0,255));	// Blue
					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_SKIPPED)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(255,255,0)); // Cyan
					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_BADCONDITION)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(255,0,255)); // Magenta
					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_UNINITIALIZED)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,255,255)); // Yellow
				}
			}
			else if((int)(freeDebugParam5+0.5f) == 5)
			{
				// Mode 5: Immature points are coloured by energy quality
				for(ImmaturePoint* ph : frameHessians[f]->immaturePoints)
				{
					if(ph==0) continue;

					if(ph->lastTraceStatus==ImmaturePointStatus::IPS_UNINITIALIZED) continue;
					float d = freeDebugParam1 * (sqrtf(ph->quality)-1);
					if(d<0) d=0;
					if(d>1) d=1;
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,d*255,(1-d)*255));
				}

			}
			else if((int)(freeDebugParam5+0.5f) == 6)
			{
				// Mode 6: Points are coloured by block level it was chosen at
				for(PointHessian* ph : frameHessians[f]->pointHessians)
				{
					if(ph==0) continue;
					if(ph->my_type==1)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,255,0)); // Green
					if(ph->my_type==2)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(255,0,0)); // Red
					if(ph->my_type==4)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,0,255)); // Blue
				}
				for(PointHessian* ph : frameHessians[f]->pointHessiansMarginalized)
				{
					if(ph==0) continue;
					if(ph->my_type==1)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,255,0)); // Green
					if(ph->my_type==2)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(255,0,0)); // Red
					if(ph->my_type==4)
						img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,0,255)); // Blue
				}

			}
			if((int)(freeDebugParam5+0.5f) == 7)
			{
				// Mode 7: Points are coloured by normalize depth, marginalized points are black
				for(PointHessian* ph : frameHessians[f]->pointHessians)
				{
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, makeJet3B((ph->idepth_scaled-minID) / ((maxID-minID))));
				}
				for(PointHessian* ph : frameHessians[f]->pointHessiansMarginalized)
				{
					if(ph==0) continue;
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,0,0));
				}
			}
		}
		if(!globalSettings.setting_disableAllDisplay) {
			IOWrap::displayImageStitch(name.c_str(), images, globalSettings.setting_maxFrames);
			IOWrap::waitKey(5);
		}

		for(unsigned int i=0;i<images.size();i++)
			delete images[i];



		if((globalSettings.setting_debugSaveImages&&false))
		{
			for(unsigned int f=0;f<frameHessians.size();f++)
			{
				MinimalImageB3* img = new MinimalImageB3(globalCalib.wG[0],globalCalib.hG[0]);
				Eigen::Vector3f* fd = frameHessians[f]->dI;

				for(int i=0;i<wh;i++)
				{
					int c = fd[i][0]*0.9f;
					if(c>255) c=255;
					img->at(i) = Vec3b(c,c,c);
				}

				for(PointHessian* ph : frameHessians[f]->pointHessians)
				{
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, makeJet3B((ph->idepth_scaled-minID) / ((maxID-minID))));
				}
				for(PointHessian* ph : frameHessians[f]->pointHessiansMarginalized)
				{
					if(ph==0) continue;
					img->setPixelCirc(ph->u+0.5f, ph->v+0.5f, Vec3b(0,0,0));
				}

				char buf[1000];
				snprintf(buf, 1000, "images_out/kf_%05d_%05d_%02u.png",
						frameHessians.back()->shell->id,  frameHessians.back()->frameID, f);
				IOWrap::writeImage(buf,img);

				delete img;
			}
		}
	}
#endif


void FullSystem::printResult(std::string file, bool onlyLogKFPoses, bool saveMetricPoses, bool useCamToTrackingRef)
{
	boost::unique_lock<boost::mutex> lock(trackMutex);
	boost::unique_lock<boost::mutex> crlock(shellPoseMutex);

	std::ofstream myfile;
	myfile.open (file.c_str());
	myfile << std::setprecision(15);

	for(FrameShell* s : allFrameHistory)
	{
		if(!s->poseValid) continue;

		if(onlyLogKFPoses && s->marginalizedAt == s->id) continue;

        // firstPose is transformFirstToWorld. We actually want camToFirst here ->
        Sophus::SE3d camToWorld = s->camToWorld;

        // Use camToTrackingReference for nonKFs and the updated camToWorld for KFs.
        if(useCamToTrackingRef && s->keyframeId == -1)
        {
            camToWorld = s->trackingRef->camToWorld * s->camToTrackingRef;
        }
        Sophus::SE3d camToFirst = firstPose.inverse() * camToWorld;

        if(saveMetricPoses)
        {
            // Transform pose to IMU frame.
            // not actually camToFirst any more...
            camToFirst = Sophus::SE3d(imuIntegration->getTransformDSOToIMU().transformPose(camToWorld.inverse().matrix()));
        }

 		myfile << s->timestamp <<
			" " << camToFirst.translation().x() <<
            " " << camToFirst.translation().y() <<
            " " << camToFirst.translation().z() <<
			" " << camToFirst.so3().unit_quaternion().x()<<
			" " << camToFirst.so3().unit_quaternion().y()<<
			" " << camToFirst.so3().unit_quaternion().z()<<
			" " << camToFirst.unit_quaternion().w() << "\n";
	}
	myfile.close();
}

void FullSystem::printPC(std::string file)
{
	boost::unique_lock<boost::mutex> lock(trackMutex);
	boost::unique_lock<boost::mutex> crlock(shellPoseMutex);

	std::ofstream myfile;
	myfile.open (file.c_str());
	myfile << std::setprecision(9);

	//unsigned long totalpcs = allMargPointsHistory.size();
	unsigned long totalpcs = allMargPointsHistory.size()+allFrameHistory.size();
	
	myfile << std::string("# .PCD v.6 - Point Cloud Data file format\n");
	myfile << std::string("FIELDS x y z rgb\n");
	myfile << std::string("SIZE 4 4 4 4\n");
	myfile << std::string("TYPE F F F F\n");
	myfile << std::string("COUNT 1 1 1 1\n");
	myfile << std::string("WIDTH ") << totalpcs << std::string("\n");
	myfile << std::string("HEIGHT 1\n");
	myfile << std::string("#VIEWPOINT 0 0 0 1 0 0 0\n");
	myfile << std::string("POINTS ") << totalpcs << std::string("\n");
	myfile << std::string("DATA ascii\n");
	
	std::unordered_map<unsigned long, PC_output>::iterator itr; 
	for (itr = allMargPointsHistory.begin(); itr != allMargPointsHistory.end(); ++itr)  
	{
		PC_output tmp_PC = itr->second;
		float rgb;
		unsigned char b[] = {tmp_PC.r, tmp_PC.g, tmp_PC.b, 0};
		memcpy(&rgb, &b, sizeof(rgb));

		myfile << tmp_PC.x << " " << tmp_PC.y << " " << tmp_PC.z << " " << rgb << "\n";
	} 

	// Show trajectory in point cloud
	for (FrameShell* s : allFrameHistory)  
	{
		Sophus::SE3d camToWorld = s->camToWorld;
		Sophus::SE3d camToFirst = firstPose.inverse() * camToWorld;
		float rgb;
		unsigned char b[] = {0, 255, 0, 0};
		memcpy(&rgb, &b, sizeof(rgb));
		
		myfile << camToFirst.translation().x() <<
            " " << camToFirst.translation().y() <<
            " " << camToFirst.translation().z() << " " << rgb << "\n";
	} 

	myfile.close();
}


/**
 * @brief For debugging
 * 
 */
void FullSystem::printLogLine()
{
    dmvio::TimeMeasurement timeMeasurementMargFrames("printLogLine");
	
	if(frameHessians.size()==0) return;

    if(!setting_debugout_runquiet && !globalSettings.no_FullSystem_debugMessage)
        printf("LOG %d: %.3f fine. Res: %d A, %d L, %d M; (%'d / %'d) forceDrop. a=%f, b=%f. Window %d (%d)\n",
                allKeyFramesHistory.back()->id,
                statistics_lastFineTrackRMSE,
                ef->resInA,
                ef->resInL,
                ef->resInM,
                (int)statistics_numForceDroppedResFwd,
                (int)statistics_numForceDroppedResBwd,
                allKeyFramesHistory.back()->aff_g2l.a,
                allKeyFramesHistory.back()->aff_g2l.b,
                frameHessians.back()->shell->id - frameHessians.front()->shell->id,
                (int)frameHessians.size());

	if(globalSettings.setting_nologStuff) return;

	if(numsLog != 0)
	{
		(*numsLog) << allKeyFramesHistory.back()->id << " "  <<
				statistics_lastFineTrackRMSE << " "  <<
				(int)statistics_numCreatedPoints << " "  <<
				(int)statistics_numActivatedPoints << " "  <<
				(int)statistics_numDroppedPoints << " "  <<
				(int)statistics_lastNumOptIts << " "  <<
				ef->resInA << " "  <<
				ef->resInL << " "  <<
				ef->resInM << " "  <<
				statistics_numMargResFwd << " "  <<
				statistics_numMargResBwd << " "  <<
				statistics_numForceDroppedResFwd << " "  <<
				statistics_numForceDroppedResBwd << " "  <<
				frameHessians.back()->aff_g2l().a << " "  <<
				frameHessians.back()->aff_g2l().b << " "  <<
				frameHessians.back()->shell->id - frameHessians.front()->shell->id << " "  <<
				(int)frameHessians.size() << " "  << "\n";
		numsLog->flush();
	}
}

/**
 * @brief For debugging the energy function
 * 
 */
void FullSystem::printEigenValLine()
{
    dmvio::TimeMeasurement timeMeasurementMargFrames("printEigenValLine");
	if(globalSettings.setting_nologStuff) return;
	if(ef->lastHS.rows() < 12) return;

	MatXX Hp = ef->lastHS.bottomRightCorner(ef->lastHS.cols()-CPARS,ef->lastHS.cols()-CPARS);
	MatXX Ha = ef->lastHS.bottomRightCorner(ef->lastHS.cols()-CPARS,ef->lastHS.cols()-CPARS);
	int n = Hp.cols()/8;
	assert(Hp.cols()%8==0);

	// sub-select
	for(int i=0;i<n;i++)
	{
		MatXX tmp6 = Hp.block(i*8,0,6,n*8);
		Hp.block(i*6,0,6,n*8) = tmp6;

		MatXX tmp2 = Ha.block(i*8+6,0,2,n*8);
		Ha.block(i*2,0,2,n*8) = tmp2;
	}
	for(int i=0;i<n;i++)
	{
		MatXX tmp6 = Hp.block(0,i*8,n*8,6);
		Hp.block(0,i*6,n*8,6) = tmp6;

		MatXX tmp2 = Ha.block(0,i*8+6,n*8,2);
		Ha.block(0,i*2,n*8,2) = tmp2;
	}

	VecX eigenvaluesAll = ef->lastHS.eigenvalues().real();
	VecX eigenP = Hp.topLeftCorner(n*6,n*6).eigenvalues().real();
	VecX eigenA = Ha.topLeftCorner(n*2,n*2).eigenvalues().real();
	VecX diagonal = ef->lastHS.diagonal();

	std::sort(eigenvaluesAll.data(), eigenvaluesAll.data()+eigenvaluesAll.size());
	std::sort(eigenP.data(), eigenP.data()+eigenP.size());
	std::sort(eigenA.data(), eigenA.data()+eigenA.size());

	int nz = std::max(100,globalSettings.setting_maxFrames*10);

	if(eigenAllLog != 0)
	{
		VecX ea = VecX::Zero(nz); ea.head(eigenvaluesAll.size()) = eigenvaluesAll;
		(*eigenAllLog) << allKeyFramesHistory.back()->id << " " <<  ea.transpose() << "\n";
		eigenAllLog->flush();
	}
	if(eigenALog != 0)
	{
		VecX ea = VecX::Zero(nz); ea.head(eigenA.size()) = eigenA;
		(*eigenALog) << allKeyFramesHistory.back()->id << " " <<  ea.transpose() << "\n";
		eigenALog->flush();
	}
	if(eigenPLog != 0)
	{
		VecX ea = VecX::Zero(nz); ea.head(eigenP.size()) = eigenP;
		(*eigenPLog) << allKeyFramesHistory.back()->id << " " <<  ea.transpose() << "\n";
		eigenPLog->flush();
	}

	if(DiagonalLog != 0)
	{
		VecX ea = VecX::Zero(nz); ea.head(diagonal.size()) = diagonal;
		(*DiagonalLog) << allKeyFramesHistory.back()->id << " " <<  ea.transpose() << "\n";
		DiagonalLog->flush();
	}

	if(variancesLog != 0)
	{
		VecX ea = VecX::Zero(nz); ea.head(diagonal.size()) = ef->lastHS.inverse().diagonal();
		(*variancesLog) << allKeyFramesHistory.back()->id << " " <<  ea.transpose() << "\n";
		variancesLog->flush();
	}

	std::vector<VecX> &nsp = ef->lastNullspaces_forLogging;
	(*nullspacesLog) << allKeyFramesHistory.back()->id << " ";
	for(unsigned int i=0;i<nsp.size();i++)
		(*nullspacesLog) << nsp[i].dot(ef->lastHS * nsp[i]) << " " << nsp[i].dot(ef->lastbS) << " " ;
	(*nullspacesLog) << "\n";
	nullspacesLog->flush();

}

/**
 * @brief For debugging the frames
 * 
 */
void FullSystem::printFrameLifetimes()
{
	if(globalSettings.setting_nologStuff) return;

	boost::unique_lock<boost::mutex> lock(trackMutex);

	std::ofstream* lg = new std::ofstream();
	lg->open("logs/lifetimeLog.txt", std::ios::trunc | std::ios::out);
	lg->precision(15);

	for(FrameShell* s : allFrameHistory)
	{
		(*lg) << s->id
			<< " " << s->marginalizedAt
			<< " " << s->statistics_goodResOnThis
			<< " " << s->statistics_outlierResOnThis
			<< " " << s->movedByOpt;

		(*lg) << "\n";
	}

	lg->close();
	delete lg;
}

void FullSystem::printEvalLine()
{
	return;
}

}
