/**
* This file is based on the file main_dso_pangolin.cpp of the project DSO written by Jakob Engel.
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

// Main file for running on datasets, based on the main file of DSO.

#include "util/MainSettings.h"
#include <thread>
#include <memory>
#include <locale.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "IOWrapper/Output3DWrapper.h"
#include "IOWrapper/ImageDisplay.h"


#include <boost/thread.hpp>
#include "dso/util/settings.h"
#include "dso/util/globalFuncs.h"
#include "dso/util/DatasetReader.h"
#include "dso/util/globalCalib.h"
#include "util/TimeMeasurement.h"

#include "dso/util/NumType.h"
#include "FullSystem/FullSystem.h"
#include "OptimizationBackend/MatrixAccumulators.h"
#include "FullSystem/PixelSelector2.h"

#include <util/SettingsUtil.h>

#include "IOWrapper/Pangolin/PangolinDSOViewer.h"
#include "IOWrapper/OutputWrapper/SampleOutputWrapper.h"

#include <cxxopts.hpp>



std::string gtFile = "";
std::string source = "";
std::string imuFile = "";
bool pauseEnd = false;

bool reverse = false;
int start = 0;
int end = 100000;
int maxPreloadImages = 0;       // If set we only preload if there are less images to be loade.
bool useSampleOutput = false;

using namespace dso;

dso::Global_Calib globalCalib;
dso::GlobalSettings globalSettings;
dmvio::MainSettings mainSettings;
dmvio::IMUCalibration imuCalibration;
dmvio::IMUSettings imuSettings;



void my_exit_handler(int s)
{
    printf("Caught signal %d\n", s);
    exit(1);
}

void exitThread()
{
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = my_exit_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    while(true) pause();
}



void run(std::shared_ptr<FullSystem> fullSystem, ImageFolderReader* reader, IOWrap::PangolinDSOViewer* viewer)
{
    // Handle Settings
    if(globalSettings.setting_photometricCalibration > 0 && reader->getPhotometricGamma() == 0)
    {
        printf("ERROR: dont't have photometric calibation. Need to use commandline options mode=1 or mode=2 ");
        exit(1);
    }

    int lstart = start;
    int lend = end;
    int linc = 1;
    if(reverse)
    {
        assert(!imuSettings.setting_useIMU);    // Reverse is not supported with IMU data at the moment
        printf("REVERSE!!!");
        lstart = end - 1;
        if(lstart >= reader->getNumImages())
            lstart = reader->getNumImages() - 1;
        lend = start;
        linc = -1;
    }

    bool linearizeOperation = (mainSettings.playbackSpeed == 0);

    if(linearizeOperation && globalSettings.setting_minFramesBetweenKeyframes < 0)
    {
        globalSettings.setting_minFramesBetweenKeyframes = -globalSettings.setting_minFramesBetweenKeyframes;
        std::cout << "Using setting_minFramesBetweenKeyframes=" << globalSettings.setting_minFramesBetweenKeyframes
                  << " because of non-realtime mode." << std::endl;
    }

    // Set GUI
    if(viewer != 0)
    {
        fullSystem->outputWrapper.push_back(viewer);
    }

    // Create Output Wrapper for GUI
    std::unique_ptr<IOWrap::SampleOutputWrapper> sampleOutPutWrapper;
    if(useSampleOutput)
    {
        sampleOutPutWrapper.reset(new IOWrap::SampleOutputWrapper());
        fullSystem->outputWrapper.push_back(sampleOutPutWrapper.get());
    }

    // Prepare images for system
    std::vector<int> idsToPlay;
    std::vector<double> timesToPlayAt;
    for(int i = lstart; i >= 0 && i < reader->getNumImages() && linc * i < linc * lend; i += linc)
    {
        idsToPlay.push_back(i);
        if(timesToPlayAt.size() == 0)
        {
            timesToPlayAt.push_back((double) 0);
        }else
        {
            double tsThis = reader->getTimestamp(idsToPlay[idsToPlay.size() - 1]);
            double tsPrev = reader->getTimestamp(idsToPlay[idsToPlay.size() - 2]);
            timesToPlayAt.push_back(timesToPlayAt.back() + fabs(tsThis - tsPrev) / mainSettings.playbackSpeed);
        }
    }

    if(mainSettings.preload && maxPreloadImages > 0)
    {
        if(reader->getNumImages() > maxPreloadImages)
        {
            printf("maxPreloadImages EXCEEDED! NOT PRELOADING!\n");
            mainSettings.preload = false;
        }
    }

    std::vector<ImageAndExposure*> preloadedImages;
    if(mainSettings.preload)
    {
        printf("LOADING ALL IMAGES!\n");
        for(int ii = 0; ii < (int) idsToPlay.size(); ii++)
        {
            int i = idsToPlay[ii];
            preloadedImages.push_back(reader->getImage(i));
        }
    }


    // Main Loop
    struct timeval tv_start;
    gettimeofday(&tv_start, NULL);
    clock_t started = clock();
    double sInitializerOffset = 0;

    bool gtDataThere = reader->loadGTData(gtFile);

    bool imuDataSkipped = false;
    dmvio::IMUData skippedIMUData;
    for(int ii = 0; ii < (int) idsToPlay.size(); ii++) // loop through all images
    {
#ifdef GRAPHICAL_DEBUG
        while (globalSettings.global_Pause && linearizeOperation)
        {
            usleep(5000);
        }
#endif

        if(!fullSystem->initialized)    // if not initialized: reset start time.
        {
            gettimeofday(&tv_start, NULL);
            started = clock();
            sInitializerOffset = timesToPlayAt[ii];
        }

        int i = idsToPlay[ii];

        // Load image
        ImageAndExposure* img;
        if(mainSettings.preload)
            img = preloadedImages[ii];
        else
            img = reader->getImage(i);

        // Set time information
        bool skipFrame = false;
        if(mainSettings.playbackSpeed != 0)
        {
            struct timeval tv_now;
            gettimeofday(&tv_now, NULL);
            double sSinceStart = sInitializerOffset + ((tv_now.tv_sec - tv_start.tv_sec) +
                                                       (tv_now.tv_usec - tv_start.tv_usec) / (1000.0f * 1000.0f));

            if(sSinceStart < timesToPlayAt[ii])
                usleep((int) ((timesToPlayAt[ii] - sSinceStart) * 1000 * 1000));
            else if(sSinceStart > timesToPlayAt[ii] + 0.5 + 0.1 * (ii % 2))
            {
                printf("SKIPFRAME %d (play at %f, now it is %f)!\n", ii, timesToPlayAt[ii], sSinceStart);
                skipFrame = true;
            }
        }

        // Get GT data if it exists
        dmvio::GTData data;
        bool found = false;
        if(gtDataThere)
        {
            data = reader->getGTData(i, found, !globalSettings.no_GT_debugMessage);
        }

        // Send image and IMU data into the system
        // THIS IS WHERE THE SYSTEM IS ACTUALLY RAN
        std::unique_ptr<dmvio::IMUData> imuData;
        if(imuSettings.setting_useIMU)
        {
            imuData = std::make_unique<dmvio::IMUData>(reader->getIMUData(i));
        }
        if(!skipFrame)
        {
            if(imuDataSkipped && imuData)
            {
                imuData->insert(imuData->begin(), skippedIMUData.begin(), skippedIMUData.end());
                skippedIMUData.clear();
                imuDataSkipped = false;
            }
            // FRAME IS ADDED TO SYSTEM!!!
            fullSystem->addActiveFrame(img, i, imuData.get(), (gtDataThere && found) ? &data : 0);
            if(gtDataThere && found && !globalSettings.setting_disableAllDisplay)
            {
                viewer->addGTCamPose(data.pose);
            }
        }else if(imuData)
        {
            imuDataSkipped = true;
            skippedIMUData.insert(skippedIMUData.end(), imuData->begin(), imuData->end());
        }

        delete img;


        // Restart system
        if(fullSystem->initFailed || setting_fullResetRequested)
        {
            if(ii < 250 || setting_fullResetRequested)
            {
                printf("RESETTING!\n");
                std::vector<IOWrap::Output3DWrapper*> wraps = fullSystem->outputWrapper;

                for(IOWrap::Output3DWrapper* ow : wraps) ow->reset();
                fullSystem->fullReset();

                setting_fullResetRequested = false;
            }
        }

        // Handle closing
        if(viewer != nullptr && viewer->shouldQuit())
        {
            std::cout << "User closed window -> Quit!" << std::endl;
            break;
        }

        if(fullSystem->isLost)
        {
            printf("LOST!!\n");
            break;
        }

    }
    // Handle end of operation messages and outputs
    fullSystem->cleanUpFrames();
    fullSystem->blockUntilMappingIsFinished();
    clock_t ended = clock();
    struct timeval tv_end;
    gettimeofday(&tv_end, NULL);


    // Output results
    fullSystem->printResult(imuSettings.resultsPrefix + "result.txt", false, false, true);
    fullSystem->printResult(imuSettings.resultsPrefix + "resultKFs.txt", true, false, false);
    fullSystem->printResult(imuSettings.resultsPrefix + "resultScaled.txt", false, true, true);
    if (globalSettings.setting_outputPC) fullSystem->printPC(imuSettings.resultsPrefix + "PC.PCD");

    dmvio::TimeMeasurement::saveResults(imuSettings.resultsPrefix + "timings.txt");


    // Print info about run to cmd line
    int numFramesProcessed = abs(idsToPlay[0] - idsToPlay.back());
    double numSecondsProcessed = fabs(reader->getTimestamp(idsToPlay[0]) - reader->getTimestamp(idsToPlay.back()));
    double MilliSecondsTakenSingle = 1000.0f * (ended - started) / (float) (CLOCKS_PER_SEC);
    double MilliSecondsTakenMT = sInitializerOffset + ((tv_end.tv_sec - tv_start.tv_sec) * 1000.0f +
                                                       (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f);
    printf("\n======================"
           "\n%d Frames (%.1f fps)"
           "\n%.2fms per frame (single core); "
           "\n%.2fms per frame (multi core); "
           "\n%.3fx (single core); "
           "\n%.3fx (multi core); "
           "\n======================\n\n",
           numFramesProcessed, numFramesProcessed / numSecondsProcessed,
           MilliSecondsTakenSingle / numFramesProcessed,
           MilliSecondsTakenMT / (float) numFramesProcessed,
           1000 / (MilliSecondsTakenSingle / numSecondsProcessed),
           1000 / (MilliSecondsTakenMT / numSecondsProcessed));
    fullSystem->printFrameLifetimes();

    if(!globalSettings.setting_nologStuff)
    {
        std::ofstream tmlog;
        tmlog.open("logs/time.txt", std::ios::trunc | std::ios::out);
        tmlog << 1000.0f * (ended - started) / (float) (CLOCKS_PER_SEC * reader->getNumImages()) << " "
              << ((tv_end.tv_sec - tv_start.tv_sec) * 1000.0f + (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f) /
                 (float) reader->getNumImages() << "\n";
        tmlog.flush();
        tmlog.close();
    }


    // Clean up
    if(!pauseEnd){
	    for(IOWrap::Output3DWrapper* ow : fullSystem->outputWrapper)
	    {
		ow->join();
	    }
    }
}


int main(int argc, char** argv)
{
    setlocale(LC_ALL, "C");

#ifdef DEBUG
    std::cout << "DEBUG MODE!" << std::endl;
#endif

    bool use16Bit = false;
    bool useColour = false;

    // Hanlde settings and command inputs
    auto settingsUtil = std::make_shared<dmvio::SettingsUtil>();

    // Dataset specific arguments
    // cmd line handling is done by cxxopt
    settingsUtil->registerArg("files", source, "f", "Inputted file directory in alphabetical order", "");
    settingsUtil->registerArg("imuFile", imuFile, "i", "Inputted IMU file directory", "");
    settingsUtil->registerArg("gtFile", gtFile, "t", "Inputted groundtruth file directory", "");

    settingsUtil->registerArg("start", start, "s", "Start frame", std::to_string(start));
    settingsUtil->registerArg("end", end, "e", "End frame", std::to_string(end));
    settingsUtil->registerArg("reverse", reverse, "r", "Reverse frames", reverse ? "1" : "0");
    settingsUtil->registerArg("use16Bit", use16Bit, "b", "16 Bit image input", use16Bit ? "1" : "0");
    settingsUtil->registerArg("useColour", useColour, "c", "Colour image input", useColour ? "1" : "0");

    settingsUtil->registerArg("maxPreloadImages", maxPreloadImages);
    settingsUtil->registerArg("sampleoutput", useSampleOutput);	

    settingsUtil->registerArg("pauseEnd", pauseEnd, "E", "Pause at end", pauseEnd ? "1" : "0");


    // Create Settings files.
    mainSettings.registerArgs(*settingsUtil, globalSettings);
    imuSettings.registerArgs(*settingsUtil);
    imuCalibration.registerArgs(*settingsUtil);
    // This is the only config setting shared between the imu and base versions
    imuSettings.imusetting_thOptIterations = globalSettings.setting_thOptIterations;

    // This call will parse all commandline arguments and potentially also read a settings yaml file if passed
    mainSettings.parseArguments(argc, argv, *settingsUtil, globalSettings);

    // Load imu calibration
    // IMU calibration is optional but recommended
    if(mainSettings.imuCalibFile != "")
    {
        imuCalibration.loadFromFile(mainSettings.imuCalibFile);
    }

    // Print settings to commandline and file.
    if(mainSettings.print_settings){
        std::cout << "Settings:\n";
        settingsUtil->printAllSettings(std::cout);
    }

    // Load imu settings
    {
        std::ofstream settingsStream;
        settingsStream.open(imuSettings.resultsPrefix + "usedSettingsdso.txt");
        settingsUtil->printAllSettings(settingsStream);
    }

    // hook crtl+C.
    boost::thread exThread = boost::thread(exitThread);

    // Create image reader
    ImageFolderReader* reader = new ImageFolderReader(source, use16Bit, useColour);
    if (reader->readCalib(mainSettings.calib, mainSettings.gammaCalib, mainSettings.vignette, globalSettings)==0){
        delete reader;
        return 0;
    }
    reader->loadIMUData(imuFile);
    reader->setGlobalCalibration(globalCalib, globalSettings.pyrLevelsUsed);

    // Create system
    std::shared_ptr<FullSystem> fullSystem;
    fullSystem = std::make_shared<FullSystem>((mainSettings.playbackSpeed == 0), globalCalib, imuCalibration, imuSettings, globalSettings);
    fullSystem->setGammaFunction(reader->getPhotometricGamma());

    // Main operationing loop
    if(!globalSettings.setting_disableAllDisplay)
    {
        IOWrap::PangolinDSOViewer* viewer = new IOWrap::PangolinDSOViewer(globalCalib.wG[0], globalCalib.hG[0], globalSettings, false, settingsUtil,
                                                                          nullptr);

        // THIS IS WHERE THE MAIN THREAD IS RUN
        boost::thread runThread = boost::thread(boost::bind(run, fullSystem, reader, viewer));

        // THIS IS WHERE THE VIEWER IS STARTED
        viewer->run();

        // Make sure that the destructor of FullSystem, etc. finishes, so all log files are properly flushed.
        runThread.join();
        if (fullSystem)
        {
            for(IOWrap::Output3DWrapper* ow : fullSystem->outputWrapper)
            {
                if(pauseEnd) ow->join();
                delete ow;
            }
        }
    }else
    {
         // THIS IS WHERE THE MAIN THREAD IS RUN
        run(fullSystem, reader, 0);
    }

    printf("DELETE FULLSYSTEM!\n");
    fullSystem.reset();

    printf("DELETE READER!\n");
    delete reader;

    printf("EXIT NOW!\n");

    return 0;
}
