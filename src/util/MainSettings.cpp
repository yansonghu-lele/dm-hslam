/**
* This file is part of DM-VIO.
*
* Copyright (c) 2022 Lukas von Stumberg <lukas dot stumberg at tum dot de>.
* for more information see <http://vision.in.tum.de/dm-vio>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* The methods parseArgument and settingsDefault are based on the file
* main_dso_pangolin.cpp of the project DSO written by Jakob Engel, but have been
* modified for the inclusion in DM-VIO. The original versions have the copyright
* Copyright 2016 Technical University of Munich and Intel.
* Developed by Jakob Engel <engelj at in dot tum dot de>,
*
* DM-VIO is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DM-VIO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DM-VIO. If not, see <http://www.gnu.org/licenses/>.
*/



#include "MainSettings.h"
#include "dso/util/settings.h"



using namespace dmvio;
using namespace dso;

/**
 * @brief Parses command line arguements
 * 
 * @param argc          cmd line number of inputs
 * @param argv          cmd line string input
 * @param settingsUtil  Class that handles settings
 */
void MainSettings::parseArguments(int argc, char** argv, SettingsUtil& settingsUtil, dso::GlobalSettings& globalSettings)
{
    cxxopts::ParseResult result = settingsUtil.cmd_options->parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << settingsUtil.cmd_options->help() << std::endl;
      exit(0);
    }

    for(const cxxopts::KeyValue &kv: result)
	{
        // Additional handling for some variables
        if(kv.key()=="quiet" && kv.value()=="1") printf("QUIET MODE!\n");
        if(kv.key()=="nolog" && kv.value()=="1") printf("DISABLE LOGGING!\n");
        if(kv.key()=="nogui" && kv.value()=="1") printf("NO GUI!\n");
        if(kv.key()=="nomt" && kv.value()=="1") printf("NO MULTITHREADING!\n");
        if(kv.key()=="save" && kv.value()=="1"){
            if(42 == system("rm -rf images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            if(42 == system("mkdir images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            if(42 == system("rm -rf images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            if(42 == system("mkdir images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            printf("SAVE IMAGES!\n");
        }

        // Process parameters
        if(kv.value() != "")
            settingsUtil.tryReadFromCommandLine(kv.key(), kv.value());
	}

    //Set photometric mode
    if(mode == 0)
    {
        printf("PHOTOMETRIC MODE WITH CALIBRATION!\n");
    }
    if(mode == 1)
    {
        printf("PHOTOMETRIC MODE WITHOUT CALIBRATION!\n");
        globalSettings.setting_photometricCalibration = 0;
        globalSettings.setting_affineOptModeA = 0; //-1: fix. >=0: optimize (with prior, if > 0).
        globalSettings.setting_affineOptModeB = 0; //-1: fix. >=0: optimize (with prior, if > 0).
    }
    if(mode == 2)
    {
        printf("PHOTOMETRIC MODE WITH PERFECT IMAGES!\n");
        globalSettings.setting_photometricCalibration = 0;
        globalSettings.setting_affineOptModeA = -1; //-1: fix. >=0: optimize (with prior, if > 0).
        globalSettings.setting_affineOptModeB = -1; //-1: fix. >=0: optimize (with prior, if > 0).
    }
    if(mode == 3)
    {
        // This mode is useful because mode 0 assumes that exposure is available (as it adds a strong prior to
        // the affine brightness change between images), and mode 1 does not use vignette at all.
        // This mode uses vignette (and response), but still fully optimizes brightness changes, hence it is
        // appropriate for sensors without exposure time but with a calibrated vignette.
        printf("PHOTOMETRIC MODE WITH CALIBRATION, BUT NO OR INACCURATE EXPOSURE!\n");
        globalSettings.setting_affineOptModeA = 0; //-1: fix. >=0: optimize (with prior, if > 0).
        globalSettings.setting_affineOptModeB = 0; //-1: fix. >=0: optimize (with prior, if > 0).
    }

    if(settingsFile!=""){
        settingsUtil.tryReadFromYaml(settingsFile);
    }

    if(result.count("print")) print_settings = true;
}

/**
 * @brief Register main settings
 * 
 * @param set 
 */
void MainSettings::registerArgs(SettingsUtil& set, dso::GlobalSettings& globalSettings)
{
    // The DM-VIO settings can also be set with commandline arguments (and also with the yaml settings file)
    set.registerArg("vignette", vignette, "V", "Photometric calibration vignette file path", "");
    set.registerArg("gamma", gammaCalib, "G", "Photometric calibration gamma file path", "");
    set.registerArg("calib", calib, "C", "Camera calibration file path", "");
    set.registerArg("imuCalib", imuCalibFile, "I", "IMU calibration file path", "");
    set.registerArg("speed", playbackSpeed, "S", "Playback speed", std::to_string(playbackSpeed));
    set.registerArg("preload", preload);

    // These are mostly the original DSO commandline arguments which also work for DM-VIO
    set.registerArg("quiet", setting_debugout_runquiet, "q", "Turn console text output off", setting_debugout_runquiet ? "1" : "0");
    set.registerArg("nolog", globalSettings.setting_nologStuff, "n", "Turn logging off", globalSettings.setting_nologStuff ? "1" : "0");
    set.registerArg("nogui", globalSettings.setting_disableAllDisplay, "g", "Turn gui output off", globalSettings.setting_disableAllDisplay ? "1" : "0");
    set.registerArg("outPC", globalSettings.setting_outputPC, "p", "Output point cloud", globalSettings.setting_outputPC ? "1" : "0");
    set.registerArg("save", globalSettings.setting_debugSaveImages, "a", "Save data", globalSettings.setting_debugSaveImages ? "1" : "0");
    set.registerArg("mode", mode, "m", "Photometric mode (1=full, 2=no calibration, 3=Synthetic, 4=none)", "0");
    set.registerArg("nomt", globalSettings.settings_no_multiThreading, "M", "Turn multithreading off", "0");
    set.registerArg("settingsFile", settingsFile, "F", "Settings file", "");

    // Register global settings
    // Mainly changed using the settings file
    set.registerArg("setting_minOptIterations", globalSettings.setting_minOptIterations);
    set.registerArg("setting_maxOptIterations", globalSettings.setting_maxOptIterations);
    set.registerArg("setting_minIdepth", globalSettings.setting_minIdepth);
    set.registerArg("setting_solverMode", globalSettings.setting_solverMode);
    set.registerArg("setting_weightZeroPriorDSOInitY", globalSettings.setting_weightZeroPriorDSOInitY);
    set.registerArg("setting_weightZeroPriorDSOInitX", globalSettings.setting_weightZeroPriorDSOInitX);
    set.registerArg("setting_forceNoKFTranslationThresh", globalSettings.setting_forceNoKFTranslationThresh);
    set.registerArg("setting_minFramesBetweenKeyframes", globalSettings.setting_minFramesBetweenKeyframes);
    set.registerArg("setting_desiredImmatureDensity", globalSettings.setting_desiredImmatureDensity);
    set.registerArg("setting_desiredPointDensity", globalSettings.setting_desiredPointDensity);
    set.registerArg("setting_minFrames", globalSettings.setting_minFrames);
    set.registerArg("setting_maxFrames", globalSettings.setting_maxFrames);
    set.registerArg("setting_GradHistCorrect", globalSettings.setting_GradHistCorrect);
    set.registerArg("setting_minGradHistAdd", globalSettings.setting_minGradHistAdd);
    set.registerArg("setting_minGradHistCut", globalSettings.setting_minGradHistCut);

    // Debug message controller
    // Turns off debug messages for functions in DSO
    set.registerArg("no_CoarseInit_debugMessage", globalSettings.no_CoarseInit_debugMessage);
    set.registerArg("no_CoarseTracker_debugMessage", globalSettings.no_CoarseTracker_debugMessage);
    set.registerArg("no_FullSystem_debugMessage", globalSettings.no_FullSystem_debugMessage);
    set.registerArg("no_Optimize_debugMessage", globalSettings.no_Optimize_debugMessage);
    set.registerArg("no_Immature_debugMessage", globalSettings.no_Immature_debugMessage);
    set.registerArg("no_Pixel_debugMessage", globalSettings.no_Pixel_debugMessage);
    set.registerArg("no_GT_debugMessage", globalSettings.no_GT_debugMessage);
}
