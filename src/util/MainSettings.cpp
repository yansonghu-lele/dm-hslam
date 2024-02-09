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

void MainSettings::parseArguments(int argc, char** argv, SettingsUtil& settingsUtil)
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
        setting_photometricCalibration = 0;
        setting_affineOptModeA = 0; //-1: fix. >=0: optimize (with prior, if > 0).
        setting_affineOptModeB = 0; //-1: fix. >=0: optimize (with prior, if > 0).
    }
    if(mode == 2)
    {
        printf("PHOTOMETRIC MODE WITH PERFECT IMAGES!\n");
        setting_photometricCalibration = 0;
        setting_affineOptModeA = -1; //-1: fix. >=0: optimize (with prior, if > 0).
        setting_affineOptModeB = -1; //-1: fix. >=0: optimize (with prior, if > 0).
        setting_minGradHistAdd = 0.005;
    }
    if(mode == 3)
    {
        // This mode is useful because mode 0 assumes that exposure is available (as it adds a strong prior to
        // the affine brightness change between images), and mode 1 does not use vignette at all.
        // This mode uses vignette (and response), but still fully optimizes brightness changes, hence it is
        // appropriate for sensors without exposure time but with a calibrated vignette.
        printf("PHOTOMETRIC MODE WITH CALIBRATION, BUT NO OR INACCURATE EXPOSURE!\n");
        setting_affineOptModeA = 0; //-1: fix. >=0: optimize (with prior, if > 0).
        setting_affineOptModeB = 0; //-1: fix. >=0: optimize (with prior, if > 0).
    }

    if(settingsFile!=""){
        settingsUtil.tryReadFromYaml(settingsFile);
    }

    if(result.count("print")) print_settings = true;
}


void MainSettings::registerArgs(SettingsUtil& set)
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
    set.registerArg("nolog", setting_logStuff, "n", "Turn logging off", setting_logStuff ? "1" : "0");
    set.registerArg("nogui", disableAllDisplay, "g", "Turn gui output off", disableAllDisplay ? "1" : "0");
    set.registerArg("outPC", outputPC, "p", "Output point cloud", outputPC ? "1" : "0");
    set.registerArg("useimu", setting_useIMU, "u", "Turn IMU on or off", setting_useIMU ? "1" : "0");
    set.registerArg("save", debugSaveImages, "a", "Save data", debugSaveImages ? "1" : "0");
    set.registerArg("mode", mode, "m", "Photometric mode (1=full, 2=no calibration, 3=Synthetic, 4=none)", "0");
    set.registerArg("settingsFile", settingsFile, "F", "Settings file", "");

    // Register global settings.
    set.registerArg("setting_minOptIterations", setting_minOptIterations);
    set.registerArg("setting_maxOptIterations", setting_maxOptIterations);
    set.registerArg("setting_minIdepth", setting_minIdepth);
    set.registerArg("setting_solverMode", setting_solverMode);
    set.registerArg("setting_weightZeroPriorDSOInitY", setting_weightZeroPriorDSOInitY);
    set.registerArg("setting_weightZeroPriorDSOInitX", setting_weightZeroPriorDSOInitX);
    set.registerArg("setting_forceNoKFTranslationThresh", setting_forceNoKFTranslationThresh);
    set.registerArg("setting_minFramesBetweenKeyframes", setting_minFramesBetweenKeyframes);
    set.registerArg("setting_desiredImmatureDensity", setting_desiredImmatureDensity);
    set.registerArg("setting_desiredPointDensity", setting_desiredPointDensity);
    set.registerArg("setting_minFrames", setting_minFrames);
    set.registerArg("setting_maxFrames", setting_maxFrames);

}
