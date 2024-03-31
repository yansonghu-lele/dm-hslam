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

#include <string.h>
#include <string>
#include <cmath>



namespace dso
{
#define SOLVER_SVD (int)1
#define SOLVER_ORTHOGONALIZE_SYSTEM (int)2
#define SOLVER_ORTHOGONALIZE_POINTMARG (int)4
#define SOLVER_ORTHOGONALIZE_FULL (int)8
#define SOLVER_SVD_CUT7 (int)16
#define SOLVER_REMOVE_POSEPRIOR (int)32
#define SOLVER_USE_GN (int)64
#define SOLVER_FIX_LAMBDA (int)128
#define SOLVER_ORTHOGONALIZE_X (int)256
#define SOLVER_MOMENTUM (int)512
#define SOLVER_STEPMOMENTUM (int)1024
#define SOLVER_ORTHOGONALIZE_X_LATER (int)2048


// ============== PARAMETERS TO BE DECIDED ON COMPILE TIME =================
#define PYR_LEVELS 6 // Max levels used. Should be between 1 and 6

struct GlobalSettings{
    int pyrLevelsUsed = PYR_LEVELS;

    // If non-zero we set a prior to the x or y direction of the translation during the coarse visual initializer (useful for car datasets).
    double setting_weightZeroPriorDSOInitY = 0.0;
    double setting_weightZeroPriorDSOInitX = 0.0;
    double setting_forceNoKFTranslationThresh = 0.0; // Force to create no KF if translation (in metric) is smaller than this.

    double setting_maxTimeBetweenKeyframes = 0;

    // If negative, the respective positive value will be used only if in non-RT mode.
    // The idea of this parameter is that in non-RT mode the systems otherwise can make successive frames keyframes, which only rarely happens in RT mode.
    // Default is -0.5 with means that the parameter is 0.5 in non-RT mode and inactive in RT mode.
    // Fractional values are also possible.
    double setting_minFramesBetweenKeyframes = -0.5;

    // minimum idepth for keeping points in the optimization window.
    float setting_minIdepth = 0.02f;

    /* Parameters controlling when KF's are taken */
    float setting_keyframesPerSecond = 0;   // if !=0, takes a fixed number of KF per second.
    bool setting_realTimeMaxKF = false;   // if true, takes as many KF's as possible (will break the system if the camera stays stationary)
    float setting_maxShiftWeightT= 0.04f * (640+480);
    float setting_maxShiftWeightR= 0.0f * (640+480);
    float setting_maxShiftWeightRT= 0.02f * (640+480);
    float setting_kfGlobalWeight = 1;   // general weight on threshold, the larger the more KF's are taken (e.g., 2 = double the amount of KF's).
    float setting_maxAffineWeight= 2;


    /* initial hessian values to fix unobservable dimensions / priors on affine lighting parameters.
    */
    float setting_idepthFixPrior = 50*50;// * 1000;
    float setting_idepthFixPriorMargFac = 600*600; // 30000*30000;
    float setting_initialRotPrior = 1e11; // 5e7;// 1e11;
    float setting_initialTransPrior = 1e10;// 1e10;
    float setting_initialAffBPrior = 1e14;
    float setting_initialAffAPrior = 1e14;
    float setting_initialCalibHessian = 5e9;



    /* some modes for solving the resulting linear system (e.g. orthogonalize wrt. unobservable dimensions) */
    //int setting_solverMode = SOLVER_FIX_LAMBDA | SOLVER_ORTHOGONALIZE_X_LATER;
    int setting_solverMode = SOLVER_ORTHOGONALIZE_X_LATER;
    double setting_solverModeDelta = 0.00001;
    bool setting_forceAceptStep = false;



    /* some thresholds on when to activate / marginalize points */
    float setting_minIdepthH_act = 100;
    float setting_minIdepthH_marg = 50;


    float setting_desiredImmatureDensity = 1500; // immature points per frame
    float setting_desiredPointDensity = 2000; // aimed total points in the active window.
    float setting_minPointsRemaining = 0.05;  // marg a frame if less than X% points remain.
    float setting_maxLogAffFacInWindow = 0.7; // marg a frame if factor between intensities to current frame is larger than 1/X or X.


    int   setting_minFrames = 5; // min frames in window.
    int   setting_maxFrames = 7; // max frames in window.
    int   setting_minFrameAge = 1;
    int   setting_maxOptIterations=6; // max GN iterations.
    int   setting_minOptIterations=1; // min GN iterations.
    float setting_thOptIterations=1.2; // factor on break threshold for GN iteration (larger = break earlier)



    /* Outlier Threshold on photometric energy */
    float setting_outlierTH = 12*12;					// higher -> less strict
    float setting_outlierTHSumComponent = 50*50; 		// higher -> less strong gradient-based reweighting .



    float setting_margWeightFac = 0.5*0.5;          // factor on hessian when marginalizing, to account for inaccurate linearization points.


    /* when to re-track a frame */
    float setting_reTrackThreshold = 1.5; // (larger = re-track more often)


    /* require some minimum number of residuals for a point to become valid */
    int   setting_minGoodActiveResForMarg=3;
    int   setting_minGoodResForMarg=4;



    // 0 = nothing.
    // 1 = apply inv. response.
    // 2 = apply inv. response & remove V.
    int setting_photometricCalibration = 2;
    bool setting_useExposure = true;
    float setting_affineOptModeA = 1e12; //-1: fix. >=0: optimize (with prior, if > 0).
    float setting_affineOptModeB = 1e8; //-1: fix. >=0: optimize (with prior, if > 0).
    float setting_affineOptModeA_huberTH = 10000;
    float setting_affineOptModeB_huberTH = 10000;
    int setting_gammaWeightsPixelSelect = 1; // 1 = use original intensity for pixel selection; 0 = use gamma-corrected intensity.



    float setting_huberTH = 9; // Huber Threshold



    // parameters controlling adaptive energy threshold computation.
    float setting_frameEnergyTHConstWeight = 0.5;
    float setting_frameEnergyTHN = 0.7f;
    float setting_frameEnergyTHFacMedian = 1.5;
    float setting_overallEnergyTHWeight = 1;
    float setting_coarseCutoffTH = 20;



    // parameters controlling pixel selection
    float setting_minGradHistCut = 0.5;
    float setting_minGradHistAdd = 0.005;
    float setting_gradDownweightPerLevel = 0.75;
    bool  setting_selectDirectionDistribution = true;



    /* settings controling initial immature point tracking */
    float setting_maxPixSearch = 0.027; // max length of the ep. line segment searched during immature point tracking. relative to image resolution.
    float setting_minTraceQuality = 3;
    int setting_minTraceTestRadius = 2;
    int setting_GNItsOnPointActivation = 3;
    float setting_trace_stepsize = 1.0;				// stepsize for initial discrete search.
    int setting_trace_GNIterations = 3;				// max # GN iterations
    float setting_trace_GNThreshold = 0.1;				// GN stop after this stepsize.
    float setting_trace_extraSlackOnTH = 1.2;			// for energy-based outlier check, be slightly more relaxed by this factor.
    float setting_trace_slackInterval = 1.5;			// if pixel-interval is smaller than this, leave it be.
    float setting_trace_minImprovementFactor = 2;		// if pixel-interval is smaller than this, leave it be.



    // for benchmarking different undistortion settings
    float benchmarkSetting_fxfyfac = 0;
    int benchmarkSetting_width = 0;
    int benchmarkSetting_height = 0;
    float benchmark_varNoise = 0;
    float benchmark_varBlurNoise = 0;
    float benchmark_initializerSlackFactor = 1;
    int benchmark_noiseGridsize = 3;



    bool setting_debugSaveImages = false;
    bool settings_no_multiThreading = false;
    bool setting_disableAllDisplay = false;
    bool setting_outputPC = false;
    bool setting_logStuff = true;



    bool setting_goStepByStep = false;


    bool setting_render_displayCoarseTrackingFull=false;
    bool setting_render_displayImmatureTracking=false;
    bool setting_render_renderWindowFrames=true;
    bool setting_render_plotTrackingFull = false;

    bool no_CoarseInit_debugMessage = false; // Coarse Initilizer Debug messages
    bool no_CoarseTracker_debugMessage = false; // Coarse Tracker Debug messages
    bool no_FullSystem_debugMessage = false; // Fullsystem Debug messages 
    bool no_Optimize_debugMessage = false; // Optimizer Debug messages
    bool no_Pixel_debugMessage = false; // Pixel Selector and Immature Point Debug messages
    bool no_GT_debugMessage = false; // GT Distance Debug messages

    bool setting_fullResetRequested = false;

    int setting_sparsityFactor = 5;	// not actually a setting, only some legacy stuff for coarse initializer.

    bool global_Pause = false;

    void handleKey(char k);
};

extern bool setting_debugout_runquiet;
extern bool setting_fullResetRequested;

extern float freeDebugParam1;
extern float freeDebugParam2;
extern float freeDebugParam3;
extern float freeDebugParam4;
extern float freeDebugParam5;


extern const int global_staticPattern[10][40][2];

#define SETTING_PATTERN 8 // point pattern used

#if SETTING_PATTERN == 0
    #define PATTERNNUM 1
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 1
#elif SETTING_PATTERN == 1
    #define PATTERNNUM 5
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 1
#elif SETTING_PATTERN == 2
    #define PATTERNNUM 5
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 1
#elif SETTING_PATTERN == 3
    #define PATTERNNUM 9
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 1
#elif SETTING_PATTERN == 4
    #define PATTERNNUM 9
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 2
#elif SETTING_PATTERN == 5
    #define PATTERNNUM 13
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 2
#elif SETTING_PATTERN == 6
    #define PATTERNNUM 25
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 2
#elif SETTING_PATTERN == 7
    #define PATTERNNUM 21
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 3
#elif SETTING_PATTERN == 8
    #define PATTERNNUM 8
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 2
#elif SETTING_PATTERN == 9
    #define PATTERNNUM 25
    #define PATTERNP global_staticPattern[SETTING_PATTERN]
    #define PATTERNPADDING 4
#else
    #error Unsupported choice setting
#endif

}
