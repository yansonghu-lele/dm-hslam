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



#include "util/globalCalib.h"
#include "stdio.h"
#include <iostream>
#include <limits>


namespace dso
{
	Global_Calib::Global_Calib()
	{
		for (size_t i = 0; i < PYR_LEVELS; i++)
		{
			wG[i] = 0;
			hG[i] = 0;

			fxG[i] = 0;
			fyG[i] = 0;
			cxG[i] = 0;
			cyG[i] = 0;

			fxiG[i] = INFINITY;
			fyiG[i] = INFINITY;
			cxiG[i] = INFINITY;
			cyiG[i] = INFINITY;

			KG[i] = Eigen::Matrix3f::Zero();
			KiG[i] = Eigen::Matrix3f::Zero();
		}
	}

	void Global_Calib::setGlobalCalib(int w, int h,const Eigen::Matrix3f &K, GlobalSettings& globalSettings)
	{
		int wlvl=w;
		int hlvl=h;
		globalSettings.pyrLevelsUsed=1;

		// Determine number of image pyramid levels to use
		while(wlvl%2==0 && hlvl%2==0 && wlvl*hlvl > 5000 && globalSettings.pyrLevelsUsed < PYR_LEVELS)
		{
			wlvl /=2;
			hlvl /=2;
			globalSettings.pyrLevelsUsed++;
		}
		printf("using pyramid levels 0 to %d. coarsest resolution: %d x %d!\n",
				globalSettings.pyrLevelsUsed-1, wlvl, hlvl);

		if(wlvl>100 && hlvl > 100)
		{
			printf("\n\n===============WARNING!===================\n "
					"using not enough pyramid levels.\n"
					"Consider scaling to a resolution that is a multiple of a power of 2.\n");
		}
		if(globalSettings.pyrLevelsUsed < 3)
		{
			printf("\n\n===============WARNING!===================\n "
					"I need higher resolution.\n"
					"I will probably segfault.\n");
		}

		// Set values for each pyramid level
		wG[0] = w;
		hG[0] = h;
		KG[0] = K;
		fxG[0] = K(0,0);
		fyG[0] = K(1,1);
		cxG[0] = K(0,2);
		cyG[0] = K(1,2);
		KiG[0] = KG[0].inverse();
		fxiG[0] = KiG[0](0,0);
		fyiG[0] = KiG[0](1,1);
		cxiG[0] = KiG[0](0,2);
		cyiG[0] = KiG[0](1,2);

		for (int level = 1; level < globalSettings.pyrLevelsUsed; ++ level)
		{
			wG[level] = w >> level;
			hG[level] = h >> level;

			fxG[level] = fxG[level-1] * 0.5;
			fyG[level] = fyG[level-1] * 0.5;
			cxG[level] = (cxG[0] + 0.5) / ((int)1<<level) - 0.5;
			cyG[level] = (cyG[0] + 0.5) / ((int)1<<level) - 0.5;

			KG[level] << fxG[level], 0.0, cxG[level], 
						 0.0, fyG[level], cyG[level], 
						 0.0, 0.0, 1.0;	// synthetic
			KiG[level] = KG[level].inverse();

			fxiG[level] = KiG[level](0,0);
			fyiG[level] = KiG[level](1,1);
			cxiG[level] = KiG[level](0,2);
			cyiG[level] = KiG[level](1,2);
		}
	}

}
