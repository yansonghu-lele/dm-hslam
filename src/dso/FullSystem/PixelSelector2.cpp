/**
* This file is part of DSO, written by Jakob Engel.
* It has been modified by Lukas von Stumberg for the inclusion in DM-VIO (http://vision.in.tum.de/dm-vio).
* It has been modified by Yan Song Hu.
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



#include "FullSystem/PixelSelector2.h"

#include "util/NumType.h"
#include "IOWrapper/ImageDisplay.h"
#include "util/globalCalib.h"
#include "FullSystem/HessianBlocks.h"
#include "util/globalFuncs.h"

# include <stdint.h>


namespace dso
{

/**
 * @brief Construct a new Pixel Selector:: Pixel Selector object
 * 
 * @param globalCalib_ 
 * @param _minGradHistCut 
 * @param _minGradHistAdd 
 * @param globalSettings_ 
 */
PixelSelector::PixelSelector(Global_Calib& globalCalib_, float _minGradHistCut, float _minGradHistAdd, GlobalSettings& globalSettings_):
globalCalib(globalCalib_), globalSettings(globalSettings_)
{
	wG0 = globalCalib.wG[0];
	hG0 = globalCalib.hG[0];
	wG1 = globalCalib.wG[1];
	wG2 = globalCalib.wG[2];
	
	minGradHistCut = _minGradHistCut;
	minGradHistAdd = _minGradHistAdd;
	
	randomPattern = new unsigned char[wG0*hG0];
	std::srand(3141592); // Want to be deterministic
	// X & 0xFF sets value between 0 and 255
	for(int i=0;i<wG0*hG0;i++) randomPattern[i] = rand() & 0xFF;

	// Grid thresholding is used to chose points

	// Block sizes for adaptive threshold grid
	// We create n blocks in width dimension, and adjust the number of blocks for the height accordingly.
    // Always use block size divisable by 16
    bW = 16;
    bH = 16;
    nbW = wG0 / bW;
    nbH = hG0 / bH;

	ths = new float[(nbW)*(nbH)];
	thsSmoothed = new float[(nbW)*(nbH)];

    if(wG0 != bW * nbW || hG0 != bH * nbH)
    {
        std::cout << "ERROR: Height or width seem to be not divisible by 16!" << std::endl;
        assert(0);
    }
    if(!setting_debugout_runquiet && !globalSettings.no_Pixel_debugMessage)
		std::cout << "PixelSelector: Using block sizes: " << bW << ", " << bH << '\n';

	gradHistFrame = 0;

	// For historgram correction
	gradHistInterceptA = (-0.5f)+sqrt(1.0f+4.0f*globalSettings.setting_GradHistCorrect)/2.0f;
	gradHistInterceptB = globalSettings.setting_GradHistCorrect/gradHistInterceptA;

	// Potential is the minimum block size for the selector
	currentPotential=3;
}

/**
 * @brief Destroy the Pixel Selector:: Pixel Selector object
 * 
 */
PixelSelector::~PixelSelector()
{
	delete[] randomPattern;
	delete[] ths;
	delete[] thsSmoothed;
}

/**
 * @brief Helper function to calcuate the approximate quantile to 1/NUM_BINS resolution
 * 
 * @param hist 		Input histogram
 * @param quantil 	Desired quantile
 * @param num_vals 	Total number of inputted values in histogram
 * 
 * @return float 	Approximate value of desired quantile
 */
inline float computeHistQuantile(const std::array <int, NUM_BINS>& hist, const float& quantile, const int& num_vals)
{
	int th = num_vals*quantile;
	for(int i=0;i<NUM_BINS;i++)
	{
		th -= hist[i];
		if(th<0){
			return (static_cast<float>(i)/NUM_BINS);
		}
	}
	return 0.99f;
}

// Fast inverse square root
float Q_rsqrt(float number)
{
	union {
		float    f;
		uint32_t i;
	} conv = { .f = number };
	conv.i  = 0x5f3759df - (conv.i >> 1);
	conv.f *= 1.5F - (number * 0.5F * conv.f * conv.f);
	return conv.f;
}

/**
 * @brief Makes thresold table to help chose points
 * 
 * The table is a grid of gradient quantiles that a gradient must be greater than to be selected as a point
 * A grid of thresholds is used instread of a global threshold to help spread out the points
 * Histogram are used to find the quantiles
 * Using Histograms allows for calculation of approximate quantiles at O(n) speed
 * Calculating actual quantile would require sorting, which is slower than O(n)
 * 
 * @param fh Input frame
 */
void PixelSelector::makeThresTable(const FrameHessian* const fh)
{
	gradHistFrame = fh;
	// absSquaredGrad contains normalized gradient values
	float * mapmax0 = fh->absSquaredGrad[0];

	int w = wG0;
	int h = hG0;

	int w32 = nbW;
	int h32 = nbH;
	thsStep = w32;

	std::array <int, NUM_BINS> gradHist;
	int num_hist_values;

	// For each grid
	for(int y=0;y<h32;y++){
		for(int x=0;x<w32;x++)
		{
			float* map0 = mapmax0+bW*x+bH*y*w;

			std::fill(std::begin(gradHist), std::end(gradHist), 0);

			num_hist_values = 0;
			// Create histogram for specific grid
			// Each histogram bin has bH*bW float intensity values between 0 and 1
			// For each values in a grid
			for(int j=0;j<bH;j++) for(int i=0;i<bW;i++)
			{
				int it = i+bW*x;
				int jt = j+bH*y;

				// Ignore border because gradients can't be calculated at the border properly
				if(it>w-2 || jt>h-2 || it<1 || jt<1) continue;

				float g = map0[i+j*w];
				// Apply pseudo histogram correction
				if(g<0) g = 0;
				g = Q_rsqrt(-globalSettings.setting_GradHistCorrect/(g+gradHistInterceptA)+gradHistInterceptB);
				if(g>0) g = 1.0f/g;
				else g = 0.0f;
				// Update histogram
				int g_int = static_cast<int>(g*NUM_BINS);
				if(g_int >= NUM_BINS) g_int = NUM_BINS-1;
				gradHist[g_int] = gradHist[g_int]+1;
				num_hist_values++;
			}

			if(!setting_debugout_runquiet && !globalSettings.no_Pixel_debugMessage){
				printf("size: %i\n", num_hist_values);
				for (size_t i = 0; i < NUM_BINS; i++)
				{
					printf("%i|", gradHist[i]);
				}
				printf("\n");
			}

			// Calculate approximate threshold
			ths[x+y*w32] = computeHistQuantile(gradHist, minGradHistCut, num_hist_values) + minGradHistAdd;
		}
	}
	

	// Smooth out the quantiles using a box kernel
	// The conditons are used to handle border conditions
	for(int y=0;y<h32;y++)
		for(int x=0;x<w32;x++)
		{
			float sum=0,num=0;
			if(x>0)
			{
				if(y>0) 	{num++; 	sum+=ths[x-1+(y-1)*w32];}
				if(y<h32-1) {num++; 	sum+=ths[x-1+(y+1)*w32];}
				num++; sum+=ths[x-1+(y)*w32];
			}

			if(x<w32-1)
			{
				if(y>0) 	{num++; 	sum+=ths[x+1+(y-1)*w32];}
				if(y<h32-1) {num++; 	sum+=ths[x+1+(y+1)*w32];}
				num++; sum+=ths[x+1+(y)*w32];
			}

			if(y>0) 	{num++; 	sum+=ths[x+(y-1)*w32];}
			if(y<h32-1) {num++; 	sum+=ths[x+(y+1)*w32];}
			num++; sum+=ths[x+y*w32];

			thsSmoothed[x+y*w32] = (sum/num) * (sum/num);
		}

	if(!setting_debugout_runquiet && !globalSettings.no_Pixel_debugMessage){
		printf("Threshold Map\n");
		for (size_t i = 0; i < (nbW)*(nbH); i++)
		{
			if(i%nbW==0) printf("\n");
			printf("%.1f ", thsSmoothed[i]*100);
		}
		printf("\n");
	}
}

/**
 * @brief Generates pixel selection map for tracking
 * 
 * Recursive function
 * Calls itself again if wrong amount of points are selected
 * 
 * @param fh 				Input frame
 * @param map_out 			Output pixels
 * @param density 			Number of points desired
 * @param recursionsLeft 	Number of recursion left
 * @param plot 				Debugging plot
 * @param thFactor 			Multiplier for minimum gradient threshold
 * @return int 				Number of output points
 */
int PixelSelector::makeMaps(
		const FrameHessian* const fh,
		float* map_out, float density, int recursionsLeft, float thFactor)
{
	float numHave = 0;
	float numWant = density;
	float quotia;
	int idealPotential = currentPotential;

	// The number of selected pixels behaves approximately as K / (pot+1)^2, where K is a scene-dependent constant.
	// We will allow sub-selecting pixels by up to a quotia of 0.25, otherwise we will re-select.
	if(fh != gradHistFrame) makeThresTable(fh);

	// Select!
	Eigen::Vector3i n = this->select(fh, map_out,currentPotential, thFactor);

	// Sub-select!
	numHave = n[0]+n[1]+n[2];
	quotia = numWant / numHave;

	// By default we want to over-sample by 40% just to be sure.
	float K = numHave * (currentPotential+1) * (currentPotential+1);
	idealPotential = sqrtf(K/numWant)-1;	// Round down.
	if(idealPotential<1) idealPotential=1;

	// Too few points if quotia > 1.25 and currentPotential>1
	if(recursionsLeft > 0 && quotia > 1.25 && currentPotential>1)
	{
		// Re-sample to get more points!
		// Potential needs to be smaller
		if(idealPotential>=currentPotential)
			idealPotential = currentPotential-1;

	if(!setting_debugout_runquiet && !globalSettings.no_Pixel_debugMessage)
		printf("PixelSelector: have %.2f%%, need %.2f%%. RESAMPLE with pot %d -> %d.\n",
				100*numHave/(float)(wG0*hG0),
				100*numWant/(float)(wG0*hG0),
				currentPotential,
				idealPotential);

		currentPotential = idealPotential;
		return makeMaps(fh,map_out, density, recursionsLeft-1, thFactor);
	}
	// Too many points if quotia < 0.25
	else if(recursionsLeft>0 && quotia < 0.25)
	{
		// Re-sample to get less points!
		// Potential needs to be greater
		if(idealPotential<=currentPotential)
			idealPotential = currentPotential+1;

	if(!setting_debugout_runquiet && !globalSettings.no_Pixel_debugMessage)
		printf("PixelSelector: have %.2f%%, need %.2f%%. RESAMPLE with pot %d -> %d.\n",
				100*numHave/(float)(wG0*hG0),
				100*numWant/(float)(wG0*hG0),
				currentPotential,
				idealPotential);

		currentPotential = idealPotential;
		return makeMaps(fh,map_out, density, recursionsLeft-1, thFactor);
	}

	// Reduce number of points randomly if number of points 5% over quotia
	int numHaveSub = numHave;
	if(quotia < 0.95)
	{
		int wh = wG0*hG0;
		int rn = 0;
		unsigned char charTH = 255*quotia;
		for(int i=0;i<wh;i++)
		{
			if(map_out[i] != 0)
			{
				if(randomPattern[rn] > charTH)
				{
					map_out[i]=0;
					numHaveSub--;
				}
				rn++;
			}
		}
	}

	if(!setting_debugout_runquiet && !globalSettings.no_Pixel_debugMessage)
		printf("PixelSelector: have %.2f%%, need %.2f%%. KEEPCURR with pot %d -> %d. Subsampled to %.2f%%\n",
				100*numHave/(float)(wG0*hG0),
				100*numWant/(float)(wG0*hG0),
				currentPotential,
				idealPotential,
				100*numHaveSub/(float)(wG0*hG0));

	currentPotential = idealPotential;

#ifdef GRAPHICAL_DEBUG
	int w = wG0;
	int h = hG0;

	if (globalSettings.setting_render_displayImmatureTracking){
		MinimalImageB3 img(w,h);

		for(int i=0;i<w*h;i++)
		{
			float c = fh->dI[i][0]*0.7;
			if(c>255) c=255;
			img.at(i) = Vec3b(c,c,c);
		}
		if(!globalSettings.setting_disableAllDisplay) IOWrap::displayImage("Selector Image", &img);

		for(int y=0; y<h;y++)
			for(int x=0;x<w;x++)
			{
				int i=x+y*w;
				if(map_out[i] == 1)
					img.setPixelCirc(x,y,Vec3b(0,255,0));
				else if(map_out[i] == 2)
					img.setPixelCirc(x,y,Vec3b(255,0,0));
				else if(map_out[i] == 4)
					img.setPixelCirc(x,y,Vec3b(0,0,255));
			}
		if(!globalSettings.setting_disableAllDisplay) IOWrap::displayImage("Selector Pixels", &img);
	}
#endif

	return numHaveSub;
}

/**
 * @brief Selects pixels for tracking
 * 
 * Algorithm attempts to chose pixels that are both spread out and easily trackable
 * 
 * @param fh 				Input frame
 * @param map_out 			Output points
 * @param pot 				Potential is relative to number of points to chose
 * @param thFactor 			Minimum gradient threshold multiplier
 * @return Eigen::Vector3i 	Number of chosen points in each pyramid level
 */
Eigen::Vector3i PixelSelector::select(const FrameHessian* const fh,
		float* map_out, int pot, float thFactor)
{
	Eigen::Vector3f const * const map0 = fh->dI;

	// Normalized gradient values between 0 and 1
	// Gradients are queried in three different pyramid levels
	float * mapmax0 = fh->absSquaredGrad[0]; // 1st pyramid level
	float * mapmax1 = fh->absSquaredGrad[1]; // 2nd pyramid level
	float * mapmax2 = fh->absSquaredGrad[2]; // 3rd pyramid level

	int w = wG0;
	int w1 = wG1;
	int w2 = wG2;
	int h = hG0;

	// Set of random directions to chose from
	// Gradient is multiplied by random direction to add variance to point selection
	const Vec2f directions[16] = {
	         Vec2f(0,    1.0000),
	         Vec2f(0.3827,    0.9239),
	         Vec2f(0.1951,    0.9808),
	         Vec2f(0.9239,    0.3827),
	         Vec2f(0.7071,    0.7071),
	         Vec2f(0.3827,   -0.9239),
	         Vec2f(0.8315,    0.5556),
	         Vec2f(0.8315,   -0.5556),
	         Vec2f(0.5556,   -0.8315),
	         Vec2f(0.9808,    0.1951),
	         Vec2f(0.9239,   -0.3827),
	         Vec2f(0.7071,   -0.7071),
	         Vec2f(0.5556,    0.8315),
	         Vec2f(0.9808,   -0.1951),
	         Vec2f(1.0000,    0.0000),
	         Vec2f(0.1951,   -0.9808)};

	// Initialize output array
	memset(map_out,0,w*h*sizeof(PixelSelectorStatus));

	// Each pyramid level has lower gradient requirements
	float dw1 = globalSettings.setting_gradDownweightPerLevel;
	float dw2 = dw1*dw1;

	// The algorithm divides the image into grids for each pyramid level
	// The highest gradient for each grid is chosen
	// Once a point is chosen in a grid, pixels cannot be chosen from that grid in higher pyramid levels
	int n3=0, n2=0, n4=0;
	for(int y4=0;y4<h;y4+=(4*pot)) for(int x4=0;x4<w;x4+=(4*pot))
	{
		int my3 = std::min((4*pot), h-y4); int mx3 = std::min((4*pot), w-x4); // Border handling

		int bestIdx4=-1; float bestVal4=0; // Reset pixel chosen flag at 3rd pyramid level
		Vec2f dir4 = directions[randomPattern[n2] & 0xF];
		for(int y3=0;y3<my3;y3+=(2*pot)) for(int x3=0;x3<mx3;x3+=(2*pot))
		{
			int x34 = x3+x4; int y34 = y3+y4; // Select index
			int my2 = std::min((2*pot), h-y34); int mx2 = std::min((2*pot), w-x34); // Border handling
			
			int bestIdx3=-1; float bestVal3=0; // Reset pixel chosen flag at 2nd pyramid level
			Vec2f dir3 = directions[randomPattern[n2] & 0xF];
			for(int y2=0;y2<my2;y2+=pot) for(int x2=0;x2<mx2;x2+=pot)
			{
				int x234 = x2+x34; int y234 = y2+y34; // Select index
				int my1 = std::min(pot, h-y234); int mx1 = std::min(pot, w-x234); // Border handling
				
				int bestIdx2=-1; float bestVal2=0; // Reset pixel chosen flag at 1st pyramid level
				Vec2f dir2 = directions[randomPattern[n2] & 0xF];
				for(int y1=0;y1<my1;y1+=1) for(int x1=0;x1<mx1;x1+=1)
				{
					assert(x1+x234 < w);
					assert(y1+y234 < h);
					
					// Select index
					int idx = x1+x234 + w*(y1+y234);
					int xf = x1+x234;
					int yf = y1+y234;

					// Don't use pixels too close to the border
					if(xf<4 || xf>=w-5 || yf<4 || yf>h-4) continue;

					// Get minimum gradient value required to be chosen
                    float pixelTH0 = thsSmoothed[xf / bW + (yf / bH) * thsStep];
					// Minimum gradient requirements decreases as pyramid level increases
					float pixelTH1 = pixelTH0*dw1;
					float pixelTH2 = pixelTH1*dw2;


					float ag0 = mapmax0[idx];
					ag0 = -globalSettings.setting_GradHistCorrect/(ag0+gradHistInterceptA)+gradHistInterceptB;

					if(ag0 > pixelTH0*thFactor)
					{
						Vec2f ag0d = map0[idx].tail<2>();
						float dirNorm = fabsf((float)(ag0d.dot(dir2)));
						if(!globalSettings.setting_selectDirectionDistribution) dirNorm = ag0;

						if(dirNorm > bestVal2)
						{ bestVal2 = dirNorm; bestIdx2 = idx; bestIdx3 = -2; bestIdx4 = -2;}
					}
					if(bestIdx3==-2) continue; // Invaliate grid for 2nd and 3rd grid levels


					float ag1 = mapmax1[(int)(xf*0.5f+0.25f) + (int)(yf*0.5f+0.25f)*w1];
					ag1 = -globalSettings.setting_GradHistCorrect/(ag1+gradHistInterceptA)+gradHistInterceptB;

					if(ag1 > pixelTH1*thFactor)
					{
						Vec2f ag0d = map0[idx].tail<2>();
						float dirNorm = fabsf((float)(ag0d.dot(dir3)));
						if(!globalSettings.setting_selectDirectionDistribution) dirNorm = ag1;

						if(dirNorm > bestVal3)
						{ bestVal3 = dirNorm; bestIdx3 = idx; bestIdx4 = -2;}
					}
					if(bestIdx4==-2) continue; // Invaliate grid for 3rd grid level


					float ag2 = mapmax2[(int)(xf*0.25f+0.125) + (int)(yf*0.25f+0.125)*w2];
					ag2 = -globalSettings.setting_GradHistCorrect/(ag2+gradHistInterceptA)+gradHistInterceptB;

					if(ag2 > pixelTH2*thFactor)
					{
						Vec2f ag0d = map0[idx].tail<2>();
						float dirNorm = fabsf((float)(ag0d.dot(dir4)));
						if(!globalSettings.setting_selectDirectionDistribution) dirNorm = ag2;

						if(dirNorm > bestVal4)
						{ bestVal4 = dirNorm; bestIdx4 = idx; }
					}
				}

				if(bestIdx2>0)
				{
					map_out[bestIdx2] = 1;
					bestVal3 = 1e10;
					n2++;
				}
			}

			if(bestIdx3>0)
			{
				map_out[bestIdx3] = 2;
				bestVal4 = 1e10;
				n3++;
			}
		}

		if(bestIdx4>0)
		{
			map_out[bestIdx4] = 4;
			n4++;
		}
	}

	return Eigen::Vector3i(n2,n3,n4);
}

}

