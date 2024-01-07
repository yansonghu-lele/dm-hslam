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


#pragma once

 
#include "util/NumType.h"

namespace dso
{
struct RawResidualJacobian
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	// ================== new structure: save independently =============
	VecNRf resF;

	// ================== Geometric Derivatives =============
	// the two rows of d[x,y]/d[xi]
	// d[x]/d[d_x], d[x]/d[d_y], d[x]/d[d_z], d[x]/d[w_i], d[x]/d[w_j], d[x]/d[w_k]
	// d[y]/d[d_x], d[y]/d[d_y], d[y]/d[d_z], d[y]/d[w_i], d[y]/d[w_j], d[y]/d[w_k]
	Vec6f Jpdxi[2];			// 2x6

	// the two rows of d[x,y]/d[C]
	// d[x]/d[f_x], d[x]/d[f_y], d[x]/d[c_x], d[x]/d[c_y]
	// d[y]/d[f_x], d[y]/d[f_y], d[y]/d[c_x], d[y]/d[c_y]
	VecCf Jpdc[2];			// 2x4

	// the two rows of d[x,y]/d[idepth]
	// d[x]/d[idepth]
	// d[y]/d[idepth]
	Vec2f Jpdd;				// 2x1

	// ================== Image Derivatives =============
	// the two columns of d[r]/d[x,y]
	// x and y pixel derivatives of the pattern points
	VecNRf JIdx[2];			// # pattern points x 2

	// ================== Photometric Derivatives =============
	// the two columns of d[r]/d[ab]
	// a and b derivatives of the pattern points
	VecNRf JabF[2];			// # pattern points x 2


	// = JIdx^T * JIdx (inner product). Only as a shorthand.
	Mat22f JIdx2;			// 2x2
	// = Jab^T * JIdx (inner product). Only as a shorthand.
	Mat22f JabJIdx;			// 2x2
	// = Jab^T * Jab (inner product). Only as a shorthand.
	Mat22f Jab2;			// 2x2

};
}

