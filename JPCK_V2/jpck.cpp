#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <filesystem>

// TO DO
// Optimize cases when lambda_d = 1. Is it possible?
// Compute static filters for all predicates
// Rewrite predicate equations, with an FMA-enabled version (use Herbie) - tried. No advantages or even slower
//

using namespace std;

string createHeadingComment()
{
	string s;
	s += "/****************************************************************************\n";
	s += "* Indirect predicates for geometric constructions					        *\n";
	s += "*                                                                           *\n";
	s += "* Consiglio Nazionale delle Ricerche                                        *\n";
	s += "* Istituto di Matematica Applicata e Tecnologie Informatiche                *\n";
	s += "* Sezione di Genova                                                         * \n";
	s += "* IMATI-GE / CNR                                                            * \n";
	s += "*                                                                           *\n";
	s += "* Authors: Marco Attene                                                     * \n";
	s += "* Copyright(C) 2019: IMATI-GE / CNR                                         * \n";
	s += "* All rights reserved.                                                      * \n";
	s += "*                                                                           *\n";
	s += "* This program is free software; you can redistribute it and/or modify      * \n";
	s += "* it under the terms of the GNU Lesser General Public License as published  * \n";
	s += "* by the Free Software Foundation; either version 3 of the License, or (at  * \n";
	s += "* your option) any later version.                                           * \n";
	s += "*                                                                           *\n";
	s += "* This program is distributed in the hope that it will be useful, but       * \n";
	s += "* WITHOUT ANY WARRANTY; without even the implied warranty of                * \n";
	s += "* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  * \n";
	s += "* General Public License for more details.                                  * \n";
	s += "*                                                                           *\n";
	s += "* You should have received a copy of the GNU Lesser General Public License  * \n";
	s += "* along with this program.  If not, see http://www.gnu.org/licenses/.       *\n";
	s += "*                                                                           *\n";
	s += "****************************************************************************/ \n";

	s += "\n/* This code was generated automatically. Do not edit unless you exactly   */\n";
	s += "/* know what you are doing!                                                */\n\n";

	return s;
}

//#define GEN_STATISTICS

#ifdef GEN_STATISTICS
const char stats_header1[] = "\
enum Stat_Type {\n\
	dotProductSign2D_num_double,\n\
	dotProductSign3D_num_double,\n\
	incircle_num_double,\n\
	inGabrielSphere_num_double,\n\
	inSphere_num_double,\n\
	dotProductSign2D_EEI_num_double,\n\
	dotProductSign2D_IEE_num_double,\n\
	dotProductSign2D_IEI_num_double,\n\
	dotProductSign2D_IIE_num_double,\n\
	dotProductSign2D_III_num_double,\n\
	dotProductSign3D_EEI_num_double,\n\
	dotProductSign3D_IEE_num_double,\n\
	dotProductSign3D_IEI_num_double,\n\
	dotProductSign3D_IIE_num_double,\n\
	dotProductSign3D_III_num_double,\n\
	incirclexy_indirect_IEEE_num_double,\n\
	incirclexy_indirect_IIEE_num_double,\n\
	incirclexy_indirect_IIIE_num_double,\n\
	incirclexy_indirect_IIII_num_double,\n\
	incircle_indirect_IEEE_num_double,\n\
	incircle_indirect_IIEE_num_double,\n\
	incircle_indirect_IIIE_num_double,\n\
	incircle_indirect_IIII_num_double,\n\
	inGabrielSphere_EIEE_num_double,\n\
	inGabrielSphere_EIIE_num_double,\n\
	inGabrielSphere_EIII_num_double,\n\
	inGabrielSphere_IEEE_num_double,\n\
	inGabrielSphere_IIEE_num_double,\n\
	inGabrielSphere_IIIE_num_double,\n\
	inGabrielSphere_IIII_num_double,\n\
	inSphere_IEEEE_num_double,\n\
	inSphere_IIEEE_num_double,\n\
	inSphere_IIIEE_num_double,\n\
	inSphere_IIIIE_num_double,\n\
	inSphere_IIIII_num_double,\n\
	lessThanOnX_IE_num_double,\n\
	lessThanOnX_II_num_double,\n\
	lessThanOnY_IE_num_double,\n\
	lessThanOnY_II_num_double,\n\
	lessThanOnZ_IE_num_double,\n\
	lessThanOnZ_II_num_double,\n\
	orient2dxy_indirect_IEE_num_double,\n\
	orient2dxy_indirect_IIE_num_double,\n\
	orient2dxy_indirect_III_num_double,\n\
	orient2dyz_indirect_IEE_num_double,\n\
	orient2dyz_indirect_IIE_num_double,\n\
	orient2dyz_indirect_III_num_double,\n\
	orient2dzx_indirect_IEE_num_double,\n\
	orient2dzx_indirect_IIE_num_double,\n\
	orient2dzx_indirect_III_num_double,\n\
	orient2d_indirect_IEE_num_double,\n\
	orient2d_indirect_IIE_num_double,\n\
	orient2d_indirect_III_num_double,\n\
	orient3d_indirect_IEEE_num_double,\n\
	orient3d_indirect_IIEE_num_double,\n\
	orient3d_indirect_IIIE_num_double,\n\
	orient3d_indirect_IIII_num_double,\n\
	dotProductSign2D_num_interval,\n\
	dotProductSign3D_num_interval,\n\
	incircle_num_interval,\n\
	inGabrielSphere_num_interval,\n\
	inSphere_num_interval,\n\
	dotProductSign2D_EEI_num_interval,\n\
	dotProductSign2D_IEE_num_interval,\n\
	dotProductSign2D_IEI_num_interval,\n\
	dotProductSign2D_IIE_num_interval,\n\
	dotProductSign2D_III_num_interval,\n\
	dotProductSign3D_EEI_num_interval,\n\
	dotProductSign3D_IEE_num_interval,\n\
	dotProductSign3D_IEI_num_interval,\n\
	dotProductSign3D_IIE_num_interval,\n\
	dotProductSign3D_III_num_interval,\n\
	incirclexy_indirect_IEEE_num_interval,\n\
	incirclexy_indirect_IIEE_num_interval,\n\
	incirclexy_indirect_IIIE_num_interval,\n\
	incirclexy_indirect_IIII_num_interval,\n\
	incircle_indirect_IEEE_num_interval,\n\
	incircle_indirect_IIEE_num_interval,\n\
	incircle_indirect_IIIE_num_interval,\n\
	incircle_indirect_IIII_num_interval,\n\
	inGabrielSphere_EIEE_num_interval,\n\
	inGabrielSphere_EIIE_num_interval,\n\
	inGabrielSphere_EIII_num_interval,\n\
	inGabrielSphere_IEEE_num_interval,\n\
	inGabrielSphere_IIEE_num_interval,\n\
	inGabrielSphere_IIIE_num_interval,\n\
	inGabrielSphere_IIII_num_interval,\n\
	inSphere_IEEEE_num_interval,\n\
	inSphere_IIEEE_num_interval,\n\
	inSphere_IIIEE_num_interval,\n\
	inSphere_IIIIE_num_interval,\n\
	inSphere_IIIII_num_interval,\n\
	lessThanOnX_IE_num_interval,\n\
	lessThanOnX_II_num_interval,\n\
	lessThanOnY_IE_num_interval,\n\
	lessThanOnY_II_num_interval,\n\
	lessThanOnZ_IE_num_interval,\n\
	lessThanOnZ_II_num_interval,\n\
	orient2dxy_indirect_IEE_num_interval,\n\
	orient2dxy_indirect_IIE_num_interval,\n\
	orient2dxy_indirect_III_num_interval,\n\
	orient2dyz_indirect_IEE_num_interval,\n\
	orient2dyz_indirect_IIE_num_interval,\n\
	orient2dyz_indirect_III_num_interval,\n\
	orient2dzx_indirect_IEE_num_interval,\n\
	orient2dzx_indirect_IIE_num_interval,\n\
	orient2dzx_indirect_III_num_interval,\n\
	orient2d_indirect_IEE_num_interval,\n\
	orient2d_indirect_IIE_num_interval,\n\
	orient2d_indirect_III_num_interval,\n\
	orient3d_indirect_IEEE_num_interval,\n\
	orient3d_indirect_IIEE_num_interval,\n\
	orient3d_indirect_IIIE_num_interval,\n\
	orient3d_indirect_IIII_num_interval,\n\
	dotProductSign2D_num_expansion,\n\
	dotProductSign3D_num_expansion,\n\
	incircle_num_expansion,\n\
	inGabrielSphere_num_expansion,\n\
	inSphere_num_expansion,\n\
	dotProductSign2D_EEI_num_expansion,\n\
	dotProductSign2D_IEE_num_expansion,\n\
	dotProductSign2D_IEI_num_expansion,\n\
	dotProductSign2D_IIE_num_expansion,\n\
	dotProductSign2D_III_num_expansion,\n\
	dotProductSign3D_EEI_num_expansion,\n\
	dotProductSign3D_IEE_num_expansion,\n\
	dotProductSign3D_IEI_num_expansion,\n\
	dotProductSign3D_IIE_num_expansion,\n\
	dotProductSign3D_III_num_expansion,\n\
	incirclexy_indirect_IEEE_num_expansion,\n\
	incirclexy_indirect_IIEE_num_expansion,\n\
	incirclexy_indirect_IIIE_num_expansion,\n\
	incirclexy_indirect_IIII_num_expansion,\n\
	incircle_indirect_IEEE_num_expansion,\n\
	incircle_indirect_IIEE_num_expansion,\n\
	incircle_indirect_IIIE_num_expansion,\n\
	incircle_indirect_IIII_num_expansion,\n\
	inGabrielSphere_EIEE_num_expansion,\n\
	inGabrielSphere_EIIE_num_expansion,\n\
	inGabrielSphere_EIII_num_expansion,\n\
	inGabrielSphere_IEEE_num_expansion,\n\
	inGabrielSphere_IIEE_num_expansion,\n\
	inGabrielSphere_IIIE_num_expansion,\n\
	inGabrielSphere_IIII_num_expansion,\n\
	inSphere_IEEEE_num_expansion,\n\
	inSphere_IIEEE_num_expansion,\n\
	inSphere_IIIEE_num_expansion,\n\
	inSphere_IIIIE_num_expansion,\n\
	inSphere_IIIII_num_expansion,\n\
	lessThanOnX_IE_num_expansion,\n\
	lessThanOnX_II_num_expansion,\n\
	lessThanOnY_IE_num_expansion,\n\
	lessThanOnY_II_num_expansion,\n\
	lessThanOnZ_IE_num_expansion,\n\
	lessThanOnZ_II_num_expansion,\n\
	orient2dxy_indirect_IEE_num_expansion,\n\
	orient2dxy_indirect_IIE_num_expansion,\n\
	orient2dxy_indirect_III_num_expansion,\n\
	orient2dyz_indirect_IEE_num_expansion,\n\
	orient2dyz_indirect_IIE_num_expansion,\n\
	orient2dyz_indirect_III_num_expansion,\n\
	orient2dzx_indirect_IEE_num_expansion,\n\
	orient2dzx_indirect_IIE_num_expansion,\n\
	orient2dzx_indirect_III_num_expansion,\n\
	orient2d_indirect_IEE_num_expansion,\n\
	orient2d_indirect_IIE_num_expansion,\n\
	orient2d_indirect_III_num_expansion,\n\
	orient3d_indirect_IEEE_num_expansion,\n\
	orient3d_indirect_IIEE_num_expansion,\n\
	orient3d_indirect_IIIE_num_expansion,\n\
	orient3d_indirect_IIII_num_expansion,\n\
	dotProductSign2D_num_bigfloat,\n\
	dotProductSign3D_num_bigfloat,\n\
	incircle_num_bigfloat,\n\
	inGabrielSphere_num_bigfloat,\n\
	inSphere_num_bigfloat,\n\
	dotProductSign2D_EEI_num_bigfloat,\n\
	dotProductSign2D_IEE_num_bigfloat,\n\
	dotProductSign2D_IEI_num_bigfloat,\n\
	dotProductSign2D_IIE_num_bigfloat,\n\
	dotProductSign2D_III_num_bigfloat,\n\
	dotProductSign3D_EEI_num_bigfloat,\n\
	dotProductSign3D_IEE_num_bigfloat,\n\
	dotProductSign3D_IEI_num_bigfloat,\n\
	dotProductSign3D_IIE_num_bigfloat,\n\
	dotProductSign3D_III_num_bigfloat,\n\
	incirclexy_indirect_IEEE_num_bigfloat,\n\
	incirclexy_indirect_IIEE_num_bigfloat,\n\
	incirclexy_indirect_IIIE_num_bigfloat,\n\
	incirclexy_indirect_IIII_num_bigfloat,\n\
	incircle_indirect_IEEE_num_bigfloat,\n\
	incircle_indirect_IIEE_num_bigfloat,\n\
	incircle_indirect_IIIE_num_bigfloat,\n\
	incircle_indirect_IIII_num_bigfloat,\n\
	inGabrielSphere_EIEE_num_bigfloat,\n\
	inGabrielSphere_EIIE_num_bigfloat,\n\
	inGabrielSphere_EIII_num_bigfloat,\n\
	inGabrielSphere_IEEE_num_bigfloat,\n\
	inGabrielSphere_IIEE_num_bigfloat,\n\
	inGabrielSphere_IIIE_num_bigfloat,\n\
	inGabrielSphere_IIII_num_bigfloat,\n\
	inSphere_IEEEE_num_bigfloat,\n\
	inSphere_IIEEE_num_bigfloat,\n\
	inSphere_IIIEE_num_bigfloat,\n\
	inSphere_IIIIE_num_bigfloat,\n\
	inSphere_IIIII_num_bigfloat,\n\
	lessThanOnX_IE_num_bigfloat,\n\
	lessThanOnX_II_num_bigfloat,\n\
	lessThanOnY_IE_num_bigfloat,\n\
	lessThanOnY_II_num_bigfloat,\n\
	lessThanOnZ_IE_num_bigfloat,\n\
	lessThanOnZ_II_num_bigfloat,\n\
	orient2dxy_indirect_IEE_num_bigfloat,\n\
	orient2dxy_indirect_IIE_num_bigfloat,\n\
	orient2dxy_indirect_III_num_bigfloat,\n\
	orient2dyz_indirect_IEE_num_bigfloat,\n\
	orient2dyz_indirect_IIE_num_bigfloat,\n\
	orient2dyz_indirect_III_num_bigfloat,\n\
	orient2dzx_indirect_IEE_num_bigfloat,\n\
	orient2dzx_indirect_IIE_num_bigfloat,\n\
	orient2dzx_indirect_III_num_bigfloat,\n\
	orient2d_indirect_IEE_num_bigfloat,\n\
	orient2d_indirect_IIE_num_bigfloat,\n\
	orient2d_indirect_III_num_bigfloat,\n\
	orient3d_indirect_IEEE_num_bigfloat,\n\
	orient3d_indirect_IIEE_num_bigfloat,\n\
	orient3d_indirect_IIIE_num_bigfloat,\n\
	orient3d_indirect_IIII_num_bigfloat,\n\
	dotProductSign2D_time_double,\n\
	dotProductSign3D_time_double,\n\
	incircle_time_double,\n\
	inGabrielSphere_time_double,\n\
	inSphere_time_double,\n\
	dotProductSign2D_EEI_time_double,\n\
	dotProductSign2D_IEE_time_double,\n\
	dotProductSign2D_IEI_time_double,\n\
	dotProductSign2D_IIE_time_double,\n\
	dotProductSign2D_III_time_double,\n\
	dotProductSign3D_EEI_time_double,\n\
	dotProductSign3D_IEE_time_double,\n\
	dotProductSign3D_IEI_time_double,\n\
	dotProductSign3D_IIE_time_double,\n\
	dotProductSign3D_III_time_double,\n\
	incirclexy_indirect_IEEE_time_double,\n\
	incirclexy_indirect_IIEE_time_double,\n\
	incirclexy_indirect_IIIE_time_double,\n\
	incirclexy_indirect_IIII_time_double,\n\
	incircle_indirect_IEEE_time_double,\n\
	incircle_indirect_IIEE_time_double,\n\
	incircle_indirect_IIIE_time_double,\n\
	incircle_indirect_IIII_time_double,\n\
	inGabrielSphere_EIEE_time_double,\n\
	inGabrielSphere_EIIE_time_double,\n\
	inGabrielSphere_EIII_time_double,\n\
	inGabrielSphere_IEEE_time_double,\n\
	inGabrielSphere_IIEE_time_double,\n\
	inGabrielSphere_IIIE_time_double,\n\
	inGabrielSphere_IIII_time_double,\n\
	inSphere_IEEEE_time_double,\n\
	inSphere_IIEEE_time_double,\n\
	inSphere_IIIEE_time_double,\n\
	inSphere_IIIIE_time_double,\n\
	inSphere_IIIII_time_double,\n\
	lessThanOnX_IE_time_double,\n\
	lessThanOnX_II_time_double,\n\
	lessThanOnY_IE_time_double,\n\
	lessThanOnY_II_time_double,\n\
	lessThanOnZ_IE_time_double,\n\
	lessThanOnZ_II_time_double,\n\
	orient2dxy_indirect_IEE_time_double,\n\
	orient2dxy_indirect_IIE_time_double,\n\
	orient2dxy_indirect_III_time_double,\n\
	orient2dyz_indirect_IEE_time_double,\n\
	orient2dyz_indirect_IIE_time_double,\n\
	orient2dyz_indirect_III_time_double,\n\
	orient2dzx_indirect_IEE_time_double,\n\
	orient2dzx_indirect_IIE_time_double,\n\
	orient2dzx_indirect_III_time_double,\n\
	orient2d_indirect_IEE_time_double,\n\
	orient2d_indirect_IIE_time_double,\n\
	orient2d_indirect_III_time_double,\n\
	orient3d_indirect_IEEE_time_double,\n\
	orient3d_indirect_IIEE_time_double,\n\
	orient3d_indirect_IIIE_time_double,\n\
	orient3d_indirect_IIII_time_double,\n\
	dotProductSign2D_time_interval,\n\
	dotProductSign3D_time_interval,\n\
	incircle_time_interval,\n\
	inGabrielSphere_time_interval,\n\
	inSphere_time_interval,\n\
	dotProductSign2D_EEI_time_interval,\n\
	dotProductSign2D_IEE_time_interval,\n\
	dotProductSign2D_IEI_time_interval,\n\
	dotProductSign2D_IIE_time_interval,\n\
	dotProductSign2D_III_time_interval,\n\
	dotProductSign3D_EEI_time_interval,\n\
	dotProductSign3D_IEE_time_interval,\n\
	dotProductSign3D_IEI_time_interval,\n\
	dotProductSign3D_IIE_time_interval,\n\
	dotProductSign3D_III_time_interval,\n\
	incirclexy_indirect_IEEE_time_interval,\n\
	incirclexy_indirect_IIEE_time_interval,\n\
	incirclexy_indirect_IIIE_time_interval,\n\
	incirclexy_indirect_IIII_time_interval,\n\
	incircle_indirect_IEEE_time_interval,\n\
	incircle_indirect_IIEE_time_interval,\n\
	incircle_indirect_IIIE_time_interval,\n\
	incircle_indirect_IIII_time_interval,\n\
	inGabrielSphere_EIEE_time_interval,\n\
	inGabrielSphere_EIIE_time_interval,\n\
	inGabrielSphere_EIII_time_interval,\n\
	inGabrielSphere_IEEE_time_interval,\n\
	inGabrielSphere_IIEE_time_interval,\n\
	inGabrielSphere_IIIE_time_interval,\n\
	inGabrielSphere_IIII_time_interval,\n\
	inSphere_IEEEE_time_interval,\n\
	inSphere_IIEEE_time_interval,\n\
	inSphere_IIIEE_time_interval,\n\
	inSphere_IIIIE_time_interval,\n\
	inSphere_IIIII_time_interval,\n\
	lessThanOnX_IE_time_interval,\n\
	lessThanOnX_II_time_interval,\n\
	lessThanOnY_IE_time_interval,\n\
	lessThanOnY_II_time_interval,\n\
	lessThanOnZ_IE_time_interval,\n\
	lessThanOnZ_II_time_interval,\n\
	orient2dxy_indirect_IEE_time_interval,\n\
	orient2dxy_indirect_IIE_time_interval,\n\
	orient2dxy_indirect_III_time_interval,\n\
	orient2dyz_indirect_IEE_time_interval,\n\
	orient2dyz_indirect_IIE_time_interval,\n\
	orient2dyz_indirect_III_time_interval,\n\
	orient2dzx_indirect_IEE_time_interval,\n\
	orient2dzx_indirect_IIE_time_interval,\n\
	orient2dzx_indirect_III_time_interval,\n\
	orient2d_indirect_IEE_time_interval,\n\
	orient2d_indirect_IIE_time_interval,\n\
	orient2d_indirect_III_time_interval,\n\
	orient3d_indirect_IEEE_time_interval,\n\
	orient3d_indirect_IIEE_time_interval,\n\
	orient3d_indirect_IIIE_time_interval,\n\
	orient3d_indirect_IIII_time_interval,\n\
	dotProductSign2D_time_expansion,\n\
	dotProductSign3D_time_expansion,\n\
	incircle_time_expansion,\n\
	inGabrielSphere_time_expansion,\n\
	inSphere_time_expansion,\n\
	dotProductSign2D_EEI_time_expansion,\n\
	dotProductSign2D_IEE_time_expansion,\n\
	dotProductSign2D_IEI_time_expansion,\n\
	dotProductSign2D_IIE_time_expansion,\n\
	dotProductSign2D_III_time_expansion,\n\
	dotProductSign3D_EEI_time_expansion,\n\
	dotProductSign3D_IEE_time_expansion,\n\
	dotProductSign3D_IEI_time_expansion,\n\
	dotProductSign3D_IIE_time_expansion,\n\
	dotProductSign3D_III_time_expansion,\n\
	incirclexy_indirect_IEEE_time_expansion,\n\
	incirclexy_indirect_IIEE_time_expansion,\n\
	incirclexy_indirect_IIIE_time_expansion,\n\
	incirclexy_indirect_IIII_time_expansion,\n\
	incircle_indirect_IEEE_time_expansion,\n\
	incircle_indirect_IIEE_time_expansion,\n\
	incircle_indirect_IIIE_time_expansion,\n\
	incircle_indirect_IIII_time_expansion,\n\
	inGabrielSphere_EIEE_time_expansion,\n\
	inGabrielSphere_EIIE_time_expansion,\n\
	inGabrielSphere_EIII_time_expansion,\n\
	inGabrielSphere_IEEE_time_expansion,\n\
	inGabrielSphere_IIEE_time_expansion,\n\
	inGabrielSphere_IIIE_time_expansion,\n\
	inGabrielSphere_IIII_time_expansion,\n\
	inSphere_IEEEE_time_expansion,\n\
	inSphere_IIEEE_time_expansion,\n\
	inSphere_IIIEE_time_expansion,\n\
	inSphere_IIIIE_time_expansion,\n\
	inSphere_IIIII_time_expansion,\n\
	lessThanOnX_IE_time_expansion,\n\
	lessThanOnX_II_time_expansion,\n\
	lessThanOnY_IE_time_expansion,\n\
	lessThanOnY_II_time_expansion,\n\
	lessThanOnZ_IE_time_expansion,\n\
	lessThanOnZ_II_time_expansion,\n\
	orient2dxy_indirect_IEE_time_expansion,\n\
	orient2dxy_indirect_IIE_time_expansion,\n\
	orient2dxy_indirect_III_time_expansion,\n\
	orient2dyz_indirect_IEE_time_expansion,\n\
	orient2dyz_indirect_IIE_time_expansion,\n\
	orient2dyz_indirect_III_time_expansion,\n\
	orient2dzx_indirect_IEE_time_expansion,\n\
	orient2dzx_indirect_IIE_time_expansion,\n\
	orient2dzx_indirect_III_time_expansion,\n\
	orient2d_indirect_IEE_time_expansion,\n\
	orient2d_indirect_IIE_time_expansion,\n\
	orient2d_indirect_III_time_expansion,\n\
	orient3d_indirect_IEEE_time_expansion,\n\
	orient3d_indirect_IIEE_time_expansion,\n\
	orient3d_indirect_IIIE_time_expansion,\n\
	orient3d_indirect_IIII_time_expansion,\n\
	dotProductSign2D_time_bigfloat,\n\
	dotProductSign3D_time_bigfloat,\n\
	incircle_time_bigfloat,\n\
	inGabrielSphere_time_bigfloat,\n\
	inSphere_time_bigfloat,\n\
	dotProductSign2D_EEI_time_bigfloat,\n\
	dotProductSign2D_IEE_time_bigfloat,\n\
	dotProductSign2D_IEI_time_bigfloat,\n\
	dotProductSign2D_IIE_time_bigfloat,\n\
	dotProductSign2D_III_time_bigfloat,\n\
	dotProductSign3D_EEI_time_bigfloat,\n\
	dotProductSign3D_IEE_time_bigfloat,\n\
	dotProductSign3D_IEI_time_bigfloat,\n\
	dotProductSign3D_IIE_time_bigfloat,\n\
	dotProductSign3D_III_time_bigfloat,\n";

const char stats_header2[] = "\
	incirclexy_indirect_IEEE_time_bigfloat,\n\
	incirclexy_indirect_IIEE_time_bigfloat,\n\
	incirclexy_indirect_IIIE_time_bigfloat,\n\
	incirclexy_indirect_IIII_time_bigfloat,\n\
	incircle_indirect_IEEE_time_bigfloat,\n\
	incircle_indirect_IIEE_time_bigfloat,\n\
	incircle_indirect_IIIE_time_bigfloat,\n\
	incircle_indirect_IIII_time_bigfloat,\n\
	inGabrielSphere_EIEE_time_bigfloat,\n\
	inGabrielSphere_EIIE_time_bigfloat,\n\
	inGabrielSphere_EIII_time_bigfloat,\n\
	inGabrielSphere_IEEE_time_bigfloat,\n\
	inGabrielSphere_IIEE_time_bigfloat,\n\
	inGabrielSphere_IIIE_time_bigfloat,\n\
	inGabrielSphere_IIII_time_bigfloat,\n\
	inSphere_IEEEE_time_bigfloat,\n\
	inSphere_IIEEE_time_bigfloat,\n\
	inSphere_IIIEE_time_bigfloat,\n\
	inSphere_IIIIE_time_bigfloat,\n\
	inSphere_IIIII_time_bigfloat,\n\
	lessThanOnX_IE_time_bigfloat,\n\
	lessThanOnX_II_time_bigfloat,\n\
	lessThanOnY_IE_time_bigfloat,\n\
	lessThanOnY_II_time_bigfloat,\n\
	lessThanOnZ_IE_time_bigfloat,\n\
	lessThanOnZ_II_time_bigfloat,\n\
	orient2dxy_indirect_IEE_time_bigfloat,\n\
	orient2dxy_indirect_IIE_time_bigfloat,\n\
	orient2dxy_indirect_III_time_bigfloat,\n\
	orient2dyz_indirect_IEE_time_bigfloat,\n\
	orient2dyz_indirect_IIE_time_bigfloat,\n\
	orient2dyz_indirect_III_time_bigfloat,\n\
	orient2dzx_indirect_IEE_time_bigfloat,\n\
	orient2dzx_indirect_IIE_time_bigfloat,\n\
	orient2dzx_indirect_III_time_bigfloat,\n\
	orient2d_indirect_IEE_time_bigfloat,\n\
	orient2d_indirect_IIE_time_bigfloat,\n\
	orient2d_indirect_III_time_bigfloat,\n\
	orient3d_indirect_IEEE_time_bigfloat,\n\
	orient3d_indirect_IIEE_time_bigfloat,\n\
	orient3d_indirect_IIIE_time_bigfloat,\n\
	orient3d_indirect_IIII_time_bigfloat,\n\
	max_expansion_size\n\
};\n\n\
std::vector<uint64_t> Pred_Stats(457, 0);\n\
";
#endif

//
//// Generic implementation of integer-exponent power functions (up to exponent 8)
//const char ipow_functions[] = "\
//template<class T> static inline T ipow2(const T& d) { return d * d; }\n\
//template<class T> static inline T ipow3(const T& d) { return d * ipow2(d); }\n\
//template<class T> static inline T ipow4(const T& d) { return ipow2(ipow2(d)); }\n\
//template<class T> static inline T ipow5(const T& d) { return d * ipow4(d); }\n\
//template<class T> static inline T ipow6(const T& d) { return ipow2(ipow3(d)); }\n\
//template<class T> static inline T ipow7(const T& d) { return d * ipow6(d); }\n\
//template<class T> static inline T ipow8(const T& d) { return ipow2(ipow2(ipow2(d))); }\n\
//";


// Generic implementation of integer-exponent power functions (up to exponent 8)
const char double_ipow_functions[] = "\
static inline double ipow2(const double& d) { return d * d; }\n\
static inline double ipow3(const double& d) { return d * ipow2(d); }\n\
";

const char interval_ipow_functions[] = "\
static inline interval_number ipow2(const interval_number& d) { return d.pow2(); }\n\
static inline interval_number ipow3(const interval_number& d) { return d.pow3(); }\n\
";

const char expansion_ipow_functions[] = "\
static inline expansion ipow2(const expansion& d) { return d.sqr(); }\n\
static inline expansion ipow3(const expansion& d) { return d.sqr()*d; }\n\
";

const char bigfloat_ipow_functions[] = "\
static inline bigfloat ipow2(const bigfloat& d) { return d * d; }\n\
static inline bigfloat ipow3(const bigfloat& d) { return d * ipow2(d); }\n\
";

const char double_sgn_function[] = "inline int sgn(const double p) { return (p>0)-(p<0); }";
const char interval_sgn_function[] = "inline int sgn(const interval_number& p) { return (p.isPositive()) ? (1) : ((p.isNegative()) ? (-1) : (0)); }";

const char double_fma_function[] = "inline double fmadd(const double& a, const double& b, const double& c) {\n\
\treturn _mm_cvtsd_f64(_mm_fmadd_sd(_mm_load_pd(&a), _mm_load_pd(&b), _mm_load_pd(&c)));\n}\n\
inline double fmsub(const double& a, const double& b, const double& c) {\n\
\treturn _mm_cvtsd_f64(_mm_fmsub_sd(_mm_load_pd(&a), _mm_load_pd(&b), _mm_load_pd(&c)));\n}\n";

const char interval_fma_function[] = "inline interval_number fmadd(const interval_number& a, const interval_number& b, const interval_number& c) {\n\
\treturn a.fmadd(b, c);\n}\n\
inline interval_number fmsub(const interval_number& a, const interval_number& b, const interval_number& c) {\n\
\treturn a.fmsub(b, c);\n}\n";

const char bigfloat_fma_function[] = "inline bigfloat fmadd(const bigfloat& a, const bigfloat& b, const bigfloat& c) { return a * b + c; }\n\
inline bigfloat fmsub(const bigfloat& a, const bigfloat& b, const bigfloat& c) { return a * b - c; }\n";

typedef double fpnumber;
#define FPN_EPSILON	DBL_EPSILON

static fpnumber ulp(fpnumber d)
{
	int exp;
	frexp(d, &exp);
	return ldexp(0.5, exp - 52);
}

void error(const char *s)
{
	cerr << s << endl;
	cerr << "Hit enter to exit\n";
	getchar();
	exit(0);
}

int myisspace(char c) { return ::isspace((unsigned char)c); }

// Bynary/ternary tree of operations on basic symbols.
// Operations supported:
// +, -, *, ^, fmadd, fmsub
// 
// Input from string accepts standard mathematical syntax
// with explicit multiplication operators (i.e. the string
// "xy" indicates a symbol called 'xy' and NOT the product
// 'x times y'.
// Multiplications are required even before/after parenthesis
// (e.g. "a(b+c)" is an invalid string. Use "a*(b+c)" instead

class polynomial {
public:
	char operation; // +, -, *, ^, f, s (the latter stand for f=fmadd, s=fmsub)
	string expression; // Non-empty for leaves, indicating a primitive symbol
	polynomial* op1, * op2, *op3; // Operands. op3!=NULL for fma only.

	fpnumber value_bound;		// Forward error analysis: bound on magnitude
	fpnumber error_bound;		// Forward error analysis: bound on error
	int error_degree;		// Forward error analysis: degree (0 means that analysis was not performed)
	bool is_a_max;			// TRUE if magnitude is relevant for error bound

	// A 'translation filter' is the difference of two leaves
	bool isTranslationFilter() const {
		return (operation == '-' && op1->operation==0 && op2->operation == 0);
	}

	// Replace each occurrence of symbol 's' with a new subtree 'p'
	void replaceSymbol(const string& s, const polynomial& p) {
		if (expression == s) {
			*this = p;
		}
		else {
			if (op1 != NULL) op1->replaceSymbol(s, p);
			if (op2 != NULL) op2->replaceSymbol(s, p);
			if (op3 != NULL) op3->replaceSymbol(s, p);
		}
	}

	// Basic constructor from string (see above for syntactical rules)
	polynomial(const string& s) : is_a_max(false), error_degree(0) {
		string e = s;
		e.erase(std::remove_if(e.begin(), e.end(), myisspace), e.end());
		init(e, 0, e.size());
	}

	polynomial(const polynomial& p) : 
		operation(p.operation), expression(p.expression), 
		value_bound(p.value_bound), error_bound(p.error_bound), 
		error_degree(p.error_degree), is_a_max(p.is_a_max) 
	{
		op1 = (p.op1 == NULL) ? (NULL) : (new polynomial(*p.op1));
		op2 = (p.op2 == NULL) ? (NULL) : (new polynomial(*p.op2));
		op3 = (p.op3 == NULL) ? (NULL) : (new polynomial(*p.op3));
	}

	polynomial& operator=(const polynomial& p) {
		operation = p.operation;
		expression = p.expression;
		value_bound = p.value_bound;
		error_bound = p.error_bound;
		error_degree = p.error_degree;
		is_a_max = p.is_a_max;
		op1 = (p.op1 == NULL) ? (NULL) : (new polynomial(*p.op1));
		op2 = (p.op2 == NULL) ? (NULL) : (new polynomial(*p.op2));
		op3 = (p.op3 == NULL) ? (NULL) : (new polynomial(*p.op3));
		return *this;
	}

	polynomial(polynomial&& p) noexcept : 
		op1(p.op1), op2(p.op2), op3(p.op3),
		value_bound(p.value_bound),
		error_bound(p.error_bound),
		error_degree(p.error_degree),
		is_a_max(p.is_a_max)
	{
		operation = p.operation;
		p.op1 = p.op2 = p.op3 = NULL;
		std::swap(expression, p.expression);
	}

	~polynomial() {
		delete op1;
		delete op2;
		delete op3;
	}

	// Converts to string
	string toString() const {
		if (operation == 0) return expression;
		else if (operation == 'f') return "fmadd(" + op1->toString() + "," + op2->toString() + "," + op3->toString() + ")";
		else if (operation == 's') return "fmsub(" + op1->toString() + "," + op2->toString() + "," + op3->toString() + ")";
		else return "(" + op1->toString() + operation + op2->toString() + ")";
	}

	// Converts to code
	string toCode() const {
		if (operation == 0) return expression;
		else if (operation == 'f') return "fmadd(" + op1->toCode() + "," + op2->toCode() + "," + op3->toCode() + ")";
		else if (operation == 's') return "fmsub(" + op1->toCode() + "," + op2->toCode() + "," + op3->toCode() + ")";
		else if (operation == '^') return "ipow" + op2->toCode() + "(" + op1->toCode() + ")";
		else return "(" + op1->toCode() + operation + op2->toCode() + ")";
	}

	// Collects all symbols
	void getSymbols(std::set<string>& symbols) const {
		if (operation != 0) {
			op1->getSymbols(symbols);
			op2->getSymbols(symbols);
			if (op3 != NULL) op3->getSymbols(symbols);
		}
		else if (!is_number(expression)) symbols.insert(expression);
	}

	// Count number of occurrences of 's' out of transaltion filters
	int countNTFSymbols(const string& s) const {
		if (isTranslationFilter()) return 0;
		else if (operation != 0) {
			int r = op1->countNTFSymbols(s);
			r += op2->countNTFSymbols(s);
			if (op3 != NULL) r += op3->countNTFSymbols(s);
			return r;
		}
		else if (s == expression) 
			return 1;
		else 
			return 0;
	}

	fpnumber getValueBound() {
		if (error_degree == 0) getErrorBound();
		return value_bound;
	}

	int getErrorDegree() {
		if (error_degree == 0) getErrorBound();
		return error_degree;
	}

	fpnumber getErrorBound() {
		if (error_degree == 0) {
			if (operation == 0) {
				value_bound = 1;
				error_bound = 0;
				error_degree = 1;
				is_a_max = true;
				return 0;
			}

			if (operation == '+' || operation == '-')
			{
				fpnumber e1 = op1->getErrorBound(), e2 = op2->getErrorBound();
				fpnumber v1 = op1->getValueBound(), v2 = op2->getValueBound();
				if (operation == '-' && e1 == 0 && e2 == 0) // Translation filter
				{
					// This is the same formula used in FPG. I am not sure it is correct.
					// Why is the value bound 1 and not 2? (1 - (-1)) = 2.
					value_bound = 1;
					error_bound = (0.5 * FPN_EPSILON);
					error_degree = 1;
					is_a_max = true;
					op1->is_a_max = op2->is_a_max = false;
				}
				else // Regular sum or subtraction
				{
					value_bound = v1 + v2;
					fpnumber u = 0.5 * ulp(value_bound);
					value_bound += u;
					error_bound = e1 + e2 + u;
					error_degree = op1->getErrorDegree();
					int td = op2->getErrorDegree();
					if (td > error_degree) error_degree = td;
				}
			}
			else if (operation == '*')
			{
				if (op1->expression == "2") {
					value_bound *= 2;
					error_bound = op2->getErrorBound() * 2;
					error_degree = op2->getErrorDegree();
				}
				else if (op2->expression == "2") {
					value_bound *= 2;
					error_bound = op1->getErrorBound() * 2;
					error_degree = op1->getErrorDegree();
				}
				else {
					fpnumber e1 = op1->getErrorBound(), e2 = op2->getErrorBound();
					fpnumber v1 = op1->getValueBound(), v2 = op2->getValueBound();
					value_bound = v1 * v2;
					fpnumber u = 0.5 * ulp(value_bound);
					value_bound += u;
					error_bound = e1 * e2 + e1 * v2 + e2 * v1 + u;
					error_degree = op1->getErrorDegree() + op2->getErrorDegree();
				}
			}
			else if (operation == '^')
			{
				int exponent = std::stoi(op2->expression);
				fpnumber v = op1->getValueBound();
				fpnumber e = op1->getErrorBound();

				if (exponent == 2) {
					value_bound = v * v;
					fpnumber u = 0.5 * ulp(value_bound);
					value_bound += u;
					error_bound = e * e + 2 * e * v + u;
				}
				else if (exponent == 3) {
					value_bound = v * v;
					fpnumber u = 0.5 * ulp(value_bound);
					value_bound = (value_bound + u) * v;
					u = 0.5 * ulp(value_bound);
					value_bound += u;

					error_bound = e * e * e + 3 * (e * e * v + e * v * v) + v * v * v + u;
				}
				else error("Exponents different than two and three are currently unsupported\n");

				error_degree = op1->getErrorDegree() * exponent;
			}
			else if (operation == 'f' || operation == 's')
			{
				fpnumber e1 = op1->getErrorBound(), e2 = op2->getErrorBound(), e3 = op3->getErrorBound();
				fpnumber v1 = op1->getValueBound(), v2 = op2->getValueBound(), v3 = op3->getValueBound();
				value_bound = v1 * v2 + v3;
				fpnumber u = 0.5 * ulp(value_bound);
				value_bound += u;
				error_bound = e1 * e2 + e1 * v2 + e2 * v1 + u + e3;
				error_degree = op1->getErrorDegree() + op2->getErrorDegree();
				int td = op3->getErrorDegree();
				if (td > error_degree) error_degree = td;
			}
		}

		return error_bound;
	}

protected:

	static bool is_number(const string& s) {
		return (!s.empty()) && (std::all_of(s.begin(), s.end(), ::isdigit));
	}

	polynomial(const string& s, size_t begin, size_t end) : is_a_max(false), error_degree(0) {
		if (s.empty()) error("Empty expression\n");
		init(s, begin, end);
	}

	void initChildren(char op, const string& s, size_t begin, size_t pos, size_t end) {
		operation = op;
		op1 = new polynomial(s, begin, pos);
		op2 = new polynomial(s, pos + 1, end);
		op3 = NULL;
	}

	void initFMAChildren(const string& s, size_t begin, size_t pos1, size_t pos2, size_t end, char aors) {
		operation = (aors == 'a') ? ('f') : ('s');
		op1 = new polynomial(s, begin, pos1);
		op2 = new polynomial(s, pos1 + 1, pos2);
		op3 = new polynomial(s, pos2 + 1, end);
	}

	void init(const string& s, size_t begin, size_t end) {
		size_t pos, pos1, pos2, pos3;
		if ((pos = findOperation(s, begin, end, '+')) != end) {
			initChildren('+', s, begin, pos, end);
		}
		else if ((pos = findOperation(s, begin, end, '-')) != end) {
			initChildren('-', s, begin, pos, end);
		}
		else if ((pos = findOperation(s, begin, end, '*')) != end) {
			initChildren('*', s, begin, pos, end);
		}
		else if ((pos = findOperation(s, begin, end, '^')) != end) {
			initChildren('^', s, begin, pos, end);
			if (!is_number(op2->expression)) error("Invalid exponent in power\n");
		}
		else if (findFMA(s, begin, end, pos, pos1, pos2, pos3, 'a')) {
			initFMAChildren(s, pos, pos1, pos2, pos3, 'a');
		}
		else if (findFMA(s, begin, end, pos, pos1, pos2, pos3, 's')) {
			initFMAChildren(s, pos, pos1, pos2, pos3, 's');
		}
		else {
			if (s[begin] == '(') {
				if (s[end - 1] != ')') error("syntax error\n");
				init(s, begin + 1, end - 1);
			}
			else {
				operation = 0;
				expression = s.substr(begin, end - begin);
				op1 = op2 = op3 = NULL;
			}
		}
	}

	static size_t findOperation(const string& s, size_t begin, size_t end, char op) {
		char c;
		int pars = 0;

		for (size_t j = begin; j < end; j++) {
			c = s[j];
			if (c == '(') pars++;
			else if (c == ')') pars--;
			else if (pars == 0 && c == op) {
				if (j + 1 == end) error("syntax error\n");
				else return j;
			}
			if (pars < 0) 
				error("Unexpeced closed parenthesis\n");
		}

		return end;
	}

	static bool findFMA(const string& s, size_t begin, size_t end, size_t& pos0, size_t& pos1, size_t& pos2, size_t& pos3, char aors) {
		char c;
		int pars = 0;
		if (end < 5) return false;

		pos0 = end;
		for (size_t j = begin; j < end-5; j++) {
			c = s[j];
			if (c == '(') pars++;
			else if (c == ')') pars--;
			else if (pars == 0 && c == 'f' && s[j + 1] == 'm' && s[j + 2] == aors && s[j + 3] == '(') {
				pos0 = j;
				break;
			}
			if (pars < 0) 
				error("Unexpeced closed parenthesis\n");
		}

		if (pos0 == end) return false;

		pos0 += 4;
		pars = 1;
		for (pos3 = pos0; pos3 < end; pos3++) {
			c = s[pos3];
			if (c == '(') pars++;
			else if (c == ')') {
				if (--pars == 0) break;
			}
		}

		if (pars != 0) 
			error("Syntax error in FMA expression\n");

		pos1 = findOperation(s, pos0, pos3, ',');
		if (pos1 == pos3) 
			error("Could not find comma ',' in FMA args\n");
		pos2 = findOperation(s, pos1 + 1, pos3, ',');
		if (pos2 == pos3) 
			error("Could not find comma ',' in FMA args\n");

		return true;
	}
};

// Basic type for a predicate equation first = second
typedef std::pair<string, polynomial> equation;

// Predicate
// A set of symbols along with a set of equations.
// The left hand side of the last equation is the
// predicate value.
// For effective filter calculation, translation filters must be
// expressed in separate equations.

class predicate {
public:
	const string predicate_name;
	std::vector<string> symbols; // Initial symbols
	std::vector<string> symbol_type; // if symbol[i] is a lambda of a point pi, then symbol_type[i]="pi", otherwise symbol_type[i]=""
	std::vector<equation> equations;

	bool isIndirect() const {
		for (auto& s : symbol_type) if (!s.empty()) return true;
		return false;
	}

	polynomial getComposedPolynomial(std::vector<string>& maxs) const {
		const equation& last = equations.back();
		const string& lastlhs = last.first;

		polynomial p(last.second);
		for (size_t i = equations.size() - 1; i > 0; i--) {
			const string& s = equations[i - 1].first;
			if (s != lastlhs) p.replaceSymbol(s, equations[i - 1].second);
		}

		// Let's detect here which symbols are relevant for the filter calculation
		// A symbol is relevant if:
		// 1) is the difference of two initial symbols (translation filter)
		// OR
		// 2) is an initial symbol used in expressions which are not (only) translation filters
		//
		// Algorithm
		// 1) create a vector SYM of the initial symbols
		// 2) create a vector SCNT with a reference count for each initial symbol
		// 3) count symbol occurrences and fill the vector of ref counts
		// 4) search translation filters, add their lhss to maxs, and decrease the count of their two symbols
		// 4) for each symbol in SYM, if the corresponding CNT is not zero, add the symbol to maxs

		std::vector<int> cnt; cnt.resize(symbols.size());
		for (size_t i = 0; i < symbols.size(); i++)	cnt[i] = p.countNTFSymbols(symbols[i]);
		for (const equation& e : equations) {
			const polynomial& ep = e.second;
			if (ep.isTranslationFilter() && hasSymbol(ep.op1->expression) && hasSymbol(ep.op2->expression)) maxs.push_back(e.first);
		}
		for (size_t i = 0; i < symbols.size(); i++) if (cnt[i] != 0) maxs.push_back(symbols[i]);

		return p;
	}

	predicate(const string& s, const string& name) : predicate_name(name) {
		stringstream ss(s);
		string line;
		size_t line_num = 0;

		while (std::getline(ss, line)) {
			if (line.empty() || line[0] == '/') continue; // Skip comments
			line_num++;
			size_t pos = line.find('=');
			if (pos == string::npos) addSymbols(line);
			else addEquation(line, pos, line_num);
		}
	}

	void parseImplicitPoint(const string& ip) {
		string ipname, last_sym;
		bool is_ip = true;
		for (size_t i = 13; i < ip.size(); i++) {
			char c = ip[i];
			if (c == ')') {
				symbols.push_back(last_sym);
				symbol_type.push_back(ipname);
				return;
			}
			else if (c == ':') is_ip = false;
			else if (c == ',') {
				symbols.push_back(last_sym);
				symbol_type.push_back(ipname);
				last_sym = "";
			}
			else if (is_ip) ipname += c;
			else last_sym += c;
		}
		error("Syntax error while parsing implicit point\n");
	}

	void addSymbols(const string& s) {
		string token;
		istringstream iss = std::istringstream(s);
		while (iss >> token) {
			if (token.substr(0, 12) == "genericPoint") {
				parseImplicitPoint(token);
			}
			else if (token[0] == '/') return;
			else {
				symbols.push_back(token);
				symbol_type.push_back("");
			}
		}
	}

	void addEquation(const string& s, size_t pos, size_t line_num) {
		istringstream iss = std::istringstream(s);
		string lhs;
		iss >> lhs;
		if (hasSymbol(lhs) || hasLHS(lhs)) {
			std::cerr << "ERROR on line " << line_num << ": Left hand side symbol '" << lhs << "' already defined\n";
			error("");
		}

		size_t nsz = 0;
		for (size_t i = pos; i < s.size(); i++) if (s[i] == '/') break; else nsz++;

		equations.push_back(equation(lhs, polynomial(s.substr(pos + 1, nsz-1))));

		std::set<string> syms;
		equations.back().second.getSymbols(syms);

		for (const string& s2 : syms)
			if (!hasSymbol(s2) && !hasLHS(s2)) {
				std::cerr << "ERROR on line " << line_num << ": Symbol '" << s2 << "' undefined\n";
				error("");
			}
	}

	bool hasSymbol(const string& s) const {
		return std::find(symbols.begin(), symbols.end(), s) != symbols.end();
	}

	bool hasLHS(const string& s) const {
		for (auto& e : equations) if (e.first == s) return true;
		return false;
	}

	string toString() const {
		stringstream ss;
		for (const string& s : symbols) ss << s << " ";
		ss << "\n\n";

		for (const equation& e : equations)
			ss << e.first << " = " << e.second.toString() << "\n";

		return ss.str();
	}

	string toCode() const {
		stringstream ss;
		const string& lhs = equations.back().first;

		// Produce templated filtered function

		ss << "template<class PT, class T> static inline int " << predicate_name << "_t(";

		if (symbols.empty()) error("Empty predicate\n");
		
		string varlist, templlist, protolist;
		string comma, last;
		std::vector<string> lambda_pts;
		std::vector<std::vector<string>> lambda_vars;

		for (size_t i = 0; i < symbols.size(); i++) {
			if (!symbol_type[i].empty()) {
				if (symbol_type[i] != last) {
					last = symbol_type[i];
					templlist += (comma + "const genericPoint& " + last);
					protolist += (comma + "const genericPoint& " + last);
					lambda_pts.push_back(last);
					lambda_vars.resize(lambda_vars.size() + 1);
					varlist += comma; varlist += last;
				}
				lambda_vars.back().push_back(symbols[i]);
			}
			else {
				templlist += (comma + "const PT& " + symbols[i]);
				protolist += (comma + "const double& " + symbols[i]);
				varlist += comma; varlist += symbols[i];
			}
			comma = ", ";
		}
		ss << templlist << ") {\n";

		ss << "\tif constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundUP();\n\n";
		ss << "\tstd::conditional_t<(std::is_same<expansion, T>::value), expansionPool, char> pool;\n\n";
		ss << "\tif constexpr (std::is_same<expansion, T>::value) {\n";
		ss << "\t\texpansion::initPool(&pool);\n";
		ss << "\t\tfeclearexcept(FE_ALL_EXCEPT);\n\t}\n\n";

		for (size_t i = 0; i < lambda_pts.size(); i++) {
			const std::vector<string>& vs = lambda_vars[i];
			string comma2;
			ss << "\tT ";
			for (const string& s : vs) {
				ss << comma2 << s;
				comma2 = ", ";
			}
			ss << ";\n";
			ss << "\tif (!" << lambda_pts[i] << ".getLambda" << vs.size()-1 << "D(";
			comma2 = "";
			for (const string& s : vs) {
				ss << comma2 << s;
				comma2 = ", ";
			}
			ss << ")) {\n\t\tif constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();\n\t\treturn 0;\n\t}\n\n";
		}

		for (const equation& e : equations)
			ss << "\tconst T " << e.first << " = " << e.second.toCode() << ";\n";

		ss << "\n\tif constexpr (std::is_same<interval_number, T>::value) setFPUModeToRoundNEAR();\n";
		ss << "\n\tif constexpr (std::is_same<expansion, T>::value) {\n";
#ifdef GEN_STATISTICS
		ss << "\t\tPred_Stats[Stat_Type::max_expansion_size] = std::max(expansion::getPoolSize(), Pred_Stats[Stat_Type::max_expansion_size]);\n";
#endif
		ss << "\t\tif (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return INT_MAX;\n\t}\n";

		std::vector<string> maxs;
		polynomial c = getComposedPolynomial(maxs);
		fpnumber epsilon = c.getErrorBound();
		int degree = c.getErrorDegree();

		if (!isIndirect()) {
			ss << "\n\tif constexpr (std::is_same<double, T>::value) {\n";

			ss << "\t\tdouble _tmp_fabs;\n";
			ss << "\t\tdouble max_var = 0.0;\n\n";

			for (const string& s : maxs)
				ss << "\t\t_tmp_fabs = fabs(" << s << "); max_var = (_tmp_fabs > max_var)?(_tmp_fabs):(max_var);\n";

			ss << "\t\tdouble epsilon = max_var;\n\n";

			int td = 1;
			int l2d = std::floor(std::log2(degree));
			for (int i = 0; i < l2d; i++) { ss << "\t\tepsilon *= epsilon;\n"; td += td; }
			for (int i = 0; i < (degree - td); i++) ss << "\t\tepsilon *= max_var;\n";
			ss << "\t\tepsilon *= " << std::setprecision(std::numeric_limits<fpnumber>::digits10 + 1) << epsilon << ";\n\n";
			ss << "\t\treturn (" << lhs << " > epsilon) - (" << lhs << " < -epsilon);\n";
			ss << "\t}\n\telse {\n";

			ss << "\t\treturn sgn(" << lhs << ");\n\t}\n}\n\n";

		}
		else {
			ss << "\treturn sgn(" << lhs << ");\n}\n\n";
		}

		// Produce multi-stage function

#ifdef GEN_STATISTICS
		ss << "inline int " << predicate_name << "(" << protolist << ") {\n";
		ss << "\tint ret;\n\tstd::chrono::steady_clock::time_point _ps_time_1, _ps_time_2;\n\n";
		ss << "\t_ps_time_1 = std::chrono::steady_clock::now();\n";
		if (!isIndirect()) {		
			ss << "\tret = " << predicate_name << "_t<double, double>(" << varlist << ");\n";
			ss << "\t_ps_time_2 = std::chrono::steady_clock::now();\n";
			ss << "\tPred_Stats[Stat_Type::" << predicate_name << "_time_double] += std::chrono::duration_cast<std::chrono::nanoseconds>(_ps_time_2 - _ps_time_1).count();\n";

			ss << "\tif (ret != 0) {\n";
			ss << "\t\tPred_Stats[Stat_Type::" << predicate_name << "_num_double]++;\n";
			ss << "\t\treturn ret;\n";
			ss << "\t}\n";
			ss << "\telse _ps_time_1 = _ps_time_2;\n";

		}
		ss << "\tret = " << predicate_name << "_t<interval_number, interval_number>(" << varlist << ");\n";
		ss << "\t_ps_time_2 = std::chrono::steady_clock::now();\n";
		ss << "\tPred_Stats[Stat_Type::" << predicate_name << "_time_interval] += std::chrono::duration_cast<std::chrono::nanoseconds>(_ps_time_2 - _ps_time_1).count();\n";

		ss << "\tif (ret != 0) {\n";
		ss << "\t\tPred_Stats[Stat_Type::" << predicate_name << "_num_interval]++;\n";
		ss << "\t\treturn ret;\n";
		ss << "\t}\n";
		ss << "\telse _ps_time_1 = _ps_time_2;\n";

		ss << "\tret = " << predicate_name << "_t<s_expansion, expansion>(" << varlist << ");\n";
		ss << "\t_ps_time_2 = std::chrono::steady_clock::now();\n";
		ss << "\tPred_Stats[Stat_Type::" << predicate_name << "_time_expansion] += std::chrono::duration_cast<std::chrono::nanoseconds>(_ps_time_2 - _ps_time_1).count();\n";

		ss << "\tif (ret != INT_MAX) {\n";
		ss << "\t\tPred_Stats[Stat_Type::" << predicate_name << "_num_expansion]++;\n";
		ss << "\t\t_ps_time_1 = _ps_time_2;\n";
		ss << "\t}\n";
		ss << "\telse _ps_time_1 = _ps_time_2;\n";

		ss << "\tret = " << predicate_name << "_t<bigfloat, bigfloat>(" << varlist << ");\n";
		ss << "\t_ps_time_2 = std::chrono::steady_clock::now();\n";

		ss << "\tPred_Stats[Stat_Type::" << predicate_name << "_num_bigfloat]++;\n";
		ss << "\tPred_Stats[Stat_Type::" << predicate_name << "_time_bigfloat] += std::chrono::duration_cast<std::chrono::nanoseconds>(_ps_time_2 - _ps_time_1).count();\n";
		ss << "\treturn ret;\n";

		ss << "}\n\n";
#else
		ss << "inline int " << predicate_name << "(" << protolist << ") {\n";
		ss << "\tint ret;\n";
		if (!isIndirect()) ss << "\tif ((ret = " << predicate_name << "_t<double, double>(" << varlist << ")) != 0) return ret;\n";
		else ss << "\tif ((ret = " << predicate_name << "_t<interval_number, interval_number>(" << varlist << ")) != 0) return ret;\n";
		if (degree < 5 && !isIndirect())
			ss << "\tif ((ret = " << predicate_name << "_t<s_expansion, expansion>(" << varlist << ")) != INT_MAX) return ret;\n";
		//ss << "\tif ((ret = " << predicate_name << "_t<s_expansion, expansion>(" << varlist << ")) != INT_MAX) {\n";
		//ss << "\t\tint ret2 = " << predicate_name << "_t<bigfloat, bigfloat>(" << varlist << ");\n";
		//ss << "\t\tassert(ret==ret2);\n";
		//ss << "\t\treturn ret;\n\t}\n";
		ss << "\treturn " << predicate_name << "_t<bigfloat, bigfloat>(" << varlist << ");\n";
		ss << "}\n\n";
#endif
		return ss.str();
	}


	static string createSingleHeader() {
		string s = createHeadingComment();
		s += "#pragma once\n\n";
		s += "#include \"numerics.h\"\n\n";
#ifdef GEN_STATISTICS
		s += "#include <chrono>\n\n";
#endif
		s += "#pragma intrinsic(fabs)\n\n";
		s += double_ipow_functions; s += "\n";
		s += interval_ipow_functions; s += "\n";
		s += expansion_ipow_functions; s += "\n";
		s += bigfloat_ipow_functions; s += "\n";

		//s += "#ifdef USE_AVX2_INSTRUCTIONS\n";
		//s += double_fma_function; s += "\n";
		//s += interval_fma_function; s += "\n";
		//s += bigfloat_fma_function; s += "\n";
		//s += "#endif\n\n";
		s += interval_sgn_function;
		s += "\n\n";

#ifdef GEN_STATISTICS
		s += stats_header1; s += stats_header2; s += "\n";
#endif
		return s;
	}
};


// In the variable declaration, the syntax genericPoint(v1:l1,l2,d) means that
// v1 is an implicit point whose lambdas are l1, l2 and d.
// fma(a,b,c) and fms(a,b,c) are accepted.

int main(int argc, char *argv[])
{
	//string prs = "\
	//	[p: px py] rx ry qx qy\n\
	//	lx = (px^3 - qx) * fma(qy, px, fma(px, px, px))\n\
	//	ly = py - qy\n\
	//	gx = rx - qx\n\
	//	gy = ry - qy\n\
	//	dx = lx * gx\n\
	//	dy = ly * gy\n\
	//	d = dx + dy\n\
	//	";

	if (argc < 2) std::cout << predicate::createSingleHeader() << "\n";
	else {
		string filename = argv[1];
		ifstream file(filename);
		if (!file.is_open()) error("Can't open file\n");
		string s((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
		file.close();
		filename = std::filesystem::path(filename).filename().string(); // Remove path
		filename = filename.substr(0, filename.size() - 4); // Remove extension
		::predicate p(s, filename);
		std::cout << p.toCode() << "\n";
	}

	return 0;
}

