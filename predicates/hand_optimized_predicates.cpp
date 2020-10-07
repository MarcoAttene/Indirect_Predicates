/****************************************************************************
* Indirect predicates for geometric constructions					        *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         * 
* IMATI-GE / CNR                                                            * 
*                                                                           *
* Authors: Marco Attene                                                     * 
* Copyright(C) 2019: IMATI-GE / CNR                                         * 
* All rights reserved.                                                      * 
*                                                                           *
* This program is free software; you can redistribute it and/or modify      * 
* it under the terms of the GNU Lesser General Public License as published  * 
* by the Free Software Foundation; either version 3 of the License, or (at  * 
* your option) any later version.                                           * 
*                                                                           *
* This program is distributed in the hope that it will be useful, but       * 
* WITHOUT ANY WARRANTY; without even the implied warranty of                * 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  * 
* General Public License for more details.                                  * 
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  * 
* along with this program.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/ 

/* Should include incircle and insphere too. */


#include "implicit_point.h"

#pragma intrinsic(fabs)


int orient2d_filtered(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
	double dl = (p2x - p1x) * (p3y - p1y);
	double dr = (p2y - p1y) * (p3x - p1x);
	double det = dl - dr;
	double eb = 3.3306690738754706e-016 * (fabs(dl) + fabs(dr));
	return ((det >= eb) - (-det >= eb));
}

int orient2d_interval(interval_number p1x, interval_number p1y, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   setFPUModeToRoundUP();
   interval_number a11(p2x - p1x);
   interval_number a12(p2y - p1y);
   interval_number a21(p3x - p1x);
   interval_number a22(p3y - p1y);
   interval_number d1(a11 * a22);
   interval_number d2(a12 * a21);
   interval_number d(d1 - d2);
   setFPUModeToRoundNEAR();

   if (!d.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return d.sign();
}

int orient2d_exact(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
   expansionObject o;
   double a11[2];
   o.two_Diff(p2x, p1x, a11);
   double a12[2];
   o.two_Diff(p2y, p1y, a12);
   double a21[2];
   o.two_Diff(p3x, p1x, a21);
   double a22[2];
   o.two_Diff(p3y, p1y, a22);
   double d1[8];
   int d1_len = o.Gen_Product(2, a11, 2, a22, d1);
   double d2[8];
   int d2_len = o.Gen_Product(2, a12, 2, a21, d2);
   double d[16];
   int d_len = o.Gen_Diff(d1_len, d1, d2_len, d2, d);

   double return_value = d[d_len - 1];

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

int orient2d(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2d_filtered(p1x, p1y, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d_interval(p1x, p1y, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d_exact(p1x, p1y, p2x, p2y, p3x, p3y);
}

int orient3d_filtered(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
	double fadx, fbdx, fcdx, fady, fbdy, fcdy, fadz, fbdz, fcdz, eb;
	double fbdxcdy, fcdxbdy, fcdxady, fadxcdy, fadxbdy, fbdxady, det;

	fadx = qx - px; fbdx = rx - px; fcdx = sx - px;
	fady = qy - py; fbdy = ry - py; fcdy = sy - py;
	fadz = qz - pz; fbdz = rz - pz; fcdz = sz - pz;

	fbdxcdy = fbdx * fcdy * fadz; fcdxbdy = fcdx * fbdy * fadz;
	fcdxady = fcdx * fady * fbdz; fadxcdy = fadx * fcdy * fbdz;
	fadxbdy = fadx * fbdy * fcdz; fbdxady = fbdx * fady * fcdz;

	det = (fbdxcdy - fcdxbdy) + (fcdxady - fadxcdy) + (fadxbdy - fbdxady);
	eb = 7.7715611723761027e-016 * (fabs(fbdxcdy) + fabs(fcdxbdy) + fabs(fcdxady) + fabs(fadxcdy) + fabs(fadxbdy) + fabs(fbdxady));
	return ((det >= eb) - (-det >= eb));
}

int orient3d_interval(interval_number px, interval_number py, interval_number pz, interval_number qx, interval_number qy, interval_number qz, interval_number rx, interval_number ry, interval_number rz, interval_number sx, interval_number sy, interval_number sz)
{
   setFPUModeToRoundUP();
   interval_number qx_px(qx - px);
   interval_number qy_py(qy - py);
   interval_number rx_px(rx - px);
   interval_number ry_py(ry - py);
   interval_number rz_pz(rz - pz);
   interval_number qz_pz(qz - pz);
   interval_number sx_px(sx - px);
   interval_number sy_py(sy - py);
   interval_number sz_pz(sz - pz);
   interval_number tmp_a(qx_px * ry_py);
   interval_number tmp_b(qy_py * rx_px);
   interval_number m01(tmp_a - tmp_b);
   interval_number tmq_a(qx_px * rz_pz);
   interval_number tmq_b(qz_pz * rx_px);
   interval_number m02(tmq_a - tmq_b);
   interval_number tmr_a(qy_py * rz_pz);
   interval_number tmr_b(qz_pz * ry_py);
   interval_number m12(tmr_a - tmr_b);
   interval_number mt1(m01 * sz_pz);
   interval_number mt2(m02 * sy_py);
   interval_number mt3(m12 * sx_px);
   interval_number mtt(mt1 - mt2);
   interval_number m012(mtt + mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return m012.sign();
}

int orient3d_exact(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   expansionObject o;
   double qx_px[2];
   o.two_Diff(qx, px, qx_px);
   double qy_py[2];
   o.two_Diff(qy, py, qy_py);
   double rx_px[2];
   o.two_Diff(rx, px, rx_px);
   double ry_py[2];
   o.two_Diff(ry, py, ry_py);
   double rz_pz[2];
   o.two_Diff(rz, pz, rz_pz);
   double qz_pz[2];
   o.two_Diff(qz, pz, qz_pz);
   double sx_px[2];
   o.two_Diff(sx, px, sx_px);
   double sy_py[2];
   o.two_Diff(sy, py, sy_py);
   double sz_pz[2];
   o.two_Diff(sz, pz, sz_pz);
   double tmp_a[8];
   int tmp_a_len = o.Gen_Product(2, qx_px, 2, ry_py, tmp_a);
   double tmp_b[8];
   int tmp_b_len = o.Gen_Product(2, qy_py, 2, rx_px, tmp_b);
   double m01[16];
   int m01_len = o.Gen_Diff(tmp_a_len, tmp_a, tmp_b_len, tmp_b, m01);
   double tmq_a[8];
   int tmq_a_len = o.Gen_Product(2, qx_px, 2, rz_pz, tmq_a);
   double tmq_b[8];
   int tmq_b_len = o.Gen_Product(2, qz_pz, 2, rx_px, tmq_b);
   double m02[16];
   int m02_len = o.Gen_Diff(tmq_a_len, tmq_a, tmq_b_len, tmq_b, m02);
   double tmr_a[8];
   int tmr_a_len = o.Gen_Product(2, qy_py, 2, rz_pz, tmr_a);
   double tmr_b[8];
   int tmr_b_len = o.Gen_Product(2, qz_pz, 2, ry_py, tmr_b);
   double m12[16];
   int m12_len = o.Gen_Diff(tmr_a_len, tmr_a, tmr_b_len, tmr_b, m12);
   double mt1[64];
   int mt1_len = o.Gen_Product(m01_len, m01, 2, sz_pz, mt1);
   double mt2[64];
   int mt2_len = o.Gen_Product(m02_len, m02, 2, sy_py, mt2);
   double mt3[64];
   int mt3_len = o.Gen_Product(m12_len, m12, 2, sx_px, mt3);
   double mtt[128];
   int mtt_len = o.Gen_Diff(mt1_len, mt1, mt2_len, mt2, mtt);
   double m012_p[128], *m012 = m012_p;
   int m012_len = o.Gen_Sum_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 128);

   double return_value = m012[m012_len - 1];
   if (m012_p != m012) free(m012);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 if (return_value == 0) return IP_Sign::ZERO;
 return IP_Sign::UNDEFINED;
}

int orient3d(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz)
{
   int ret;
   ret = orient3d_filtered(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient3d_interval(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_exact(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
}
