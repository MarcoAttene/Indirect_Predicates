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

/* This code was generated automatically. Do not edit unless you exactly   */
/* know what you are doing!                                                */

#include "implicit_point.h"

#pragma intrinsic(fabs)

int incircle_filtered(double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   double adx = pax - pdx;
   double ady = pay - pdy;
   double bdx = pbx - pdx;
   double bdy = pby - pdy;
   double cdx = pcx - pdx;
   double cdy = pcy - pdy;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double alift = alifta + aliftb;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double blift = blifta + bliftb;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double clift = clifta + cliftb;
   double la = alift * bcdet;
   double lb = blift * cadet;
   double lc = clift * abdet;
   double lab = la + lb;
   double L = lab + lc;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(adx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ady)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.376676550535194e-14;

   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_interval(interval_number pax, interval_number pay, interval_number pbx, interval_number pby, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   setFPUModeToRoundUP();
   interval_number adx(pax - pdx);
   interval_number ady(pay - pdy);
   interval_number bdx(pbx - pdx);
   interval_number bdy(pby - pdy);
   interval_number cdx(pcx - pdx);
   interval_number cdy(pcy - pdy);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number alift(alifta + aliftb);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number blift(blifta + bliftb);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number clift(clifta + cliftb);
   interval_number la(alift * bcdet);
   interval_number lb(blift * cadet);
   interval_number lc(clift * abdet);
   interval_number lab(la + lb);
   interval_number L(lab + lc);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_exact(double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   expansionObject o;
   double adx[2];
   o.two_Diff(pax, pdx, adx);
   double ady[2];
   o.two_Diff(pay, pdy, ady);
   double bdx[2];
   o.two_Diff(pbx, pdx, bdx);
   double bdy[2];
   o.two_Diff(pby, pdy, bdy);
   double cdx[2];
   o.two_Diff(pcx, pdx, cdx);
   double cdy[2];
   o.two_Diff(pcy, pdy, cdy);
   double abdeta[8];
   int abdeta_len = o.Gen_Product(2, adx, 2, bdy, abdeta);
   double abdetb[8];
   int abdetb_len = o.Gen_Product(2, bdx, 2, ady, abdetb);
   double abdet[16];
   int abdet_len = o.Gen_Diff(abdeta_len, abdeta, abdetb_len, abdetb, abdet);
   double bcdeta[8];
   int bcdeta_len = o.Gen_Product(2, bdx, 2, cdy, bcdeta);
   double bcdetb[8];
   int bcdetb_len = o.Gen_Product(2, cdx, 2, bdy, bcdetb);
   double bcdet[16];
   int bcdet_len = o.Gen_Diff(bcdeta_len, bcdeta, bcdetb_len, bcdetb, bcdet);
   double cadeta[8];
   int cadeta_len = o.Gen_Product(2, cdx, 2, ady, cadeta);
   double cadetb[8];
   int cadetb_len = o.Gen_Product(2, adx, 2, cdy, cadetb);
   double cadet[16];
   int cadet_len = o.Gen_Diff(cadeta_len, cadeta, cadetb_len, cadetb, cadet);
   double alifta[8];
   int alifta_len = o.Gen_Product(2, adx, 2, adx, alifta);
   double aliftb[8];
   int aliftb_len = o.Gen_Product(2, ady, 2, ady, aliftb);
   double alift[16];
   int alift_len = o.Gen_Sum(alifta_len, alifta, aliftb_len, aliftb, alift);
   double blifta[8];
   int blifta_len = o.Gen_Product(2, bdx, 2, bdx, blifta);
   double bliftb[8];
   int bliftb_len = o.Gen_Product(2, bdy, 2, bdy, bliftb);
   double blift[16];
   int blift_len = o.Gen_Sum(blifta_len, blifta, bliftb_len, bliftb, blift);
   double clifta[8];
   int clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
   double cliftb[8];
   int cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
   double clift[16];
   int clift_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, clift);
   double la_p[128], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 128);
   double lb_p[128], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 128);
   double lc_p[128], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 128);
   double lab_p[128], *lab = lab_p;
   int lab_len = o.Gen_Sum_With_PreAlloc(la_len, la, lb_len, lb, &lab, 128);
   double L_p[128], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, lc_len, lc, &L, 128);

   double return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab_p != lab) free(lab);
   if (lc_p != lc) free(lc);
   if (lb_p != lb) free(lb);
   if (la_p != la) free(la);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle(double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_filtered(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_interval(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_exact(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
}

int orient2d_filtered(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
   double a11 = p2x - p1x;
   double a12 = p2y - p1y;
   double a21 = p3x - p1x;
   double a22 = p3y - p1y;
   double d1 = a11 * a22;
   double d2 = a12 * a21;
   double d = d1 - d2;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(a11)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a12)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a21)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a22)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= 8.881784197001252e-16;
   if (d > epsilon) return IP_Sign::POSITIVE;
   if (-d > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
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
 return IP_Sign::ZERO;
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
   double qx_px = qx - px;
   double qy_py = qy - py;
   double rx_px = rx - px;
   double ry_py = ry - py;
   double rz_pz = rz - pz;
   double qz_pz = qz - pz;
   double sx_px = sx - px;
   double sy_py = sy - py;
   double sz_pz = sz - pz;
   double tmp_a = qx_px * ry_py;
   double tmp_b = qy_py * rx_px;
   double m01 = tmp_a - tmp_b;
   double tmq_a = qx_px * rz_pz;
   double tmq_b = qz_pz * rx_px;
   double m02 = tmq_a - tmq_b;
   double tmr_a = qy_py * rz_pz;
   double tmr_b = qz_pz * ry_py;
   double m12 = tmr_a - tmr_b;
   double mt1 = m01 * sz_pz;
   double mt2 = m02 * sy_py;
   double mt3 = m12 * sx_px;
   double mtt = mt1 - mt2;
   double m012 = mtt + mt3;

   double _tmp_fabs;

   double max_var = 0.0;
   if ((_tmp_fabs = fabs(qx_px)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qy_py)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(rx_px)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ry_py)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(rz_pz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(qz_pz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sx_px)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sy_py)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(sz_pz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= 4.884981308350689e-15;
   if (m012 > epsilon) return IP_Sign::POSITIVE;
   if (-m012 > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
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
 return IP_Sign::ZERO;
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

int incircle_indirect_LEEE_filtered(implicitPoint3D_LPI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double pdxt = pdx * d1;
   double pdyt = pdy * d1;
   double adx = l1x - pdxt;
   double ady = l1y - pdyt;
   double bdx = pbx - pdx;
   double bdy = pby - pdy;
   double cdx = pcx - pdx;
   double cdy = pcy - pdy;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double alift = alifta + aliftb;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double blift = blifta + bliftb;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double clift = clifta + cliftb;
   double la = alift * bcdet;
   double lbt = blift * cadet;
   double lb = lbt * d1;
   double lct = clift * abdet;
   double lc = lct * d1;
   double lab = la + lb;
   double L = lab + lc;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(pdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 6.463567528425267e-12;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_LEEE_interval(implicitPoint3D_LPI& p1, interval_number pbx, interval_number pby, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number pdxt(pdx * d1);
   interval_number pdyt(pdy * d1);
   interval_number adx(l1x - pdxt);
   interval_number ady(l1y - pdyt);
   interval_number bdx(pbx - pdx);
   interval_number bdy(pby - pdy);
   interval_number cdx(pcx - pdx);
   interval_number cdy(pcy - pdy);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number alift(alifta + aliftb);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number blift(blifta + bliftb);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number clift(clifta + cliftb);
   interval_number la(alift * bcdet);
   interval_number lbt(blift * cadet);
   interval_number lb(lbt * d1);
   interval_number lct(clift * abdet);
   interval_number lc(lct * d1);
   interval_number lab(la + lb);
   interval_number L(lab + lc);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_LEEE_exact(implicitPoint3D_LPI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double pdxt_p[64], *pdxt = pdxt_p;
   int pdxt_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdxt, 64);
   double pdyt_p[64], *pdyt = pdyt_p;
   int pdyt_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdyt, 64);
   double adx_p[64], *adx = adx_p;
   int adx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdxt_len, pdxt, &adx, 64);
   double ady_p[64], *ady = ady_p;
   int ady_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdyt_len, pdyt, &ady, 64);
   double bdx[2];
   o.two_Diff(pbx, pdx, bdx);
   double bdy[2];
   o.two_Diff(pby, pdy, bdy);
   double cdx[2];
   o.two_Diff(pcx, pdx, cdx);
   double cdy[2];
   o.two_Diff(pcy, pdy, cdy);
   double abdeta_p[64], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, 2, bdy, &abdeta, 64);
   double abdetb_p[64], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(2, bdx, ady_len, ady, &abdetb, 64);
   double abdet_p[64], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 64);
   double bcdeta[8];
   int bcdeta_len = o.Gen_Product(2, bdx, 2, cdy, bcdeta);
   double bcdetb[8];
   int bcdetb_len = o.Gen_Product(2, cdx, 2, bdy, bcdetb);
   double bcdet[16];
   int bcdet_len = o.Gen_Diff(bcdeta_len, bcdeta, bcdetb_len, bcdetb, bcdet);
   double cadeta_p[64], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 64);
   double cadetb_p[64], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 64);
   double cadet_p[64], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 64);
   double alifta_p[64], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 64);
   double aliftb_p[64], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 64);
   double alift_p[64], *alift = alift_p;
   int alift_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &alift, 64);
   double blifta[8];
   int blifta_len = o.Gen_Product(2, bdx, 2, bdx, blifta);
   double bliftb[8];
   int bliftb_len = o.Gen_Product(2, bdy, 2, bdy, bliftb);
   double blift[16];
   int blift_len = o.Gen_Sum(blifta_len, blifta, bliftb_len, bliftb, blift);
   double clifta[8];
   int clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
   double cliftb[8];
   int cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
   double clift[16];
   int clift_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, clift);
   double la_p[64], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 64);
   double lbt_p[64], *lbt = lbt_p;
   int lbt_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lbt, 64);
   double lb_p[64], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(lbt_len, lbt, d1_len, d1, &lb, 64);
   double lct_p[64], *lct = lct_p;
   int lct_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lct, 64);
   double lc_p[64], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(lct_len, lct, d1_len, d1, &lc, 64);
   double lab_p[64], *lab = lab_p;
   int lab_len = o.Gen_Sum_With_PreAlloc(la_len, la, lb_len, lb, &lab, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, lc_len, lc, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab_p != lab) free(lab);
   if (lc_p != lc) free(lc);
   if (lct_p != lct) free(lct);
   if (lb_p != lb) free(lb);
   if (lbt_p != lbt) free(lbt);
   if (la_p != la) free(la);
   if (alift_p != alift) free(alift);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
   if (ady_p != ady) free(ady);
   if (adx_p != adx) free(adx);
   if (pdyt_p != pdyt) free(pdyt);
   if (pdxt_p != pdxt) free(pdxt);
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_LEEE(implicitPoint3D_LPI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_LEEE_filtered(p1, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_LEEE_interval(p1, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_LEEE_exact(p1, pbx, pby, pcx, pcy, pdx, pdy);
}

int incircle_indirect_LLEE_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double pcx, double pcy, double pdx, double pdy)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double pdx1 = pdx * d1;
   double pdy1 = pdy * d1;
   double adx = l1x - pdx1;
   double ady = l1y - pdy1;
   double pdx2 = pdx * d2;
   double pdy2 = pdy * d2;
   double bdx = l2x - pdx2;
   double bdy = l2y - pdy2;
   double cdx = pcx - pdx;
   double cdy = pcy - pdy;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double aliftt = alifta + aliftb;
   double alift = aliftt * d2;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double blift = blifta + bliftb;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double cliftt = clifta + cliftb;
   double clift = cliftt * d2;
   double la = alift * bcdet;
   double lb = blift * cadet;
   double lc = clift * abdet;
   double lab = lc + lb;
   double lab2 = lab * d1;
   double L = lab2 + la;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(pdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.564426935218861e-09;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_LLEE_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number pdx1(pdx * d1);
   interval_number pdy1(pdy * d1);
   interval_number adx(l1x - pdx1);
   interval_number ady(l1y - pdy1);
   interval_number pdx2(pdx * d2);
   interval_number pdy2(pdy * d2);
   interval_number bdx(l2x - pdx2);
   interval_number bdy(l2y - pdy2);
   interval_number cdx(pcx - pdx);
   interval_number cdy(pcy - pdy);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number aliftt(alifta + aliftb);
   interval_number alift(aliftt * d2);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number blift(blifta + bliftb);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number cliftt(clifta + cliftb);
   interval_number clift(cliftt * d2);
   interval_number la(alift * bcdet);
   interval_number lb(blift * cadet);
   interval_number lc(clift * abdet);
   interval_number lab(lc + lb);
   interval_number lab2(lab * d1);
   interval_number L(lab2 + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_LLEE_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double pcx, double pcy, double pdx, double pdy)
{
 double return_value = 0.0;
 double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32], *l1z = l1z_p, d1_p[32], *d1 = d1_p, l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p, l2z_p[32], *l2z = l2z_p, d2_p[32], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double pdx1_p[32], *pdx1 = pdx1_p;
   int pdx1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdx1, 32);
   double pdy1_p[32], *pdy1 = pdy1_p;
   int pdy1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdy1, 32);
   double adx_p[32], *adx = adx_p;
   int adx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
   double ady_p[32], *ady = ady_p;
   int ady_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
   double pdx2_p[32], *pdx2 = pdx2_p;
   int pdx2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdx, &pdx2, 32);
   double pdy2_p[32], *pdy2 = pdy2_p;
   int pdy2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdy, &pdy2, 32);
   double bdx_p[32], *bdx = bdx_p;
   int bdx_len = o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
   double bdy_p[32], *bdy = bdy_p;
   int bdy_len = o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
   double cdx[2];
   o.two_Diff(pcx, pdx, cdx);
   double cdy[2];
   o.two_Diff(pcy, pdy, cdy);
   double abdeta_p[32], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
   double abdetb_p[32], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
   double abdet_p[32], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 32);
   double bcdeta_p[32], *bcdeta = bcdeta_p;
   int bcdeta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, 2, cdy, &bcdeta, 32);
   double bcdetb_p[32], *bcdetb = bcdetb_p;
   int bcdetb_len = o.Gen_Product_With_PreAlloc(2, cdx, bdy_len, bdy, &bcdetb, 32);
   double bcdet_p[32], *bcdet = bcdet_p;
   int bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len, bcdetb, &bcdet, 32);
   double cadeta_p[32], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 32);
   double cadetb_p[32], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 32);
   double cadet_p[32], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 32);
   double alifta_p[32], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
   double aliftb_p[32], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
   double aliftt_p[32], *aliftt = aliftt_p;
   int aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &aliftt, 32);
   double alift_p[32], *alift = alift_p;
   int alift_len = o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift, 32);
   double blifta_p[32], *blifta = blifta_p;
   int blifta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
   double bliftb_p[32], *bliftb = bliftb_p;
   int bliftb_len = o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
   double blift_p[32], *blift = blift_p;
   int blift_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len, bliftb, &blift, 32);
   double clifta[8];
   int clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
   double cliftb[8];
   int cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
   double cliftt[16];
   int cliftt_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, cliftt);
   double clift_p[32], *clift = clift_p;
   int clift_len = o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
   double la_p[32], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
   double lb_p[32], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
   double lc_p[32], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
   double lab_p[32], *lab = lab_p;
   int lab_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab, 32);
   double lab2_p[32], *lab2 = lab2_p;
   int lab2_len = o.Gen_Product_With_PreAlloc(lab_len, lab, d1_len, d1, &lab2, 32);
   double L_p[32], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab2_len, lab2, la_len, la, &L, 32);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab2_p != lab2) free(lab2);
   if (lab_p != lab) free(lab);
   if (lc_p != lc) free(lc);
   if (lb_p != lb) free(lb);
   if (la_p != la) free(la);
   if (clift_p != clift) free(clift);
   if (blift_p != blift) free(blift);
   if (bliftb_p != bliftb) free(bliftb);
   if (blifta_p != blifta) free(blifta);
   if (alift_p != alift) free(alift);
   if (aliftt_p != aliftt) free(aliftt);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (bcdet_p != bcdet) free(bcdet);
   if (bcdetb_p != bcdetb) free(bcdetb);
   if (bcdeta_p != bcdeta) free(bcdeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
   if (bdy_p != bdy) free(bdy);
   if (bdx_p != bdx) free(bdx);
   if (pdy2_p != pdy2) free(pdy2);
   if (pdx2_p != pdx2) free(pdx2);
   if (ady_p != ady) free(ady);
   if (adx_p != adx) free(adx);
   if (pdy1_p != pdy1) free(pdy1);
   if (pdx1_p != pdx1) free(pdx1);
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_LLEE(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_LLEE_filtered(p1, p2, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_LLEE_interval(p1, p2, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_LLEE_exact(p1, p2, pcx, pcy, pdx, pdy);
}

int incircle_indirect_LLLE_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, double pdx, double pdy)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double pdx1 = pdx * d1;
   double pdy1 = pdy * d1;
   double adx = l1x - pdx1;
   double ady = l1y - pdy1;
   double pdx2 = pdx * d2;
   double pdy2 = pdy * d2;
   double bdx = l2x - pdx2;
   double bdy = l2y - pdy2;
   double pdx3 = pdx * d3;
   double pdy3 = pdy * d3;
   double cdx = l3x - pdx3;
   double cdy = l3y - pdy3;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double aliftt = alifta + aliftb;
   double alift2 = aliftt * d2;
   double alift = alift2 * d3;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double bliftt = blifta + bliftb;
   double blift = bliftt * d3;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double cliftt = clifta + cliftb;
   double clift = cliftt * d2;
   double la = alift * bcdet;
   double lb = blift * cadet;
   double lc = clift * abdet;
   double lab2 = lc + lb;
   double lab = lab2 * d1;
   double L = lab + la;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(pdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.004469273093523e-07;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_LLLE_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number pdx1(pdx * d1);
   interval_number pdy1(pdy * d1);
   interval_number adx(l1x - pdx1);
   interval_number ady(l1y - pdy1);
   interval_number pdx2(pdx * d2);
   interval_number pdy2(pdy * d2);
   interval_number bdx(l2x - pdx2);
   interval_number bdy(l2y - pdy2);
   interval_number pdx3(pdx * d3);
   interval_number pdy3(pdy * d3);
   interval_number cdx(l3x - pdx3);
   interval_number cdy(l3y - pdy3);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number aliftt(alifta + aliftb);
   interval_number alift2(aliftt * d2);
   interval_number alift(alift2 * d3);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number bliftt(blifta + bliftb);
   interval_number blift(bliftt * d3);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number cliftt(clifta + cliftb);
   interval_number clift(cliftt * d2);
   interval_number la(alift * bcdet);
   interval_number lb(blift * cadet);
   interval_number lc(clift * abdet);
   interval_number lab2(lc + lb);
   interval_number lab(lab2 * d1);
   interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_LLLE_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, double pdx, double pdy)
{
 double return_value = 0.0;
 double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32], *l1z = l1z_p, d1_p[32], *d1 = d1_p, l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p, l2z_p[32], *l2z = l2z_p, d2_p[32], *d2 = d2_p, l3x_p[32], *l3x = l3x_p, l3y_p[32], *l3y = l3y_p, l3z_p[32], *l3z = l3z_p, d3_p[32], *d3 = d3_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
 {
   expansionObject o;
   double pdx1_p[32], *pdx1 = pdx1_p;
   int pdx1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdx1, 32);
   double pdy1_p[32], *pdy1 = pdy1_p;
   int pdy1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdy1, 32);
   double adx_p[32], *adx = adx_p;
   int adx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
   double ady_p[32], *ady = ady_p;
   int ady_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
   double pdx2_p[32], *pdx2 = pdx2_p;
   int pdx2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdx, &pdx2, 32);
   double pdy2_p[32], *pdy2 = pdy2_p;
   int pdy2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdy, &pdy2, 32);
   double bdx_p[32], *bdx = bdx_p;
   int bdx_len = o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
   double bdy_p[32], *bdy = bdy_p;
   int bdy_len = o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
   double pdx3_p[32], *pdx3 = pdx3_p;
   int pdx3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pdx, &pdx3, 32);
   double pdy3_p[32], *pdy3 = pdy3_p;
   int pdy3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pdy, &pdy3, 32);
   double cdx_p[32], *cdx = cdx_p;
   int cdx_len = o.Gen_Diff_With_PreAlloc(l3x_len, l3x, pdx3_len, pdx3, &cdx, 32);
   double cdy_p[32], *cdy = cdy_p;
   int cdy_len = o.Gen_Diff_With_PreAlloc(l3y_len, l3y, pdy3_len, pdy3, &cdy, 32);
   double abdeta_p[32], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
   double abdetb_p[32], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
   double abdet_p[32], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 32);
   double bcdeta_p[32], *bcdeta = bcdeta_p;
   int bcdeta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 32);
   double bcdetb_p[32], *bcdetb = bcdetb_p;
   int bcdetb_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 32);
   double bcdet_p[32], *bcdet = bcdet_p;
   int bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len, bcdetb, &bcdet, 32);
   double cadeta_p[32], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 32);
   double cadetb_p[32], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 32);
   double cadet_p[32], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 32);
   double alifta_p[32], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
   double aliftb_p[32], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
   double aliftt_p[32], *aliftt = aliftt_p;
   int aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &aliftt, 32);
   double alift2_p[32], *alift2 = alift2_p;
   int alift2_len = o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 32);
   double alift_p[32], *alift = alift_p;
   int alift_len = o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 32);
   double blifta_p[32], *blifta = blifta_p;
   int blifta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
   double bliftb_p[32], *bliftb = bliftb_p;
   int bliftb_len = o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
   double bliftt_p[32], *bliftt = bliftt_p;
   int bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len, bliftb, &bliftt, 32);
   double blift_p[32], *blift = blift_p;
   int blift_len = o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 32);
   double clifta_p[32], *clifta = clifta_p;
   int clifta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 32);
   double cliftb_p[32], *cliftb = cliftb_p;
   int cliftb_len = o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 32);
   double cliftt_p[32], *cliftt = cliftt_p;
   int cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len, cliftb, &cliftt, 32);
   double clift_p[32], *clift = clift_p;
   int clift_len = o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
   double la_p[32], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
   double lb_p[32], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
   double lc_p[32], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
   double lab2_p[32], *lab2 = lab2_p;
   int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 32);
   double lab_p[32], *lab = lab_p;
   int lab_len = o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 32);
   double L_p[32], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 32);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab_p != lab) free(lab);
   if (lab2_p != lab2) free(lab2);
   if (lc_p != lc) free(lc);
   if (lb_p != lb) free(lb);
   if (la_p != la) free(la);
   if (clift_p != clift) free(clift);
   if (cliftt_p != cliftt) free(cliftt);
   if (cliftb_p != cliftb) free(cliftb);
   if (clifta_p != clifta) free(clifta);
   if (blift_p != blift) free(blift);
   if (bliftt_p != bliftt) free(bliftt);
   if (bliftb_p != bliftb) free(bliftb);
   if (blifta_p != blifta) free(blifta);
   if (alift_p != alift) free(alift);
   if (alift2_p != alift2) free(alift2);
   if (aliftt_p != aliftt) free(aliftt);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (bcdet_p != bcdet) free(bcdet);
   if (bcdetb_p != bcdetb) free(bcdetb);
   if (bcdeta_p != bcdeta) free(bcdeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
   if (cdy_p != cdy) free(cdy);
   if (cdx_p != cdx) free(cdx);
   if (pdy3_p != pdy3) free(pdy3);
   if (pdx3_p != pdx3) free(pdx3);
   if (bdy_p != bdy) free(bdy);
   if (bdx_p != bdx) free(bdx);
   if (pdy2_p != pdy2) free(pdy2);
   if (pdx2_p != pdx2) free(pdx2);
   if (ady_p != ady) free(ady);
   if (adx_p != adx) free(adx);
   if (pdy1_p != pdy1) free(pdy1);
   if (pdx1_p != pdx1) free(pdx1);
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);
 if (l3x_p != l3x) free(l3x);
 if (l3y_p != l3y) free(l3y);
 if (l3z_p != l3z) free(l3z);
 if (d3_p != d3) free(d3);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_LLLE(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, double pdx, double pdy)
{
   int ret;
//   ret = incircle_indirect_LLLE_filtered(p1, p2, p3, pdx, pdy);
//   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_LLLE_interval(p1, p2, p3, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_LLLE_exact(p1, p2, p3, pdx, pdy);
}

int incircle_indirect_LLLL_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, implicitPoint3D_LPI& p4)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var)
       || !p4.getFilteredLambda(l4x, l4y, l4z, d4, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double l1xt = l1x * d4;
   double l1yt = l1y * d4;
   double l2xt = l2x * d4;
   double l2yt = l2y * d4;
   double l3xt = l3x * d4;
   double l3yt = l3y * d4;
   double l4x1 = l4x * d1;
   double l4y1 = l4y * d1;
   double adx = l1xt - l4x1;
   double ady = l1yt - l4y1;
   double l4x2 = l4x * d2;
   double l4y2 = l4y * d2;
   double bdx = l2xt - l4x2;
   double bdy = l2yt - l4y2;
   double l4x3 = l4x * d3;
   double l4y3 = l4y * d3;
   double cdx = l3xt - l4x3;
   double cdy = l3yt - l4y3;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double aliftt = alifta + aliftb;
   double alift2 = aliftt * d2;
   double alift = alift2 * d3;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double bliftt = blifta + bliftb;
   double blift = bliftt * d3;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double cliftt = clifta + cliftb;
   double clift = cliftt * d2;
   double la = alift * bcdet;
   double lb = blift * cadet;
   double lc = clift * abdet;
   double lab2 = lc + lb;
   double lab = lab2 * d1;
   double L = lab + la;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 0.001941157917940494;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_LLLL_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, implicitPoint3D_LPI& p4)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   || !p4.getIntervalLambda(l4x, l4y, l4z, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number l1xt(l1x * d4);
   interval_number l1yt(l1y * d4);
   interval_number l2xt(l2x * d4);
   interval_number l2yt(l2y * d4);
   interval_number l3xt(l3x * d4);
   interval_number l3yt(l3y * d4);
   interval_number l4x1(l4x * d1);
   interval_number l4y1(l4y * d1);
   interval_number adx(l1xt - l4x1);
   interval_number ady(l1yt - l4y1);
   interval_number l4x2(l4x * d2);
   interval_number l4y2(l4y * d2);
   interval_number bdx(l2xt - l4x2);
   interval_number bdy(l2yt - l4y2);
   interval_number l4x3(l4x * d3);
   interval_number l4y3(l4y * d3);
   interval_number cdx(l3xt - l4x3);
   interval_number cdy(l3yt - l4y3);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number aliftt(alifta + aliftb);
   interval_number alift2(aliftt * d2);
   interval_number alift(alift2 * d3);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number bliftt(blifta + bliftb);
   interval_number blift(bliftt * d3);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number cliftt(clifta + cliftb);
   interval_number clift(cliftt * d2);
   interval_number la(alift * bcdet);
   interval_number lb(blift * cadet);
   interval_number lc(clift * abdet);
   interval_number lab2(lc + lb);
   interval_number lab(lab2 * d1);
   interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_LLLL_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, implicitPoint3D_LPI& p4)
{
 double return_value = 0.0;
 double l1x_p[16], *l1x = l1x_p, l1y_p[16], *l1y = l1y_p, l1z_p[16], *l1z = l1z_p, d1_p[16], *d1 = d1_p, l2x_p[16], *l2x = l2x_p, l2y_p[16], *l2y = l2y_p, l2z_p[16], *l2z = l2z_p, d2_p[16], *d2 = d2_p, l3x_p[16], *l3x = l3x_p, l3y_p[16], *l3y = l3y_p, l3z_p[16], *l3z = l3z_p, d3_p[16], *d3 = d3_p, l4x_p[16], *l4x = l4x_p, l4y_p[16], *l4y = l4y_p, l4z_p[16], *l4z = l4z_p, d4_p[16], *d4 = d4_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len, l4x_len, l4y_len, l4z_len, d4_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
 p4.getExactLambda(l4x, l4x_len, l4y, l4y_len, l4z, l4z_len, d4, d4_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) && (d4[d4_len - 1] != 0))
 {
   expansionObject o;
   double l1xt_p[16], *l1xt = l1xt_p;
   int l1xt_len = o.Gen_Product_With_PreAlloc(l1x_len, l1x, d4_len, d4, &l1xt, 16);
   double l1yt_p[16], *l1yt = l1yt_p;
   int l1yt_len = o.Gen_Product_With_PreAlloc(l1y_len, l1y, d4_len, d4, &l1yt, 16);
   double l2xt_p[16], *l2xt = l2xt_p;
   int l2xt_len = o.Gen_Product_With_PreAlloc(l2x_len, l2x, d4_len, d4, &l2xt, 16);
   double l2yt_p[16], *l2yt = l2yt_p;
   int l2yt_len = o.Gen_Product_With_PreAlloc(l2y_len, l2y, d4_len, d4, &l2yt, 16);
   double l3xt_p[16], *l3xt = l3xt_p;
   int l3xt_len = o.Gen_Product_With_PreAlloc(l3x_len, l3x, d4_len, d4, &l3xt, 16);
   double l3yt_p[16], *l3yt = l3yt_p;
   int l3yt_len = o.Gen_Product_With_PreAlloc(l3y_len, l3y, d4_len, d4, &l3yt, 16);
   double l4x1_p[16], *l4x1 = l4x1_p;
   int l4x1_len = o.Gen_Product_With_PreAlloc(l4x_len, l4x, d1_len, d1, &l4x1, 16);
   double l4y1_p[16], *l4y1 = l4y1_p;
   int l4y1_len = o.Gen_Product_With_PreAlloc(l4y_len, l4y, d1_len, d1, &l4y1, 16);
   double adx_p[16], *adx = adx_p;
   int adx_len = o.Gen_Diff_With_PreAlloc(l1xt_len, l1xt, l4x1_len, l4x1, &adx, 16);
   double ady_p[16], *ady = ady_p;
   int ady_len = o.Gen_Diff_With_PreAlloc(l1yt_len, l1yt, l4y1_len, l4y1, &ady, 16);
   double l4x2_p[16], *l4x2 = l4x2_p;
   int l4x2_len = o.Gen_Product_With_PreAlloc(l4x_len, l4x, d2_len, d2, &l4x2, 16);
   double l4y2_p[16], *l4y2 = l4y2_p;
   int l4y2_len = o.Gen_Product_With_PreAlloc(l4y_len, l4y, d2_len, d2, &l4y2, 16);
   double bdx_p[16], *bdx = bdx_p;
   int bdx_len = o.Gen_Diff_With_PreAlloc(l2xt_len, l2xt, l4x2_len, l4x2, &bdx, 16);
   double bdy_p[16], *bdy = bdy_p;
   int bdy_len = o.Gen_Diff_With_PreAlloc(l2yt_len, l2yt, l4y2_len, l4y2, &bdy, 16);
   double l4x3_p[16], *l4x3 = l4x3_p;
   int l4x3_len = o.Gen_Product_With_PreAlloc(l4x_len, l4x, d3_len, d3, &l4x3, 16);
   double l4y3_p[16], *l4y3 = l4y3_p;
   int l4y3_len = o.Gen_Product_With_PreAlloc(l4y_len, l4y, d3_len, d3, &l4y3, 16);
   double cdx_p[16], *cdx = cdx_p;
   int cdx_len = o.Gen_Diff_With_PreAlloc(l3xt_len, l3xt, l4x3_len, l4x3, &cdx, 16);
   double cdy_p[16], *cdy = cdy_p;
   int cdy_len = o.Gen_Diff_With_PreAlloc(l3yt_len, l3yt, l4y3_len, l4y3, &cdy, 16);
   double abdeta_p[16], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 16);
   double abdetb_p[16], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 16);
   double abdet_p[16], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 16);
   double bcdeta_p[16], *bcdeta = bcdeta_p;
   int bcdeta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 16);
   double bcdetb_p[16], *bcdetb = bcdetb_p;
   int bcdetb_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 16);
   double bcdet_p[16], *bcdet = bcdet_p;
   int bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len, bcdetb, &bcdet, 16);
   double cadeta_p[16], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 16);
   double cadetb_p[16], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 16);
   double cadet_p[16], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 16);
   double alifta_p[16], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 16);
   double aliftb_p[16], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 16);
   double aliftt_p[16], *aliftt = aliftt_p;
   int aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &aliftt, 16);
   double alift2_p[16], *alift2 = alift2_p;
   int alift2_len = o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 16);
   double alift_p[16], *alift = alift_p;
   int alift_len = o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 16);
   double blifta_p[16], *blifta = blifta_p;
   int blifta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 16);
   double bliftb_p[16], *bliftb = bliftb_p;
   int bliftb_len = o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 16);
   double bliftt_p[16], *bliftt = bliftt_p;
   int bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len, bliftb, &bliftt, 16);
   double blift_p[16], *blift = blift_p;
   int blift_len = o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 16);
   double clifta_p[16], *clifta = clifta_p;
   int clifta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 16);
   double cliftb_p[16], *cliftb = cliftb_p;
   int cliftb_len = o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 16);
   double cliftt_p[16], *cliftt = cliftt_p;
   int cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len, cliftb, &cliftt, 16);
   double clift_p[16], *clift = clift_p;
   int clift_len = o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 16);
   double la_p[16], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 16);
   double lb_p[16], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 16);
   double lc_p[16], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 16);
   double lab2_p[16], *lab2 = lab2_p;
   int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 16);
   double lab_p[16], *lab = lab_p;
   int lab_len = o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 16);
   double L_p[16], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 16);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab_p != lab) free(lab);
   if (lab2_p != lab2) free(lab2);
   if (lc_p != lc) free(lc);
   if (lb_p != lb) free(lb);
   if (la_p != la) free(la);
   if (clift_p != clift) free(clift);
   if (cliftt_p != cliftt) free(cliftt);
   if (cliftb_p != cliftb) free(cliftb);
   if (clifta_p != clifta) free(clifta);
   if (blift_p != blift) free(blift);
   if (bliftt_p != bliftt) free(bliftt);
   if (bliftb_p != bliftb) free(bliftb);
   if (blifta_p != blifta) free(blifta);
   if (alift_p != alift) free(alift);
   if (alift2_p != alift2) free(alift2);
   if (aliftt_p != aliftt) free(aliftt);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (bcdet_p != bcdet) free(bcdet);
   if (bcdetb_p != bcdetb) free(bcdetb);
   if (bcdeta_p != bcdeta) free(bcdeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
   if (cdy_p != cdy) free(cdy);
   if (cdx_p != cdx) free(cdx);
   if (l4y3_p != l4y3) free(l4y3);
   if (l4x3_p != l4x3) free(l4x3);
   if (bdy_p != bdy) free(bdy);
   if (bdx_p != bdx) free(bdx);
   if (l4y2_p != l4y2) free(l4y2);
   if (l4x2_p != l4x2) free(l4x2);
   if (ady_p != ady) free(ady);
   if (adx_p != adx) free(adx);
   if (l4y1_p != l4y1) free(l4y1);
   if (l4x1_p != l4x1) free(l4x1);
   if (l3yt_p != l3yt) free(l3yt);
   if (l3xt_p != l3xt) free(l3xt);
   if (l2yt_p != l2yt) free(l2yt);
   if (l2xt_p != l2xt) free(l2xt);
   if (l1yt_p != l1yt) free(l1yt);
   if (l1xt_p != l1xt) free(l1xt);
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);
 if (l3x_p != l3x) free(l3x);
 if (l3y_p != l3y) free(l3y);
 if (l3z_p != l3z) free(l3z);
 if (d3_p != d3) free(d3);
 if (l4x_p != l4x) free(l4x);
 if (l4y_p != l4y) free(l4y);
 if (l4z_p != l4z) free(l4z);
 if (d4_p != d4) free(d4);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_LLLL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, implicitPoint3D_LPI& p4)
{
   int ret;
//   ret = incircle_indirect_LLLL_filtered(p1, p2, p3, p4);
//   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_LLLL_interval(p1, p2, p3, p4);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_LLLL_exact(p1, p2, p3, p4);
}

int incircle_indirect_SEEE_filtered(implicitPoint2D_SSI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   double l1x, l1y, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double pdxt = pdx * d1;
   double pdyt = pdy * d1;
   double adx = l1x - pdxt;
   double ady = l1y - pdyt;
   double bdx = pbx - pdx;
   double bdy = pby - pdy;
   double cdx = pcx - pdx;
   double cdy = pcy - pdy;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double alift = alifta + aliftb;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double blift = blifta + bliftb;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double clift = clifta + cliftb;
   double la = alift * bcdet;
   double lbt = blift * cadet;
   double lb = lbt * d1;
   double lct = clift * abdet;
   double lc = lct * d1;
   double lab = la + lb;
   double L = lab + lc;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(pdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 5.010124259907664e-13;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_SEEE_interval(implicitPoint2D_SSI& p1, interval_number pbx, interval_number pby, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number pdxt(pdx * d1);
   interval_number pdyt(pdy * d1);
   interval_number adx(l1x - pdxt);
   interval_number ady(l1y - pdyt);
   interval_number bdx(pbx - pdx);
   interval_number bdy(pby - pdy);
   interval_number cdx(pcx - pdx);
   interval_number cdy(pcy - pdy);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number alift(alifta + aliftb);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number blift(blifta + bliftb);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number clift(clifta + cliftb);
   interval_number la(alift * bcdet);
   interval_number lbt(blift * cadet);
   interval_number lb(lbt * d1);
   interval_number lct(clift * abdet);
   interval_number lc(lct * d1);
   interval_number lab(la + lb);
   interval_number L(lab + lc);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_SEEE_exact(implicitPoint2D_SSI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
 double return_value = 0.0;
 double l1x[32], l1y[32], d1[16];
 int l1x_len, l1y_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double pdxt[32];
   int pdxt_len = o.Gen_Scale(d1_len, d1, pdx, pdxt);
   double pdyt[32];
   int pdyt_len = o.Gen_Scale(d1_len, d1, pdy, pdyt);
   double adx[64];
   int adx_len = o.Gen_Diff(l1x_len, l1x, pdxt_len, pdxt, adx);
   double ady[64];
   int ady_len = o.Gen_Diff(l1y_len, l1y, pdyt_len, pdyt, ady);
   double bdx[2];
   o.two_Diff(pbx, pdx, bdx);
   double bdy[2];
   o.two_Diff(pby, pdy, bdy);
   double cdx[2];
   o.two_Diff(pcx, pdx, cdx);
   double cdy[2];
   o.two_Diff(pcy, pdy, cdy);
   double abdeta_p[64], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, 2, bdy, &abdeta, 64);
   double abdetb_p[64], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(2, bdx, ady_len, ady, &abdetb, 64);
   double abdet_p[64], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 64);
   double bcdeta[8];
   int bcdeta_len = o.Gen_Product(2, bdx, 2, cdy, bcdeta);
   double bcdetb[8];
   int bcdetb_len = o.Gen_Product(2, cdx, 2, bdy, bcdetb);
   double bcdet[16];
   int bcdet_len = o.Gen_Diff(bcdeta_len, bcdeta, bcdetb_len, bcdetb, bcdet);
   double cadeta_p[64], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 64);
   double cadetb_p[64], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 64);
   double cadet_p[64], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 64);
   double alifta_p[64], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 64);
   double aliftb_p[64], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 64);
   double alift_p[64], *alift = alift_p;
   int alift_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &alift, 64);
   double blifta[8];
   int blifta_len = o.Gen_Product(2, bdx, 2, bdx, blifta);
   double bliftb[8];
   int bliftb_len = o.Gen_Product(2, bdy, 2, bdy, bliftb);
   double blift[16];
   int blift_len = o.Gen_Sum(blifta_len, blifta, bliftb_len, bliftb, blift);
   double clifta[8];
   int clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
   double cliftb[8];
   int cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
   double clift[16];
   int clift_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, clift);
   double la_p[64], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 64);
   double lbt_p[64], *lbt = lbt_p;
   int lbt_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lbt, 64);
   double lb_p[64], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(lbt_len, lbt, d1_len, d1, &lb, 64);
   double lct_p[64], *lct = lct_p;
   int lct_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lct, 64);
   double lc_p[64], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(lct_len, lct, d1_len, d1, &lc, 64);
   double lab_p[64], *lab = lab_p;
   int lab_len = o.Gen_Sum_With_PreAlloc(la_len, la, lb_len, lb, &lab, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, lc_len, lc, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab_p != lab) free(lab);
   if (lc_p != lc) free(lc);
   if (lct_p != lct) free(lct);
   if (lb_p != lb) free(lb);
   if (lbt_p != lbt) free(lbt);
   if (la_p != la) free(la);
   if (alift_p != alift) free(alift);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
 }


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_SEEE(implicitPoint2D_SSI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_SEEE_filtered(p1, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_SEEE_interval(p1, pbx, pby, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_SEEE_exact(p1, pbx, pby, pcx, pcy, pdx, pdy);
}

int incircle_indirect_SSEE_filtered(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double pcx, double pcy, double pdx, double pdy)
{
   double l1x, l1y, d1, l2x, l2y, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double pdx1 = pdx * d1;
   double pdy1 = pdy * d1;
   double adx = l1x - pdx1;
   double ady = l1y - pdy1;
   double pdx2 = pdx * d2;
   double pdy2 = pdy * d2;
   double bdx = l2x - pdx2;
   double bdy = l2y - pdy2;
   double cdx = pcx - pdx;
   double cdy = pcy - pdy;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double aliftt = alifta + aliftb;
   double alift = aliftt * d2;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double blift = blifta + bliftb;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double cliftt = clifta + cliftb;
   double clift = cliftt * d2;
   double la = alift * bcdet;
   double lb = blift * cadet;
   double lc = clift * abdet;
   double lab = lc + lb;
   double lab2 = lab * d1;
   double L = lab2 + la;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(pdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pdy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 1.244418457169165e-11;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_SSEE_interval(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, interval_number pcx, interval_number pcy, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number pdx1(pdx * d1);
   interval_number pdy1(pdy * d1);
   interval_number adx(l1x - pdx1);
   interval_number ady(l1y - pdy1);
   interval_number pdx2(pdx * d2);
   interval_number pdy2(pdy * d2);
   interval_number bdx(l2x - pdx2);
   interval_number bdy(l2y - pdy2);
   interval_number cdx(pcx - pdx);
   interval_number cdy(pcy - pdy);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number aliftt(alifta + aliftb);
   interval_number alift(aliftt * d2);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number blift(blifta + bliftb);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number cliftt(clifta + cliftb);
   interval_number clift(cliftt * d2);
   interval_number la(alift * bcdet);
   interval_number lb(blift * cadet);
   interval_number lc(clift * abdet);
   interval_number lab(lc + lb);
   interval_number lab2(lab * d1);
   interval_number L(lab2 + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_SSEE_exact(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double pcx, double pcy, double pdx, double pdy)
{
 double return_value = 0.0;
 double l1x[32], l1y[32], d1[16], l2x[32], l2y[32], d2[16];
 int l1x_len, l1y_len, d1_len, l2x_len, l2y_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double pdx1[32];
   int pdx1_len = o.Gen_Scale(d1_len, d1, pdx, pdx1);
   double pdy1[32];
   int pdy1_len = o.Gen_Scale(d1_len, d1, pdy, pdy1);
   double adx_p[32], *adx = adx_p;
   int adx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
   double ady_p[32], *ady = ady_p;
   int ady_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
   double pdx2[32];
   int pdx2_len = o.Gen_Scale(d2_len, d2, pdx, pdx2);
   double pdy2[32];
   int pdy2_len = o.Gen_Scale(d2_len, d2, pdy, pdy2);
   double bdx_p[32], *bdx = bdx_p;
   int bdx_len = o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
   double bdy_p[32], *bdy = bdy_p;
   int bdy_len = o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
   double cdx[2];
   o.two_Diff(pcx, pdx, cdx);
   double cdy[2];
   o.two_Diff(pcy, pdy, cdy);
   double abdeta_p[32], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
   double abdetb_p[32], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
   double abdet_p[32], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 32);
   double bcdeta_p[32], *bcdeta = bcdeta_p;
   int bcdeta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, 2, cdy, &bcdeta, 32);
   double bcdetb_p[32], *bcdetb = bcdetb_p;
   int bcdetb_len = o.Gen_Product_With_PreAlloc(2, cdx, bdy_len, bdy, &bcdetb, 32);
   double bcdet_p[32], *bcdet = bcdet_p;
   int bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len, bcdetb, &bcdet, 32);
   double cadeta_p[32], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 32);
   double cadetb_p[32], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 32);
   double cadet_p[32], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 32);
   double alifta_p[32], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
   double aliftb_p[32], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
   double aliftt_p[32], *aliftt = aliftt_p;
   int aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &aliftt, 32);
   double alift_p[32], *alift = alift_p;
   int alift_len = o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift, 32);
   double blifta_p[32], *blifta = blifta_p;
   int blifta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
   double bliftb_p[32], *bliftb = bliftb_p;
   int bliftb_len = o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
   double blift_p[32], *blift = blift_p;
   int blift_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len, bliftb, &blift, 32);
   double clifta[8];
   int clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
   double cliftb[8];
   int cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
   double cliftt[16];
   int cliftt_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, cliftt);
   double clift_p[32], *clift = clift_p;
   int clift_len = o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
   double la_p[32], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
   double lb_p[32], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
   double lc_p[32], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
   double lab_p[32], *lab = lab_p;
   int lab_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab, 32);
   double lab2_p[32], *lab2 = lab2_p;
   int lab2_len = o.Gen_Product_With_PreAlloc(lab_len, lab, d1_len, d1, &lab2, 32);
   double L_p[32], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab2_len, lab2, la_len, la, &L, 32);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab2_p != lab2) free(lab2);
   if (lab_p != lab) free(lab);
   if (lc_p != lc) free(lc);
   if (lb_p != lb) free(lb);
   if (la_p != la) free(la);
   if (clift_p != clift) free(clift);
   if (blift_p != blift) free(blift);
   if (bliftb_p != bliftb) free(bliftb);
   if (blifta_p != blifta) free(blifta);
   if (alift_p != alift) free(alift);
   if (aliftt_p != aliftt) free(aliftt);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (bcdet_p != bcdet) free(bcdet);
   if (bcdetb_p != bcdetb) free(bcdetb);
   if (bcdeta_p != bcdeta) free(bcdeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
   if (bdy_p != bdy) free(bdy);
   if (bdx_p != bdx) free(bdx);
   if (ady_p != ady) free(ady);
   if (adx_p != adx) free(adx);
 }


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_SSEE(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double pcx, double pcy, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_SSEE_filtered(p1, p2, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_SSEE_interval(p1, p2, pcx, pcy, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_SSEE_exact(p1, p2, pcx, pcy, pdx, pdy);
}

int incircle_indirect_SSSE_filtered(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, double pdx, double pdy)
{
   double l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, d3, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double pdx1 = pdx * d1;
   double pdy1 = pdy * d1;
   double adx = l1x - pdx1;
   double ady = l1y - pdy1;
   double pdx2 = pdx * d2;
   double pdy2 = pdy * d2;
   double bdx = l2x - pdx2;
   double bdy = l2y - pdy2;
   double pdx3 = pdx * d3;
   double pdy3 = pdy * d3;
   double cdx = l3x - pdx3;
   double cdy = l3y - pdy3;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double aliftt = alifta + aliftb;
   double alift2 = aliftt * d2;
   double alift = alift2 * d3;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double bliftt = blifta + bliftb;
   double blift = bliftt * d3;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double cliftt = clifta + cliftb;
   double clift = cliftt * d2;
   double la = alift * bcdet;
   double lb = blift * cadet;
   double lc = clift * abdet;
   double lab2 = lc + lb;
   double lab = lab2 * d1;
   double L = lab + la;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(pdx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pdy)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 2.609697458133566e-10;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_SSSE_interval(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, interval_number pdx, interval_number pdy)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   || !p3.getIntervalLambda(l3x, l3y, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number pdx1(pdx * d1);
   interval_number pdy1(pdy * d1);
   interval_number adx(l1x - pdx1);
   interval_number ady(l1y - pdy1);
   interval_number pdx2(pdx * d2);
   interval_number pdy2(pdy * d2);
   interval_number bdx(l2x - pdx2);
   interval_number bdy(l2y - pdy2);
   interval_number pdx3(pdx * d3);
   interval_number pdy3(pdy * d3);
   interval_number cdx(l3x - pdx3);
   interval_number cdy(l3y - pdy3);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number aliftt(alifta + aliftb);
   interval_number alift2(aliftt * d2);
   interval_number alift(alift2 * d3);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number bliftt(blifta + bliftb);
   interval_number blift(bliftt * d3);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number cliftt(clifta + cliftb);
   interval_number clift(cliftt * d2);
   interval_number la(alift * bcdet);
   interval_number lb(blift * cadet);
   interval_number lc(clift * abdet);
   interval_number lab2(lc + lb);
   interval_number lab(lab2 * d1);
   interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_SSSE_exact(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, double pdx, double pdy)
{
 double return_value = 0.0;
 double l1x[32], l1y[32], d1[16], l2x[32], l2y[32], d2[16], l3x[32], l3y[32], d3[16];
 int l1x_len, l1y_len, d1_len, l2x_len, l2y_len, d2_len, l3x_len, l3y_len, d3_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, d3, d3_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
 {
   expansionObject o;
   double pdx1[32];
   int pdx1_len = o.Gen_Scale(d1_len, d1, pdx, pdx1);
   double pdy1[32];
   int pdy1_len = o.Gen_Scale(d1_len, d1, pdy, pdy1);
   double adx_p[32], *adx = adx_p;
   int adx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
   double ady_p[32], *ady = ady_p;
   int ady_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
   double pdx2[32];
   int pdx2_len = o.Gen_Scale(d2_len, d2, pdx, pdx2);
   double pdy2[32];
   int pdy2_len = o.Gen_Scale(d2_len, d2, pdy, pdy2);
   double bdx_p[32], *bdx = bdx_p;
   int bdx_len = o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
   double bdy_p[32], *bdy = bdy_p;
   int bdy_len = o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
   double pdx3[32];
   int pdx3_len = o.Gen_Scale(d3_len, d3, pdx, pdx3);
   double pdy3[32];
   int pdy3_len = o.Gen_Scale(d3_len, d3, pdy, pdy3);
   double cdx_p[32], *cdx = cdx_p;
   int cdx_len = o.Gen_Diff_With_PreAlloc(l3x_len, l3x, pdx3_len, pdx3, &cdx, 32);
   double cdy_p[32], *cdy = cdy_p;
   int cdy_len = o.Gen_Diff_With_PreAlloc(l3y_len, l3y, pdy3_len, pdy3, &cdy, 32);
   double abdeta_p[32], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
   double abdetb_p[32], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
   double abdet_p[32], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 32);
   double bcdeta_p[32], *bcdeta = bcdeta_p;
   int bcdeta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 32);
   double bcdetb_p[32], *bcdetb = bcdetb_p;
   int bcdetb_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 32);
   double bcdet_p[32], *bcdet = bcdet_p;
   int bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len, bcdetb, &bcdet, 32);
   double cadeta_p[32], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 32);
   double cadetb_p[32], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 32);
   double cadet_p[32], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 32);
   double alifta_p[32], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
   double aliftb_p[32], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
   double aliftt_p[32], *aliftt = aliftt_p;
   int aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &aliftt, 32);
   double alift2_p[32], *alift2 = alift2_p;
   int alift2_len = o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 32);
   double alift_p[32], *alift = alift_p;
   int alift_len = o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 32);
   double blifta_p[32], *blifta = blifta_p;
   int blifta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
   double bliftb_p[32], *bliftb = bliftb_p;
   int bliftb_len = o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
   double bliftt_p[32], *bliftt = bliftt_p;
   int bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len, bliftb, &bliftt, 32);
   double blift_p[32], *blift = blift_p;
   int blift_len = o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 32);
   double clifta_p[32], *clifta = clifta_p;
   int clifta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 32);
   double cliftb_p[32], *cliftb = cliftb_p;
   int cliftb_len = o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 32);
   double cliftt_p[32], *cliftt = cliftt_p;
   int cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len, cliftb, &cliftt, 32);
   double clift_p[32], *clift = clift_p;
   int clift_len = o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
   double la_p[32], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
   double lb_p[32], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
   double lc_p[32], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
   double lab2_p[32], *lab2 = lab2_p;
   int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 32);
   double lab_p[32], *lab = lab_p;
   int lab_len = o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 32);
   double L_p[32], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 32);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab_p != lab) free(lab);
   if (lab2_p != lab2) free(lab2);
   if (lc_p != lc) free(lc);
   if (lb_p != lb) free(lb);
   if (la_p != la) free(la);
   if (clift_p != clift) free(clift);
   if (cliftt_p != cliftt) free(cliftt);
   if (cliftb_p != cliftb) free(cliftb);
   if (clifta_p != clifta) free(clifta);
   if (blift_p != blift) free(blift);
   if (bliftt_p != bliftt) free(bliftt);
   if (bliftb_p != bliftb) free(bliftb);
   if (blifta_p != blifta) free(blifta);
   if (alift_p != alift) free(alift);
   if (alift2_p != alift2) free(alift2);
   if (aliftt_p != aliftt) free(aliftt);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (bcdet_p != bcdet) free(bcdet);
   if (bcdetb_p != bcdetb) free(bcdetb);
   if (bcdeta_p != bcdeta) free(bcdeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
   if (cdy_p != cdy) free(cdy);
   if (cdx_p != cdx) free(cdx);
   if (bdy_p != bdy) free(bdy);
   if (bdx_p != bdx) free(bdx);
   if (ady_p != ady) free(ady);
   if (adx_p != adx) free(adx);
 }


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_SSSE(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, double pdx, double pdy)
{
   int ret;
   ret = incircle_indirect_SSSE_filtered(p1, p2, p3, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_SSSE_interval(p1, p2, p3, pdx, pdy);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_SSSE_exact(p1, p2, p3, pdx, pdy);
}

int incircle_indirect_SSSS_filtered(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, implicitPoint2D_SSI& p4)
{
   double l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, l4x, l4y, d4, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, d3, max_var)
       || !p4.getFilteredLambda(l4x, l4y, d4, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double l1xt = l1x * d4;
   double l1yt = l1y * d4;
   double l2xt = l2x * d4;
   double l2yt = l2y * d4;
   double l3xt = l3x * d4;
   double l3yt = l3y * d4;
   double l4x1 = l4x * d1;
   double l4y1 = l4y * d1;
   double adx = l1xt - l4x1;
   double ady = l1yt - l4y1;
   double l4x2 = l4x * d2;
   double l4y2 = l4y * d2;
   double bdx = l2xt - l4x2;
   double bdy = l2yt - l4y2;
   double l4x3 = l4x * d3;
   double l4y3 = l4y * d3;
   double cdx = l3xt - l4x3;
   double cdy = l3yt - l4y3;
   double abdeta = adx * bdy;
   double abdetb = bdx * ady;
   double abdet = abdeta - abdetb;
   double bcdeta = bdx * cdy;
   double bcdetb = cdx * bdy;
   double bcdet = bcdeta - bcdetb;
   double cadeta = cdx * ady;
   double cadetb = adx * cdy;
   double cadet = cadeta - cadetb;
   double alifta = adx * adx;
   double aliftb = ady * ady;
   double aliftt = alifta + aliftb;
   double alift2 = aliftt * d2;
   double alift = alift2 * d3;
   double blifta = bdx * bdx;
   double bliftb = bdy * bdy;
   double bliftt = blifta + bliftb;
   double blift = bliftt * d3;
   double clifta = cdx * cdx;
   double cliftb = cdy * cdy;
   double cliftt = clifta + cliftb;
   double clift = cliftt * d2;
   double la = alift * bcdet;
   double lb = blift * cadet;
   double lc = clift * abdet;
   double lab2 = lc + lb;
   double lab = lab2 * d1;
   double L = lab + la;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 2.107594809785969e-08;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int incircle_indirect_SSSS_interval(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, implicitPoint2D_SSI& p4)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, l4x, l4y, d4;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   || !p3.getIntervalLambda(l3x, l3y, d3)
   || !p4.getIntervalLambda(l4x, l4y, d4)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number l1xt(l1x * d4);
   interval_number l1yt(l1y * d4);
   interval_number l2xt(l2x * d4);
   interval_number l2yt(l2y * d4);
   interval_number l3xt(l3x * d4);
   interval_number l3yt(l3y * d4);
   interval_number l4x1(l4x * d1);
   interval_number l4y1(l4y * d1);
   interval_number adx(l1xt - l4x1);
   interval_number ady(l1yt - l4y1);
   interval_number l4x2(l4x * d2);
   interval_number l4y2(l4y * d2);
   interval_number bdx(l2xt - l4x2);
   interval_number bdy(l2yt - l4y2);
   interval_number l4x3(l4x * d3);
   interval_number l4y3(l4y * d3);
   interval_number cdx(l3xt - l4x3);
   interval_number cdy(l3yt - l4y3);
   interval_number abdeta(adx * bdy);
   interval_number abdetb(bdx * ady);
   interval_number abdet(abdeta - abdetb);
   interval_number bcdeta(bdx * cdy);
   interval_number bcdetb(cdx * bdy);
   interval_number bcdet(bcdeta - bcdetb);
   interval_number cadeta(cdx * ady);
   interval_number cadetb(adx * cdy);
   interval_number cadet(cadeta - cadetb);
   interval_number alifta(adx * adx);
   interval_number aliftb(ady * ady);
   interval_number aliftt(alifta + aliftb);
   interval_number alift2(aliftt * d2);
   interval_number alift(alift2 * d3);
   interval_number blifta(bdx * bdx);
   interval_number bliftb(bdy * bdy);
   interval_number bliftt(blifta + bliftb);
   interval_number blift(bliftt * d3);
   interval_number clifta(cdx * cdx);
   interval_number cliftb(cdy * cdy);
   interval_number cliftt(clifta + cliftb);
   interval_number clift(cliftt * d2);
   interval_number la(alift * bcdet);
   interval_number lb(blift * cadet);
   interval_number lc(clift * abdet);
   interval_number lab2(lc + lb);
   interval_number lab(lab2 * d1);
   interval_number L(lab + la);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int incircle_indirect_SSSS_exact(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, implicitPoint2D_SSI& p4)
{
 double return_value = 0.0;
 double l1x[32], l1y[32], d1[16], l2x[32], l2y[32], d2[16], l3x[32], l3y[32], d3[16], l4x[32], l4y[32], d4[16];
 int l1x_len, l1y_len, d1_len, l2x_len, l2y_len, d2_len, l3x_len, l3y_len, d3_len, l4x_len, l4y_len, d4_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, d3, d3_len);
 p4.getExactLambda(l4x, l4x_len, l4y, l4y_len, d4, d4_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) && (d4[d4_len - 1] != 0))
 {
   expansionObject o;
   double l1xt_p[32], *l1xt = l1xt_p;
   int l1xt_len = o.Gen_Product_With_PreAlloc(l1x_len, l1x, d4_len, d4, &l1xt, 32);
   double l1yt_p[32], *l1yt = l1yt_p;
   int l1yt_len = o.Gen_Product_With_PreAlloc(l1y_len, l1y, d4_len, d4, &l1yt, 32);
   double l2xt_p[32], *l2xt = l2xt_p;
   int l2xt_len = o.Gen_Product_With_PreAlloc(l2x_len, l2x, d4_len, d4, &l2xt, 32);
   double l2yt_p[32], *l2yt = l2yt_p;
   int l2yt_len = o.Gen_Product_With_PreAlloc(l2y_len, l2y, d4_len, d4, &l2yt, 32);
   double l3xt_p[32], *l3xt = l3xt_p;
   int l3xt_len = o.Gen_Product_With_PreAlloc(l3x_len, l3x, d4_len, d4, &l3xt, 32);
   double l3yt_p[32], *l3yt = l3yt_p;
   int l3yt_len = o.Gen_Product_With_PreAlloc(l3y_len, l3y, d4_len, d4, &l3yt, 32);
   double l4x1_p[32], *l4x1 = l4x1_p;
   int l4x1_len = o.Gen_Product_With_PreAlloc(l4x_len, l4x, d1_len, d1, &l4x1, 32);
   double l4y1_p[32], *l4y1 = l4y1_p;
   int l4y1_len = o.Gen_Product_With_PreAlloc(l4y_len, l4y, d1_len, d1, &l4y1, 32);
   double adx_p[32], *adx = adx_p;
   int adx_len = o.Gen_Diff_With_PreAlloc(l1xt_len, l1xt, l4x1_len, l4x1, &adx, 32);
   double ady_p[32], *ady = ady_p;
   int ady_len = o.Gen_Diff_With_PreAlloc(l1yt_len, l1yt, l4y1_len, l4y1, &ady, 32);
   double l4x2_p[32], *l4x2 = l4x2_p;
   int l4x2_len = o.Gen_Product_With_PreAlloc(l4x_len, l4x, d2_len, d2, &l4x2, 32);
   double l4y2_p[32], *l4y2 = l4y2_p;
   int l4y2_len = o.Gen_Product_With_PreAlloc(l4y_len, l4y, d2_len, d2, &l4y2, 32);
   double bdx_p[32], *bdx = bdx_p;
   int bdx_len = o.Gen_Diff_With_PreAlloc(l2xt_len, l2xt, l4x2_len, l4x2, &bdx, 32);
   double bdy_p[32], *bdy = bdy_p;
   int bdy_len = o.Gen_Diff_With_PreAlloc(l2yt_len, l2yt, l4y2_len, l4y2, &bdy, 32);
   double l4x3_p[32], *l4x3 = l4x3_p;
   int l4x3_len = o.Gen_Product_With_PreAlloc(l4x_len, l4x, d3_len, d3, &l4x3, 32);
   double l4y3_p[32], *l4y3 = l4y3_p;
   int l4y3_len = o.Gen_Product_With_PreAlloc(l4y_len, l4y, d3_len, d3, &l4y3, 32);
   double cdx_p[32], *cdx = cdx_p;
   int cdx_len = o.Gen_Diff_With_PreAlloc(l3xt_len, l3xt, l4x3_len, l4x3, &cdx, 32);
   double cdy_p[32], *cdy = cdy_p;
   int cdy_len = o.Gen_Diff_With_PreAlloc(l3yt_len, l3yt, l4y3_len, l4y3, &cdy, 32);
   double abdeta_p[32], *abdeta = abdeta_p;
   int abdeta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
   double abdetb_p[32], *abdetb = abdetb_p;
   int abdetb_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
   double abdet_p[32], *abdet = abdet_p;
   int abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len, abdetb, &abdet, 32);
   double bcdeta_p[32], *bcdeta = bcdeta_p;
   int bcdeta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 32);
   double bcdetb_p[32], *bcdetb = bcdetb_p;
   int bcdetb_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 32);
   double bcdet_p[32], *bcdet = bcdet_p;
   int bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len, bcdetb, &bcdet, 32);
   double cadeta_p[32], *cadeta = cadeta_p;
   int cadeta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 32);
   double cadetb_p[32], *cadetb = cadetb_p;
   int cadetb_len = o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 32);
   double cadet_p[32], *cadet = cadet_p;
   int cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len, cadetb, &cadet, 32);
   double alifta_p[32], *alifta = alifta_p;
   int alifta_len = o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
   double aliftb_p[32], *aliftb = aliftb_p;
   int aliftb_len = o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
   double aliftt_p[32], *aliftt = aliftt_p;
   int aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len, aliftb, &aliftt, 32);
   double alift2_p[32], *alift2 = alift2_p;
   int alift2_len = o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 32);
   double alift_p[32], *alift = alift_p;
   int alift_len = o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 32);
   double blifta_p[32], *blifta = blifta_p;
   int blifta_len = o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
   double bliftb_p[32], *bliftb = bliftb_p;
   int bliftb_len = o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
   double bliftt_p[32], *bliftt = bliftt_p;
   int bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len, bliftb, &bliftt, 32);
   double blift_p[32], *blift = blift_p;
   int blift_len = o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 32);
   double clifta_p[32], *clifta = clifta_p;
   int clifta_len = o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 32);
   double cliftb_p[32], *cliftb = cliftb_p;
   int cliftb_len = o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 32);
   double cliftt_p[32], *cliftt = cliftt_p;
   int cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len, cliftb, &cliftt, 32);
   double clift_p[32], *clift = clift_p;
   int clift_len = o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
   double la_p[32], *la = la_p;
   int la_len = o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
   double lb_p[32], *lb = lb_p;
   int lb_len = o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
   double lc_p[32], *lc = lc_p;
   int lc_len = o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
   double lab2_p[32], *lab2 = lab2_p;
   int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 32);
   double lab_p[32], *lab = lab_p;
   int lab_len = o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 32);
   double L_p[32], *L = L_p;
   int L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 32);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (lab_p != lab) free(lab);
   if (lab2_p != lab2) free(lab2);
   if (lc_p != lc) free(lc);
   if (lb_p != lb) free(lb);
   if (la_p != la) free(la);
   if (clift_p != clift) free(clift);
   if (cliftt_p != cliftt) free(cliftt);
   if (cliftb_p != cliftb) free(cliftb);
   if (clifta_p != clifta) free(clifta);
   if (blift_p != blift) free(blift);
   if (bliftt_p != bliftt) free(bliftt);
   if (bliftb_p != bliftb) free(bliftb);
   if (blifta_p != blifta) free(blifta);
   if (alift_p != alift) free(alift);
   if (alift2_p != alift2) free(alift2);
   if (aliftt_p != aliftt) free(aliftt);
   if (aliftb_p != aliftb) free(aliftb);
   if (alifta_p != alifta) free(alifta);
   if (cadet_p != cadet) free(cadet);
   if (cadetb_p != cadetb) free(cadetb);
   if (cadeta_p != cadeta) free(cadeta);
   if (bcdet_p != bcdet) free(bcdet);
   if (bcdetb_p != bcdetb) free(bcdetb);
   if (bcdeta_p != bcdeta) free(bcdeta);
   if (abdet_p != abdet) free(abdet);
   if (abdetb_p != abdetb) free(abdetb);
   if (abdeta_p != abdeta) free(abdeta);
   if (cdy_p != cdy) free(cdy);
   if (cdx_p != cdx) free(cdx);
   if (l4y3_p != l4y3) free(l4y3);
   if (l4x3_p != l4x3) free(l4x3);
   if (bdy_p != bdy) free(bdy);
   if (bdx_p != bdx) free(bdx);
   if (l4y2_p != l4y2) free(l4y2);
   if (l4x2_p != l4x2) free(l4x2);
   if (ady_p != ady) free(ady);
   if (adx_p != adx) free(adx);
   if (l4y1_p != l4y1) free(l4y1);
   if (l4x1_p != l4x1) free(l4x1);
   if (l3yt_p != l3yt) free(l3yt);
   if (l3xt_p != l3xt) free(l3xt);
   if (l2yt_p != l2yt) free(l2yt);
   if (l2xt_p != l2xt) free(l2xt);
   if (l1yt_p != l1yt) free(l1yt);
   if (l1xt_p != l1xt) free(l1xt);
 }


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int incircle_indirect_SSSS(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, implicitPoint2D_SSI& p4)
{
   int ret;
//   ret = incircle_indirect_SSSS_filtered(p1, p2, p3, p4);
//   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = incircle_indirect_SSSS_interval(p1, p2, p3, p4);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return incircle_indirect_SSSS_exact(p1, p2, p3, p4);
}

bool lambda2d_SSI_filtered(double ea1x, double ea1y, double ea2x, double ea2y, double eb1x, double eb1y, double eb2x, double eb2y, double& lambda_x, double& lambda_y, double& lambda_det, double& max_var)
{
   double t1a = ea1x * ea2y;
   double t1b = ea2x * ea1y;
   double t1 = t1a - t1b;
   double tx2 = eb1x - eb2x;
   double t3a = eb1x * eb2y;
   double t3b = eb2x * eb1y;
   double t3 = t3a - t3b;
   double tx4 = ea1x - ea2x;
   double ty2 = eb1y - eb2y;
   double ty4 = ea1y - ea2y;
   double lxa = t1 * tx2;
   double lxb = t3 * tx4;
   lambda_x = lxa - lxb;
   double lya = t1 * ty2;
   double lyb = t3 * ty4;
   lambda_y = lya - lyb;
   double deta = tx4 * ty2;
   double detb = tx2 * ty4;
   lambda_det = deta - detb;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(ea1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ea1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ea2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ea2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(eb1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(eb1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(eb2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(eb2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(tx2)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(tx4)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ty2)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ty4)) > max_var) max_var = _tmp_fabs;
   double lambda_det_eps = max_var;
   lambda_det_eps *= lambda_det_eps;
   lambda_det_eps *= 8.881784197001252e-16;

   return ( (lambda_det > lambda_det_eps || lambda_det < -lambda_det_eps) );
}

bool lambda2d_SSI_interval(interval_number ea1x, interval_number ea1y, interval_number ea2x, interval_number ea2y, interval_number eb1x, interval_number eb1y, interval_number eb2x, interval_number eb2y, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_det)
{
   setFPUModeToRoundUP();
   interval_number t1a(ea1x * ea2y);
   interval_number t1b(ea2x * ea1y);
   interval_number t1(t1a - t1b);
   interval_number tx2(eb1x - eb2x);
   interval_number t3a(eb1x * eb2y);
   interval_number t3b(eb2x * eb1y);
   interval_number t3(t3a - t3b);
   interval_number tx4(ea1x - ea2x);
   interval_number ty2(eb1y - eb2y);
   interval_number ty4(ea1y - ea2y);
   interval_number lxa(t1 * tx2);
   interval_number lxb(t3 * tx4);
   lambda_x = lxa - lxb;
   interval_number lya(t1 * ty2);
   interval_number lyb(t3 * ty4);
   lambda_y = lya - lyb;
   interval_number deta(tx4 * ty2);
   interval_number detb(tx2 * ty4);
   lambda_det = deta - detb;
   setFPUModeToRoundNEAR();


   return (
      lambda_det.signIsReliable()
   );
}

void lambda2d_SSI_exact(double ea1x, double ea1y, double ea2x, double ea2y, double eb1x, double eb1y, double eb2x, double eb2y, double *lambda_x, int& lambda_x_len, double *lambda_y, int& lambda_y_len, double *lambda_det, int& lambda_det_len)
{
   expansionObject o;
   double t1a[2];
   o.Two_Prod(ea1x, ea2y, t1a);
   double t1b[2];
   o.Two_Prod(ea2x, ea1y, t1b);
   double t1[4];
   o.Two_Two_Diff(t1a, t1b, t1);
   double tx2[2];
   o.two_Diff(eb1x, eb2x, tx2);
   double t3a[2];
   o.Two_Prod(eb1x, eb2y, t3a);
   double t3b[2];
   o.Two_Prod(eb2x, eb1y, t3b);
   double t3[4];
   o.Two_Two_Diff(t3a, t3b, t3);
   double tx4[2];
   o.two_Diff(ea1x, ea2x, tx4);
   double ty2[2];
   o.two_Diff(eb1y, eb2y, ty2);
   double ty4[2];
   o.two_Diff(ea1y, ea2y, ty4);
   double lxa[16];
   int lxa_len = o.Gen_Product(4, t1, 2, tx2, lxa);
   double lxb[16];
   int lxb_len = o.Gen_Product(4, t3, 2, tx4, lxb);
   lambda_x_len = o.Gen_Diff(lxa_len, lxa, lxb_len, lxb, lambda_x);
   double lya[16];
   int lya_len = o.Gen_Product(4, t1, 2, ty2, lya);
   double lyb[16];
   int lyb_len = o.Gen_Product(4, t3, 2, ty4, lyb);
   lambda_y_len = o.Gen_Diff(lya_len, lya, lyb_len, lyb, lambda_y);
   double deta[8];
   int deta_len = o.Gen_Product(2, tx4, 2, ty2, deta);
   double detb[8];
   int detb_len = o.Gen_Product(2, tx2, 2, ty4, detb);
   lambda_det_len = o.Gen_Diff(deta_len, deta, detb_len, detb, lambda_det);

}

bool lambda3d_LPI_filtered(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz, double& lambda_d, double& lambda_x, double& lambda_y, double& lambda_z, double& max_var)
{
   double a11 = px - qx;
   double a12 = py - qy;
   double a13 = pz - qz;
   double a21 = sx - rx;
   double a22 = sy - ry;
   double a23 = sz - rz;
   double a31 = tx - rx;
   double a32 = ty - ry;
   double a33 = tz - rz;
   double tv1 = a22 * a33;
   double tv2 = a23 * a32;
   double a2233 = tv1 - tv2;
   double tv3 = a21 * a33;
   double tv4 = a23 * a31;
   double a2133 = tv3 - tv4;
   double tv5 = a21 * a32;
   double tv6 = a22 * a31;
   double a2132 = tv5 - tv6;
   double tv7 = a11 * a2233;
   double tv8 = a12 * a2133;
   double tv9 = a13 * a2132;
   double tt1 = tv7 - tv8;
   lambda_d = tt1 + tv9;
   double px_rx = px - rx;
   double py_ry = py - ry;
   double pz_rz = pz - rz;
   double tt2 = py_ry * a2133;
   double tt3 = px_rx * a2233;
   double tt4 = pz_rz * a2132;
   double tt5 = tt3 + tt4;
   double n = tt5 - tt2;
   double ax = a11 * n;
   double ay = a12 * n;
   double az = a13 * n;
   double dpx = lambda_d * px;
   double dpy = lambda_d * py;
   double dpz = lambda_d * pz;
   lambda_x = dpx - ax;
   lambda_y = dpy - ay;
   lambda_z = dpz - az;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(px)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(py)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a11)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a12)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a13)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a21)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a22)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a23)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a31)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a32)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(a33)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(px_rx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(py_ry)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(pz_rz)) > max_var) max_var = _tmp_fabs;
   double lambda_d_eps = max_var;
   lambda_d_eps *= lambda_d_eps;
   lambda_d_eps *= max_var;
   lambda_d_eps *= 4.884981308350689e-15;

   return ( (lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps) );
}

bool lambda3d_LPI_interval(interval_number px, interval_number py, interval_number pz, interval_number qx, interval_number qy, interval_number qz, interval_number rx, interval_number ry, interval_number rz, interval_number sx, interval_number sy, interval_number sz, interval_number tx, interval_number ty, interval_number tz, interval_number& lambda_d, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z)
{
   setFPUModeToRoundUP();
   interval_number a11(px - qx);
   interval_number a12(py - qy);
   interval_number a13(pz - qz);
   interval_number a21(sx - rx);
   interval_number a22(sy - ry);
   interval_number a23(sz - rz);
   interval_number a31(tx - rx);
   interval_number a32(ty - ry);
   interval_number a33(tz - rz);
   interval_number tv1(a22 * a33);
   interval_number tv2(a23 * a32);
   interval_number a2233(tv1 - tv2);
   interval_number tv3(a21 * a33);
   interval_number tv4(a23 * a31);
   interval_number a2133(tv3 - tv4);
   interval_number tv5(a21 * a32);
   interval_number tv6(a22 * a31);
   interval_number a2132(tv5 - tv6);
   interval_number tv7(a11 * a2233);
   interval_number tv8(a12 * a2133);
   interval_number tv9(a13 * a2132);
   interval_number tt1(tv7 - tv8);
   lambda_d = tt1 + tv9;
   interval_number px_rx(px - rx);
   interval_number py_ry(py - ry);
   interval_number pz_rz(pz - rz);
   interval_number tt2(py_ry * a2133);
   interval_number tt3(px_rx * a2233);
   interval_number tt4(pz_rz * a2132);
   interval_number tt5(tt3 + tt4);
   interval_number n(tt5 - tt2);
   interval_number ax(a11 * n);
   interval_number ay(a12 * n);
   interval_number az(a13 * n);
   interval_number dpx(lambda_d * px);
   interval_number dpy(lambda_d * py);
   interval_number dpz(lambda_d * pz);
   lambda_x = dpx - ax;
   lambda_y = dpy - ay;
   lambda_z = dpz - az;
   setFPUModeToRoundNEAR();


   return (
      lambda_d.signIsReliable()
   );
}

void lambda3d_LPI_exact(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz, double *lambda_d, int& lambda_d_len, double *lambda_x, int& lambda_x_len, double *lambda_y, int& lambda_y_len, double *lambda_z, int& lambda_z_len)
{
   expansionObject o;
   double a11[2];
   o.two_Diff(px, qx, a11);
   double a12[2];
   o.two_Diff(py, qy, a12);
   double a13[2];
   o.two_Diff(pz, qz, a13);
   double a21[2];
   o.two_Diff(sx, rx, a21);
   double a22[2];
   o.two_Diff(sy, ry, a22);
   double a23[2];
   o.two_Diff(sz, rz, a23);
   double a31[2];
   o.two_Diff(tx, rx, a31);
   double a32[2];
   o.two_Diff(ty, ry, a32);
   double a33[2];
   o.two_Diff(tz, rz, a33);
   double tv1[8];
   int tv1_len = o.Gen_Product(2, a22, 2, a33, tv1);
   double tv2[8];
   int tv2_len = o.Gen_Product(2, a23, 2, a32, tv2);
   double a2233[16];
   int a2233_len = o.Gen_Diff(tv1_len, tv1, tv2_len, tv2, a2233);
   double tv3[8];
   int tv3_len = o.Gen_Product(2, a21, 2, a33, tv3);
   double tv4[8];
   int tv4_len = o.Gen_Product(2, a23, 2, a31, tv4);
   double a2133[16];
   int a2133_len = o.Gen_Diff(tv3_len, tv3, tv4_len, tv4, a2133);
   double tv5[8];
   int tv5_len = o.Gen_Product(2, a21, 2, a32, tv5);
   double tv6[8];
   int tv6_len = o.Gen_Product(2, a22, 2, a31, tv6);
   double a2132[16];
   int a2132_len = o.Gen_Diff(tv5_len, tv5, tv6_len, tv6, a2132);
   double tv7[64];
   int tv7_len = o.Gen_Product(2, a11, a2233_len, a2233, tv7);
   double tv8[64];
   int tv8_len = o.Gen_Product(2, a12, a2133_len, a2133, tv8);
   double tv9[64];
   int tv9_len = o.Gen_Product(2, a13, a2132_len, a2132, tv9);
   double tt1[128];
   int tt1_len = o.Gen_Diff(tv7_len, tv7, tv8_len, tv8, tt1);
   lambda_d_len = o.Gen_Sum_With_PreAlloc(tt1_len, tt1, tv9_len, tv9, &lambda_d, 128);
   double px_rx[2];
   o.two_Diff(px, rx, px_rx);
   double py_ry[2];
   o.two_Diff(py, ry, py_ry);
   double pz_rz[2];
   o.two_Diff(pz, rz, pz_rz);
   double tt2[64];
   int tt2_len = o.Gen_Product(2, py_ry, a2133_len, a2133, tt2);
   double tt3[64];
   int tt3_len = o.Gen_Product(2, px_rx, a2233_len, a2233, tt3);
   double tt4[64];
   int tt4_len = o.Gen_Product(2, pz_rz, a2132_len, a2132, tt4);
   double tt5[128];
   int tt5_len = o.Gen_Sum(tt3_len, tt3, tt4_len, tt4, tt5);
   double n_p[128], *n = n_p;
   int n_len = o.Gen_Diff_With_PreAlloc(tt5_len, tt5, tt2_len, tt2, &n, 128);
   double ax_p[128], *ax = ax_p;
   int ax_len = o.Gen_Product_With_PreAlloc(2, a11, n_len, n, &ax, 128);
   double ay_p[128], *ay = ay_p;
   int ay_len = o.Gen_Product_With_PreAlloc(2, a12, n_len, n, &ay, 128);
   double az_p[128], *az = az_p;
   int az_len = o.Gen_Product_With_PreAlloc(2, a13, n_len, n, &az, 128);
   double dpx_p[128], *dpx = dpx_p;
   int dpx_len = o.Gen_Scale_With_PreAlloc(lambda_d_len, lambda_d, px, &dpx, 128);
   double dpy_p[128], *dpy = dpy_p;
   int dpy_len = o.Gen_Scale_With_PreAlloc(lambda_d_len, lambda_d, py, &dpy, 128);
   double dpz_p[128], *dpz = dpz_p;
   int dpz_len = o.Gen_Scale_With_PreAlloc(lambda_d_len, lambda_d, pz, &dpz, 128);
   lambda_x_len = o.Gen_Diff_With_PreAlloc(dpx_len, dpx, ax_len, ax, &lambda_x, 128);
   lambda_y_len = o.Gen_Diff_With_PreAlloc(dpy_len, dpy, ay_len, ay, &lambda_y, 128);
   lambda_z_len = o.Gen_Diff_With_PreAlloc(dpz_len, dpz, az_len, az, &lambda_z, 128);

   if (dpz_p != dpz) free(dpz);
   if (dpy_p != dpy) free(dpy);
   if (dpx_p != dpx) free(dpx);
   if (az_p != az) free(az);
   if (ay_p != ay) free(ay);
   if (ax_p != ax) free(ax);
   if (n_p != n) free(n);
}

bool lambda3d_TPI_filtered(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z, double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z, double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z, double ou3x, double ou3y, double ou3z, double& lambda_x, double& lambda_y, double& lambda_z, double& lambda_d, double& max_var)
{
   double v3x = ov3x - ov2x;
   double v3y = ov3y - ov2y;
   double v3z = ov3z - ov2z;
   double v2x = ov2x - ov1x;
   double v2y = ov2y - ov1y;
   double v2z = ov2z - ov1z;
   double w3x = ow3x - ow2x;
   double w3y = ow3y - ow2y;
   double w3z = ow3z - ow2z;
   double w2x = ow2x - ow1x;
   double w2y = ow2y - ow1y;
   double w2z = ow2z - ow1z;
   double u3x = ou3x - ou2x;
   double u3y = ou3y - ou2y;
   double u3z = ou3z - ou2z;
   double u2x = ou2x - ou1x;
   double u2y = ou2y - ou1y;
   double u2z = ou2z - ou1z;
   double nvx1 = v2y * v3z;
   double nvx2 = v2z * v3y;
   double nvx = nvx1 - nvx2;
   double nvy1 = v3x * v2z;
   double nvy2 = v3z * v2x;
   double nvy = nvy1 - nvy2;
   double nvz1 = v2x * v3y;
   double nvz2 = v2y * v3x;
   double nvz = nvz1 - nvz2;
   double nwx1 = w2y * w3z;
   double nwx2 = w2z * w3y;
   double nwx = nwx1 - nwx2;
   double nwy1 = w3x * w2z;
   double nwy2 = w3z * w2x;
   double nwy = nwy1 - nwy2;
   double nwz1 = w2x * w3y;
   double nwz2 = w2y * w3x;
   double nwz = nwz1 - nwz2;
   double nux1 = u2y * u3z;
   double nux2 = u2z * u3y;
   double nux = nux1 - nux2;
   double nuy1 = u3x * u2z;
   double nuy2 = u3z * u2x;
   double nuy = nuy1 - nuy2;
   double nuz1 = u2x * u3y;
   double nuz2 = u2y * u3x;
   double nuz = nuz1 - nuz2;
   double nwyuz1 = nwy * nuz;
   double nwyuz2 = nwz * nuy;
   double nwyuz = nwyuz1 - nwyuz2;
   double nwxuz1 = nwx * nuz;
   double nwxuz2 = nwz * nux;
   double nwxuz = nwxuz1 - nwxuz2;
   double nwxuy1 = nwx * nuy;
   double nwxuy2 = nwy * nux;
   double nwxuy = nwxuy1 - nwxuy2;
   double nvyuz1 = nvy * nuz;
   double nvyuz2 = nvz * nuy;
   double nvyuz = nvyuz1 - nvyuz2;
   double nvxuz1 = nvx * nuz;
   double nvxuz2 = nvz * nux;
   double nvxuz = nvxuz1 - nvxuz2;
   double nvxuy1 = nvx * nuy;
   double nvxuy2 = nvy * nux;
   double nvxuy = nvxuy1 - nvxuy2;
   double nvywz1 = nvy * nwz;
   double nvywz2 = nvz * nwy;
   double nvywz = nvywz1 - nvywz2;
   double nvxwz1 = nvx * nwz;
   double nvxwz2 = nvz * nwx;
   double nvxwz = nvxwz1 - nvxwz2;
   double nvxwy1 = nvx * nwy;
   double nvxwy2 = nvy * nwx;
   double nvxwy = nvxwy1 - nvxwy2;
   double p1a = nvx * ov1x;
   double p1b = nvy * ov1y;
   double p1c = nvz * ov1z;
   double p1ab = p1a + p1b;
   double p1 = p1ab + p1c;
   double p2a = nwx * ow1x;
   double p2b = nwy * ow1y;
   double p2c = nwz * ow1z;
   double p2ab = p2a + p2b;
   double p2 = p2ab + p2c;
   double p3a = nux * ou1x;
   double p3b = nuy * ou1y;
   double p3c = nuz * ou1z;
   double p3ab = p3a + p3b;
   double p3 = p3ab + p3c;
   double lxa = p1 * nwyuz;
   double lxb = p3 * nvywz;
   double lxc = p2 * nvyuz;
   double lxab = lxa + lxb;
   lambda_x = lxab - lxc;
   double lya = p2 * nvxuz;
   double lyb = p3 * nvxwz;
   double lyc = p1 * nwxuz;
   double lybc = lyc + lyb;
   lambda_y = lya - lybc;
   double lza = p3 * nvxwy;
   double lzb = p1 * nwxuy;
   double lzc = p2 * nvxuy;
   double lzab = lza + lzb;
   lambda_z = lzab - lzc;
   double da = nvx * nwyuz;
   double db = nvz * nwxuy;
   double dc = nvy * nwxuz;
   double dab = da + db;
   lambda_d = dab - dc;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(ov1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ov1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ov1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ow1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ow1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ow1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ou1x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ou1y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ou1z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(v3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(v3y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(v3z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(v2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(v2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(v2z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(w3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(w3y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(w3z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(w2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(w2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(w2z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(u3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(u3y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(u3z)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(u2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(u2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(u2z)) > max_var) max_var = _tmp_fabs;
   double lambda_d_eps = max_var;
   lambda_d_eps *= lambda_d_eps;
   lambda_d_eps *= lambda_d_eps;
   lambda_d_eps *= max_var;
   lambda_d_eps *= max_var;
   lambda_d_eps *= 8.704148513061234e-14;

   return ( (lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps) );
}

bool lambda3d_TPI_interval(interval_number ov1x, interval_number ov1y, interval_number ov1z, interval_number ov2x, interval_number ov2y, interval_number ov2z, interval_number ov3x, interval_number ov3y, interval_number ov3z, interval_number ow1x, interval_number ow1y, interval_number ow1z, interval_number ow2x, interval_number ow2y, interval_number ow2z, interval_number ow3x, interval_number ow3y, interval_number ow3z, interval_number ou1x, interval_number ou1y, interval_number ou1z, interval_number ou2x, interval_number ou2y, interval_number ou2z, interval_number ou3x, interval_number ou3y, interval_number ou3z, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z, interval_number& lambda_d)
{
   setFPUModeToRoundUP();
   interval_number v3x(ov3x - ov2x);
   interval_number v3y(ov3y - ov2y);
   interval_number v3z(ov3z - ov2z);
   interval_number v2x(ov2x - ov1x);
   interval_number v2y(ov2y - ov1y);
   interval_number v2z(ov2z - ov1z);
   interval_number w3x(ow3x - ow2x);
   interval_number w3y(ow3y - ow2y);
   interval_number w3z(ow3z - ow2z);
   interval_number w2x(ow2x - ow1x);
   interval_number w2y(ow2y - ow1y);
   interval_number w2z(ow2z - ow1z);
   interval_number u3x(ou3x - ou2x);
   interval_number u3y(ou3y - ou2y);
   interval_number u3z(ou3z - ou2z);
   interval_number u2x(ou2x - ou1x);
   interval_number u2y(ou2y - ou1y);
   interval_number u2z(ou2z - ou1z);
   interval_number nvx1(v2y * v3z);
   interval_number nvx2(v2z * v3y);
   interval_number nvx(nvx1 - nvx2);
   interval_number nvy1(v3x * v2z);
   interval_number nvy2(v3z * v2x);
   interval_number nvy(nvy1 - nvy2);
   interval_number nvz1(v2x * v3y);
   interval_number nvz2(v2y * v3x);
   interval_number nvz(nvz1 - nvz2);
   interval_number nwx1(w2y * w3z);
   interval_number nwx2(w2z * w3y);
   interval_number nwx(nwx1 - nwx2);
   interval_number nwy1(w3x * w2z);
   interval_number nwy2(w3z * w2x);
   interval_number nwy(nwy1 - nwy2);
   interval_number nwz1(w2x * w3y);
   interval_number nwz2(w2y * w3x);
   interval_number nwz(nwz1 - nwz2);
   interval_number nux1(u2y * u3z);
   interval_number nux2(u2z * u3y);
   interval_number nux(nux1 - nux2);
   interval_number nuy1(u3x * u2z);
   interval_number nuy2(u3z * u2x);
   interval_number nuy(nuy1 - nuy2);
   interval_number nuz1(u2x * u3y);
   interval_number nuz2(u2y * u3x);
   interval_number nuz(nuz1 - nuz2);
   interval_number nwyuz1(nwy * nuz);
   interval_number nwyuz2(nwz * nuy);
   interval_number nwyuz(nwyuz1 - nwyuz2);
   interval_number nwxuz1(nwx * nuz);
   interval_number nwxuz2(nwz * nux);
   interval_number nwxuz(nwxuz1 - nwxuz2);
   interval_number nwxuy1(nwx * nuy);
   interval_number nwxuy2(nwy * nux);
   interval_number nwxuy(nwxuy1 - nwxuy2);
   interval_number nvyuz1(nvy * nuz);
   interval_number nvyuz2(nvz * nuy);
   interval_number nvyuz(nvyuz1 - nvyuz2);
   interval_number nvxuz1(nvx * nuz);
   interval_number nvxuz2(nvz * nux);
   interval_number nvxuz(nvxuz1 - nvxuz2);
   interval_number nvxuy1(nvx * nuy);
   interval_number nvxuy2(nvy * nux);
   interval_number nvxuy(nvxuy1 - nvxuy2);
   interval_number nvywz1(nvy * nwz);
   interval_number nvywz2(nvz * nwy);
   interval_number nvywz(nvywz1 - nvywz2);
   interval_number nvxwz1(nvx * nwz);
   interval_number nvxwz2(nvz * nwx);
   interval_number nvxwz(nvxwz1 - nvxwz2);
   interval_number nvxwy1(nvx * nwy);
   interval_number nvxwy2(nvy * nwx);
   interval_number nvxwy(nvxwy1 - nvxwy2);
   interval_number p1a(nvx * ov1x);
   interval_number p1b(nvy * ov1y);
   interval_number p1c(nvz * ov1z);
   interval_number p1ab(p1a + p1b);
   interval_number p1(p1ab + p1c);
   interval_number p2a(nwx * ow1x);
   interval_number p2b(nwy * ow1y);
   interval_number p2c(nwz * ow1z);
   interval_number p2ab(p2a + p2b);
   interval_number p2(p2ab + p2c);
   interval_number p3a(nux * ou1x);
   interval_number p3b(nuy * ou1y);
   interval_number p3c(nuz * ou1z);
   interval_number p3ab(p3a + p3b);
   interval_number p3(p3ab + p3c);
   interval_number lxa(p1 * nwyuz);
   interval_number lxb(p3 * nvywz);
   interval_number lxc(p2 * nvyuz);
   interval_number lxab(lxa + lxb);
   lambda_x = lxab - lxc;
   interval_number lya(p2 * nvxuz);
   interval_number lyb(p3 * nvxwz);
   interval_number lyc(p1 * nwxuz);
   interval_number lybc(lyc + lyb);
   lambda_y = lya - lybc;
   interval_number lza(p3 * nvxwy);
   interval_number lzb(p1 * nwxuy);
   interval_number lzc(p2 * nvxuy);
   interval_number lzab(lza + lzb);
   lambda_z = lzab - lzc;
   interval_number da(nvx * nwyuz);
   interval_number db(nvz * nwxuy);
   interval_number dc(nvy * nwxuz);
   interval_number dab(da + db);
   lambda_d = dab - dc;
   setFPUModeToRoundNEAR();


   return (
      lambda_d.signIsReliable()
   );
}

void lambda3d_TPI_exact(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z, double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z, double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z, double ou3x, double ou3y, double ou3z, double *lambda_x, int& lambda_x_len, double *lambda_y, int& lambda_y_len, double *lambda_z, int& lambda_z_len, double *lambda_d, int& lambda_d_len)
{
   expansionObject o;
   double v3x[2];
   o.two_Diff(ov3x, ov2x, v3x);
   double v3y[2];
   o.two_Diff(ov3y, ov2y, v3y);
   double v3z[2];
   o.two_Diff(ov3z, ov2z, v3z);
   double v2x[2];
   o.two_Diff(ov2x, ov1x, v2x);
   double v2y[2];
   o.two_Diff(ov2y, ov1y, v2y);
   double v2z[2];
   o.two_Diff(ov2z, ov1z, v2z);
   double w3x[2];
   o.two_Diff(ow3x, ow2x, w3x);
   double w3y[2];
   o.two_Diff(ow3y, ow2y, w3y);
   double w3z[2];
   o.two_Diff(ow3z, ow2z, w3z);
   double w2x[2];
   o.two_Diff(ow2x, ow1x, w2x);
   double w2y[2];
   o.two_Diff(ow2y, ow1y, w2y);
   double w2z[2];
   o.two_Diff(ow2z, ow1z, w2z);
   double u3x[2];
   o.two_Diff(ou3x, ou2x, u3x);
   double u3y[2];
   o.two_Diff(ou3y, ou2y, u3y);
   double u3z[2];
   o.two_Diff(ou3z, ou2z, u3z);
   double u2x[2];
   o.two_Diff(ou2x, ou1x, u2x);
   double u2y[2];
   o.two_Diff(ou2y, ou1y, u2y);
   double u2z[2];
   o.two_Diff(ou2z, ou1z, u2z);
   double nvx1[8];
   int nvx1_len = o.Gen_Product(2, v2y, 2, v3z, nvx1);
   double nvx2[8];
   int nvx2_len = o.Gen_Product(2, v2z, 2, v3y, nvx2);
   double nvx[16];
   int nvx_len = o.Gen_Diff(nvx1_len, nvx1, nvx2_len, nvx2, nvx);
   double nvy1[8];
   int nvy1_len = o.Gen_Product(2, v3x, 2, v2z, nvy1);
   double nvy2[8];
   int nvy2_len = o.Gen_Product(2, v3z, 2, v2x, nvy2);
   double nvy[16];
   int nvy_len = o.Gen_Diff(nvy1_len, nvy1, nvy2_len, nvy2, nvy);
   double nvz1[8];
   int nvz1_len = o.Gen_Product(2, v2x, 2, v3y, nvz1);
   double nvz2[8];
   int nvz2_len = o.Gen_Product(2, v2y, 2, v3x, nvz2);
   double nvz[16];
   int nvz_len = o.Gen_Diff(nvz1_len, nvz1, nvz2_len, nvz2, nvz);
   double nwx1[8];
   int nwx1_len = o.Gen_Product(2, w2y, 2, w3z, nwx1);
   double nwx2[8];
   int nwx2_len = o.Gen_Product(2, w2z, 2, w3y, nwx2);
   double nwx[16];
   int nwx_len = o.Gen_Diff(nwx1_len, nwx1, nwx2_len, nwx2, nwx);
   double nwy1[8];
   int nwy1_len = o.Gen_Product(2, w3x, 2, w2z, nwy1);
   double nwy2[8];
   int nwy2_len = o.Gen_Product(2, w3z, 2, w2x, nwy2);
   double nwy[16];
   int nwy_len = o.Gen_Diff(nwy1_len, nwy1, nwy2_len, nwy2, nwy);
   double nwz1[8];
   int nwz1_len = o.Gen_Product(2, w2x, 2, w3y, nwz1);
   double nwz2[8];
   int nwz2_len = o.Gen_Product(2, w2y, 2, w3x, nwz2);
   double nwz[16];
   int nwz_len = o.Gen_Diff(nwz1_len, nwz1, nwz2_len, nwz2, nwz);
   double nux1[8];
   int nux1_len = o.Gen_Product(2, u2y, 2, u3z, nux1);
   double nux2[8];
   int nux2_len = o.Gen_Product(2, u2z, 2, u3y, nux2);
   double nux[16];
   int nux_len = o.Gen_Diff(nux1_len, nux1, nux2_len, nux2, nux);
   double nuy1[8];
   int nuy1_len = o.Gen_Product(2, u3x, 2, u2z, nuy1);
   double nuy2[8];
   int nuy2_len = o.Gen_Product(2, u3z, 2, u2x, nuy2);
   double nuy[16];
   int nuy_len = o.Gen_Diff(nuy1_len, nuy1, nuy2_len, nuy2, nuy);
   double nuz1[8];
   int nuz1_len = o.Gen_Product(2, u2x, 2, u3y, nuz1);
   double nuz2[8];
   int nuz2_len = o.Gen_Product(2, u2y, 2, u3x, nuz2);
   double nuz[16];
   int nuz_len = o.Gen_Diff(nuz1_len, nuz1, nuz2_len, nuz2, nuz);
   double nwyuz1_p[16], *nwyuz1 = nwyuz1_p;
   int nwyuz1_len = o.Gen_Product_With_PreAlloc(nwy_len, nwy, nuz_len, nuz, &nwyuz1, 16);
   double nwyuz2_p[16], *nwyuz2 = nwyuz2_p;
   int nwyuz2_len = o.Gen_Product_With_PreAlloc(nwz_len, nwz, nuy_len, nuy, &nwyuz2, 16);
   double nwyuz_p[16], *nwyuz = nwyuz_p;
   int nwyuz_len = o.Gen_Diff_With_PreAlloc(nwyuz1_len, nwyuz1, nwyuz2_len, nwyuz2, &nwyuz, 16);
   double nwxuz1_p[16], *nwxuz1 = nwxuz1_p;
   int nwxuz1_len = o.Gen_Product_With_PreAlloc(nwx_len, nwx, nuz_len, nuz, &nwxuz1, 16);
   double nwxuz2_p[16], *nwxuz2 = nwxuz2_p;
   int nwxuz2_len = o.Gen_Product_With_PreAlloc(nwz_len, nwz, nux_len, nux, &nwxuz2, 16);
   double nwxuz_p[16], *nwxuz = nwxuz_p;
   int nwxuz_len = o.Gen_Diff_With_PreAlloc(nwxuz1_len, nwxuz1, nwxuz2_len, nwxuz2, &nwxuz, 16);
   double nwxuy1_p[16], *nwxuy1 = nwxuy1_p;
   int nwxuy1_len = o.Gen_Product_With_PreAlloc(nwx_len, nwx, nuy_len, nuy, &nwxuy1, 16);
   double nwxuy2_p[16], *nwxuy2 = nwxuy2_p;
   int nwxuy2_len = o.Gen_Product_With_PreAlloc(nwy_len, nwy, nux_len, nux, &nwxuy2, 16);
   double nwxuy_p[16], *nwxuy = nwxuy_p;
   int nwxuy_len = o.Gen_Diff_With_PreAlloc(nwxuy1_len, nwxuy1, nwxuy2_len, nwxuy2, &nwxuy, 16);
   double nvyuz1_p[16], *nvyuz1 = nvyuz1_p;
   int nvyuz1_len = o.Gen_Product_With_PreAlloc(nvy_len, nvy, nuz_len, nuz, &nvyuz1, 16);
   double nvyuz2_p[16], *nvyuz2 = nvyuz2_p;
   int nvyuz2_len = o.Gen_Product_With_PreAlloc(nvz_len, nvz, nuy_len, nuy, &nvyuz2, 16);
   double nvyuz_p[16], *nvyuz = nvyuz_p;
   int nvyuz_len = o.Gen_Diff_With_PreAlloc(nvyuz1_len, nvyuz1, nvyuz2_len, nvyuz2, &nvyuz, 16);
   double nvxuz1_p[16], *nvxuz1 = nvxuz1_p;
   int nvxuz1_len = o.Gen_Product_With_PreAlloc(nvx_len, nvx, nuz_len, nuz, &nvxuz1, 16);
   double nvxuz2_p[16], *nvxuz2 = nvxuz2_p;
   int nvxuz2_len = o.Gen_Product_With_PreAlloc(nvz_len, nvz, nux_len, nux, &nvxuz2, 16);
   double nvxuz_p[16], *nvxuz = nvxuz_p;
   int nvxuz_len = o.Gen_Diff_With_PreAlloc(nvxuz1_len, nvxuz1, nvxuz2_len, nvxuz2, &nvxuz, 16);
   double nvxuy1_p[16], *nvxuy1 = nvxuy1_p;
   int nvxuy1_len = o.Gen_Product_With_PreAlloc(nvx_len, nvx, nuy_len, nuy, &nvxuy1, 16);
   double nvxuy2_p[16], *nvxuy2 = nvxuy2_p;
   int nvxuy2_len = o.Gen_Product_With_PreAlloc(nvy_len, nvy, nux_len, nux, &nvxuy2, 16);
   double nvxuy_p[16], *nvxuy = nvxuy_p;
   int nvxuy_len = o.Gen_Diff_With_PreAlloc(nvxuy1_len, nvxuy1, nvxuy2_len, nvxuy2, &nvxuy, 16);
   double nvywz1_p[16], *nvywz1 = nvywz1_p;
   int nvywz1_len = o.Gen_Product_With_PreAlloc(nvy_len, nvy, nwz_len, nwz, &nvywz1, 16);
   double nvywz2_p[16], *nvywz2 = nvywz2_p;
   int nvywz2_len = o.Gen_Product_With_PreAlloc(nvz_len, nvz, nwy_len, nwy, &nvywz2, 16);
   double nvywz_p[16], *nvywz = nvywz_p;
   int nvywz_len = o.Gen_Diff_With_PreAlloc(nvywz1_len, nvywz1, nvywz2_len, nvywz2, &nvywz, 16);
   double nvxwz1_p[16], *nvxwz1 = nvxwz1_p;
   int nvxwz1_len = o.Gen_Product_With_PreAlloc(nvx_len, nvx, nwz_len, nwz, &nvxwz1, 16);
   double nvxwz2_p[16], *nvxwz2 = nvxwz2_p;
   int nvxwz2_len = o.Gen_Product_With_PreAlloc(nvz_len, nvz, nwx_len, nwx, &nvxwz2, 16);
   double nvxwz_p[16], *nvxwz = nvxwz_p;
   int nvxwz_len = o.Gen_Diff_With_PreAlloc(nvxwz1_len, nvxwz1, nvxwz2_len, nvxwz2, &nvxwz, 16);
   double nvxwy1_p[16], *nvxwy1 = nvxwy1_p;
   int nvxwy1_len = o.Gen_Product_With_PreAlloc(nvx_len, nvx, nwy_len, nwy, &nvxwy1, 16);
   double nvxwy2_p[16], *nvxwy2 = nvxwy2_p;
   int nvxwy2_len = o.Gen_Product_With_PreAlloc(nvy_len, nvy, nwx_len, nwx, &nvxwy2, 16);
   double nvxwy_p[16], *nvxwy = nvxwy_p;
   int nvxwy_len = o.Gen_Diff_With_PreAlloc(nvxwy1_len, nvxwy1, nvxwy2_len, nvxwy2, &nvxwy, 16);
   double p1a_p[16], *p1a = p1a_p;
   int p1a_len = o.Gen_Scale_With_PreAlloc(nvx_len, nvx, ov1x, &p1a, 16);
   double p1b_p[16], *p1b = p1b_p;
   int p1b_len = o.Gen_Scale_With_PreAlloc(nvy_len, nvy, ov1y, &p1b, 16);
   double p1c_p[16], *p1c = p1c_p;
   int p1c_len = o.Gen_Scale_With_PreAlloc(nvz_len, nvz, ov1z, &p1c, 16);
   double p1ab_p[16], *p1ab = p1ab_p;
   int p1ab_len = o.Gen_Sum_With_PreAlloc(p1a_len, p1a, p1b_len, p1b, &p1ab, 16);
   double p1_p[16], *p1 = p1_p;
   int p1_len = o.Gen_Sum_With_PreAlloc(p1ab_len, p1ab, p1c_len, p1c, &p1, 16);
   double p2a_p[16], *p2a = p2a_p;
   int p2a_len = o.Gen_Scale_With_PreAlloc(nwx_len, nwx, ow1x, &p2a, 16);
   double p2b_p[16], *p2b = p2b_p;
   int p2b_len = o.Gen_Scale_With_PreAlloc(nwy_len, nwy, ow1y, &p2b, 16);
   double p2c_p[16], *p2c = p2c_p;
   int p2c_len = o.Gen_Scale_With_PreAlloc(nwz_len, nwz, ow1z, &p2c, 16);
   double p2ab_p[16], *p2ab = p2ab_p;
   int p2ab_len = o.Gen_Sum_With_PreAlloc(p2a_len, p2a, p2b_len, p2b, &p2ab, 16);
   double p2_p[16], *p2 = p2_p;
   int p2_len = o.Gen_Sum_With_PreAlloc(p2ab_len, p2ab, p2c_len, p2c, &p2, 16);
   double p3a_p[16], *p3a = p3a_p;
   int p3a_len = o.Gen_Scale_With_PreAlloc(nux_len, nux, ou1x, &p3a, 16);
   double p3b_p[16], *p3b = p3b_p;
   int p3b_len = o.Gen_Scale_With_PreAlloc(nuy_len, nuy, ou1y, &p3b, 16);
   double p3c_p[16], *p3c = p3c_p;
   int p3c_len = o.Gen_Scale_With_PreAlloc(nuz_len, nuz, ou1z, &p3c, 16);
   double p3ab_p[16], *p3ab = p3ab_p;
   int p3ab_len = o.Gen_Sum_With_PreAlloc(p3a_len, p3a, p3b_len, p3b, &p3ab, 16);
   double p3_p[16], *p3 = p3_p;
   int p3_len = o.Gen_Sum_With_PreAlloc(p3ab_len, p3ab, p3c_len, p3c, &p3, 16);
   double lxa_p[16], *lxa = lxa_p;
   int lxa_len = o.Gen_Product_With_PreAlloc(p1_len, p1, nwyuz_len, nwyuz, &lxa, 16);
   double lxb_p[16], *lxb = lxb_p;
   int lxb_len = o.Gen_Product_With_PreAlloc(p3_len, p3, nvywz_len, nvywz, &lxb, 16);
   double lxc_p[16], *lxc = lxc_p;
   int lxc_len = o.Gen_Product_With_PreAlloc(p2_len, p2, nvyuz_len, nvyuz, &lxc, 16);
   double lxab_p[16], *lxab = lxab_p;
   int lxab_len = o.Gen_Sum_With_PreAlloc(lxa_len, lxa, lxb_len, lxb, &lxab, 16);
   lambda_x_len = o.Gen_Diff_With_PreAlloc(lxab_len, lxab, lxc_len, lxc, &lambda_x, 128);
   double lya_p[16], *lya = lya_p;
   int lya_len = o.Gen_Product_With_PreAlloc(p2_len, p2, nvxuz_len, nvxuz, &lya, 16);
   double lyb_p[16], *lyb = lyb_p;
   int lyb_len = o.Gen_Product_With_PreAlloc(p3_len, p3, nvxwz_len, nvxwz, &lyb, 16);
   double lyc_p[16], *lyc = lyc_p;
   int lyc_len = o.Gen_Product_With_PreAlloc(p1_len, p1, nwxuz_len, nwxuz, &lyc, 16);
   double lybc_p[16], *lybc = lybc_p;
   int lybc_len = o.Gen_Sum_With_PreAlloc(lyc_len, lyc, lyb_len, lyb, &lybc, 16);
   lambda_y_len = o.Gen_Diff_With_PreAlloc(lya_len, lya, lybc_len, lybc, &lambda_y, 128);
   double lza_p[16], *lza = lza_p;
   int lza_len = o.Gen_Product_With_PreAlloc(p3_len, p3, nvxwy_len, nvxwy, &lza, 16);
   double lzb_p[16], *lzb = lzb_p;
   int lzb_len = o.Gen_Product_With_PreAlloc(p1_len, p1, nwxuy_len, nwxuy, &lzb, 16);
   double lzc_p[16], *lzc = lzc_p;
   int lzc_len = o.Gen_Product_With_PreAlloc(p2_len, p2, nvxuy_len, nvxuy, &lzc, 16);
   double lzab_p[16], *lzab = lzab_p;
   int lzab_len = o.Gen_Sum_With_PreAlloc(lza_len, lza, lzb_len, lzb, &lzab, 16);
   lambda_z_len = o.Gen_Diff_With_PreAlloc(lzab_len, lzab, lzc_len, lzc, &lambda_z, 128);
   double da_p[16], *da = da_p;
   int da_len = o.Gen_Product_With_PreAlloc(nvx_len, nvx, nwyuz_len, nwyuz, &da, 16);
   double db_p[16], *db = db_p;
   int db_len = o.Gen_Product_With_PreAlloc(nvz_len, nvz, nwxuy_len, nwxuy, &db, 16);
   double dc_p[16], *dc = dc_p;
   int dc_len = o.Gen_Product_With_PreAlloc(nvy_len, nvy, nwxuz_len, nwxuz, &dc, 16);
   double dab_p[16], *dab = dab_p;
   int dab_len = o.Gen_Sum_With_PreAlloc(da_len, da, db_len, db, &dab, 16);
   lambda_d_len = o.Gen_Diff_With_PreAlloc(dab_len, dab, dc_len, dc, &lambda_d, 128);

   if (dab_p != dab) free(dab);
   if (dc_p != dc) free(dc);
   if (db_p != db) free(db);
   if (da_p != da) free(da);
   if (lzab_p != lzab) free(lzab);
   if (lzc_p != lzc) free(lzc);
   if (lzb_p != lzb) free(lzb);
   if (lza_p != lza) free(lza);
   if (lybc_p != lybc) free(lybc);
   if (lyc_p != lyc) free(lyc);
   if (lyb_p != lyb) free(lyb);
   if (lya_p != lya) free(lya);
   if (lxab_p != lxab) free(lxab);
   if (lxc_p != lxc) free(lxc);
   if (lxb_p != lxb) free(lxb);
   if (lxa_p != lxa) free(lxa);
   if (p3_p != p3) free(p3);
   if (p3ab_p != p3ab) free(p3ab);
   if (p3c_p != p3c) free(p3c);
   if (p3b_p != p3b) free(p3b);
   if (p3a_p != p3a) free(p3a);
   if (p2_p != p2) free(p2);
   if (p2ab_p != p2ab) free(p2ab);
   if (p2c_p != p2c) free(p2c);
   if (p2b_p != p2b) free(p2b);
   if (p2a_p != p2a) free(p2a);
   if (p1_p != p1) free(p1);
   if (p1ab_p != p1ab) free(p1ab);
   if (p1c_p != p1c) free(p1c);
   if (p1b_p != p1b) free(p1b);
   if (p1a_p != p1a) free(p1a);
   if (nvxwy_p != nvxwy) free(nvxwy);
   if (nvxwy2_p != nvxwy2) free(nvxwy2);
   if (nvxwy1_p != nvxwy1) free(nvxwy1);
   if (nvxwz_p != nvxwz) free(nvxwz);
   if (nvxwz2_p != nvxwz2) free(nvxwz2);
   if (nvxwz1_p != nvxwz1) free(nvxwz1);
   if (nvywz_p != nvywz) free(nvywz);
   if (nvywz2_p != nvywz2) free(nvywz2);
   if (nvywz1_p != nvywz1) free(nvywz1);
   if (nvxuy_p != nvxuy) free(nvxuy);
   if (nvxuy2_p != nvxuy2) free(nvxuy2);
   if (nvxuy1_p != nvxuy1) free(nvxuy1);
   if (nvxuz_p != nvxuz) free(nvxuz);
   if (nvxuz2_p != nvxuz2) free(nvxuz2);
   if (nvxuz1_p != nvxuz1) free(nvxuz1);
   if (nvyuz_p != nvyuz) free(nvyuz);
   if (nvyuz2_p != nvyuz2) free(nvyuz2);
   if (nvyuz1_p != nvyuz1) free(nvyuz1);
   if (nwxuy_p != nwxuy) free(nwxuy);
   if (nwxuy2_p != nwxuy2) free(nwxuy2);
   if (nwxuy1_p != nwxuy1) free(nwxuy1);
   if (nwxuz_p != nwxuz) free(nwxuz);
   if (nwxuz2_p != nwxuz2) free(nwxuz2);
   if (nwxuz1_p != nwxuz1) free(nwxuz1);
   if (nwyuz_p != nwyuz) free(nwyuz);
   if (nwyuz2_p != nwyuz2) free(nwyuz2);
   if (nwyuz1_p != nwyuz1) free(nwyuz1);
}

int lessThanOnX_LE_filtered(implicitPoint3D_LPI& p1, double bx)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dbx = bx * d1;
   double kx = l1x - dbx;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(bx)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.932297637868842e-14;
   if (((d1 < 0))) kx = -kx;
   if (kx > epsilon) return IP_Sign::POSITIVE;
   if (-kx > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnX_LE_interval(implicitPoint3D_LPI& p1, interval_number bx)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dbx(bx * d1);
   interval_number kx(l1x - dbx);
   setFPUModeToRoundNEAR();

   if (!kx.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -kx.sign();
   else   return kx.sign();
}

int lessThanOnX_LE_exact(implicitPoint3D_LPI& p1, double bx)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dbx_p[128], *dbx = dbx_p;
   int dbx_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, bx, &dbx, 128);
   double kx_p[128], *kx = kx_p;
   int kx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, dbx_len, dbx, &kx, 128);

   return_value = kx[kx_len - 1];
   if (kx_p != kx) free(kx);
   if (dbx_p != dbx) free(dbx);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnX_LE(implicitPoint3D_LPI& p1, double bx)
{
   int ret;
   ret = lessThanOnX_LE_filtered(p1, bx);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnX_LE_interval(p1, bx);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnX_LE_exact(p1, bx);
}

int lessThanOnX_LL_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1x;
   double k2 = d1 * l2x;
   double kx = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 2.922887626377606e-13;
   if (((d1 < 0) + (d2 < 0)) & 1) kx = -kx;
   if (kx > epsilon) return IP_Sign::POSITIVE;
   if (-kx > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnX_LL_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1x);
   interval_number k2(d1 * l2x);
   interval_number kx(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kx.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -kx.sign();
   else   return kx.sign();
}

int lessThanOnX_LL_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &k2, 128);
   double kx_p[128], *kx = kx_p;
   int kx_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kx, 128);

   return_value = kx[kx_len - 1];
   if (kx_p != kx) free(kx);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnX_LL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   int ret;
   ret = lessThanOnX_LL_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnX_LL_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnX_LL_exact(p1, p2);
}

int lessThanOnX_LT_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1x;
   double k2 = d1 * l2x;
   double kx = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 4.321380059346694e-12;
   if (((d1 < 0) + (d2 < 0)) & 1) kx = -kx;
   if (kx > epsilon) return IP_Sign::POSITIVE;
   if (-kx > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnX_LT_interval(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1x);
   interval_number k2(d1 * l2x);
   interval_number kx(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kx.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -kx.sign();
   else   return kx.sign();
}

int lessThanOnX_LT_exact(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &k2, 128);
   double kx_p[128], *kx = kx_p;
   int kx_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kx, 128);

   return_value = kx[kx_len - 1];
   if (kx_p != kx) free(kx);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnX_LT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   int ret;
   ret = lessThanOnX_LT_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnX_LT_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnX_LT_exact(p1, p2);
}

int lessThanOnX_TE_filtered(implicitPoint3D_TPI& p1, double bx)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dbx = bx * d1;
   double kx = l1x - dbx;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(bx)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.980270973924514e-13;
   if (((d1 < 0))) kx = -kx;
   if (kx > epsilon) return IP_Sign::POSITIVE;
   if (-kx > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnX_TE_interval(implicitPoint3D_TPI& p1, interval_number bx)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dbx(bx * d1);
   interval_number kx(l1x - dbx);
   setFPUModeToRoundNEAR();

   if (!kx.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -kx.sign();
   else   return kx.sign();
}

int lessThanOnX_TE_exact(implicitPoint3D_TPI& p1, double bx)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dbx_p[128], *dbx = dbx_p;
   int dbx_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, bx, &dbx, 128);
   double kx_p[128], *kx = kx_p;
   int kx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, dbx_len, dbx, &kx, 128);

   return_value = kx[kx_len - 1];
   if (kx_p != kx) free(kx);
   if (dbx_p != dbx) free(dbx);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnX_TE(implicitPoint3D_TPI& p1, double bx)
{
   int ret;
   ret = lessThanOnX_TE_filtered(p1, bx);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnX_TE_interval(p1, bx);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnX_TE_exact(p1, bx);
}

int lessThanOnX_TT_filtered(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1x;
   double k2 = d1 * l2x;
   double kx = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 5.504141586953918e-11;
   if (((d1 < 0) + (d2 < 0)) & 1) kx = -kx;
   if (kx > epsilon) return IP_Sign::POSITIVE;
   if (-kx > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnX_TT_interval(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1x);
   interval_number k2(d1 * l2x);
   interval_number kx(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kx.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -kx.sign();
   else   return kx.sign();
}

int lessThanOnX_TT_exact(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &k2, 128);
   double kx_p[128], *kx = kx_p;
   int kx_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kx, 128);

   return_value = kx[kx_len - 1];
   if (kx_p != kx) free(kx);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnX_TT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   int ret;
   ret = lessThanOnX_TT_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnX_TT_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnX_TT_exact(p1, p2);
}

int lessThanOnY_LE_filtered(implicitPoint3D_LPI& p1, double by)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dby = by * d1;
   double ky = l1y - dby;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(by)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.932297637868842e-14;
   if (((d1 < 0))) ky = -ky;
   if (ky > epsilon) return IP_Sign::POSITIVE;
   if (-ky > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnY_LE_interval(implicitPoint3D_LPI& p1, interval_number by)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dby(by * d1);
   interval_number ky(l1y - dby);
   setFPUModeToRoundNEAR();

   if (!ky.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -ky.sign();
   else   return ky.sign();
}

int lessThanOnY_LE_exact(implicitPoint3D_LPI& p1, double by)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dby_p[128], *dby = dby_p;
   int dby_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, by, &dby, 128);
   double ky_p[128], *ky = ky_p;
   int ky_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, dby_len, dby, &ky, 128);

   return_value = ky[ky_len - 1];
   if (ky_p != ky) free(ky);
   if (dby_p != dby) free(dby);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnY_LE(implicitPoint3D_LPI& p1, double by)
{
   int ret;
   ret = lessThanOnY_LE_filtered(p1, by);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnY_LE_interval(p1, by);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnY_LE_exact(p1, by);
}

int lessThanOnY_LL_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1y;
   double k2 = d1 * l2y;
   double ky = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 2.922887626377606e-13;
   if (((d1 < 0) + (d2 < 0)) & 1) ky = -ky;
   if (ky > epsilon) return IP_Sign::POSITIVE;
   if (-ky > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnY_LL_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1y);
   interval_number k2(d1 * l2y);
   interval_number ky(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!ky.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -ky.sign();
   else   return ky.sign();
}

int lessThanOnY_LL_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &k2, 128);
   double ky_p[128], *ky = ky_p;
   int ky_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &ky, 128);

   return_value = ky[ky_len - 1];
   if (ky_p != ky) free(ky);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnY_LL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   int ret;
   ret = lessThanOnY_LL_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnY_LL_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnY_LL_exact(p1, p2);
}

int lessThanOnY_LT_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1y;
   double k2 = d1 * l2y;
   double ky = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 4.321380059346694e-12;
   if (((d1 < 0) + (d2 < 0)) & 1) ky = -ky;
   if (ky > epsilon) return IP_Sign::POSITIVE;
   if (-ky > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnY_LT_interval(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1y);
   interval_number k2(d1 * l2y);
   interval_number ky(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!ky.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -ky.sign();
   else   return ky.sign();
}

int lessThanOnY_LT_exact(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &k2, 128);
   double ky_p[128], *ky = ky_p;
   int ky_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &ky, 128);

   return_value = ky[ky_len - 1];
   if (ky_p != ky) free(ky);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnY_LT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   int ret;
   ret = lessThanOnY_LT_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnY_LT_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnY_LT_exact(p1, p2);
}

int lessThanOnY_TE_filtered(implicitPoint3D_TPI& p1, double by)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dby = by * d1;
   double ky = l1y - dby;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(by)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.980270973924514e-13;
   if (((d1 < 0))) ky = -ky;
   if (ky > epsilon) return IP_Sign::POSITIVE;
   if (-ky > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnY_TE_interval(implicitPoint3D_TPI& p1, interval_number by)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dby(by * d1);
   interval_number ky(l1y - dby);
   setFPUModeToRoundNEAR();

   if (!ky.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -ky.sign();
   else   return ky.sign();
}

int lessThanOnY_TE_exact(implicitPoint3D_TPI& p1, double by)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dby_p[128], *dby = dby_p;
   int dby_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, by, &dby, 128);
   double ky_p[128], *ky = ky_p;
   int ky_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, dby_len, dby, &ky, 128);

   return_value = ky[ky_len - 1];
   if (ky_p != ky) free(ky);
   if (dby_p != dby) free(dby);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnY_TE(implicitPoint3D_TPI& p1, double by)
{
   int ret;
   ret = lessThanOnY_TE_filtered(p1, by);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnY_TE_interval(p1, by);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnY_TE_exact(p1, by);
}

int lessThanOnY_TT_filtered(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1y;
   double k2 = d1 * l2y;
   double ky = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 5.504141586953918e-11;
   if (((d1 < 0) + (d2 < 0)) & 1) ky = -ky;
   if (ky > epsilon) return IP_Sign::POSITIVE;
   if (-ky > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnY_TT_interval(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1y);
   interval_number k2(d1 * l2y);
   interval_number ky(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!ky.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -ky.sign();
   else   return ky.sign();
}

int lessThanOnY_TT_exact(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &k2, 128);
   double ky_p[128], *ky = ky_p;
   int ky_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &ky, 128);

   return_value = ky[ky_len - 1];
   if (ky_p != ky) free(ky);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnY_TT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   int ret;
   ret = lessThanOnY_TT_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnY_TT_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnY_TT_exact(p1, p2);
}

int lessThanOnZ_LE_filtered(implicitPoint3D_LPI& p1, double bz)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dbz = bz * d1;
   double kz = l1z - dbz;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(bz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.932297637868842e-14;
   if (((d1 < 0))) kz = -kz;
   if (kz > epsilon) return IP_Sign::POSITIVE;
   if (-kz > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnZ_LE_interval(implicitPoint3D_LPI& p1, interval_number bz)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dbz(bz * d1);
   interval_number kz(l1z - dbz);
   setFPUModeToRoundNEAR();

   if (!kz.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -kz.sign();
   else   return kz.sign();
}

int lessThanOnZ_LE_exact(implicitPoint3D_LPI& p1, double bz)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dbz_p[128], *dbz = dbz_p;
   int dbz_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, bz, &dbz, 128);
   double kz_p[128], *kz = kz_p;
   int kz_len = o.Gen_Diff_With_PreAlloc(l1z_len, l1z, dbz_len, dbz, &kz, 128);

   return_value = kz[kz_len - 1];
   if (kz_p != kz) free(kz);
   if (dbz_p != dbz) free(dbz);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnZ_LE(implicitPoint3D_LPI& p1, double bz)
{
   int ret;
   ret = lessThanOnZ_LE_filtered(p1, bz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnZ_LE_interval(p1, bz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnZ_LE_exact(p1, bz);
}

int lessThanOnZ_LL_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1z;
   double k2 = d1 * l2z;
   double kz = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 2.922887626377606e-13;
   if (((d1 < 0) + (d2 < 0)) & 1) kz = -kz;
   if (kz > epsilon) return IP_Sign::POSITIVE;
   if (-kz > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnZ_LL_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1z);
   interval_number k2(d1 * l2z);
   interval_number kz(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kz.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -kz.sign();
   else   return kz.sign();
}

int lessThanOnZ_LL_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &k2, 128);
   double kz_p[128], *kz = kz_p;
   int kz_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kz, 128);

   return_value = kz[kz_len - 1];
   if (kz_p != kz) free(kz);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnZ_LL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
   int ret;
   ret = lessThanOnZ_LL_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnZ_LL_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnZ_LL_exact(p1, p2);
}

int lessThanOnZ_LT_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1z;
   double k2 = d1 * l2z;
   double kz = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 4.321380059346694e-12;
   if (((d1 < 0) + (d2 < 0)) & 1) kz = -kz;
   if (kz > epsilon) return IP_Sign::POSITIVE;
   if (-kz > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnZ_LT_interval(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1z);
   interval_number k2(d1 * l2z);
   interval_number kz(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kz.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -kz.sign();
   else   return kz.sign();
}

int lessThanOnZ_LT_exact(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &k2, 128);
   double kz_p[128], *kz = kz_p;
   int kz_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kz, 128);

   return_value = kz[kz_len - 1];
   if (kz_p != kz) free(kz);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnZ_LT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
   int ret;
   ret = lessThanOnZ_LT_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnZ_LT_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnZ_LT_exact(p1, p2);
}

int lessThanOnZ_TE_filtered(implicitPoint3D_TPI& p1, double bz)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dbz = bz * d1;
   double kz = l1z - dbz;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(bz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.980270973924514e-13;
   if (((d1 < 0))) kz = -kz;
   if (kz > epsilon) return IP_Sign::POSITIVE;
   if (-kz > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnZ_TE_interval(implicitPoint3D_TPI& p1, interval_number bz)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dbz(bz * d1);
   interval_number kz(l1z - dbz);
   setFPUModeToRoundNEAR();

   if (!kz.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -kz.sign();
   else   return kz.sign();
}

int lessThanOnZ_TE_exact(implicitPoint3D_TPI& p1, double bz)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dbz_p[128], *dbz = dbz_p;
   int dbz_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, bz, &dbz, 128);
   double kz_p[128], *kz = kz_p;
   int kz_len = o.Gen_Diff_With_PreAlloc(l1z_len, l1z, dbz_len, dbz, &kz, 128);

   return_value = kz[kz_len - 1];
   if (kz_p != kz) free(kz);
   if (dbz_p != dbz) free(dbz);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnZ_TE(implicitPoint3D_TPI& p1, double bz)
{
   int ret;
   ret = lessThanOnZ_TE_filtered(p1, bz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnZ_TE_interval(p1, bz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnZ_TE_exact(p1, bz);
}

int lessThanOnZ_TT_filtered(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double k1 = d2 * l1z;
   double k2 = d1 * l2z;
   double kz = k1 - k2;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 5.504141586953918e-11;
   if (((d1 < 0) + (d2 < 0)) & 1) kz = -kz;
   if (kz > epsilon) return IP_Sign::POSITIVE;
   if (-kz > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int lessThanOnZ_TT_interval(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number k1(d2 * l1z);
   interval_number k2(d1 * l2z);
   interval_number kz(k1 - k2);
   setFPUModeToRoundNEAR();

   if (!kz.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0) + (d2 < 0)) & 1) return -kz.sign();
   else   return kz.sign();
}

int lessThanOnZ_TT_exact(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double k1_p[128], *k1 = k1_p;
   int k1_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &k1, 128);
   double k2_p[128], *k2 = k2_p;
   int k2_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &k2, 128);
   double kz_p[128], *kz = kz_p;
   int kz_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kz, 128);

   return_value = kz[kz_len - 1];
   if (kz_p != kz) free(kz);
   if (k2_p != k2) free(k2);
   if (k1_p != k1) free(k1);
   if (( (d1[d1_len -1] < 0) + (d2[d2_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int lessThanOnZ_TT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
   int ret;
   ret = lessThanOnZ_TT_filtered(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = lessThanOnZ_TT_interval(p1, p2);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return lessThanOnZ_TT_exact(p1, p2);
}

int orient2d3d_indirect_LEE_filtered(implicitPoint3D_LPI& p1, double p2x, double p2y, double p3x, double p3y)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * p2x;
   double c = d1 * p3y;
   double e = d1 * p2y;
   double g = d1 * p3x;
   double ab = a - l1x;
   double cd = c - l1y;
   double ef = e - l1y;
   double gh = g - l1x;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(p2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3y)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 1.504941136987184e-12;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_LEE_interval(implicitPoint3D_LPI& p1, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * p2x);
   interval_number c(d1 * p3y);
   interval_number e(d1 * p2y);
   interval_number g(d1 * p3x);
   interval_number ab(a - l1x);
   interval_number cd(c - l1y);
   interval_number ef(e - l1y);
   interval_number gh(g - l1x);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int orient2d3d_indirect_LEE_exact(implicitPoint3D_LPI& p1, double p2x, double p2y, double p3x, double p3y)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double a_p[128], *a = a_p;
   int a_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p2x, &a, 128);
   double c_p[128], *c = c_p;
   int c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3y, &c, 128);
   double e_p[128], *e = e_p;
   int e_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p2y, &e, 128);
   double g_p[128], *g = g_p;
   int g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3x, &g, 128);
   double ab_p[128], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, l1x_len, l1x, &ab, 128);
   double cd_p[128], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1y_len, l1y, &cd, 128);
   double ef_p[128], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, l1y_len, l1y, &ef, 128);
   double gh_p[128], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1x_len, l1x, &gh, 128);
   double abcd_p[128], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 128);
   double efgh_p[128], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 128);
   double L_p[128], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 128);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (g_p != g) free(g);
   if (e_p != e) free(e);
   if (c_p != c) free(c);
   if (a_p != a) free(a);
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_LEE(implicitPoint3D_LPI& p1, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2d3d_indirect_LEE_filtered(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_LEE_interval(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_LEE_exact(p1, p2x, p2y, p3x, p3y);
}

int orient2d3d_indirect_LLE_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double op3x, double op3y)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * op3y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * op3x;
   double ab = a - b;
   double cd = c - l1y;
   double ef = e - f;
   double gh = g - l1x;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(op3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(op3y)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 1.699690735379461e-11;
   if (((d2 < 0))) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_LLE_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, interval_number op3x, interval_number op3y)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * op3y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * op3x);
   interval_number ab(a - b);
   interval_number cd(c - l1y);
   interval_number ef(e - f);
   interval_number gh(g - l1x);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0))) return -L.sign();
   else   return L.sign();
}

int orient2d3d_indirect_LLE_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double op3x, double op3y)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3y, &c, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3x, &g, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1y_len, l1y, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1x_len, l1x, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_LLE(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double op3x, double op3y)
{
   int ret;
   ret = orient2d3d_indirect_LLE_filtered(p1, p2, op3x, op3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_LLE_interval(p1, p2, op3x, op3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_LLE_exact(p1, p2, op3x, op3y);
}

int orient2d3d_indirect_LLL_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * l3y;
   double d = d3 * l1y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * l3x;
   double h = d3 * l1x;
   double ab = a - b;
   double cd = c - d;
   double ef = e - f;
   double gh = g - h;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 1.75634284893534e-10;
   if (((d2 < 0) + (d3 < 0)) & 1) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_LLL_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * l3y);
   interval_number d(d3 * l1y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * l3x);
   interval_number h(d3 * l1x);
   interval_number ab(a - b);
   interval_number cd(c - d);
   interval_number ef(e - f);
   interval_number gh(g - h);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0) + (d3 < 0)) & 1) return -L.sign();
   else   return L.sign();
}

int orient2d3d_indirect_LLL_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64], *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64], *l3y = l3y_p, l3z_p[64], *l3z = l3z_p, d3_p[64], *d3 = d3_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3y_len, l3y, &c, 64);
   double d_p[64], *d = d_p;
   int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1y_len, l1y, &d, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3x_len, l3x, &g, 64);
   double h_p[64], *h = h_p;
   int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1x_len, l1x, &h, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (h_p != h) free(h);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (d_p != d) free(d);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0) + (d3[d3_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);
 if (l3x_p != l3x) free(l3x);
 if (l3y_p != l3y) free(l3y);
 if (l3z_p != l3z) free(l3z);
 if (d3_p != d3) free(d3);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_LLL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3)
{
   int ret;
   ret = orient2d3d_indirect_LLL_filtered(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_LLL_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_LLL_exact(p1, p2, p3);
}

int orient2d3d_indirect_LLT_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_TPI& p3)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * l3y;
   double d = d3 * l1y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * l3x;
   double h = d3 * l1x;
   double ab = a - b;
   double cd = c - d;
   double ef = e - f;
   double gh = g - h;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= 2.144556754402072e-09;
   if (((d2 < 0) + (d3 < 0)) & 1) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_LLT_interval(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_TPI& p3)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * l3y);
   interval_number d(d3 * l1y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * l3x);
   interval_number h(d3 * l1x);
   interval_number ab(a - b);
   interval_number cd(c - d);
   interval_number ef(e - f);
   interval_number gh(g - h);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0) + (d3 < 0)) & 1) return -L.sign();
   else   return L.sign();
}

int orient2d3d_indirect_LLT_exact(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_TPI& p3)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64], *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64], *l3y = l3y_p, l3z_p[64], *l3z = l3z_p, d3_p[64], *d3 = d3_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3y_len, l3y, &c, 64);
   double d_p[64], *d = d_p;
   int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1y_len, l1y, &d, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3x_len, l3x, &g, 64);
   double h_p[64], *h = h_p;
   int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1x_len, l1x, &h, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (h_p != h) free(h);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (d_p != d) free(d);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0) + (d3[d3_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);
 if (l3x_p != l3x) free(l3x);
 if (l3y_p != l3y) free(l3y);
 if (l3z_p != l3z) free(l3z);
 if (d3_p != d3) free(d3);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_LLT(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_TPI& p3)
{
   int ret;
   ret = orient2d3d_indirect_LLT_filtered(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_LLT_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_LLT_exact(p1, p2, p3);
}

int orient2d3d_indirect_LTE_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * p3y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * p3x;
   double ab = a - b;
   double cd = c - l1y;
   double ef = e - f;
   double gh = g - l1x;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(p3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3y)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 2.184958117212875e-10;
   if (((d2 < 0))) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_LTE_interval(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * p3y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * p3x);
   interval_number ab(a - b);
   interval_number cd(c - l1y);
   interval_number ef(e - f);
   interval_number gh(g - l1x);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0))) return -L.sign();
   else   return L.sign();
}

int orient2d3d_indirect_LTE_exact(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3y, &c, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3x, &g, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1y_len, l1y, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1x_len, l1x, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_LTE(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y)
{
   int ret;
   ret = orient2d3d_indirect_LTE_filtered(p1, p2, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_LTE_interval(p1, p2, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_LTE_exact(p1, p2, p3x, p3y);
}

int orient2d3d_indirect_LTT_filtered(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * l3y;
   double d = d3 * l1y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * l3x;
   double h = d3 * l1x;
   double ab = a - b;
   double cd = c - d;
   double ef = e - f;
   double gh = g - h;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 2.535681042914479e-08;
   if (((d2 < 0) + (d3 < 0)) & 1) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_LTT_interval(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * l3y);
   interval_number d(d3 * l1y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * l3x);
   interval_number h(d3 * l1x);
   interval_number ab(a - b);
   interval_number cd(c - d);
   interval_number ef(e - f);
   interval_number gh(g - h);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0) + (d3 < 0)) & 1) return -L.sign();
   else   return L.sign();
}

int orient2d3d_indirect_LTT_exact(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64], *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64], *l3y = l3y_p, l3z_p[64], *l3z = l3z_p, d3_p[64], *d3 = d3_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3y_len, l3y, &c, 64);
   double d_p[64], *d = d_p;
   int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1y_len, l1y, &d, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3x_len, l3x, &g, 64);
   double h_p[64], *h = h_p;
   int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1x_len, l1x, &h, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (h_p != h) free(h);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (d_p != d) free(d);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0) + (d3[d3_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);
 if (l3x_p != l3x) free(l3x);
 if (l3y_p != l3y) free(l3y);
 if (l3z_p != l3z) free(l3z);
 if (d3_p != d3) free(d3);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_LTT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
   int ret;
   ret = orient2d3d_indirect_LTT_filtered(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_LTT_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_LTT_exact(p1, p2, p3);
}

int orient2d3d_indirect_TEE_filtered(implicitPoint3D_TPI& p1, double p2x, double p2y, double p3x, double p3y)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * p2x;
   double c = d1 * p3y;
   double e = d1 * p2y;
   double g = d1 * p3x;
   double ab = a - l1x;
   double cd = c - l1y;
   double ef = e - l1y;
   double gh = g - l1x;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(p2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3y)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.202367260257704e-10;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_TEE_interval(implicitPoint3D_TPI& p1, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * p2x);
   interval_number c(d1 * p3y);
   interval_number e(d1 * p2y);
   interval_number g(d1 * p3x);
   interval_number ab(a - l1x);
   interval_number cd(c - l1y);
   interval_number ef(e - l1y);
   interval_number gh(g - l1x);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int orient2d3d_indirect_TEE_exact(implicitPoint3D_TPI& p1, double p2x, double p2y, double p3x, double p3y)
{
 double return_value = 0.0;
 double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128], *l1z = l1z_p, d1_p[128], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double a_p[128], *a = a_p;
   int a_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p2x, &a, 128);
   double c_p[128], *c = c_p;
   int c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3y, &c, 128);
   double e_p[128], *e = e_p;
   int e_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p2y, &e, 128);
   double g_p[128], *g = g_p;
   int g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3x, &g, 128);
   double ab_p[128], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, l1x_len, l1x, &ab, 128);
   double cd_p[128], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1y_len, l1y, &cd, 128);
   double ef_p[128], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, l1y_len, l1y, &ef, 128);
   double gh_p[128], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1x_len, l1x, &gh, 128);
   double abcd_p[128], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 128);
   double efgh_p[128], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 128);
   double L_p[128], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 128);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (g_p != g) free(g);
   if (e_p != e) free(e);
   if (c_p != c) free(c);
   if (a_p != a) free(a);
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_TEE(implicitPoint3D_TPI& p1, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2d3d_indirect_TEE_filtered(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_TEE_interval(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_TEE_exact(p1, p2x, p2y, p3x, p3y);
}

int orient2d3d_indirect_TTE_filtered(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * p3y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * p3x;
   double ab = a - b;
   double cd = c - l1y;
   double ef = e - f;
   double gh = g - l1x;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(p3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3y)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.307187945722513e-08;
   if (((d2 < 0))) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_TTE_interval(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * p3y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * p3x);
   interval_number ab(a - b);
   interval_number cd(c - l1y);
   interval_number ef(e - f);
   interval_number gh(g - l1x);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0))) return -L.sign();
   else   return L.sign();
}

int orient2d3d_indirect_TTE_exact(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64], *d2 = d2_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3y, &c, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3x, &g, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1y_len, l1y, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1x_len, l1x, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_TTE(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y)
{
   int ret;
   ret = orient2d3d_indirect_TTE_filtered(p1, p2, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_TTE_interval(p1, p2, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_TTE_exact(p1, p2, p3x, p3y);
}

int orient2d3d_indirect_TTT_filtered(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
   double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * l3y;
   double d = d3 * l1y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * l3x;
   double h = d3 * l1x;
   double ab = a - b;
   double cd = c - d;
   double ef = e - f;
   double gh = g - h;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 3.103174776697444e-06;
   if (((d2 < 0) + (d3 < 0)) & 1) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d3d_indirect_TTT_interval(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
   interval_number l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   || !p2.getIntervalLambda(l2x, l2y, l2z, d2)
   || !p3.getIntervalLambda(l3x, l3y, l3z, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * l3y);
   interval_number d(d3 * l1y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * l3x);
   interval_number h(d3 * l1x);
   interval_number ab(a - b);
   interval_number cd(c - d);
   interval_number ef(e - f);
   interval_number gh(g - h);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0) + (d3 < 0)) & 1) return -L.sign();
   else   return L.sign();
}

int orient2d3d_indirect_TTT_exact(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64], *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64], *l3y = l3y_p, l3z_p[64], *l3z = l3z_p, d3_p[64], *d3 = d3_p;
 int l1x_len, l1y_len, l1z_len, d1_len, l2x_len, l2y_len, l2z_len, d2_len, l3x_len, l3y_len, l3z_len, d3_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, l2z, l2z_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, l3z, l3z_len, d3, d3_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3y_len, l3y, &c, 64);
   double d_p[64], *d = d_p;
   int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1y_len, l1y, &d, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3x_len, l3x, &g, 64);
   double h_p[64], *h = h_p;
   int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1x_len, l1x, &h, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (h_p != h) free(h);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (d_p != d) free(d);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0) + (d3[d3_len -1] < 0)) & 1) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);
 if (l2x_p != l2x) free(l2x);
 if (l2y_p != l2y) free(l2y);
 if (l2z_p != l2z) free(l2z);
 if (d2_p != d2) free(d2);
 if (l3x_p != l3x) free(l3x);
 if (l3y_p != l3y) free(l3y);
 if (l3z_p != l3z) free(l3z);
 if (d3_p != d3) free(d3);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d3d_indirect_TTT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3)
{
   int ret;
//   ret = orient2d3d_indirect_TTT_filtered(p1, p2, p3);
//   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d3d_indirect_TTT_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d3d_indirect_TTT_exact(p1, p2, p3);
}

int orient2d_indirect_SEE_filtered(implicitPoint2D_SSI& p1, double p2x, double p2y, double p3x, double p3y)
{
   double l1x, l1y, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * p2x;
   double c = d1 * p3y;
   double e = d1 * p2y;
   double g = d1 * p3x;
   double ab = a - l1x;
   double cd = c - l1y;
   double ef = e - l1y;
   double gh = g - l1x;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(p2x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p2y)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3y)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 1.048310743767546e-13;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d_indirect_SEE_interval(implicitPoint2D_SSI& p1, interval_number p2x, interval_number p2y, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * p2x);
   interval_number c(d1 * p3y);
   interval_number e(d1 * p2y);
   interval_number g(d1 * p3x);
   interval_number ab(a - l1x);
   interval_number cd(c - l1y);
   interval_number ef(e - l1y);
   interval_number gh(g - l1x);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   return L.sign();
}

int orient2d_indirect_SEE_exact(implicitPoint2D_SSI& p1, double p2x, double p2y, double p3x, double p3y)
{
 double return_value = 0.0;
 double l1x[32], l1y[32], d1[16];
 int l1x_len, l1y_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double a[32];
   int a_len = o.Gen_Scale(d1_len, d1, p2x, a);
   double c[32];
   int c_len = o.Gen_Scale(d1_len, d1, p3y, c);
   double e[32];
   int e_len = o.Gen_Scale(d1_len, d1, p2y, e);
   double g[32];
   int g_len = o.Gen_Scale(d1_len, d1, p3x, g);
   double ab[64];
   int ab_len = o.Gen_Diff(a_len, a, l1x_len, l1x, ab);
   double cd[64];
   int cd_len = o.Gen_Diff(c_len, c, l1y_len, l1y, cd);
   double ef[64];
   int ef_len = o.Gen_Diff(e_len, e, l1y_len, l1y, ef);
   double gh[64];
   int gh_len = o.Gen_Diff(g_len, g, l1x_len, l1x, gh);
   double abcd_p[128], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 128);
   double efgh_p[128], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 128);
   double L_p[128], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 128);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
 }


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d_indirect_SEE(implicitPoint2D_SSI& p1, double p2x, double p2y, double p3x, double p3y)
{
   int ret;
   ret = orient2d_indirect_SEE_filtered(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d_indirect_SEE_interval(p1, p2x, p2y, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d_indirect_SEE_exact(p1, p2x, p2y, p3x, p3y);
}

int orient2d_indirect_SSE_filtered(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double p3x, double p3y)
{
   double l1x, l1y, d1, l2x, l2y, d2, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, d2, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * p3y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * p3x;
   double ab = a - b;
   double cd = c - l1y;
   double ef = e - f;
   double gh = g - l1x;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(p3x)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(p3y)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= 3.837902218251096e-13;
   if (((d2 < 0))) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d_indirect_SSE_interval(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, interval_number p3x, interval_number p3y)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * p3y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * p3x);
   interval_number ab(a - b);
   interval_number cd(c - l1y);
   interval_number ef(e - f);
   interval_number gh(g - l1x);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0))) return -L.sign();
   else   return L.sign();
}

int orient2d_indirect_SSE_exact(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double p3x, double p3y)
{
 double return_value = 0.0;
 double l1x[32], l1y[32], d1[16], l2x[32], l2y[32], d2[16];
 int l1x_len, l1y_len, d1_len, l2x_len, l2y_len, d2_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, d2, d2_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
 {
   expansionObject o;
   double a_p[128], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 128);
   double b_p[128], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 128);
   double c[32];
   int c_len = o.Gen_Scale(d1_len, d1, p3y, c);
   double e_p[128], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 128);
   double f_p[128], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 128);
   double g[32];
   int g_len = o.Gen_Scale(d1_len, d1, p3x, g);
   double ab_p[128], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 128);
   double cd[64];
   int cd_len = o.Gen_Diff(c_len, c, l1y_len, l1y, cd);
   double ef_p[128], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 128);
   double gh[64];
   int gh_len = o.Gen_Diff(g_len, g, l1x_len, l1x, gh);
   double abcd_p[128], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 128);
   double efgh_p[128], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 128);
   double L_p[128], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 128);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (ef_p != ef) free(ef);
   if (ab_p != ab) free(ab);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0))) return_value = -return_value;
 }


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d_indirect_SSE(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double p3x, double p3y)
{
   int ret;
   ret = orient2d_indirect_SSE_filtered(p1, p2, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d_indirect_SSE_interval(p1, p2, p3x, p3y);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d_indirect_SSE_exact(p1, p2, p3x, p3y);
}

int orient2d_indirect_SSS_filtered(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3)
{
   double l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, d1, max_var)
       || !p2.getFilteredLambda(l2x, l2y, d2, max_var)
       || !p3.getFilteredLambda(l3x, l3y, d3, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double a = d1 * l2x;
   double b = d2 * l1x;
   double c = d1 * l3y;
   double d = d3 * l1y;
   double e = d1 * l2y;
   double f = d2 * l1y;
   double g = d1 * l3x;
   double h = d3 * l1x;
   double ab = a - b;
   double cd = c - d;
   double ef = e - f;
   double gh = g - h;
   double abcd = ab * cd;
   double efgh = ef * gh;
   double L = abcd - efgh;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 1.364575119566784e-12;
   if (((d2 < 0) + (d3 < 0)) & 1) L = -L;
   if (L > epsilon) return IP_Sign::POSITIVE;
   if (-L > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient2d_indirect_SSS_interval(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3)
{
   interval_number l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
   if (
   !p1.getIntervalLambda(l1x, l1y, d1)
   || !p2.getIntervalLambda(l2x, l2y, d2)
   || !p3.getIntervalLambda(l3x, l3y, d3)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number a(d1 * l2x);
   interval_number b(d2 * l1x);
   interval_number c(d1 * l3y);
   interval_number d(d3 * l1y);
   interval_number e(d1 * l2y);
   interval_number f(d2 * l1y);
   interval_number g(d1 * l3x);
   interval_number h(d3 * l1x);
   interval_number ab(a - b);
   interval_number cd(c - d);
   interval_number ef(e - f);
   interval_number gh(g - h);
   interval_number abcd(ab * cd);
   interval_number efgh(ef * gh);
   interval_number L(abcd - efgh);
   setFPUModeToRoundNEAR();

   if (!L.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d2 < 0) + (d3 < 0)) & 1) return -L.sign();
   else   return L.sign();
}

int orient2d_indirect_SSS_exact(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3)
{
 double return_value = 0.0;
 double l1x[32], l1y[32], d1[16], l2x[32], l2y[32], d2[16], l3x[32], l3y[32], d3[16];
 int l1x_len, l1y_len, d1_len, l2x_len, l2y_len, d2_len, l3x_len, l3y_len, d3_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, d1, d1_len);
 p2.getExactLambda(l2x, l2x_len, l2y, l2y_len, d2, d2_len);
 p3.getExactLambda(l3x, l3x_len, l3y, l3y_len, d3, d3_len);
 if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
 {
   expansionObject o;
   double a_p[64], *a = a_p;
   int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
   double b_p[64], *b = b_p;
   int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
   double c_p[64], *c = c_p;
   int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3y_len, l3y, &c, 64);
   double d_p[64], *d = d_p;
   int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1y_len, l1y, &d, 64);
   double e_p[64], *e = e_p;
   int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
   double f_p[64], *f = f_p;
   int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
   double g_p[64], *g = g_p;
   int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3x_len, l3x, &g, 64);
   double h_p[64], *h = h_p;
   int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1x_len, l1x, &h, 64);
   double ab_p[64], *ab = ab_p;
   int ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
   double cd_p[64], *cd = cd_p;
   int cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
   double ef_p[64], *ef = ef_p;
   int ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
   double gh_p[64], *gh = gh_p;
   int gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);
   double abcd_p[64], *abcd = abcd_p;
   int abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
   double efgh_p[64], *efgh = efgh_p;
   int efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
   double L_p[64], *L = L_p;
   int L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

   return_value = L[L_len - 1];
   if (L_p != L) free(L);
   if (efgh_p != efgh) free(efgh);
   if (abcd_p != abcd) free(abcd);
   if (gh_p != gh) free(gh);
   if (ef_p != ef) free(ef);
   if (cd_p != cd) free(cd);
   if (ab_p != ab) free(ab);
   if (h_p != h) free(h);
   if (g_p != g) free(g);
   if (f_p != f) free(f);
   if (e_p != e) free(e);
   if (d_p != d) free(d);
   if (c_p != c) free(c);
   if (b_p != b) free(b);
   if (a_p != a) free(a);
   if (( (d2[d2_len -1] < 0) + (d3[d3_len -1] < 0)) & 1) return_value = -return_value;
 }


 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient2d_indirect_SSS(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3)
{
   int ret;
   ret = orient2d_indirect_SSS_filtered(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient2d_indirect_SSS_interval(p1, p2, p3);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient2d_indirect_SSS_exact(p1, p2, p3);
}

int orient3d_indirect_LEEE_filtered(implicitPoint3D_LPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dcx = d1 * cx;
   double dcy = d1 * cy;
   double dcz = d1 * cz;
   double ix_cx = l1x - dcx;
   double iy_cy = l1y - dcy;
   double ax_cx = ax - cx;
   double ay_cy = ay - cy;
   double az_cz = az - cz;
   double iz_cz = l1z - dcz;
   double bx_cx = bx - cx;
   double by_cy = by - cy;
   double bz_cz = bz - cz;
   double tmc_a = ix_cx * ay_cy;
   double tmc_b = iy_cy * ax_cx;
   double m01 = tmc_a - tmc_b;
   double tmi_a = ix_cx * az_cz;
   double tmi_b = iz_cz * ax_cx;
   double m02 = tmi_a - tmi_b;
   double tma_a = iy_cy * az_cz;
   double tma_b = iz_cz * ay_cy;
   double m12 = tma_a - tma_b;
   double mt1 = m01 * bz_cz;
   double mt2 = m02 * by_cy;
   double mt3 = m12 * bx_cx;
   double mtt = mt2 - mt1;
   double m012 = mtt - mt3;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(cx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ax_cx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ay_cy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(az_cz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bx_cx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(by_cy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bz_cz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= max_var;
   epsilon *= 1.861039534284405e-13;
   if (((d1 < 0))) m012 = -m012;
   if (m012 > epsilon) return IP_Sign::POSITIVE;
   if (-m012 > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient3d_indirect_LEEE_interval(implicitPoint3D_LPI& p1, interval_number ax, interval_number ay, interval_number az, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dcx(d1 * cx);
   interval_number dcy(d1 * cy);
   interval_number dcz(d1 * cz);
   interval_number ix_cx(l1x - dcx);
   interval_number iy_cy(l1y - dcy);
   interval_number ax_cx(ax - cx);
   interval_number ay_cy(ay - cy);
   interval_number az_cz(az - cz);
   interval_number iz_cz(l1z - dcz);
   interval_number bx_cx(bx - cx);
   interval_number by_cy(by - cy);
   interval_number bz_cz(bz - cz);
   interval_number tmc_a(ix_cx * ay_cy);
   interval_number tmc_b(iy_cy * ax_cx);
   interval_number m01(tmc_a - tmc_b);
   interval_number tmi_a(ix_cx * az_cz);
   interval_number tmi_b(iz_cz * ax_cx);
   interval_number m02(tmi_a - tmi_b);
   interval_number tma_a(iy_cy * az_cz);
   interval_number tma_b(iz_cz * ay_cy);
   interval_number m12(tma_a - tma_b);
   interval_number mt1(m01 * bz_cz);
   interval_number mt2(m02 * by_cy);
   interval_number mt3(m12 * bx_cx);
   interval_number mtt(mt2 - mt1);
   interval_number m012(mtt - mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -m012.sign();
   else   return m012.sign();
}

int orient3d_indirect_LEEE_exact(implicitPoint3D_LPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dcx_p[64], *dcx = dcx_p;
   int dcx_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cx, &dcx, 64);
   double dcy_p[64], *dcy = dcy_p;
   int dcy_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cy, &dcy, 64);
   double dcz_p[64], *dcz = dcz_p;
   int dcz_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cz, &dcz, 64);
   double ix_cx_p[64], *ix_cx = ix_cx_p;
   int ix_cx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, dcx_len, dcx, &ix_cx, 64);
   double iy_cy_p[64], *iy_cy = iy_cy_p;
   int iy_cy_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, dcy_len, dcy, &iy_cy, 64);
   double ax_cx[2];
   o.two_Diff(ax, cx, ax_cx);
   double ay_cy[2];
   o.two_Diff(ay, cy, ay_cy);
   double az_cz[2];
   o.two_Diff(az, cz, az_cz);
   double iz_cz_p[64], *iz_cz = iz_cz_p;
   int iz_cz_len = o.Gen_Diff_With_PreAlloc(l1z_len, l1z, dcz_len, dcz, &iz_cz, 64);
   double bx_cx[2];
   o.two_Diff(bx, cx, bx_cx);
   double by_cy[2];
   o.two_Diff(by, cy, by_cy);
   double bz_cz[2];
   o.two_Diff(bz, cz, bz_cz);
   double tmc_a_p[64], *tmc_a = tmc_a_p;
   int tmc_a_len = o.Gen_Product_With_PreAlloc(ix_cx_len, ix_cx, 2, ay_cy, &tmc_a, 64);
   double tmc_b_p[64], *tmc_b = tmc_b_p;
   int tmc_b_len = o.Gen_Product_With_PreAlloc(iy_cy_len, iy_cy, 2, ax_cx, &tmc_b, 64);
   double m01_p[64], *m01 = m01_p;
   int m01_len = o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 64);
   double tmi_a_p[64], *tmi_a = tmi_a_p;
   int tmi_a_len = o.Gen_Product_With_PreAlloc(ix_cx_len, ix_cx, 2, az_cz, &tmi_a, 64);
   double tmi_b_p[64], *tmi_b = tmi_b_p;
   int tmi_b_len = o.Gen_Product_With_PreAlloc(iz_cz_len, iz_cz, 2, ax_cx, &tmi_b, 64);
   double m02_p[64], *m02 = m02_p;
   int m02_len = o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 64);
   double tma_a_p[64], *tma_a = tma_a_p;
   int tma_a_len = o.Gen_Product_With_PreAlloc(iy_cy_len, iy_cy, 2, az_cz, &tma_a, 64);
   double tma_b_p[64], *tma_b = tma_b_p;
   int tma_b_len = o.Gen_Product_With_PreAlloc(iz_cz_len, iz_cz, 2, ay_cy, &tma_b, 64);
   double m12_p[64], *m12 = m12_p;
   int m12_len = o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 64);
   double mt1_p[64], *mt1 = mt1_p;
   int mt1_len = o.Gen_Product_With_PreAlloc(m01_len, m01, 2, bz_cz, &mt1, 64);
   double mt2_p[64], *mt2 = mt2_p;
   int mt2_len = o.Gen_Product_With_PreAlloc(m02_len, m02, 2, by_cy, &mt2, 64);
   double mt3_p[64], *mt3 = mt3_p;
   int mt3_len = o.Gen_Product_With_PreAlloc(m12_len, m12, 2, bx_cx, &mt3, 64);
   double mtt_p[64], *mtt = mtt_p;
   int mtt_len = o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 64);
   double m012_p[64], *m012 = m012_p;
   int m012_len = o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 64);

   return_value = m012[m012_len - 1];
   if (m012_p != m012) free(m012);
   if (mtt_p != mtt) free(mtt);
   if (mt3_p != mt3) free(mt3);
   if (mt2_p != mt2) free(mt2);
   if (mt1_p != mt1) free(mt1);
   if (m12_p != m12) free(m12);
   if (tma_b_p != tma_b) free(tma_b);
   if (tma_a_p != tma_a) free(tma_a);
   if (m02_p != m02) free(m02);
   if (tmi_b_p != tmi_b) free(tmi_b);
   if (tmi_a_p != tmi_a) free(tmi_a);
   if (m01_p != m01) free(m01);
   if (tmc_b_p != tmc_b) free(tmc_b);
   if (tmc_a_p != tmc_a) free(tmc_a);
   if (iz_cz_p != iz_cz) free(iz_cz);
   if (iy_cy_p != iy_cy) free(iy_cy);
   if (ix_cx_p != ix_cx) free(ix_cx);
   if (dcz_p != dcz) free(dcz);
   if (dcy_p != dcy) free(dcy);
   if (dcx_p != dcx) free(dcx);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient3d_indirect_LEEE(implicitPoint3D_LPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   int ret;
   ret = orient3d_indirect_LEEE_filtered(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient3d_indirect_LEEE_interval(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_indirect_LEEE_exact(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
}

int orient3d_indirect_TEEE_filtered(implicitPoint3D_TPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   double l1x, l1y, l1z, d1, max_var = 0;
    if (
       !p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var)
   ) return Filtered_Sign::UNCERTAIN;

   double dcx = d1 * cx;
   double dcy = d1 * cy;
   double dcz = d1 * cz;
   double ix_cx = l1x - dcx;
   double iy_cy = l1y - dcy;
   double ax_cx = ax - cx;
   double ay_cy = ay - cy;
   double az_cz = az - cz;
   double iz_cz = l1z - dcz;
   double bx_cx = bx - cx;
   double by_cy = by - cy;
   double bz_cz = bz - cz;
   double tmc_a = ix_cx * ay_cy;
   double tmc_b = iy_cy * ax_cx;
   double m01 = tmc_a - tmc_b;
   double tmi_a = ix_cx * az_cz;
   double tmi_b = iz_cz * ax_cx;
   double m02 = tmi_a - tmi_b;
   double tma_a = iy_cy * az_cz;
   double tma_b = iz_cz * ay_cy;
   double m12 = tma_a - tma_b;
   double mt1 = m01 * bz_cz;
   double mt2 = m02 * by_cy;
   double mt3 = m12 * bx_cx;
   double mtt = mt2 - mt1;
   double m012 = mtt - mt3;

   double _tmp_fabs;
   if ((_tmp_fabs = fabs(cx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(cz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ax_cx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(ay_cy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(az_cz)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bx_cx)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(by_cy)) > max_var) max_var = _tmp_fabs;
   if ((_tmp_fabs = fabs(bz_cz)) > max_var) max_var = _tmp_fabs;
   double epsilon = max_var;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= epsilon;
   epsilon *= max_var;
   epsilon *= 3.070283610684406e-12;
   if (((d1 < 0))) m012 = -m012;
   if (m012 > epsilon) return IP_Sign::POSITIVE;
   if (-m012 > epsilon) return IP_Sign::NEGATIVE;
   return Filtered_Sign::UNCERTAIN;
}

int orient3d_indirect_TEEE_interval(implicitPoint3D_TPI& p1, interval_number ax, interval_number ay, interval_number az, interval_number bx, interval_number by, interval_number bz, interval_number cx, interval_number cy, interval_number cz)
{
   interval_number l1x, l1y, l1z, d1;
   if (
   !p1.getIntervalLambda(l1x, l1y, l1z, d1)
   ) return Filtered_Sign::UNCERTAIN;

   setFPUModeToRoundUP();
   interval_number dcx(d1 * cx);
   interval_number dcy(d1 * cy);
   interval_number dcz(d1 * cz);
   interval_number ix_cx(l1x - dcx);
   interval_number iy_cy(l1y - dcy);
   interval_number ax_cx(ax - cx);
   interval_number ay_cy(ay - cy);
   interval_number az_cz(az - cz);
   interval_number iz_cz(l1z - dcz);
   interval_number bx_cx(bx - cx);
   interval_number by_cy(by - cy);
   interval_number bz_cz(bz - cz);
   interval_number tmc_a(ix_cx * ay_cy);
   interval_number tmc_b(iy_cy * ax_cx);
   interval_number m01(tmc_a - tmc_b);
   interval_number tmi_a(ix_cx * az_cz);
   interval_number tmi_b(iz_cz * ax_cx);
   interval_number m02(tmi_a - tmi_b);
   interval_number tma_a(iy_cy * az_cz);
   interval_number tma_b(iz_cz * ay_cy);
   interval_number m12(tma_a - tma_b);
   interval_number mt1(m01 * bz_cz);
   interval_number mt2(m02 * by_cy);
   interval_number mt3(m12 * bx_cx);
   interval_number mtt(mt2 - mt1);
   interval_number m012(mtt - mt3);
   setFPUModeToRoundNEAR();

   if (!m012.signIsReliable()) return Filtered_Sign::UNCERTAIN;
   if (((d1 < 0))) return -m012.sign();
   else   return m012.sign();
}

int orient3d_indirect_TEEE_exact(implicitPoint3D_TPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
 double return_value = 0.0;
 double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64], *l1z = l1z_p, d1_p[64], *d1 = d1_p;
 int l1x_len, l1y_len, l1z_len, d1_len;
 p1.getExactLambda(l1x, l1x_len, l1y, l1y_len, l1z, l1z_len, d1, d1_len);
 if ((d1[d1_len - 1] != 0))
 {
   expansionObject o;
   double dcx_p[64], *dcx = dcx_p;
   int dcx_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cx, &dcx, 64);
   double dcy_p[64], *dcy = dcy_p;
   int dcy_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cy, &dcy, 64);
   double dcz_p[64], *dcz = dcz_p;
   int dcz_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cz, &dcz, 64);
   double ix_cx_p[64], *ix_cx = ix_cx_p;
   int ix_cx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, dcx_len, dcx, &ix_cx, 64);
   double iy_cy_p[64], *iy_cy = iy_cy_p;
   int iy_cy_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, dcy_len, dcy, &iy_cy, 64);
   double ax_cx[2];
   o.two_Diff(ax, cx, ax_cx);
   double ay_cy[2];
   o.two_Diff(ay, cy, ay_cy);
   double az_cz[2];
   o.two_Diff(az, cz, az_cz);
   double iz_cz_p[64], *iz_cz = iz_cz_p;
   int iz_cz_len = o.Gen_Diff_With_PreAlloc(l1z_len, l1z, dcz_len, dcz, &iz_cz, 64);
   double bx_cx[2];
   o.two_Diff(bx, cx, bx_cx);
   double by_cy[2];
   o.two_Diff(by, cy, by_cy);
   double bz_cz[2];
   o.two_Diff(bz, cz, bz_cz);
   double tmc_a_p[64], *tmc_a = tmc_a_p;
   int tmc_a_len = o.Gen_Product_With_PreAlloc(ix_cx_len, ix_cx, 2, ay_cy, &tmc_a, 64);
   double tmc_b_p[64], *tmc_b = tmc_b_p;
   int tmc_b_len = o.Gen_Product_With_PreAlloc(iy_cy_len, iy_cy, 2, ax_cx, &tmc_b, 64);
   double m01_p[64], *m01 = m01_p;
   int m01_len = o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 64);
   double tmi_a_p[64], *tmi_a = tmi_a_p;
   int tmi_a_len = o.Gen_Product_With_PreAlloc(ix_cx_len, ix_cx, 2, az_cz, &tmi_a, 64);
   double tmi_b_p[64], *tmi_b = tmi_b_p;
   int tmi_b_len = o.Gen_Product_With_PreAlloc(iz_cz_len, iz_cz, 2, ax_cx, &tmi_b, 64);
   double m02_p[64], *m02 = m02_p;
   int m02_len = o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 64);
   double tma_a_p[64], *tma_a = tma_a_p;
   int tma_a_len = o.Gen_Product_With_PreAlloc(iy_cy_len, iy_cy, 2, az_cz, &tma_a, 64);
   double tma_b_p[64], *tma_b = tma_b_p;
   int tma_b_len = o.Gen_Product_With_PreAlloc(iz_cz_len, iz_cz, 2, ay_cy, &tma_b, 64);
   double m12_p[64], *m12 = m12_p;
   int m12_len = o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 64);
   double mt1_p[64], *mt1 = mt1_p;
   int mt1_len = o.Gen_Product_With_PreAlloc(m01_len, m01, 2, bz_cz, &mt1, 64);
   double mt2_p[64], *mt2 = mt2_p;
   int mt2_len = o.Gen_Product_With_PreAlloc(m02_len, m02, 2, by_cy, &mt2, 64);
   double mt3_p[64], *mt3 = mt3_p;
   int mt3_len = o.Gen_Product_With_PreAlloc(m12_len, m12, 2, bx_cx, &mt3, 64);
   double mtt_p[64], *mtt = mtt_p;
   int mtt_len = o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 64);
   double m012_p[64], *m012 = m012_p;
   int m012_len = o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 64);

   return_value = m012[m012_len - 1];
   if (m012_p != m012) free(m012);
   if (mtt_p != mtt) free(mtt);
   if (mt3_p != mt3) free(mt3);
   if (mt2_p != mt2) free(mt2);
   if (mt1_p != mt1) free(mt1);
   if (m12_p != m12) free(m12);
   if (tma_b_p != tma_b) free(tma_b);
   if (tma_a_p != tma_a) free(tma_a);
   if (m02_p != m02) free(m02);
   if (tmi_b_p != tmi_b) free(tmi_b);
   if (tmi_a_p != tmi_a) free(tmi_a);
   if (m01_p != m01) free(m01);
   if (tmc_b_p != tmc_b) free(tmc_b);
   if (tmc_a_p != tmc_a) free(tmc_a);
   if (iz_cz_p != iz_cz) free(iz_cz);
   if (iy_cy_p != iy_cy) free(iy_cy);
   if (ix_cx_p != ix_cx) free(ix_cx);
   if (dcz_p != dcz) free(dcz);
   if (dcy_p != dcy) free(dcy);
   if (dcx_p != dcx) free(dcx);
   if (( (d1[d1_len -1] < 0))) return_value = -return_value;
 }

 if (l1x_p != l1x) free(l1x);
 if (l1y_p != l1y) free(l1y);
 if (l1z_p != l1z) free(l1z);
 if (d1_p != d1) free(d1);

 if (return_value > 0) return IP_Sign::POSITIVE;
 if (return_value < 0) return IP_Sign::NEGATIVE;
 return IP_Sign::ZERO;
}

int orient3d_indirect_TEEE(implicitPoint3D_TPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
{
   int ret;
   ret = orient3d_indirect_TEEE_filtered(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   ret = orient3d_indirect_TEEE_interval(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
   if (ret != Filtered_Sign::UNCERTAIN) return ret;
   return orient3d_indirect_TEEE_exact(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
}

