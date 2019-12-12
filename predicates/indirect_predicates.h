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

int incircle(double pax, double pay, double pbx, double pby, double pcx, double pcy, double pdx, double pdy);
int orient2d(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y);
int orient3d(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz);
int incircle_indirect_LEEE(implicitPoint3D_LPI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy);
int incircle_indirect_LLEE(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double pcx, double pcy, double pdx, double pdy);
int incircle_indirect_LLLE(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, double pdx, double pdy);
int incircle_indirect_LLLL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3, implicitPoint3D_LPI& p4);
int incircle_indirect_SEEE(implicitPoint2D_SSI& p1, double pbx, double pby, double pcx, double pcy, double pdx, double pdy);
int incircle_indirect_SSEE(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double pcx, double pcy, double pdx, double pdy);
int incircle_indirect_SSSE(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, double pdx, double pdy);
int incircle_indirect_SSSS(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3, implicitPoint2D_SSI& p4);
bool lambda2d_SSI_filtered(double ea1x, double ea1y, double ea2x, double ea2y, double eb1x, double eb1y, double eb2x, double eb2y, double& lambda_x, double& lambda_y, double& lambda_det, double& max_var);
bool lambda2d_SSI_interval(interval_number ea1x, interval_number ea1y, interval_number ea2x, interval_number ea2y, interval_number eb1x, interval_number eb1y, interval_number eb2x, interval_number eb2y, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_det);
void lambda2d_SSI_exact(double ea1x, double ea1y, double ea2x, double ea2y, double eb1x, double eb1y, double eb2x, double eb2y, double *lambda_x, int& lambda_x_len, double *lambda_y, int& lambda_y_len, double *lambda_det, int& lambda_det_len);
bool lambda3d_LPI_filtered(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz, double& lambda_d, double& lambda_x, double& lambda_y, double& lambda_z, double& max_var);
bool lambda3d_LPI_interval(interval_number px, interval_number py, interval_number pz, interval_number qx, interval_number qy, interval_number qz, interval_number rx, interval_number ry, interval_number rz, interval_number sx, interval_number sy, interval_number sz, interval_number tx, interval_number ty, interval_number tz, interval_number& lambda_d, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z);
void lambda3d_LPI_exact(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz, double tx, double ty, double tz, double *lambda_d, int& lambda_d_len, double *lambda_x, int& lambda_x_len, double *lambda_y, int& lambda_y_len, double *lambda_z, int& lambda_z_len);
bool lambda3d_TPI_filtered(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z, double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z, double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z, double ou3x, double ou3y, double ou3z, double& lambda_x, double& lambda_y, double& lambda_z, double& lambda_d, double& max_var);
bool lambda3d_TPI_interval(interval_number ov1x, interval_number ov1y, interval_number ov1z, interval_number ov2x, interval_number ov2y, interval_number ov2z, interval_number ov3x, interval_number ov3y, interval_number ov3z, interval_number ow1x, interval_number ow1y, interval_number ow1z, interval_number ow2x, interval_number ow2y, interval_number ow2z, interval_number ow3x, interval_number ow3y, interval_number ow3z, interval_number ou1x, interval_number ou1y, interval_number ou1z, interval_number ou2x, interval_number ou2y, interval_number ou2z, interval_number ou3x, interval_number ou3y, interval_number ou3z, interval_number& lambda_x, interval_number& lambda_y, interval_number& lambda_z, interval_number& lambda_d);
void lambda3d_TPI_exact(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z, double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z, double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z, double ou3x, double ou3y, double ou3z, double *lambda_x, int& lambda_x_len, double *lambda_y, int& lambda_y_len, double *lambda_z, int& lambda_z_len, double *lambda_d, int& lambda_d_len);
int lessThanOnX_LE(implicitPoint3D_LPI& p1, double bx);
int lessThanOnX_LL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2);
int lessThanOnX_LT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2);
int lessThanOnX_TE(implicitPoint3D_TPI& p1, double bx);
int lessThanOnX_TT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2);
int lessThanOnY_LE(implicitPoint3D_LPI& p1, double by);
int lessThanOnY_LL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2);
int lessThanOnY_LT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2);
int lessThanOnY_TE(implicitPoint3D_TPI& p1, double by);
int lessThanOnY_TT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2);
int lessThanOnZ_LE(implicitPoint3D_LPI& p1, double bz);
int lessThanOnZ_LL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2);
int lessThanOnZ_LT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2);
int lessThanOnZ_TE(implicitPoint3D_TPI& p1, double bz);
int lessThanOnZ_TT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2);
int orient2d3d_indirect_LEE(implicitPoint3D_LPI& p1, double p2x, double p2y, double p3x, double p3y);
int orient2d3d_indirect_LLE(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, double op3x, double op3y);
int orient2d3d_indirect_LLL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_LPI& p3);
int orient2d3d_indirect_LLT(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2, implicitPoint3D_TPI& p3);
int orient2d3d_indirect_LTE(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y);
int orient2d3d_indirect_LTT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3);
int orient2d3d_indirect_TEE(implicitPoint3D_TPI& p1, double p2x, double p2y, double p3x, double p3y);
int orient2d3d_indirect_TTE(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, double p3x, double p3y);
int orient2d3d_indirect_TTT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2, implicitPoint3D_TPI& p3);
int orient2d_indirect_SEE(implicitPoint2D_SSI& p1, double p2x, double p2y, double p3x, double p3y);
int orient2d_indirect_SSE(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, double p3x, double p3y);
int orient2d_indirect_SSS(implicitPoint2D_SSI& p1, implicitPoint2D_SSI& p2, implicitPoint2D_SSI& p3);
int orient3d_indirect_LEEE(implicitPoint3D_LPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz);
int orient3d_indirect_TEEE(implicitPoint3D_TPI& p1, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz);
