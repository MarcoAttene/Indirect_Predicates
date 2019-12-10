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

#include "implicit_point.h"
#include "indirect_predicates.h"

int lessThan_EE(double x1, double y1, double z1, double x2, double y2, double z2)
{
	int ret;
	if ((ret = ((x2 > x1) - (x2 < x1)))) return ret;
	if ((ret = ((y2 > y1) - (y2 < y1)))) return ret;
	return ((z2 > z1) - (z2 < z1));
}

int lessThan_LE(implicitPoint3D_LPI& p1, double x2, double y2, double z2)
{
	int ret;
	if ((ret = lessThanOnX_LE(p1, x2))) return ret;
	if ((ret = lessThanOnY_LE(p1, y2))) return ret;
	return lessThanOnZ_LE(p1, z2);
}

int lessThan_TE(implicitPoint3D_TPI& p1, double x2, double y2, double z2)
{
	int ret;
	if ((ret = lessThanOnX_TE(p1, x2))) return ret;
	if ((ret = lessThanOnY_TE(p1, y2))) return ret;
	return lessThanOnZ_TE(p1, z2);
}

int lessThan_LL(implicitPoint3D_LPI& p1, implicitPoint3D_LPI& p2)
{
	int ret;
	if ((ret = lessThanOnX_LL(p1, p2))) return ret;
	if ((ret = lessThanOnY_LL(p1, p2))) return ret;
	return lessThanOnZ_LL(p1, p2);
}

int lessThan_TT(implicitPoint3D_TPI& p1, implicitPoint3D_TPI& p2)
{
	int ret;
	if ((ret = lessThanOnX_TT(p1, p2))) return ret;
	if ((ret = lessThanOnY_TT(p1, p2))) return ret;
	return lessThanOnZ_TT(p1, p2);
}

int lessThan_LT(implicitPoint3D_LPI& p1, implicitPoint3D_TPI& p2)
{
	int ret;
	if ((ret = lessThanOnX_LT(p1, p2))) return ret;
	if ((ret = lessThanOnY_LT(p1, p2))) return ret;
	return lessThanOnZ_LT(p1, p2);
}

int genericPoint::lessThan(genericPoint& a, genericPoint& b)
{
	if (a.is3D()) // Here we may want to check that b is also 3D, but for now we just assume it is
	{
		int config = (a.getType() - Point_Type::EXPLICIT3D) + ((b.getType() - Point_Type::EXPLICIT3D) << 2);
		switch (config)
		{
		case 0: // EE
			return lessThan_EE(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z());
		case 1: // LE
			return lessThan_LE(a.toLPI(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z());
		case 2: // TE
			return lessThan_TE(a.toTPI(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z());
		case 4: // EL
			return -lessThan_LE(b.toLPI(), a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z());
		case 5: // LL
			return lessThan_LL(a.toLPI(), b.toLPI());
		case 6: // TL
			return -lessThan_LT(b.toLPI(), a.toTPI());
		case 8:
			return -lessThan_TE(b.toTPI(), a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z());
		case 9: // LT
			return lessThan_LT(a.toLPI(), b.toTPI());
		case 10: // TT
			return lessThan_TT(a.toTPI(), b.toTPI());
		}
	}
	ip_error("genericPoint::lessThan - points must be 3D\n");
}

int genericPoint::orient2D(genericPoint& a, genericPoint& b, genericPoint& c)
{
	if (a.is2D()) // Here we may want to check that b and c are also 2D, but for now we just assume they are
	{	
		if (a.isExplicit2D() && b.isExplicit2D() && c.isExplicit2D())
			return orient2d(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y());
		if (a.isSSI() && b.isExplicit2D() && c.isExplicit2D())
			return orient2d_indirect_SEE(a.toSSI(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y());
		if (a.isSSI() && b.isSSI() && c.isExplicit2D())
			return orient2d_indirect_SSE(a.toSSI(), b.toSSI(), c.toExplicit2D().X(), c.toExplicit2D().Y());
		if (a.isSSI() && b.isSSI() && c.isSSI())
			return orient2d_indirect_SSS(a.toSSI(), b.toSSI(), c.toSSI());
		if (a.isExplicit2D() && b.isSSI() && c.isExplicit2D()) //return orient2D(b, c, a);
			return orient2d_indirect_SEE(b.toSSI(), c.toExplicit2D().X(), c.toExplicit2D().Y(), a.toExplicit2D().X(), a.toExplicit2D().Y());
		if (a.isExplicit2D() && b.isSSI() && c.isSSI()) //return orient2D(b, c, a);
			return orient2d_indirect_SSE(b.toSSI(), c.toSSI(), a.toExplicit2D().X(), a.toExplicit2D().Y());
		if (a.isExplicit2D() && b.isExplicit2D() && c.isSSI()) //return orient2D(c, a, b);
			return orient2d_indirect_SEE(c.toSSI(), a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y());
		if (a.isSSI() && b.isExplicit2D() && c.isSSI()) //return orient2D(c, a, b);
			return orient2d_indirect_SSE(c.toSSI(), a.toSSI(), b.toExplicit2D().X(), b.toExplicit2D().Y());
	}
	
	if (a.is3D()) // Here we may want to check that b and c are also 3D, but for now we just assume they are
	{		
		if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D())
			return orient2d(a.toExplicit3D().X(), a.toExplicit3D().Y(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y());
		if (a.isLPI() && b.isExplicit3D() && c.isExplicit3D())
			return orient2d3d_indirect_LEE(a.toLPI(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y());
		if (a.isLPI() && b.isLPI() && c.isExplicit3D())
			return orient2d3d_indirect_LLE(a.toLPI(), b.toLPI(), c.toExplicit3D().X(), c.toExplicit3D().Y());
		if (a.isLPI() && b.isLPI() && c.isLPI())
			return orient2d3d_indirect_LLL(a.toLPI(), b.toLPI(), c.toLPI());

		if (a.isExplicit3D() && b.isLPI() && c.isExplicit3D()) return orient2D(b, c, a);
		if (a.isExplicit3D() && b.isLPI() && c.isLPI()) return orient2D(b, c, a);
		if (a.isExplicit3D() && b.isExplicit3D() && c.isLPI()) return orient2D(c, a, b);
		if (a.isLPI() && b.isExplicit3D() && c.isLPI()) return orient2D(c, a, b);

		if (a.isTPI() && b.isExplicit3D() && c.isExplicit3D())
			return orient2d3d_indirect_TEE(a.toTPI(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y());
		if (a.isTPI() && b.isTPI() && c.isExplicit3D())
			return orient2d3d_indirect_TTE(a.toTPI(), b.toTPI(), c.toExplicit3D().X(), c.toExplicit3D().Y());
		if (a.isTPI() && b.isTPI() && c.isTPI())
			return orient2d3d_indirect_TTT(a.toTPI(), b.toTPI(), c.toTPI());

		if (a.isExplicit3D() && b.isTPI() && c.isExplicit3D()) return orient2D(b, c, a);
		if (a.isExplicit3D() && b.isTPI() && c.isTPI()) return orient2D(b, c, a);
		if (a.isExplicit3D() && b.isExplicit3D() && c.isTPI()) return orient2D(c, a, b);
		if (a.isTPI() && b.isExplicit3D() && c.isTPI()) return orient2D(c, a, b);

		if (a.isLPI() && b.isTPI() && c.isExplicit3D())
			return orient2d3d_indirect_LTE(a.toLPI(), b.toTPI(), c.toExplicit3D().X(), c.toExplicit3D().Y());
		if (a.isLPI() && b.isLPI() && c.isTPI())
			return orient2d3d_indirect_LLT(a.toLPI(), b.toLPI(), c.toTPI());
		if (a.isLPI() && b.isTPI() && c.isTPI())
			return orient2d3d_indirect_LTT(a.toLPI(), b.toTPI(), c.toTPI());

		if (a.isExplicit3D() && b.isLPI() && c.isTPI()) return orient2D(b, c, a);
		if (a.isExplicit3D() && b.isTPI() && c.isLPI()) return -orient2D(c, b, a);
		if (a.isLPI() && b.isExplicit3D() && c.isTPI()) return -orient2D(a, c, b);
		if (a.isTPI() && b.isExplicit3D() && c.isLPI()) return orient2D(c, a, b);
		if (a.isTPI() && b.isLPI() && c.isExplicit3D()) return -orient2D(b, a, c);

		if (a.isLPI() && b.isTPI() && c.isLPI()) return orient2D(c, a, b);
		if (a.isTPI() && b.isLPI() && c.isLPI()) return orient2D(b, c, a);

		if (a.isTPI() && b.isLPI() && c.isTPI()) return orient2D(b, c, a);
		if (a.isTPI() && b.isTPI() && c.isLPI()) return orient2D(c, a, b);

		ip_error("genericPoint::orient2d - should not happen\n");
	}

	ip_error("genericPoint::orient2d - points must be 2D or 3D\n");
	return 0;
}

int genericPoint::orient3D(genericPoint& a, genericPoint& b, genericPoint& c, genericPoint& d)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D() && d.isExplicit3D())
		return orient3d(
			a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(),
			b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
			c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
			d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
	if (a.isLPI() && b.isExplicit3D() && c.isExplicit3D() && d.isExplicit3D())
		return orient3d_indirect_LEEE(a.toLPI(),
			b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
			c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
			d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
	if (a.isTPI() && b.isExplicit3D() && c.isExplicit3D() && d.isExplicit3D())
		return orient3d_indirect_TEEE(a.toTPI(),
			b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
			c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
			d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D() && d.isExplicit3D()) return -orient3D(b, c, d, a);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D() && d.isExplicit3D()) return orient3D(c, d, a, b);
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D() && !d.isExplicit3D()) return -orient3D(d, a, b, c);

	ip_error("genericPoint::orient3d - should not happen/unsupported\n");
	return 0;
}

int genericPoint::incircle(genericPoint& a, genericPoint& b, genericPoint& c, genericPoint& d)
{
	if (a.is2D())
	{
		if (a.isExplicit2D() && b.isExplicit2D() && c.isExplicit2D() && d.isExplicit2D())
			return ::incircle(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
		if (a.isSSI() && b.isExplicit2D() && c.isExplicit2D() && d.isExplicit2D())
			return incircle_indirect_SEEE(a.toSSI(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
		if (a.isSSI() && b.isSSI() && c.isExplicit2D() && d.isExplicit2D())
			return incircle_indirect_SSEE(a.toSSI(), b.toSSI(), c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
		if (a.isSSI() && b.isSSI() && c.isSSI() && d.isExplicit2D())
			return incircle_indirect_SSSE(a.toSSI(), b.toSSI(), c.toSSI(), d.toExplicit2D().X(), d.toExplicit2D().Y());
		if (a.isSSI() && b.isSSI() && c.isSSI() && d.isSSI())
			return incircle_indirect_SSSS(a.toSSI(), b.toSSI(), c.toSSI(), d.toSSI());

		if (a.isExplicit2D() && b.isExplicit2D() && c.isExplicit2D() && d.isSSI()) return -incircle(d, a, b, c);
		if (a.isExplicit2D() && b.isExplicit2D() && c.isSSI() && d.isExplicit2D()) return incircle(c, a, b, d);
		if (a.isExplicit2D() && b.isExplicit2D() && c.isSSI() && d.isSSI()) return incircle(c, d, a, b);
		if (a.isExplicit2D() && b.isSSI() && c.isExplicit2D() && d.isExplicit2D()) return -incircle(b, a, c, d);
		if (a.isExplicit2D() && b.isSSI() && c.isExplicit2D() && d.isSSI()) return -incircle(b, d, a, c);
		if (a.isExplicit2D() && b.isSSI() && c.isSSI() && d.isExplicit2D()) return incircle(b, c, a, d);
		if (a.isExplicit2D() && b.isSSI() && c.isSSI() && d.isSSI()) return -incircle(b, c, d, a);
		if (a.isSSI() && b.isExplicit2D() && c.isExplicit2D() && d.isSSI()) return incircle(a, d, b, c);
		if (a.isSSI() && b.isExplicit2D() && c.isSSI() && d.isExplicit2D()) return -incircle(a, c, b, d);
		if (a.isSSI() && b.isExplicit2D() && c.isSSI() && d.isSSI()) return -incircle(a, c, b, d);
		if (a.isSSI() && b.isSSI() && c.isExplicit2D() && d.isSSI()) return -incircle(a, b, d, c);

		ip_error("genericPoint::incircle - should not happen\n");
	}

	if (a.is3D())
	{
		if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D() && d.isExplicit3D())
			return ::incircle(a.toExplicit3D().X(), a.toExplicit3D().Y(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
		if (a.isLPI() && b.isExplicit3D() && c.isExplicit3D() && d.isExplicit3D())
			return incircle_indirect_LEEE(a.toLPI(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
		if (a.isLPI() && b.isLPI() && c.isExplicit3D() && d.isExplicit3D())
			return incircle_indirect_LLEE(a.toLPI(), b.toLPI(), c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
		if (a.isLPI() && b.isLPI() && c.isLPI() && d.isExplicit3D())
			return incircle_indirect_LLLE(a.toLPI(), b.toLPI(), c.toLPI(), d.toExplicit3D().X(), d.toExplicit3D().Y());
		if (a.isLPI() && b.isLPI() && c.isLPI() && d.isLPI())
			return incircle_indirect_LLLL(a.toLPI(), b.toLPI(), c.toLPI(), d.toLPI());

		if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D() && d.isLPI()) return -incircle(d, a, b, c);
		if (a.isExplicit3D() && b.isExplicit3D() && c.isLPI() && d.isExplicit3D()) return incircle(c, a, b, d);
		if (a.isExplicit3D() && b.isExplicit3D() && c.isLPI() && d.isLPI()) return incircle(c, d, a, b);
		if (a.isExplicit3D() && b.isLPI() && c.isExplicit3D() && d.isExplicit3D()) return -incircle(b, a, c, d);
		if (a.isExplicit3D() && b.isLPI() && c.isExplicit3D() && d.isLPI()) return -incircle(b, d, a, c);
		if (a.isExplicit3D() && b.isLPI() && c.isLPI() && d.isExplicit3D()) return incircle(b, c, a, d);
		if (a.isExplicit3D() && b.isLPI() && c.isLPI() && d.isLPI()) return -incircle(b, c, d, a);
		if (a.isLPI() && b.isExplicit3D() && c.isExplicit3D() && d.isLPI()) return incircle(a, d, b, c);
		if (a.isLPI() && b.isExplicit3D() && c.isLPI() && d.isExplicit3D()) return -incircle(a, c, b, d);
		if (a.isLPI() && b.isExplicit3D() && c.isLPI() && d.isLPI()) return -incircle(a, c, b, d);
		if (a.isLPI() && b.isLPI() && c.isExplicit3D() && d.isLPI()) return -incircle(a, b, d, c);

		ip_error("genericPoint::incircle - unsupported configuration (TPIs not implemented yet)\n");
	}

	ip_error("genericPoint::incircle - should not happen\n");
	return 0;
}


#ifdef USE_CACHED_VALUES

bool implicitPoint2D_SSI::getFilteredLambda(double& lx, double& ly, double &d, double& mv)
{
	if (needsFilteredLambda())
		lambda2d_SSI_filtered(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_denominator, ssfilter_max_val);

	lx = ssfilter_lambda_x;
	ly = ssfilter_lambda_y;
	d = ssfilter_denominator;
	if (ssfilter_max_val > mv) mv = ssfilter_max_val;
	return (ssfilter_denominator != 0);
}

bool implicitPoint2D_SSI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number &d)
{
	if (needsIntervalLambda())
		lambda2d_SSI_interval(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), dfilter_lambda_x, dfilter_lambda_y, dfilter_denominator);

	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

bool implicitPoint3D_LPI::getFilteredLambda(double& lx, double& ly, double& lz, double &d, double& mv)
{
	if (needsFilteredLambda())
		lambda3d_LPI_filtered(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), ssfilter_denominator, ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_lambda_z, ssfilter_max_val);

	lx = ssfilter_lambda_x;
	ly = ssfilter_lambda_y;
	lz = ssfilter_lambda_z;
	d = ssfilter_denominator;
	if (ssfilter_max_val > mv) mv = ssfilter_max_val;
	return (ssfilter_denominator != 0);
}

bool implicitPoint3D_LPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number &d)
{
	if (needsIntervalLambda())
		lambda3d_LPI_interval(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), dfilter_denominator, dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z);

	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

bool implicitPoint3D_TPI::getFilteredLambda(double& lx, double& ly, double& lz, double &d, double& mv)
{
	if (needsFilteredLambda())
		lambda3d_TPI_filtered(
		V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(),
		ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_lambda_z, ssfilter_denominator, ssfilter_max_val);

	lx = ssfilter_lambda_x;
	ly = ssfilter_lambda_y;
	lz = ssfilter_lambda_z;
	d = ssfilter_denominator;
	if (ssfilter_max_val > mv) mv = ssfilter_max_val;
	return (ssfilter_denominator != 0);
}

bool implicitPoint3D_TPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number &d)
{
	if (needsIntervalLambda())
		lambda3d_TPI_interval(
		V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(),
		dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);

	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

#else
bool implicitPoint2D_SSI::getFilteredLambda(double& lx, double& ly, double &d, double& mv)
{
	return lambda2d_SSI_filtered(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), lx, ly, d, mv);
}

bool implicitPoint2D_SSI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number &d)
{
	return lambda2d_SSI_interval(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), lx, ly, d);
}

bool implicitPoint3D_LPI::getFilteredLambda(double& lx, double& ly, double& lz, double& d, double& mv)
{
	return
		lambda3d_LPI_filtered(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), d, lx, ly, lz, mv);
}

bool implicitPoint3D_LPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d)
{
	return
		lambda3d_LPI_interval(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), d, lx, ly, lz);
}

bool implicitPoint3D_TPI::getFilteredLambda(double& lx, double& ly, double& lz, double& d, double& mv)
{
	return lambda3d_TPI_filtered(
		V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(),
		lx, ly, lz, d, mv);
}

bool implicitPoint3D_TPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d)
{
	return lambda3d_TPI_interval(
		V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(),
		lx, ly, lz, d);
}

#endif

void implicitPoint2D_SSI::getExactLambda(double *lx, int& lxl, double *ly, int& lyl, double *d, int& dl)
{
	lambda2d_SSI_exact(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), lx, lxl, ly, lyl, d, dl);
}

void implicitPoint3D_LPI::getExactLambda(double *lx, int& lxl, double *ly, int& lyl, double *lz, int& lzl, double *d, int& dl)
{
	lambda3d_LPI_exact(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), d, dl, lx, lxl, ly, lyl, lz, lzl);
}

void implicitPoint3D_TPI::getExactLambda(double *lx, int& lxl, double *ly, int& lyl, double *lz, int& lzl, double *d, int& dl)
{
	lambda3d_TPI_exact(V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(), lx, lxl, ly, lyl, lz, lzl, d, dl);
}


bool implicitPoint2D_SSI::approxExplicit(explicitPoint2D& e) const
{
	double lambda_x, lambda_y, lambda_d, max_var;
	if (!lambda2d_SSI_filtered(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), lambda_x, lambda_y, lambda_d, max_var))
		return false;
	e = explicitPoint2D(lambda_x / lambda_d, lambda_y / lambda_d);
	return true;
}

bool implicitPoint3D_LPI::approxExplicit(explicitPoint3D& e) const
{
	double lambda_x, lambda_y, lambda_z, lambda_d, max_var;
	if (!lambda3d_LPI_filtered(ip.X(), ip.Y(), ip.Z(), iq.X(), iq.Y(), iq.Z(), ir.X(), ir.Y(), ir.Z(), is.X(), is.Y(), is.Z(), it.X(), it.Y(), it.Z(), lambda_d, lambda_x, lambda_y, lambda_z, max_var))
		return false;
	e = explicitPoint3D(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
	return true;
}

bool implicitPoint3D_TPI::approxExplicit(explicitPoint3D& e) const
{
	double lambda_x, lambda_y, lambda_z, lambda_d, max_var;
	if (!lambda3d_TPI_filtered(iv1.X(), iv1.Y(), iv1.Z(), iv2.X(), iv2.Y(), iv2.Z(), iv3.X(), iv3.Y(), iv3.Z(), iw1.X(), iw1.Y(), iw1.Z(), iw2.X(), iw2.Y(), iw2.Z(), iw3.X(), iw3.Y(), iw3.Z(), iu1.X(), iu1.Y(), iu1.Z(), iu2.X(), iu2.Y(), iu2.Z(), iu3.X(), iu3.Y(), iu3.Z(), lambda_x, lambda_y, lambda_z, lambda_d, max_var))
		return false;
	e = explicitPoint3D(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
	return true;
}

bool genericPoint::getApproxXYCoordinates(double& x, double& y) const
{
	if (is2D())
	{
		explicitPoint2D op;
		const explicitPoint2D* p = &op;
		if (isExplicit2D()) p = &toExplicit2D();
		else if (isSSI()) { if (!toSSI().approxExplicit(op)) return false; }
		x = p->X(); y = p->Y();
		return true;
	}
	if (is3D())
	{
		explicitPoint3D op;
		const explicitPoint3D* p = &op;
		if (isExplicit3D()) p = &toExplicit3D();
		else if (isLPI()) { if (!toLPI().approxExplicit(op)) return false; }
		else if (isTPI()) { if (!toTPI().approxExplicit(op)) return false; }
		x = p->X(); y = p->Y();
		return true;
	}
	ip_error("genericPoint::getApproxXYCoordinates - should not happen\n");
	return false;
}

ostream& operator<<(ostream& os, const genericPoint& p)
{
	if (p.isExplicit2D()) return os << p.toExplicit2D();
	else if (p.isExplicit3D()) return os << p.toExplicit3D();
	else if (p.isSSI()) return os << p.toSSI();
	else if (p.isLPI()) return os << p.toLPI();
	else if (p.isTPI()) return os << p.toTPI();
	else ip_error("genericPoint::operator<< - should not happen\n");
	return os;
}
