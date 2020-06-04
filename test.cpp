#include "implicit_point.h"

void main(int argc, char *argv[])
{
	initFPU();

	explicitPoint2D a(1, 1), b(3, 3), c(3, 1), d(1, 3);
	implicitPoint2D_SSI i(a, b, c, d);

	if (genericPoint::orient2D(a, i, b) == 0) std::cout << "Collinear\n";
	else std::cout << "Not collinear\n";
}
