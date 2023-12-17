#include "implicit_point.h"

int main(int argc, char *argv[])
{
	initFPU();

	explicitPoint2D a(1, 1), b(3, 3), c(2, 1), d(1, 2);
	implicitPoint2D_SSI i(a, b, c, d);

	if (genericPoint::orient2D(a, i, b) == 0) std::cout << "Collinear\n";
	else std::cout << "Not collinear\n";

	return 0;
}
