#include "implicit_point.h"

int main(int argc, char *argv[])
{
	explicitPoint3D p1(0, 0, 0), p2(123, 122, 124);
	explicitPoint3D p3(1, 0, 0), p4(0, 1, 0), p5(0, 0, 1);
	implicitPoint3D_LPI l(p1, p2, p3, p4, p5);
	implicitPoint3D_BPT b(p3, p4, l, 0.1, 0.3);

	if (genericPoint::orient3D(l, b, p4, p5) == 0) std::cout << "Coplanar - Test succeeded.\n";
	else std::cout << "Not coplanar - Test failed.\n";

	//explicitPoint2D a(1, 1), b(3, 3), c(2, 1), d(1, 2);
	//implicitPoint2D_SSI i(a, b, c, d);

	//if (genericPoint::orient2D(a, i, b) == 0) std::cout << "Collinear - Test succeeded.\n";
	//else std::cout << "Not collinear - Test failed.\n";

	return 0;
}
