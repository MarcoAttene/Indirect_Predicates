#include "implicit_point.h"

using namespace IPs;

int main(int argc, char *argv[])
{
	explicitPoint3D p1(0, 0, 0), p2(123, 122, 124);
	explicitPoint3D p3(1, 0, 0), p4(0, 1, 0), p5(0, 0, 1);
	implicitPoint3D_LPI l(p1, p2, p3, p4, p5);
	implicitPoint3D_BPT b(p3, p4, l, 0.1, 0.3);

	if (genericPoint::orient3D(l, b, p4, p5) == 0) std::cout << "Coplanar - Test succeeded.\n";
	else std::cout << "Not coplanar - Test failed.\n";

	return 0;
}
