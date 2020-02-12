#include <iostream>
#include "../../inc/plot3d.h"

using namespace GridTool;

void test(const std::string &dir_name)
{
	std::cout << "\tCase \"" << dir_name << "\" ..." << std::endl;
	PLOT3D::GRID p3d(dir_name + "/xyz.fmt");
	p3d.writeToFile(dir_name + "/xyz_blessed.fmt");
}

int main(int argc, char *argv[])
{
	std::cout << "Testing I/O of \"PLOT3D\" grid ..." << std::endl;

	test("Planar1"); // 2D grid in 2D form
	test("Shell1"); // 2D grid in 3D form
	test("Cube1"); // 3D single-block grid

	std::cout << "Done!" << std::endl;
	return 0;
}
