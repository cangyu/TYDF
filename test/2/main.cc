#include "plot3d.hpp"

// 2D grid in 3D form
void test0()
{
	PLOT3D::GRID p3d;
	int ret = 0;

	ret = p3d.readFromFile("xyz1.fmt");
	if (ret != 0)
		throw std::runtime_error("Failed to read grid from file.");

	ret = p3d.writeToFile("xyz1_blessed.fmt");
	if (ret != 0)
		throw std::runtime_error("Failed to write grid to file.");
}

// 2D grid in 2D form
void test1()
{
	PLOT3D::GRID p3d;
	int ret = 0;

	ret = p3d.readFromFile("xyz2.fmt");
	if (ret != 0)
		throw std::runtime_error("Failed to read grid from file.");

	ret = p3d.writeToFile("xyz2_blessed.fmt");
	if (ret != 0)
		throw std::runtime_error("Failed to write grid to file.");
}

int main(int argc, char *argv[])
{
	test0();
	test1();

	return 0;
}

