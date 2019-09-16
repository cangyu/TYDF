#include "plot3d.h"

// 3D single-block grid
void test1()
{
	PLOT3D p3d;
	int ret = 0;

	ret = p3d.readFromFile("xyz.fmt");
	if (ret != 0)
		throw std::runtime_error("Failed to read grid from file.");

	ret = p3d.writeToFile("xyz_blessed.fmt");
	if (ret != 0)
		throw std::runtime_error("Failed to write grid to file.");
}

int main(int argc, char *argv[])
{
	test1();

	return 0;
}

