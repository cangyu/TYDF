#include <iostream>
#include "../../inc/xf.h"

using namespace GridTool;

void test(const std::string &dir_name)
{
	std::cout << "Case \"" << dir_name << "\" ..." << std::endl;
	std::ofstream fout(dir_name + "/report.txt");
	if (fout.fail())
		throw std::runtime_error("Failed to open report file.");

	XF::MESH msh(dir_name + "/fluent.msh", fout);
	fout.close();

	msh.writeToFile(dir_name + "/blessed.msh");
}

int main(int argc, char *argv[])
{
	std::cout << "Testing I/O of \"FLUENT\" unstructured mesh ..." << std::endl;

	test("Cube1"); // 3D tet

	return 0;
}
