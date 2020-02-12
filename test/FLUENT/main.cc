#include <iostream>
#include "../../inc/xf.h"

using namespace GridTool;

void test(const std::string &case_name, const std::string &case_dir)
{
	std::cout << "Case \"" << case_name << "\" ..." << std::endl;
	
	std::ofstream fout(case_dir + "report.txt");
	if (fout.fail())
		throw std::runtime_error("Failed to open report file.");

	XF::MESH msh(case_dir + "fluent.msh", fout);
	fout.close();

	msh.writeToFile(case_dir + "blessed.msh");
}

int main(int argc, char *argv[])
{
	std::cout << "Testing I/O of \"FLUENT\" unstructured mesh ..." << std::endl;

	test("Cavity64", "../../case/Cavity/FLUENT/");
	test("Structure1 from LiuSha tutorial", "../../case/LS1/FLUENT/");

	return 0;
}
