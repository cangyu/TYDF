#include <iostream>
#include <fstream>
#include <vector>
#include "../../inc/nmf.h"

using namespace GridTool;

void test(const std::string &dir_name)
{
	std::cout << "Case \"" << dir_name << "\" ..." << std::endl;
	const std::string dir_path("./" + dir_name + "/");
	NMF::Mapping3D mapping(dir_path + "map.nmf");
	std::ofstream fout(dir_path + "report.txt");
	mapping.summary(fout);
	fout.close();
	mapping.numbering();
	mapping.writeToFile(dir_path + "map_blessed.nmf");
}

int main(int argc, char *argv[])
{
	std::cout << "Test the \"Neutral Map File\" utilities." << std::endl;

	test("Cavity"); // Single block.
	test("Sky1"); // 2 blocks connected through 1 surface.
	test("Langley"); // 4 blocks in 2x2 form.

	return 0;
}
