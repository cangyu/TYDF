#include <iostream>
#include <fstream>
#include <vector>
#include "../../inc/nmf.h"

using namespace GridTool;

const std::vector<std::string> test
{
	"Cavity", // Single block.
	"Sky1", // 2 blocks connected through 1 surface.
	"Langley", // 4 blocks in 2x2 form.
};

int main(int argc, char *argv[])
{
	int cnt = 0;
	int failure = 0;
	std::cout << "Testing cases for the Neutral Map File utilities." << std::endl;
	for (const auto &dir_name : test)
	{
		std::cout << "Case " << ++cnt << " ... " << std::endl;
		try
		{
			const std::string dir_path("./" + dir_name + "/");
			NMF::Mapping3D mapping(dir_path + "map.nmf");
			std::ofstream fout(dir_path + "report.txt");
			mapping.summary(fout);
			fout.close();
			mapping.numbering();
			mapping.writeToFile(dir_path + "map_blessed.nmf");
		}
		catch (std::exception &e)
		{
			std::cout << e.what() << std::endl;
			++failure;
		}
	}
	std::cout << "Done with " << failure << " failed." << std::endl;
	return 0;
}
