#include <iostream>
#include <fstream>
#include <vector>
#include "nmf.hpp"

// 4 blocks in 2x2 form.
void case1()
{
	const std::string dir("./Langley GeoLab/");
	NMF::Mapping3D mapping(dir + "map.nmf");
	std::ofstream fout(dir + "report.txt");
	mapping.summary(fout);
	fout.close();
	mapping.writeToFile(dir + "map_blessed.nmf");
}

// Simple 2 blocks connected through 1 face.
void case2()
{
	const std::string dir("./Sky1/");
	NMF::Mapping3D mapping(dir + "map.nmf");
	std::ofstream fout(dir + "report.txt");
	mapping.summary(fout);
	fout.close();
	mapping.writeToFile(dir + "map_blessed.nmf");
}

typedef void(*pTestFunction)(void);
const std::vector<pTestFunction> func{
	case1,
	case2
};

int main(int argc, char *argv[])
{
	int cnt = 0;
	int failure = 0;
	std::cout << "Testing cases for the Neutral Map File utilities." << std::endl;
	for (auto f : func)
	{
		std::cout << "Case " << ++cnt << " ... " << std::endl;
		try {
			f();
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
