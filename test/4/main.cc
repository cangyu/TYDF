#include <iostream>
#include <vector>
#include <fstream>
#include "nmf.hpp"

void LangleyExample()
{
	// 4 blocks in 2x2 form.
	NMF::Mapping3D mapping("0/map.nmf");
	std::ofstream fout("0/report.txt");
	mapping.summary(fout);
	fout.close();
	mapping.writeToFile("0/map_blessed.nmf");
}

void TwoBlocks()
{
	NMF::Mapping3D mapping("1/map.nmf");
	std::ofstream fout("1/report.txt");
	mapping.summary(fout);
	fout.close();
	mapping.writeToFile("1/map_blessed.nmf");
}

typedef void(*pTestFunction)(void);
const std::vector<pTestFunction> func{
	LangleyExample,
	TwoBlocks
};

int main(int argc, char *argv[])
{
	int cnt = 0;
	int failure = 0;
	std::cout << "Testing cases for the Neutral Map File(NMF) utilities." << std::endl;
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
