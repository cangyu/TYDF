#include <iostream>
#include <vector>
#include "nmf.hpp"

void LangleyExample()
{
	// 4 blocks in 2x2 form.
	NMF::Mapping3D mapping("0/map.nmf");
	mapping.summary();
	mapping.writeToFile("0/map_blessed.nmf");
}

void TwoBlocks()
{
	NMF::Mapping3D mapping("1/map.nmf");
	mapping.summary();
	mapping.writeToFile("1/map_blessed.nmf");
}

typedef void(*pTestFunction)(void);
const std::vector<pTestFunction> func = {
	LangleyExample,
	TwoBlocks
};

int main(int argc, char *argv[])
{
	int cnt = 0;
	for (auto f : func)
	{
		std::cout << "\nCase " << ++cnt << " ..." << std::endl;
		try {
			f();
		}
		catch (std::exception &e)
		{
			std::cout << e.what() << std::endl;
		}
		catch (...)
		{
			throw std::runtime_error("Unexpected.");
		}
	}
	return 0;
}
