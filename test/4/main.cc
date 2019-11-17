#include <iostream>
#include "nmf.hpp"

int test0()
{
	NMF::Mapping3D mapping;
	mapping.readFromFile("mapping0.nmf");
	mapping.writeToFile("mapping0_blessed.nmf");

	mapping.compute_topology();

	return 0;
}

int test1()
{
	NMF::Mapping3D mapping;
	mapping.readFromFile("mapping1.nmf");
	mapping.writeToFile("mapping1_blessed.nmf");
	return 0;
}

int main(int argc, char *argv[])
{
	int ret = 0;

    std::cout << "test0 ... " << std::endl;
    ret = test0();
	if (ret)
		throw std::runtime_error("Error: " + std::to_string(ret));

    std::cout << "test1 ... " << std::endl;
    ret = test1();
	if (ret)
		throw std::runtime_error("Error: " + std::to_string(ret));

	return ret;
}
