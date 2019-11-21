#include <iostream>
#include "nmf.hpp"

void test0()
{
	NMF::Mapping3D mapping;
	mapping.readFromFile("mapping0.nmf");
	mapping.compute_topology();
	mapping.summary();
	mapping.writeToFile("mapping0_blessed.nmf");
}

void test1()
{
	NMF::Mapping3D mapping("mapping1.nmf");
	mapping.summary();
	mapping.writeToFile("mapping1_blessed.nmf");
}

int main(int argc, char *argv[])
{
	std::cout << "test0 ... " << std::endl;
	try {
		test0();
	}
	catch (std::exception &e)
	{
		std::cout << e.what() << std::endl;
	}
	catch (...)
	{
		throw std::runtime_error("Unexpected.");
	}

	std::cout << "test1 ... " << std::endl;
	try {
		test1();
	}
	catch (std::exception &e)
	{
		std::cout << e.what() << std::endl;
	}
	catch (...)
	{
		throw std::runtime_error("Unexpected.");
	}

	return 0;
}
