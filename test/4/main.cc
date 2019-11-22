#include <iostream>
#include "nmf.hpp"

void test0()
{
	NMF::Mapping3D mapping("0/map.nmf");
	mapping.summary();
	mapping.writeToFile("0/map_blessed.nmf");
}

void test1()
{
	NMF::Mapping3D mapping("1/map.nmf");
	mapping.summary();
	mapping.writeToFile("1/map_blessed.nmf");
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
