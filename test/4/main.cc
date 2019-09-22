#include "nmf.h"

int test0()
{
    std::cout << "test0" << std::endl;

    NMF mapping;
    mapping.readFromFile("mapping0.nmf");
    mapping.writeToFile("mapping0_blessed.nmf");
    return 0;
}

int test1()
{
    std::cout << "test1" << std::endl;

    NMF mapping;
    mapping.readFromFile("mapping1.nmf");
    mapping.writeToFile("mapping1_blessed.nmf");
    return 0;
}

int main(int argc, char *argv[])
{
    int ret = 0;

    ret = test0();
    if(ret)
        throw std::runtime_error("Error: " + std::to_string(ret));

    ret = test1();
    if(ret)
        throw std::runtime_error("Error: " + std::to_string(ret));

	return ret;
}
