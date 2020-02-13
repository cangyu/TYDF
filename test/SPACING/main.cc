#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include "../../inc/spacing.h"

using namespace GridTool;

void writeToFile(const DIST_ARR &src, const std::string &dst)
{
    std::ofstream fout(dst);
    if (fout.fail())
        throw std::runtime_error("Failed to open target output file.");

    for (const auto &e : src)
        fout << std::to_string(e) << std::endl;

    fout.close();
}

int main(int argc, char *argv[])
{
    std::cout << "Test the \"SPACING\" utilities." << std::endl;

    DIST_ARR c1;
    SPACING::uniform(100, c1);
    writeToFile(c1, "uniform100.txt");

    return 0;
}
