#include <iostream>
#include "../../inc/plot3d.h"

using namespace GridTool;

static const std::string CASTE_SEP = "  ";

void test(const std::string &case_name, const std::string &case_desc, const std::string &file_dir, const std::string &file_name)
{
    const std::string GRID_PATH = file_dir + file_name + ".fmt";
    const std::string TRANSCRIPT_PATH = file_dir + file_name + "_blessed.fmt";

    std::cout << "Case \"" << case_name << "\"," << case_desc << " ..." << std::endl;

    std::cout << CASTE_SEP << "Reading ..." << std::endl;
    PLOT3D::GRID p3d(GRID_PATH);

    std::cout << CASTE_SEP << "Transcribing ..." << std::endl;
    p3d.writeToFile(TRANSCRIPT_PATH);

    std::cout << CASTE_SEP << "Done!" << std::endl;
}

int main(int argc, char *argv[])
{
    std::cout << "Testing I/O of \"PLOT3D\" grid ..." << std::endl;

    test("Planar1", "a 2D grid in 2D form", "../../case/PLOT3D/", "xyz");
    test("Shell1", "a 2D grid in 3D form", "../../case/PLOT3D/", "xyz");
    test("Cube1", "a 3D single-block grid", "../../case/PLOT3D/", "xyz");

    return 0;
}
