#include <iostream>
#include "../../inc/xf.h"

using namespace GridTool;

static const std::string CASTE_SEP = "  ";

void test(const std::string &case_name, const std::string &case_desc, const std::string &MAP_PATH, const std::string &GRID_PATH, const std::string &MESH_DIR, const std::string &MESH_NAME)
{
    const std::string MESH_PATH = MESH_DIR + MESH_NAME + ".msh";
    const std::string REPORT_PATH = MESH_DIR + MESH_NAME + "_report.txt";

    try
    {
        std::cout << "Case \"" << case_name << "\"," << case_desc << " ..." << std::endl;
        std::ofstream frpt(REPORT_PATH);
        if (frpt.fail())
            throw std::runtime_error("Failed to open target report file.");

        std::cout << CASTE_SEP << "Combining ..." << std::endl;
        const XF::MESH mesh(MAP_PATH, GRID_PATH, frpt);
        frpt.close();

        std::cout << CASTE_SEP << "Writing ..." << std::endl;
        mesh.writeToFile(MESH_PATH);

        std::cout << CASTE_SEP << "Done!" << std::endl;
    }
    catch (std::exception &e)
    {
        std::cout << e.what() << std::endl;
    }
};

int main(int argc, char *argv[])
{
    std::cout << "Test the \"Block-Glue\" utilities." << std::endl;

    return 0;
}
