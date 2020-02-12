#include <iostream>
#include "../../inc/xf.h"

using namespace GridTool;

static const std::string CASTE_SEP = "  ";

void test(const std::string &case_name, const std::string &case_desc, const std::string &file_dir, const std::string &file_name)
{
    const std::string REPORT_PATH = file_dir + file_name + "_report.txt";
    const std::string MESH_PATH = file_dir + file_name + ".msh";
    const std::string TRANSCRIPT_PATH = file_dir + file_name + "_blessed.msh";

    std::cout << "Case \"" << case_name << "\"," << case_desc << " ..." << std::endl;
    std::ofstream fout(REPORT_PATH);
    if (fout.fail())
        throw std::runtime_error("Failed to open report file.");

    std::cout << CASTE_SEP << "Reading ..." << std::endl;
    XF::MESH msh(MESH_PATH, fout);
    fout.close();

    std::cout << CASTE_SEP << "Transcribing ..." << std::endl;
    msh.writeToFile(TRANSCRIPT_PATH);

    std::cout << CASTE_SEP << "Done!" << std::endl;
}

int main(int argc, char *argv[])
{
    std::cout << "Testing I/O of \"FLUENT\" unstructured mesh ..." << std::endl;

    test("Cavity1", "a 32 x 32 x 32 cube", "../../case/Cavity/FLUENT/", "grid32");
    test("Cavity2", "a 64 x 64 x 64 cube", "../../case/Cavity/FLUENT/", "grid64");
    test("Cavity3", "a 128 x 128 x 128 cube", "../../case/Cavity/FLUENT/", "grid128");
    test("Cavity4", "a 256 x 256 x 256 cube", "../../case/Cavity/FLUENT/", "grid256");

    test("Structure1", "an example from LiuSha's tutorial", "../../case/LS1/FLUENT/", "fluent");

    return 0;
}
