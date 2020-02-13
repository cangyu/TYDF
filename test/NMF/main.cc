#include <iostream>
#include "../../inc/nmf.h"

using namespace GridTool;

static const std::string CASTE_SEP = "  ";

void test(const std::string &case_name, const std::string &case_desc, const std::string &file_dir, const std::string &file_name)
{
    const std::string REPORT_PATH = file_dir + file_name + "_report.txt";
    const std::string MAP_PATH = file_dir + file_name + ".nmf";
    const std::string TRANSCRIPT_PATH = file_dir + file_name + "_blessed.nmf";

    std::cout << "Case \"" << case_name << "\"," << case_desc << " ..." << std::endl;

    std::cout << CASTE_SEP << "Reading ..." << std::endl;
    NMF::Mapping3D mapping(MAP_PATH);

    std::cout << CASTE_SEP << "Diagnosing ..." << std::endl;
    std::ofstream fout(REPORT_PATH);
    mapping.summary(fout);
    fout.close();

    std::cout << CASTE_SEP << "Numbering ..." << std::endl;
    mapping.numbering();

    std::cout << CASTE_SEP << "Transcribing ..." << std::endl;
    mapping.writeToFile(TRANSCRIPT_PATH);

    std::cout << CASTE_SEP << "Done!" << std::endl;
}

int main(int argc, char *argv[])
{
    std::cout << "Test the \"Neutral Map File\" utilities." << std::endl;

    test("Cavity", "a single block", "../../case/Cavity/NMF/", "map");
    test("Sky1", "2 blocks connected through 1 surface", "../../case/Sky1/NMF/", "map");
    test("Langley", "4 blocks in 2 x 2 form", "../../case/Langley/NMF/", "map");

    return 0;
}
