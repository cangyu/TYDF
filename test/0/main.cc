#include <iostream>
#include <cstddef>
#include "xf_msh.h"

using namespace std;

const string INPUT_MESH_PATH("./fluent.msh");
const string OUTPUT_MESH_PATH("./blessed.msh");

XF_MSH msh;

int main(int argc, char *argv[])
{
	int ret = 0;

	cout << "Reading mesh: \"" << INPUT_MESH_PATH << "\" ..." << endl;
	ret = msh.readFromFile(INPUT_MESH_PATH);
	if (!ret)
		cout << "\"" << INPUT_MESH_PATH << "\" reading successfully!" << endl;
	else
		cout << "Failure: " << ret << endl;

	cout << "Writing mesh: \"" << OUTPUT_MESH_PATH << "\" ..." << endl;
	ret = msh.writeToFile(OUTPUT_MESH_PATH);
	if (!ret)
		cout << "\"" << OUTPUT_MESH_PATH << "\" writing successfully!" << endl;
	else
		cout << "Failure: " << ret << endl;

	return ret;
}
