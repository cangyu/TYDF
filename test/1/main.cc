#include "xf_msh.hpp"

using namespace std;

int test0()
{
	const string MESH_PATH("fluent.msh");
	const string OUTPUT_PATH("blessed.msh");

	int ret = 0;
	XF_MSH msh;

	cout << "Reading mesh: \"" << MESH_PATH << "\" ..." << endl;

	ret = msh.readFromFile(MESH_PATH);
	if (!ret)
		cout << "\"" << MESH_PATH << "\" read successfully!" << endl;
	else
		cout << "Failure: " << ret << endl;

	ret = msh.writeToFile(OUTPUT_PATH);

	return ret;
}

int main(int argc, char *argv[])
{
	test0();

	return 0;
}
