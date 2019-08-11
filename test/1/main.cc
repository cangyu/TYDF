#include "topology.h"
#include "xf_msh.h"


using namespace std;

const string MESH_PATH("./fluent.msh");
const string OUTPUT_PATH("./blessed.msh");

XF_MSH msh;

int main(int argc, char *argv[])
{
	cout << "Reading mesh: \"" << MESH_PATH << "\" ..." << endl;

	int ret = msh.readFromFile(MESH_PATH);
	
	if (!ret)
		cout << "\"" << MESH_PATH << "\" read successfully!" << endl;
	else
		cout << "Failure: " << ret << endl;

	ret = msh.writeToFile(OUTPUT_PATH);

    return 0;
}
