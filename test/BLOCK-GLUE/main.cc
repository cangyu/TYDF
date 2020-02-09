#include <iostream>
#include <fstream>
#include <vector>
#include "../../inc/xf.h"

const std::vector<std::string> test
{

};

int main(int argc, char *argv[])
{
	int cnt = 0;
	int failure = 0;
	std::cout << "Testing cases for the Grid-Glue utilities." << std::endl;
	for (const auto &dir_name : test)
	{
		std::cout << "Case " << ++cnt << " ... " << std::endl;
		try
		{

		}
		catch (std::exception &e)
		{
			std::cout << e.what() << std::endl;
			++failure;
		}
	}
	std::cout << "Done with " << failure << " failed." << std::endl;
	return 0;
}
