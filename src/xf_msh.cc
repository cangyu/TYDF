#include "xf_msh.hpp"

const std::map<int, std::string> XF_FACE::MAPPING_Idx2Str{
	std::pair<int, std::string>(XF_FACE::MIXED, "MIXED"),
	std::pair<int, std::string>(XF_FACE::LINEAR, "LINEAR"),
	std::pair<int, std::string>(XF_FACE::TRIANGULAR, "TRIANGULAR"),
	std::pair<int, std::string>(XF_FACE::QUADRILATERAL, "QUADRILATERAL"),
	std::pair<int, std::string>(XF_FACE::POLYGONAL, "POLYGONAL")
};

const std::map<std::string, int> XF_FACE::MAPPING_Str2Idx{
	// MIXED
	std::pair<std::string, int>("MIXED", XF_FACE::MIXED),
	std::pair<std::string, int>("Mixed", XF_FACE::MIXED),
	std::pair<std::string, int>("mixed", XF_FACE::MIXED),
	// LINEAR
	std::pair<std::string, int>("LINEAR", XF_FACE::LINEAR),
	std::pair<std::string, int>("Linear", XF_FACE::LINEAR),
	std::pair<std::string, int>("linear", XF_FACE::LINEAR),
	std::pair<std::string, int>("Line", XF_FACE::LINEAR),
	std::pair<std::string, int>("line", XF_FACE::LINEAR),
	// TRIANGULAR
	std::pair<std::string, int>("TRIANGULAR", XF_FACE::TRIANGULAR),
	std::pair<std::string, int>("Triangular", XF_FACE::TRIANGULAR),
	std::pair<std::string, int>("triangular", XF_FACE::TRIANGULAR),
	std::pair<std::string, int>("TRI", XF_FACE::TRIANGULAR),
	std::pair<std::string, int>("Tri", XF_FACE::TRIANGULAR),
	std::pair<std::string, int>("tri", XF_FACE::TRIANGULAR),
	// QUADRILATERAL
	std::pair<std::string, int>("QUADRILATERAL", XF_FACE::QUADRILATERAL),
	std::pair<std::string, int>("Quadrilateral", XF_FACE::QUADRILATERAL),
	std::pair<std::string, int>("quadrilateral", XF_FACE::QUADRILATERAL),
	std::pair<std::string, int>("QUAD", XF_FACE::QUADRILATERAL),
	std::pair<std::string, int>("Quad", XF_FACE::QUADRILATERAL),
	std::pair<std::string, int>("quad", XF_FACE::QUADRILATERAL),
	// POLYGONAL
	std::pair<std::string, int>("POLYGONAL", XF_FACE::POLYGONAL),
	std::pair<std::string, int>("Polygonal", XF_FACE::POLYGONAL),
	std::pair<std::string, int>("polygonal", XF_FACE::POLYGONAL)
};
