#include "xf_msh.h"

const std::map<int, std::string> XF_BC::MAPPING_Idx2Str{
	std::pair<int, std::string>(XF_BC::INTERIOR, "INTERIOR"),
	std::pair<int, std::string>(XF_BC::WALL, "WALL"),
	std::pair<int, std::string>(XF_BC::PRESSURE_INLET, "PRESSURE_INLET"),
	std::pair<int, std::string>(XF_BC::PRESSURE_OUTLET, "PRESSURE_OUTLET"),
	std::pair<int, std::string>(XF_BC::SYMMETRY, "SYMMETRY"),
	std::pair<int, std::string>(XF_BC::PERIODIC_SHADOW, "PERIODIC_SHADOW"),
	std::pair<int, std::string>(XF_BC::PRESSURE_FAR_FIELD, "PRESSURE_FAR_FIELD"),
	std::pair<int, std::string>(XF_BC::VELOCITY_INLET, "VELOCITY_INLET"),
	std::pair<int, std::string>(XF_BC::PERIODIC, "PERIODIC"),
	std::pair<int, std::string>(XF_BC::FAN, "FAN"),
	std::pair<int, std::string>(XF_BC::MASS_FLOW_INLET, "MASS_FLOW_INLET"),
	std::pair<int, std::string>(XF_BC::INTERFACE, "INTERFACE"),
	std::pair<int, std::string>(XF_BC::PARENT, "PARENT"),
	std::pair<int, std::string>(XF_BC::OUTFLOW, "OUTFLOW"),
	std::pair<int, std::string>(XF_BC::AXIS, "AXIS")
};

const std::map<std::string, int> XF_BC::MAPPING_Str2Idx{
	// INTERIOR
	std::pair<std::string, int>("INTERIOR", XF_BC::INTERIOR),
	std::pair<std::string, int>("Interior", XF_BC::INTERIOR),
	std::pair<std::string, int>("interior", XF_BC::INTERIOR),
	// WALL
	std::pair<std::string, int>("WALL", XF_BC::WALL),
	std::pair<std::string, int>("Wall", XF_BC::WALL),
	std::pair<std::string, int>("wall", XF_BC::WALL),
	// PRESSURE_INLET
	std::pair<std::string, int>("PRESSURE_INLET", XF_BC::PRESSURE_INLET),
	std::pair<std::string, int>("Pressure_Inlet", XF_BC::PRESSURE_INLET),
	std::pair<std::string, int>("pressure_inlet", XF_BC::PRESSURE_INLET),
	// PRESSURE_OUTLET
	std::pair<std::string, int>("PRESSURE_OUTLET", XF_BC::PRESSURE_OUTLET),
	std::pair<std::string, int>("Pressure_Outlet", XF_BC::PRESSURE_OUTLET),
	std::pair<std::string, int>("pressure_outlet", XF_BC::PRESSURE_OUTLET),
	// SYMMETRY
	std::pair<std::string, int>("SYMMETRY", XF_BC::SYMMETRY),
	std::pair<std::string, int>("Symmetry", XF_BC::SYMMETRY),
	std::pair<std::string, int>("symmetry", XF_BC::SYMMETRY),
	// PERIODIC_SHADOW
	std::pair<std::string, int>("PERIODIC_SHADOW", XF_BC::PERIODIC_SHADOW),
	std::pair<std::string, int>("Periodic_Shadow", XF_BC::PERIODIC_SHADOW),
	std::pair<std::string, int>("periodic_shadow", XF_BC::PERIODIC_SHADOW),
	// PRESSURE_FAR_FIELD
	std::pair<std::string, int>("PRESSURE_FAR_FIELD", XF_BC::PRESSURE_FAR_FIELD),
	std::pair<std::string, int>("Pressure_Far_Field", XF_BC::PRESSURE_FAR_FIELD),
	std::pair<std::string, int>("pressure_far_field", XF_BC::PRESSURE_FAR_FIELD),
	// VELOCITY_INLET
	std::pair<std::string, int>("VELOCITY_INLET", XF_BC::VELOCITY_INLET),
	std::pair<std::string, int>("Velocity_Inlet", XF_BC::VELOCITY_INLET),
	std::pair<std::string, int>("velocity_inlet", XF_BC::VELOCITY_INLET),
	// PERIODIC
	std::pair<std::string, int>("PERIODIC", XF_BC::PERIODIC),
	std::pair<std::string, int>("Periodic", XF_BC::PERIODIC),
	std::pair<std::string, int>("periodic", XF_BC::PERIODIC),
	// FAN
	std::pair<std::string, int>("FAN", XF_BC::FAN),
	std::pair<std::string, int>("Fan", XF_BC::FAN),
	std::pair<std::string, int>("fan", XF_BC::FAN),
	// MASS_FLOW_INLET
	std::pair<std::string, int>("MASS_FLOW_INLET", XF_BC::MASS_FLOW_INLET),
	std::pair<std::string, int>("Mass_Flow_Inlet", XF_BC::MASS_FLOW_INLET),
	std::pair<std::string, int>("mass_flow_inlet", XF_BC::MASS_FLOW_INLET),
	// INTERFACE
	std::pair<std::string, int>("INTERFACE", XF_BC::INTERFACE),
	std::pair<std::string, int>("Interface", XF_BC::INTERFACE),
	std::pair<std::string, int>("interface", XF_BC::INTERFACE),
	// PARENT
	std::pair<std::string, int>("PARENT", XF_BC::PARENT),
	std::pair<std::string, int>("Parent", XF_BC::PARENT),
	std::pair<std::string, int>("parent", XF_BC::PARENT),
	// OUTFLOW
	std::pair<std::string, int>("OUTFLOW", XF_BC::OUTFLOW),
	std::pair<std::string, int>("Outflow", XF_BC::OUTFLOW),
	std::pair<std::string, int>("outflow", XF_BC::OUTFLOW),
	// AXIS
	std::pair<std::string, int>("AXIS", XF_BC::AXIS),
	std::pair<std::string, int>("Axis", XF_BC::AXIS),
	std::pair<std::string, int>("axis", XF_BC::AXIS)
};

const std::map<int, std::string> XF_CELL::TYPE_MAPPING_Idx2Str{
	std::pair<int, std::string>(XF_CELL::FLUID, "FLUID"),
	std::pair<int, std::string>(XF_CELL::SOLID, "SOLID"),
	std::pair<int, std::string>(XF_CELL::DEAD, "DEAD")
};

const std::map<std::string, int> XF_CELL::TYPE_MAPPING_Str2Idx{
	// FLUID
	std::pair<std::string, int>("FLUID", XF_CELL::FLUID),
	std::pair<std::string, int>("Fluid", XF_CELL::FLUID),
	std::pair<std::string, int>("fluid", XF_CELL::FLUID),
	// SOLID
	std::pair<std::string, int>("SOLID", XF_CELL::SOLID),
	std::pair<std::string, int>("Solid", XF_CELL::SOLID),
	std::pair<std::string, int>("solid", XF_CELL::SOLID),
	// DEAD
	std::pair<std::string, int>("DEAD", XF_CELL::DEAD),
	std::pair<std::string, int>("Dead", XF_CELL::DEAD),
	std::pair<std::string, int>("dead", XF_CELL::DEAD)
};

const std::map<int, std::string> XF_CELL::ELEM_MAPPING_Idx2Str{
	std::pair<int, std::string>(XF_CELL::MIXED, "MIXED"),
	std::pair<int, std::string>(XF_CELL::TRIANGULAR, "TRIANGULAR"),
	std::pair<int, std::string>(XF_CELL::TETRAHEDRAL, "TETRAHEDRAL"),
	std::pair<int, std::string>(XF_CELL::QUADRILATERAL, "QUADRILATERAL"),
	std::pair<int, std::string>(XF_CELL::HEXAHEDRAL, "HEXAHEDRAL"),
	std::pair<int, std::string>(XF_CELL::PYRAMID, "PYRAMID"),
	std::pair<int, std::string>(XF_CELL::WEDGE, "WEDGE"),
	std::pair<int, std::string>(XF_CELL::POLYHEDRAL, "POLYHEDRAL")
};

const std::map<std::string, int> XF_CELL::ELEM_MAPPING_Str2Idx{
	// MIXED
	std::pair<std::string, int>("MIXED", XF_CELL::MIXED),
	std::pair<std::string, int>("Mixed", XF_CELL::MIXED),
	std::pair<std::string, int>("mixed", XF_CELL::MIXED),
	// TRIANGULAR
	std::pair<std::string, int>("TRIANGULAR", XF_CELL::TRIANGULAR),
	std::pair<std::string, int>("Triangular", XF_CELL::TRIANGULAR),
	std::pair<std::string, int>("triangular", XF_CELL::TRIANGULAR),
	std::pair<std::string, int>("TRI", XF_CELL::TRIANGULAR),
	std::pair<std::string, int>("Tri", XF_CELL::TRIANGULAR),
	std::pair<std::string, int>("tri", XF_CELL::TRIANGULAR),
	// TETRAHEDRAL
	std::pair<std::string, int>("TETRAHEDRAL", XF_CELL::TETRAHEDRAL),
	std::pair<std::string, int>("Tetrahedral", XF_CELL::TETRAHEDRAL),
	std::pair<std::string, int>("tetrahedral", XF_CELL::TETRAHEDRAL),
	std::pair<std::string, int>("TET", XF_CELL::TETRAHEDRAL),
	std::pair<std::string, int>("Tet", XF_CELL::TETRAHEDRAL),
	std::pair<std::string, int>("tet", XF_CELL::TETRAHEDRAL),
	// QUADRILATERAL
	std::pair<std::string, int>("QUADRILATERAL", XF_CELL::QUADRILATERAL),
	std::pair<std::string, int>("Quadrilateral", XF_CELL::QUADRILATERAL),
	std::pair<std::string, int>("quadrilateral", XF_CELL::QUADRILATERAL),
	std::pair<std::string, int>("QUAD", XF_CELL::QUADRILATERAL),
	std::pair<std::string, int>("Quad", XF_CELL::QUADRILATERAL),
	std::pair<std::string, int>("quad", XF_CELL::QUADRILATERAL),
	// HEXAHEDRAL
	std::pair<std::string, int>("HEXAHEDRAL", XF_CELL::HEXAHEDRAL),
	std::pair<std::string, int>("Hexahedral", XF_CELL::HEXAHEDRAL),
	std::pair<std::string, int>("hexahedral", XF_CELL::HEXAHEDRAL),
	std::pair<std::string, int>("HEX", XF_CELL::HEXAHEDRAL),
	std::pair<std::string, int>("Hex", XF_CELL::HEXAHEDRAL),
	std::pair<std::string, int>("hex", XF_CELL::HEXAHEDRAL),
	// PYRAMID
	std::pair<std::string, int>("PYRAMID", XF_CELL::PYRAMID),
	std::pair<std::string, int>("Pyramid", XF_CELL::PYRAMID),
	std::pair<std::string, int>("pyramid", XF_CELL::PYRAMID),
	// WEDGE
	std::pair<std::string, int>("WEDGE", XF_CELL::WEDGE),
	std::pair<std::string, int>("Wedge", XF_CELL::WEDGE),
	std::pair<std::string, int>("wedge", XF_CELL::WEDGE),
	// POLYHEDRAL
	std::pair<std::string, int>("POLYHEDRAL", XF_CELL::POLYHEDRAL),
	std::pair<std::string, int>("Polyhedral", XF_CELL::POLYHEDRAL),
	std::pair<std::string, int>("polyhedral", XF_CELL::POLYHEDRAL)
};

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
