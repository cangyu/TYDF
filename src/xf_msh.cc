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

static inline void eat(std::istream &in, char c)
{
	char tmp;

	do {
		in >> tmp;
	} while (tmp != c);
}

static inline void skip_white(std::istream &in)
{
	char tmp;

	do {
		in >> tmp;
	} while (tmp == ' ' || tmp == '\t' || tmp == '\n');

	if (!in.eof())
		in.unget();
}

typedef std::vector<double> Array1D;

static inline double dot_product(const Array1D &na, const Array1D &nb)
{
	const size_t ND = na.size();
	double ret = 0.0;
	for (size_t i = 0; i < ND; ++i)
		ret += na[i] * nb[i];
	return ret;
}

static inline void cross_product(const Array1D &a, const Array1D &b, Array1D &dst)
{
	dst[0] = a[1] * b[2] - a[2] * b[1];
	dst[1] = a[2] * b[0] - a[0] * b[2];
	dst[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void cross_product(double *a, double *b, double *dst)
{
	dst[0] = a[1] * b[2] - a[2] * b[1];
	dst[1] = a[2] * b[0] - a[0] * b[2];
	dst[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void delta(const Array1D &na, const Array1D &nb, Array1D &dst)
{
	const size_t ND = dst.size();
	for (size_t i = 0; i < ND; ++i)
		dst[i] = nb[i] - na[i];
}

static inline void delta(double *na, double *nb, double *dst)
{
	for (size_t i = 0; i < 3; ++i)
		dst[i] = nb[i] - na[i];
}

static inline void normalize(const Array1D &src, Array1D &dst)
{
	const size_t ND = dst.size();
	double L = 0.0;
	for (size_t i = 0; i < ND; ++i)
		L += src[i] * src[i];
	L = std::sqrt(L);
	for (size_t i = 0; i < ND; ++i)
		dst[i] = src[i] / L;
}

static inline void normalize(double *src, double *dst)
{
	double L = 0.0;
	for (size_t i = 0; i < 3; ++i)
		L += src[i] * src[i];
	L = std::sqrt(L);
	for (size_t i = 0; i < 3; ++i)
		dst[i] = src[i] / L;
}

static inline double distance(const Array1D &na, const Array1D &nb)
{
	const size_t ND = na.size();
	double L = 0.0;
	for (size_t i = 0; i < ND; ++i)
	{
		double di = nb[i] - na[i];
		L += di * di;
	}
	return std::sqrt(L);
}

static inline double distance(double *na, double *nb)
{
	double L = 0.0;
	for (size_t i = 0; i < 3; ++i)
	{
		double di = nb[i] - na[i];
		L += di * di;
	}
	return std::sqrt(L);
}

static inline void line_center(const Array1D &na, const Array1D &nb, Array1D &dst)
{
	const size_t ND = dst.size();
	for (size_t i = 0; i < ND; ++i)
		dst[i] = 0.5*(na[i] + nb[i]);
}

static inline void line_center(double *na, double *nb, double *dst)
{
	for (size_t i = 0; i < 3; ++i)
		dst[i] = 0.5*(na[i] + nb[i]);
}

static inline void line_normal(double *na, double *nb, double *dst, double *dst_r)
{
	delta(na, nb, dst);
	std::swap(dst[0], dst[1]);
	dst[1] = -dst[1];
	normalize(dst, dst);

	for (size_t i = 0; i < 3; ++i)
		dst_r[i] = -dst[i];
}

static inline void triangle_center(const Array1D &na, const Array1D &nb, const Array1D &nc, Array1D &dst)
{
	const size_t ND = dst.size();
	for (size_t i = 0; i < ND; ++i)
		dst[i] = (na[i] + nb[i] + nc[i]) / 3.0;
}

static inline void triangle_center(double *na, double *nb, double *nc, double *dst)
{
	for (size_t i = 0; i < 3; ++i)
		dst[i] = (na[i] + nb[i] + nc[i]) / 3.0;
}

static inline double triangle_area(const Array1D &na, const Array1D &nb, const Array1D &nc)
{
	double c = distance(na, nb);
	double a = distance(nb, nc);
	double b = distance(nc, na);
	double p = 0.5*(a + b + c);

	// Heron's formula
	return std::sqrt(p*(p - a)*(p - b)*(p - c));
}

static inline double triangle_area(double *na, double *nb, double *nc)
{
	double c = distance(na, nb);
	double a = distance(nb, nc);
	double b = distance(nc, na);
	double p = 0.5*(a + b + c);

	// Heron's formula
	return std::sqrt(p*(p - a)*(p - b)*(p - c));
}

static inline void triangle_normal(double *na, double *nb, double *nc, double *dst, double *dst_r)
{
	// Order of nodes follows the right-hand convention.
	double rab[3], rac[3];
	delta(na, nb, rab);
	delta(na, nc, rac);
	cross_product(rac, rab, dst); // Take cross product to find normal direction
	normalize(dst, dst); // Normalize

	for (size_t i = 0; i < 3; ++i)
		dst_r[i] = -dst[i];
}

static inline void quadrilateral_center(const Array1D &n1, const Array1D &n2, const Array1D &n3, const Array1D &n4, Array1D &dst)
{
	// 1, 2, 3, 4 are in anti-clockwise direction.
	const double S123 = triangle_area(n1, n2, n3);
	const double S134 = triangle_area(n1, n3, n4);

	auto rc123(dst), rc134(dst);
	triangle_center(n1, n2, n3, rc123);
	triangle_center(n1, n3, n4, rc134);

	const double alpha = S123 / (S123 + S134);
	const double beta = 1.0 - alpha;

	const size_t ND = dst.size();
	for (size_t i = 0; i < ND; ++i)
		dst[i] = alpha * rc123[i] + beta * rc134[i];
}

static inline void quadrilateral_center(double *n1, double *n2, double *n3, double *n4, double *dst)
{
	// 1, 2, 3, 4 are in anti-clockwise direction.
	const double S123 = triangle_area(n1, n2, n3);
	const double S134 = triangle_area(n1, n3, n4);

	double rc123[3], rc134[3];
	triangle_center(n1, n2, n3, rc123);
	triangle_center(n1, n3, n4, rc134);

	const double alpha = S123 / (S123 + S134);
	const double beta = 1.0 - alpha;

	for (size_t i = 0; i < 3; ++i)
		dst[i] = alpha * rc123[i] + beta * rc134[i];
}

static inline double quadrilateral_area(const Array1D &n1, const Array1D &n2, const Array1D &n3, const Array1D &n4)
{
	// 1, 2, 3, 4 are in anti-clockwise direction.
	const double S123 = triangle_area(n1, n2, n3);
	const double S134 = triangle_area(n1, n3, n4);
	return S123 + S134;
}

static inline double quadrilateral_area(double *n1, double *n2, double *n3, double *n4)
{
	// 1, 2, 3, 4 are in anti-clockwise direction.
	const double S123 = triangle_area(n1, n2, n3);
	const double S134 = triangle_area(n1, n3, n4);
	return S123 + S134;
}

static inline void quadrilateral_normal(double *n1, double *n2, double *n3, double *n4, double *dst, double *dst_r)
{
	// See (5.12) of Jiri Blazek's CFD book.
	const double dxa = n4[0] - n2[0], dxb = n3[0] - n1[0];
	const double dya = n4[1] - n2[1], dyb = n3[1] - n1[1];
	const double dza = n4[2] - n2[2], dzb = n3[2] - n1[2];

	// See (5.13) of Jiri Blazek's CFD book.
	dst[0] = 0.5*(dza * dyb - dya * dzb);
	dst[1] = 0.5*(dxa * dzb - dza * dxb);
	dst[2] = 0.5*(dya * dxb - dxa * dyb);

	// Normalize
	normalize(dst, dst);

	for (size_t i = 0; i < 3; ++i)
		dst_r[i] = -dst[i];
}

int XF_MSH::readFromFile(const std::string &src)
{
	// Open grid file
	std::ifstream fin(src);
	if (fin.fail())
		throw std::runtime_error("Failed to open input grid file: " + src);

	// Clear existing records if any.
	clear_entry();

	// Read contents
	while (!fin.eof())
	{
		skip_white(fin);
		eat(fin, '(');
		int ti;
		fin >> std::dec >> ti;
		if (ti == XF_SECTION::COMMENT)
		{
			eat(fin, '\"');
			std::string ts;
			char tc;
			while ((tc = fin.get()) != '\"')
				ts.push_back(tc);
			eat(fin, ')');
			add_entry(new XF_COMMENT(ts));
			skip_white(fin);
		}
		else if (ti == XF_SECTION::HEADER)
		{
			eat(fin, '\"');
			std::string ts;
			char tc;
			while ((tc = fin.get()) != '\"')
				ts.push_back(tc);
			eat(fin, ')');
			add_entry(new XF_HEADER(ts));
			skip_white(fin);
		}
		else if (ti == XF_SECTION::DIMENSION)
		{
			int nd = 0;
			fin >> std::dec >> nd;
			eat(fin, ')');
			add_entry(new XF_DIMENSION(nd));
			skip_white(fin);
			m_dim = nd;
			m_is3D = (nd == 3);
		}
		else if (ti == XF_SECTION::NODE)
		{
			eat(fin, '(');
			int zone;
			fin >> std::hex >> zone;
			if (zone == 0)
			{
				// If zone-id is 0, indicating total number of nodes in the mesh.
				int tmp;
				fin >> std::hex;
				fin >> tmp;
				if (tmp != 1)
					throw std::runtime_error("Invalid \"first-index\" in NODE declaration!");
				fin >> m_totalNodeNum;
				std::cout << "Total number of nodes: " << m_totalNodeNum << std::endl;
				fin >> tmp;
				if (tmp != 0)
					throw std::runtime_error("Invalid \"type\" in NODE declaration!");
				char ndc = fin.get();
				if (ndc != ')')
				{
					fin >> tmp;
					eat(fin, ')');
				}
				eat(fin, ')');
			}
			else
			{
				// If zone-id is positive, it indicates the zone to which the nodes belong.
				int first, last, tp, nd;
				fin >> std::hex;
				fin >> first >> last;
				fin >> tp >> nd;
				auto e = new XF_NODE(zone, first, last, tp, nd);
				eat(fin, ')');
				eat(fin, '(');
				std::cout << "Reading " << e->num() << " nodes in zone " << zone << " (from " << first << " to " << last << ") ..." << std::endl;

				if (nd != dimension())
					throw std::runtime_error("Inconsistent with previous DIMENSION declaration!");

				if (nd == 3)
				{
					double x, y, z;
					for (int i = first; i <= last; ++i)
					{
						size_t i_loc = i - first;
						fin >> x >> y >> z;
						e->set_coordinate(i_loc, x, y, z);
					}
				}
				else
				{
					double x, y;
					for (int i = first; i <= last; ++i)
					{
						size_t i_loc = i - first;
						fin >> x >> y;
						e->set_coordinate(i_loc, x, y);
					}
				}
				eat(fin, ')');
				eat(fin, ')');
				std::cout << "Done!" << std::endl;
				add_entry(e);
			}
			skip_white(fin);
		}
		else if (ti == XF_SECTION::CELL)
		{
			eat(fin, '(');
			int zone;
			fin >> std::hex >> zone;
			if (zone == 0)
			{
				// If zone-id is 0, indicating total number of cells in the mesh.
				int tmp;
				fin >> std::hex;
				fin >> tmp;
				if (tmp != 1)
					throw std::runtime_error("Invalid \"first-index\" in CELL declaration!");
				fin >> m_totalCellNum;
				std::cout << "Total number of cells: " << m_totalCellNum << std::endl;
				fin >> tmp;
				if (tmp != 0)
					throw std::runtime_error("Invalid \"type\" in CELL declaration!");
				char ndc = fin.get();
				if (ndc != ')')
				{
					fin >> tmp;
					eat(fin, ')');
				}
				eat(fin, ')');
			}
			else
			{
				// If zone-id is positive, it indicates the zone to which the cells belong.
				int first, last, tp, elem;
				fin >> std::hex;
				fin >> first >> last;
				fin >> tp >> elem;
				auto e = new XF_CELL(zone, first, last, tp, elem);
				eat(fin, ')');

				if (elem == 0)
				{
					std::cout << "Reading " << e->num() << " mixed cells in zone " << zone << " (from " << first << " to " << last << ") ..." << std::endl;
					eat(fin, '(');
					for (int i = first; i <= last; ++i)
					{
						fin >> elem;
						size_t i_loc = i - first;
						switch (elem)
						{
						case 1:
							e->elem(i_loc) = XF_CELL::TRIANGULAR;
							break;
						case 2:
							e->elem(i_loc) = XF_CELL::TETRAHEDRAL;
							break;
						case 3:
							e->elem(i_loc) = XF_CELL::QUADRILATERAL;
							break;
						case 4:
							e->elem(i_loc) = XF_CELL::HEXAHEDRAL;
							break;
						case 5:
							e->elem(i_loc) = XF_CELL::PYRAMID;
							break;
						case 6:
							e->elem(i_loc) = XF_CELL::WEDGE;
							break;
						case 7:
							e->elem(i_loc) = XF_CELL::POLYHEDRAL;
							break;
						default:
							throw std::runtime_error("Invalid cell type!");
						}
					}
					eat(fin, ')');
					std::cout << "Done!" << std::endl;
				}
				else
					std::cout << e->num() << " " << XF_CELL::ELEM_MAPPING_Idx2Str.at(elem) << " in zone " << zone << " (from " << first << " to " << last << ")" << std::endl;

				eat(fin, ')');
				add_entry(e);
			}
			skip_white(fin);
		}
		else if (ti == XF_SECTION::FACE)
		{
			eat(fin, '(');
			int zone;
			fin >> std::hex >> zone;
			if (zone == 0)
			{
				// If zone-id is 0, indicating total number of faces in the mesh.
				int tmp;
				fin >> tmp;
				if (tmp != 1)
					throw std::runtime_error("Invalid \"first-index\" in FACE declaration!");
				fin >> m_totalFaceNum;
				std::cout << "Total number of faces: " << m_totalFaceNum << std::endl;
				fin >> tmp;
				char ndc = fin.get();
				if (ndc != ')')
				{
					fin >> tmp;
					eat(fin, ')');
				}
				eat(fin, ')');
			}
			else
			{
				// If zone-id is positive, it indicates a regular face section and will be
				// followed by a body containing information about the grid connectivity.
				size_t first, last;
				int bc, face;
				fin >> first >> last;
				fin >> bc >> face;
				auto e = new XF_FACE(zone, first, last, bc, face);
				eat(fin, ')');
				eat(fin, '(');

				std::cout << "Reading " << e->num() << " faces in zone " << zone << " (from " << first << " to " << last << ") ..." << std::endl;

				size_t tmp_n[4];
				size_t tmp_c[2];
				for (size_t i = first; i <= last; ++i)
				{
					// Local index
					size_t i_loc = i - first;

					// Read connectivity record
					for (int j = 0; j < face; ++j)
						fin >> tmp_n[j];
					fin >> tmp_c[0] >> tmp_c[1];

					// Store current connectivity info
					e->connectivity(i_loc).set(face, tmp_n, tmp_c);
				}
				eat(fin, ')');
				eat(fin, ')');
				std::cout << "Done!" << std::endl;
				add_entry(e);
			}
			skip_white(fin);
		}
		else if (ti == XF_SECTION::ZONE)
		{
			eat(fin, '(');
			int zone;
			fin >> std::dec >> zone;
			std::string ztp;
			fin >> ztp;
			skip_white(fin);
			std::string zname;
			char t0;
			while ((t0 = fin.get()) != ')')
				zname.push_back(t0);
			eat(fin, '(');
			eat(fin, ')');
			eat(fin, ')');
			auto e = new XF_ZONE(zone, ztp, zname);
			add_entry(e);
			skip_white(fin);
			std::cout << "ZONE " << e->zone() << ", named " << R"(")" << e->zone_name() << R"(", )" << "is " << R"(")" << e->zone_type() << R"(")" << std::endl;
		}
		else
			throw std::runtime_error("Unsupported section index: " + std::to_string(ti));
	}

	// Close grid file
	fin.close();

	// Re-orginize grid connectivities in a much easier way,
	// and compute some derived quantities.
	std::cout << "Converting into high-level representation ..." << std::endl;
	raw2derived();
	std::cout << "Done!" << std::endl;

	// Finalize
	return 0;
}

int XF_MSH::writeToFile(const std::string &dst) const
{
	// Open grid file
	std::ofstream fout(dst);
	if (fout.fail())
		throw std::runtime_error("Failed to open output grid file: " + dst);

	if (!numOfCell() || !numOfFace() || !numOfNode())
		return -2;

	const size_t N = m_content.size();
	if (!N)
		return -3;

	size_t i = 0;

	// Write until dimension declaration
	while (true)
	{
		m_content[i]->repr(fout);
		bool flag = dynamic_cast<XF_DIMENSION*>(m_content[i]) != nullptr;
		++i;
		if (flag)
			break;
	}

	// Declaration of NODE, FACE, CELL
	fout << "(" << std::dec << XF_SECTION::NODE << " (";
	fout << std::hex << 0 << " " << 1 << " " << m_totalNodeNum << " ";
	fout << std::dec << 0 << " " << (m_is3D ? 3 : 2) << "))" << std::endl;
	fout << "(" << std::dec << XF_SECTION::CELL << " (";
	fout << std::hex << 0 << " " << 1 << " " << m_totalCellNum << " ";
	fout << std::dec << 0 << " " << 0 << "))" << std::endl;
	fout << "(" << std::dec << XF_SECTION::FACE << " (";
	fout << std::hex << 0 << " " << 1 << " " << m_totalFaceNum << " ";
	fout << std::dec << 0 << " " << 0 << "))" << std::endl;

	// Contents
	for (; i < N; ++i)
		m_content[i]->repr(fout);

	// Close grid file
	fout.close();

	// Finalize
	return 0;
}

void XF_MSH::raw2derived()
{
	// Allocate storage
	m_node.resize(numOfNode());
	m_face.resize(numOfFace());
	m_cell.resize(numOfCell());

	// Init
	for (auto &e : m_node)
	{
		e.z() = 0.0;
		e.atBdry = false;
	}
	for (auto &e : m_face)
	{
		e.center[2] = 0.0;
		e.n_LR[2] = 0.0;
		e.n_RL[2] = 0.0;
	}
	for (auto &e : m_cell)
	{
		e.center[2] = 0.0;
	}

	// Parse
	for (auto curPtr : m_content)
	{
		// Node
		if (curPtr->identity() == XF_SECTION::NODE)
		{
			auto curObj = dynamic_cast<XF_NODE*>(curPtr);

			// Node type within this zone
			const bool flag = curObj->is_boundary_node();

			// 1-based global node index
			const int cur_first = curObj->first_index();
			const int cur_last = curObj->last_index();

			for (int i = cur_first; i <= cur_last; ++i)
			{
				// Node Coordinates
				curObj->get_coordinate(i - cur_first, node(i).coordinate);

				// Node on boundary or not
				node(i).atBdry = flag;
			}
		}

		// Face
		if (curPtr->identity() == XF_SECTION::FACE)
		{
			auto curObj = dynamic_cast<XF_FACE*>(curPtr);

			// 1-based global face index
			const int cur_first = curObj->first_index();
			const int cur_last = curObj->last_index();

			// Face type of this zone
			const int ft = curObj->face_type();

			for (int i = cur_first; i <= cur_last; ++i)
			{
				const auto &cnct = curObj->connectivity(i - cur_first);

				// Face type
				if (cnct.x != ft)
					throw std::runtime_error("Internal error!");
				else
					face(i).type = ft;

				// Nodes within this face, 1-based node index are stored, right-hand convention is preserved.
				face(i).node.assign(cnct.n, cnct.n + cnct.x);

				// Adjacent cells, 1-based cell index are stored, 0 stands for boundary, right-hand convention is preserved.
				size_t lc = cnct.cl(), rc = cnct.cr();
				face(i).leftCell = lc;
				face(i).rightCell = rc;
				if (lc != 0)
					cell(lc).face.push_back(i);
				if (rc != 0)
					cell(rc).face.push_back(i);

				// Face on boundary or not
				face(i).atBdry = (cnct.c0() == 0 || cnct.c1() == 0);

				//Face area & center
				if (cnct.x == XF_FACE::LINEAR)
				{
					auto na = cnct.n[0], nb = cnct.n[1];
					auto p1 = node(na).coordinate, p2 = node(nb).coordinate;

					face(i).area = distance(p1, p2);
					line_center(p1, p2, face(i).center);
					line_normal(p1, p2, face(i).n_LR, face(i).n_RL);
				}
				else if (cnct.x == XF_FACE::TRIANGULAR)
				{
					auto na = cnct.n[0], nb = cnct.n[1], nc = cnct.n[2];
					auto p1 = node(na).coordinate, p2 = node(nb).coordinate, p3 = node(nc).coordinate;

					face(i).area = triangle_area(p1, p2, p3);
					triangle_center(p1, p2, p3, face(i).center);
					triangle_normal(p1, p2, p3, face(i).n_LR, face(i).n_RL);
				}
				else if (cnct.x == XF_FACE::QUADRILATERAL)
				{
					auto na = cnct.n[0], nb = cnct.n[1], nc = cnct.n[2], nd = cnct.n[3];
					auto p1 = node(na).coordinate, p2 = node(nb).coordinate, p3 = node(nc).coordinate, p4 = node(nd).coordinate;

					face(i).area = quadrilateral_area(p1, p2, p3, p4);
					quadrilateral_center(p1, p2, p3, p4, face(i).center);
					quadrilateral_normal(p1, p2, p3, p4, face(i).n_LR, face(i).n_RL);
				}
				else if (cnct.x == XF_FACE::POLYGONAL)
					throw std::runtime_error("Not supported currently!");
				else
					throw std::runtime_error("Internal error!");
			}
		}

		// Cell
		if (curPtr->identity() == XF_SECTION::CELL)
		{
			auto curObj = dynamic_cast<XF_CELL*>(curPtr);

			// 1-based global face index
			const int cur_first = curObj->first_index();
			const int cur_last = curObj->last_index();

			for (int i = cur_first; i <= cur_last; ++i)
			{
				auto &curCell = cell(i);

				// Element type of cells in this zone
				curCell.type = curObj->elem(i - cur_first);

				// Organize order of included nodes and faces
				switch (curCell.type)
				{
				case XF_CELL::TETRAHEDRAL:
					tet_standardization(curCell);
					break;
				case XF_CELL::HEXAHEDRAL:
					hex_standardization(curCell);
					break;
				case XF_CELL::PYRAMID:
					pyramid_standardization(curCell);
					break;
				case XF_CELL::WEDGE:
					prism_standardization(curCell);
					break;
				default:
					throw std::runtime_error("Invalid cell element type: " + std::to_string(cell(i).type) + ", maybe caused by internal error.");
				}
			}
		}
	}
}

void XF_MSH::derived2raw()
{
	// TODO
}
