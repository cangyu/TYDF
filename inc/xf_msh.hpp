#ifndef __XF_MSH_HPP__
#define __XF_MSH_HPP__

#include <cstddef>
#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <stdexcept>

typedef double XF_Scalar;

class XF_Vector
{
private:
	XF_Scalar m_data[3];

	void exchange(XF_Vector &rhs)
	{
		for (int i = 0; i < 3; ++i)
			std::swap(m_data[i], rhs.m_data[i]);
	}

public:
	XF_Vector() : m_data{ 0.0, 0.0, 0.0 } {}

	XF_Vector(XF_Scalar val) : m_data{ val, val, val } {}

	XF_Vector(XF_Scalar v1, XF_Scalar v2, XF_Scalar v3) : m_data{ v1, v2, v3 } {}

	XF_Vector(const XF_Vector &rhs) = default;

	~XF_Vector() = default;

	XF_Vector &operator=(XF_Vector rhs)
	{
		exchange(rhs);
		return *this;
	}

	XF_Scalar *data() { return m_data; }

	// 0-based indexing
	XF_Scalar at(size_t idx) const { return m_data[idx]; }
	XF_Scalar &at(size_t idx) { return m_data[idx]; }

	// 1-based indexing
	XF_Scalar operator()(int idx) const { return m_data[idx - 1]; }
	XF_Scalar &operator()(int idx) { return m_data[idx - 1]; }

	// Access through component
	XF_Scalar x() const { return m_data[0]; }
	XF_Scalar y() const { return m_data[1]; }
	XF_Scalar z() const { return m_data[2]; }

	XF_Scalar &x() { return m_data[0]; }
	XF_Scalar &y() { return m_data[1]; }
	XF_Scalar &z() { return m_data[2]; }
};

template<typename T>
class XF_Array1D : public std::vector<T>
{
public:
	XF_Array1D(size_t n = 0) : std::vector<T>(n) {}

	XF_Array1D(size_t n, const T &val) : std::vector<T>(n, val) {}

	~XF_Array1D() = default;

	// 1-based indexing
	T &operator()(size_t i) { return this->at(i - 1); }
	T operator()(size_t i) const { return this->at(i - 1); }

	// Check includances
	bool contains(const T &x) const
	{
		const int N = this->size();
		for (int i = 0; i < N; ++i)
			if (x == this->at(i))
				return true;

		return false;
	}

	bool contains(const T &a, const T &b) const
	{
		bool flag_a = false, flag_b = false;
		const int N = this->size();
		for (int i = 0; i < N; ++i)
		{
			const T &x = this->at(i);
			if (!flag_a && a == x)
				flag_a = true;
			if (!flag_b && b == x)
				flag_b = true;

			if (flag_a && flag_b)
				return true;
		}
		return false;
	}
};

class XF_BC
{
public:
	enum {
		INTERIOR = 2,
		WALL = 3,
		PRESSURE_INLET = 4,
		INLET_VENT = 4,
		INTAKE_FAN = 4,
		PRESSURE_OUTLET = 5,
		EXHAUST_FAN = 5,
		OUTLET_VENT = 5,
		SYMMETRY = 7,
		PERIODIC_SHADOW = 8,
		PRESSURE_FAR_FIELD = 9,
		VELOCITY_INLET = 10,
		PERIODIC = 12,
		FAN = 14,
		POROUS_JUMP = 14,
		RADIATOR = 14,
		MASS_FLOW_INLET = 20,
		INTERFACE = 24,
		PARENT = 31,
		OUTFLOW = 36,
		AXIS = 37
	};

	static const std::map<int, std::string> MAPPING_Idx2Str;

	static const std::map<std::string, int> MAPPING_Str2Idx;

	XF_BC() = delete;

	~XF_BC() = default;
};

class XF_SECTION
{
private:
	int m_identity;

public:
	enum { COMMENT = 0, HEADER = 1, DIMENSION = 2, NODE = 10, CELL = 12, FACE = 13, EDGE = 11, ZONE = 39 };

	XF_SECTION(int id) : m_identity(id) {}

	virtual ~XF_SECTION() = default;

	virtual void repr(std::ostream &out) = 0;

	int identity() const { return m_identity; }
};

class XF_STR
{
private:
	std::string m_msg;

public:
	XF_STR(const std::string &msg) : m_msg(msg) {}

	virtual ~XF_STR() = default;

	const std::string &str() const { return m_msg; }
};

class XF_COMMENT : public XF_SECTION, public XF_STR
{
public:
	XF_COMMENT(const std::string &info) : XF_SECTION(XF_SECTION::COMMENT), XF_STR(info) {}

	~XF_COMMENT() = default;

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " \"" << str() << "\")" << std::endl;
	}
};

class XF_HEADER : public XF_SECTION, public XF_STR
{
public:
	XF_HEADER(const std::string &info) : XF_SECTION(XF_SECTION::HEADER), XF_STR(info) {}

	~XF_HEADER() = default;

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " \"" << str() << "\")" << std::endl;
	}
};

class XF_DIM
{
protected:
	bool m_is3D;
	int m_dim;

public:
	XF_DIM(int dim)
	{
		if (dim == 2)
			m_is3D = false;
		else if (dim == 3)
			m_is3D = true;
		else
			throw std::runtime_error("Invalid dimension: " + std::to_string(dim));

		m_dim = dim;
	}

	XF_DIM(bool is3d) : m_is3D(is3d), m_dim(is3d ? 3 : 2) {}

	virtual ~XF_DIM() = default;

	bool is3D() const { return m_is3D; }

	int dimension() const { return m_dim; }
};

class XF_DIMENSION :public XF_SECTION, public XF_DIM
{
public:
	XF_DIMENSION(int dim) : XF_SECTION(XF_SECTION::DIMENSION), XF_DIM(dim) {}

	~XF_DIMENSION() = default;

	int ND() const { return dimension(); }

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " " << ND() << ")" << std::endl;
	}
};

class XF_RANGE : public XF_SECTION
{
protected:
	size_t m_zone;
	size_t m_first, m_last;

public:
	XF_RANGE(int id, size_t zone, size_t first, size_t last) : XF_SECTION(id), m_zone(zone), m_first(first), m_last(last)
	{
		if (first > last)
			throw std::runtime_error("Invalid node index!");
	}

	virtual ~XF_RANGE() = default;

	size_t zone() const { return m_zone; }

	size_t first_index() const { return m_first; }

	size_t last_index() const { return m_last; }

	size_t num() const { return (last_index() - first_index() + 1); }
};

class XF_NODE : public XF_RANGE, public XF_DIM
{
private:
	int m_type;
	XF_Array1D<XF_Vector> m_node;

public:
	enum { VIRTUAL = 0, ANY = 1, BOUNDARY = 2 };

	XF_NODE(int zone, int first, int last, int type, int ND) :
		XF_RANGE(XF_SECTION::NODE, zone, first, last),
		XF_DIM(ND),
		m_type(type),
		m_node(num())
	{
		if (type != VIRTUAL && type != ANY && type != BOUNDARY)
			throw std::runtime_error("Invalid description of node type!");
	}

	~XF_NODE() = default;

	int type() const { return m_type; }

	int ND() const { return dimension(); }

	void get_coordinate(size_t loc_idx, std::vector<double> &dst) const
	{
		const auto &node = m_node.at(loc_idx); // 0-based indexing
		for (int i = 0; i < m_dim; ++i)
			dst[i] = node.at(i);
	}

	void get_coordinate(size_t loc_idx, double *dst) const
	{
		const auto &node = m_node.at(loc_idx); // 0-based indexing
		for (int i = 0; i < m_dim; ++i)
			dst[i] = node.at(i);
	}

	void set_coordinate(size_t loc_idx, double x0, double x1, double x2 = 0.0)
	{
		auto &node = m_node.at(loc_idx);
		node.x() = x0;
		node.y() = x1;
		node.z() = x2;
	}

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity();
		out << " (" << std::hex << zone() << " " << first_index() << " " << last_index() << " ";
		out << std::dec << type() << " " << ND() << ")(" << std::endl;

		out.precision(12);
		const int N = num();
		for (int i = 0; i < N; ++i)
		{
			const auto &node = m_node.at(i);
			for (int k = 0; k < m_dim; ++k)
				out << " " << node.at(k);
			out << std::endl;
		}
		out << "))" << std::endl;
	}

	bool is_virtual_node() const { return m_type == XF_NODE::VIRTUAL; }
	bool is_boundary_node() const { return m_type == XF_NODE::BOUNDARY; }
	bool is_internal_node() const { return m_type == XF_NODE::ANY; }
};

class XF_CELL : public XF_RANGE
{
private:
	int m_type;
	int m_elem;
	std::vector<int> m_mixedElemDesc; // Only effective when 'm_elem == MIXED'.

public:
	enum { DEAD = 0, FLUID = 1, SOLID = 17 }; // Cell type.

	enum { MIXED = 0, TRIANGULAR = 1, TETRAHEDRAL = 2, QUADRILATERAL = 3, HEXAHEDRAL = 4, PYRAMID = 5, WEDGE = 6, POLYHEDRAL = 7 }; // Cell element type.

	static const std::map<int, std::string> TYPE_MAPPING_Idx2Str;

	static const std::map<std::string, int> TYPE_MAPPING_Str2Idx;

	static const std::map<int, std::string> ELEM_MAPPING_Idx2Str;

	static const std::map<std::string, int> ELEM_MAPPING_Str2Idx;

	XF_CELL(int zone, int first, int last, int type, int elem_type) : XF_RANGE(XF_SECTION::CELL, zone, first, last)
	{
		// Check cell type before assign
		auto it1 = XF_CELL::TYPE_MAPPING_Idx2Str.find(type);
		if (it1 == XF_CELL::TYPE_MAPPING_Idx2Str.end())
			throw std::runtime_error("Invalid cell type: " + std::to_string(type));
		else
			m_type = type;

		// Check cell elem before assign
		auto it2 = XF_CELL::ELEM_MAPPING_Idx2Str.find(elem_type);
		if (it2 == XF_CELL::ELEM_MAPPING_Idx2Str.end())
			throw std::runtime_error("Invalid cell element type: " + std::to_string(elem_type));
		else
			m_elem = elem_type;

		// Special treatment for mixed cell
		if (elem_type == XF_CELL::MIXED)
		{
			m_mixedElemDesc.resize(num());
			std::fill(m_mixedElemDesc.begin(), m_mixedElemDesc.end(), XF_CELL::MIXED);
		}
		else
			m_mixedElemDesc.resize(0);
	}

	~XF_CELL() = default;

	int type() const { return m_type; }
	int &type() { return m_type; }

	int element_type() const { return m_elem; }
	int &element_type() { return m_elem; }

	int elem(size_t loc_idx) const
	{
		int et = element_type();
		if (et == XF_CELL::MIXED)
			return m_mixedElemDesc[loc_idx];
		else
			return et;
	}

	int &elem(size_t loc_idx)
	{
		auto &et = element_type();
		if (et == XF_CELL::MIXED)
			return m_mixedElemDesc[loc_idx];
		else
			return et;
	}

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " (";
		out << std::hex;
		out << zone() << " " << first_index() << " " << last_index() << " ";
		out << m_type << " " << m_elem << ")";

		if (m_elem != XF_CELL::MIXED)
			out << ")" << std::endl;
		else
		{
			out << "(";

			const size_t N = num();
			for (size_t i = 0; i < N; ++i)
			{
				if (i % 40 == 0)
					out << std::endl;
				out << " " << elem(i);
			}
			out << std::endl << "))" << std::endl;
		}
	}
};

class XF_CONNECTIVITY
{
public:
	int x; // Num of nodes
	size_t n[4];
	size_t c[2];

public:
	XF_CONNECTIVITY() : x(1), n{ 0, 0, 0, 0 }, c{ 0, 0 } {}

	~XF_CONNECTIVITY() = default;

	size_t cl() const { return c[0]; }

	size_t cr() const { return c[1]; }

	size_t c0() const { return c[0]; }

	size_t c1() const { return c[1]; }

	void set(int x_, size_t *n_, size_t *c_)
	{
		x = x_;
		c[0] = c_[0];
		c[1] = c_[1];

		int i;
		for (i = 0; i < x_; ++i)
			n[i] = n_[i];
		while (i < 4)
		{
			n[i] = 0;
			++i;
		}
	}

	size_t leftAdj(int loc_idx) const
	{
		if (loc_idx == 0)
			return n[x - 1];
		else
			return n[loc_idx - 1];
	}

	size_t rightAdj(int loc_idx) const
	{
		if (loc_idx == x - 1)
			return n[0];
		else
			return n[loc_idx + 1];
	}
};

class XF_FACE : public XF_RANGE
{
private:
	int m_bc;
	int m_face;
	std::vector<XF_CONNECTIVITY> m_connectivity;

public:
	enum { MIXED = 0, LINEAR = 2, TRIANGULAR = 3, QUADRILATERAL = 4, POLYGONAL = 5 };

	static const std::map<int, std::string> MAPPING_Idx2Str;

	static const std::map<std::string, int> MAPPING_Str2Idx;

	XF_FACE(int zone, int first, int last, int bc, int face) : XF_RANGE(XF_SECTION::FACE, zone, first, last)
	{
		// Check B.C. before assign
		auto it1 = XF_BC::MAPPING_Idx2Str.find(bc);
		if (it1 == XF_BC::MAPPING_Idx2Str.end())
			throw std::runtime_error("Invalid B.C. type: " + std::to_string(bc));
		else
			m_bc = bc;

		// Check face type before assign
		auto it2 = XF_FACE::MAPPING_Idx2Str.find(face);
		if (it2 == XF_FACE::MAPPING_Idx2Str.end())
			throw std::runtime_error("Invalid face type: " + std::to_string(face));
		else
			m_face = face;

		// Exception currently
		if (m_face == XF_FACE::POLYGONAL)
			throw std::runtime_error("Not supported face type: " + XF_FACE::MAPPING_Idx2Str.at(m_face));

		// Resize local storage
		m_connectivity.resize(num());
	}

	~XF_FACE() = default;

	int bc_type() const { return m_bc; }

	int face_type() const { return m_face; }

	// 0-based local indexing
	const XF_CONNECTIVITY &connectivity(size_t loc_idx) const { return m_connectivity[loc_idx]; }
	XF_CONNECTIVITY &connectivity(size_t loc_idx) { return m_connectivity[loc_idx]; }

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " (";
		out << std::hex;
		out << zone() << " " << first_index() << " " << last_index() << " ";
		out << bc_type() << " " << face_type() << ")(" << std::endl;

		const int N = num();
		if (m_face == XF_FACE::MIXED)
		{
			for (int i = 0; i < N; ++i)
			{
				const auto &loc_cnect = m_connectivity[i];
				out << " " << loc_cnect.x;
				for (int j = 0; j < loc_cnect.x; ++j)
					out << " " << loc_cnect.n[j];
				out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
			}
		}
		else
		{
			for (int i = 0; i < N; ++i)
			{
				const auto &loc_cnect = m_connectivity[i];
				for (int j = 0; j < loc_cnect.x; ++j)
					out << " " << loc_cnect.n[j];
				out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
			}
		}

		out << "))" << std::endl;
	}
};

class XF_ZONE :public XF_SECTION
{
private:
	int m_zoneID;
	std::string m_zoneType, m_zoneName;
	int m_domainID;

public:
	XF_ZONE(int zone, const std::string &type, const std::string &name) : XF_SECTION(XF_SECTION::ZONE), m_zoneID(zone), m_zoneType(type), m_zoneName(name), m_domainID(0) {}

	~XF_ZONE() = default;

	int zone() const { return m_zoneID; }

	const std::string &type() const { return m_zoneType; }

	const std::string &name() const { return m_zoneName; }

	int domain() const { return m_domainID; }

	void repr(std::ostream &out)
	{
		out << std::dec << "(" << identity() << " (" << zone() << " " << type() << " " << name() << ")())" << std::endl;
	}
};

class XF_MSH : public XF_DIM
{
private:
	// Index of node, face, cell and zone all start from 1 and increase continuously.
	struct NODE_ELEM
	{
		XF_Vector coordinate;
		bool atBdry;
	};

	struct FACE_ELEM
	{
		int type;
		XF_Vector center;
		double area;
		XF_Array1D<size_t> node;
		size_t leftCell, rightCell;
		bool atBdry;
		XF_Vector n_LR, n_RL; // Surface unit normal
	};

	struct CELL_ELEM
	{
		int type;
		XF_Vector center;
		double volume;
		XF_Array1D<size_t> face;
		XF_Array1D<size_t> node;
		XF_Array1D<size_t> adjCell;
		XF_Array1D<XF_Vector> n, S;
	};

	struct ZONE_ELEM
	{
		std::string type;
		std::string name;
		XF_RANGE *obj;
	};

	// Raw
	std::vector<XF_SECTION*> m_content;
	size_t m_totalNodeNum, m_totalCellNum, m_totalFaceNum;

	// Derived
	XF_Array1D<NODE_ELEM> m_node;
	XF_Array1D<FACE_ELEM> m_face;
	XF_Array1D<CELL_ELEM> m_cell;

	size_t m_totalZoneNum;
	XF_Array1D<ZONE_ELEM> m_zone; // May be the group of nodes, the group of faces or the group of cells

public:
	XF_MSH() : XF_DIM(3), m_totalNodeNum(0), m_totalCellNum(0), m_totalFaceNum(0), m_totalZoneNum(0) {}

	XF_MSH(const std::string &inp) : XF_DIM(3), m_totalNodeNum(0), m_totalCellNum(0), m_totalFaceNum(0), m_totalZoneNum(0)
	{
		readFromFile(inp);
	}

	~XF_MSH() { clear_entry(); }

	int readFromFile(const std::string &src)
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
			int ti = -1;
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
					std::cout << "Reading " << e->num() << " nodes in zone " << zone << " (from " << first << " to " << last << ") ... ";

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
						std::cout << "Reading " << e->num() << " mixed cells in zone " << zone << " (from " << first << " to " << last << ") ... ";
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

					std::cout << "Reading " << e->num() << " faces in zone " << zone << " (from " << first << " to " << last << ") ... ";

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
				std::cout << "ZONE " << e->zone() << ", named " << R"(")" << e->name() << R"(", )" << "is " << R"(")" << e->type() << R"(")" << std::endl;
			}
			else
				throw std::runtime_error("Unsupported section index: " + std::to_string(ti));
		}

		// Close grid file
		fin.close();

		// Re-orginize grid connectivities in a much easier way,
		// and compute some derived quantities.
		std::cout << "Converting into high-level representation ... ";
		raw2derived();
		std::cout << "Done!" << std::endl;

		// Finalize
		return 0;
	}

	int writeToFile(const std::string &dst) const
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

	size_t numOfNode() const { return m_totalNodeNum; }
	size_t numOfFace() const { return m_totalFaceNum; }
	size_t numOfCell() const { return m_totalCellNum; }
	size_t numOfZone() const { return m_totalZoneNum; }

	// 1-based indexing
	NODE_ELEM &node(size_t idx) { return m_node(idx); }
	FACE_ELEM &face(size_t idx) { return m_face(idx); }
	CELL_ELEM &cell(size_t idx) { return m_cell(idx); }
	ZONE_ELEM &zone(size_t idx) { return m_zone(idx); }

private:
	static void eat(std::istream &in, char c)
	{
		char tmp = 0;
		do {
			in >> tmp;
		} while (tmp != c);
	}

	static void skip_white(std::istream &in)
	{
		char tmp = 0;
		do {
			in >> tmp;
		} while (tmp == ' ' || tmp == '\t' || tmp == '\n');

		if (!in.eof())
			in.unget();
	}

	static double dot_product(const XF_Vector &na, const XF_Vector &nb)
	{
		double ret = 0.0;
		for (int i = 1; i <= 3; ++i)
			ret += na(i) * nb(i);
		return ret;
	}

	static void cross_product(const double *a, const double *b, double *dst)
	{
		dst[0] = a[1] * b[2] - a[2] * b[1];
		dst[1] = a[2] * b[0] - a[0] * b[2];
		dst[2] = a[0] * b[1] - a[1] * b[0];
	}

	static void delta(double *na, double *nb, double *dst)
	{
		for (size_t i = 0; i < 3; ++i)
			dst[i] = nb[i] - na[i];
	}

	static void normalize(double *src, double *dst)
	{
		double L = 0.0;
		for (size_t i = 0; i < 3; ++i)
			L += src[i] * src[i];
		L = std::sqrt(L);
		for (size_t i = 0; i < 3; ++i)
			dst[i] = src[i] / L;
	}

	static double distance(double *na, double *nb)
	{
		double L = 0.0;
		for (size_t i = 0; i < 3; ++i)
		{
			double di = nb[i] - na[i];
			L += di * di;
		}
		return std::sqrt(L);
	}

	static void line_center(double *na, double *nb, double *dst)
	{
		for (size_t i = 0; i < 3; ++i)
			dst[i] = 0.5*(na[i] + nb[i]);
	}

	static void line_normal(double *na, double *nb, double *dst, double *dst_r)
	{
		// dst: unit normal vector from left cell to right cell.
		// dst_r: unit normal vector from right cell to left cell.

		delta(na, nb, dst);
		// Rotate 90 deg in clockwise direction
		std::swap(dst[0], dst[1]);
		dst[1] = -dst[1];
		normalize(dst, dst);

		for (size_t i = 0; i < 3; ++i)
			dst_r[i] = -dst[i];
	}

	static void triangle_center(double *na, double *nb, double *nc, double *dst)
	{
		for (size_t i = 0; i < 3; ++i)
			dst[i] = (na[i] + nb[i] + nc[i]) / 3.0;
	}

	static double triangle_area(double *na, double *nb, double *nc)
	{
		const double c = distance(na, nb);
		const double a = distance(nb, nc);
		const double b = distance(nc, na);
		const double p = 0.5*(a + b + c);
		return std::sqrt(p*(p - a)*(p - b)*(p - c)); // Heron's formula
	}

	static void triangle_normal(double *na, double *nb, double *nc, double *dst, double *dst_r)
	{
		// dst: unit normal vector from left cell to right cell.
		// dst_r: unit normal vector from right cell to left cell.
		// Order of "na, nb, nc" follows the right-hand convention.

		double rab[3], rac[3];
		delta(na, nb, rab);
		delta(na, nc, rac);
		cross_product(rac, rab, dst); // Take cross product to find normal direction
		normalize(dst, dst); // Normalize

		for (size_t i = 0; i < 3; ++i)
			dst_r[i] = -dst[i];
	}

	static void quadrilateral_center(double *n1, double *n2, double *n3, double *n4, double *dst)
	{
		// Order of "n1, n2, n3, n4" follows the right-hand convention.
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

	static double quadrilateral_area(double *n1, double *n2, double *n3, double *n4)
	{
		// Order of "n1, n2, n3, n4" follows the right-hand convention.
		const double S123 = triangle_area(n1, n2, n3);
		const double S134 = triangle_area(n1, n3, n4);
		return S123 + S134;
	}

	static void quadrilateral_normal(double *n1, double *n2, double *n3, double *n4, double *dst, double *dst_r)
	{
		// dst: unit normal vector from left cell to right cell.
		// dst_r: unit normal vector from right cell to left cell.
		// Order of "n1, n2, n3, n4" follows the right-hand convention.

		double ra[3] = { 0 }, rb[3] = { 0 };
		delta(n2, n4, ra);
		delta(n1, n3, rb);
		cross_product(ra, rb, dst);
		normalize(dst, dst);
		for (size_t i = 0; i < 3; ++i)
			dst_r[i] = -dst[i];
	}

	void add_entry(XF_SECTION *e) { m_content.push_back(e); }

	void clear_entry()
	{
		// Release previous contents
		for (auto ptr : m_content)
			if (ptr)
				delete ptr;

		// Clear container
		m_content.clear();
	}

	void raw2derived()
	{
		// Allocate storage
		m_node.resize(numOfNode());
		m_face.resize(numOfFace());
		m_cell.resize(numOfCell());

		// Set initial values
		for (auto &e : m_node)
		{
			e.coordinate.z() = 0.0;
			e.atBdry = false;
		}
		for (auto &e : m_face)
		{
			e.center.z() = 0.0;
			e.n_LR.z() = 0.0;
			e.n_RL.z() = 0.0;
		}
		for (auto &e : m_cell)
			e.center.z() = 0.0;

		// Parse records of node and face
		for (auto curPtr : m_content)
		{
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
					curObj->get_coordinate(i - cur_first, node(i).coordinate.data());

					// Node on boundary or not
					node(i).atBdry = flag;
				}
			}

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
						auto p1 = node(na).coordinate.data(), p2 = node(nb).coordinate.data();

						face(i).area = distance(p1, p2);
						line_center(p1, p2, face(i).center.data());
						line_normal(p1, p2, face(i).n_LR.data(), face(i).n_RL.data());
					}
					else if (cnct.x == XF_FACE::TRIANGULAR)
					{
						auto na = cnct.n[0], nb = cnct.n[1], nc = cnct.n[2];
						auto p1 = node(na).coordinate.data(), p2 = node(nb).coordinate.data(), p3 = node(nc).coordinate.data();

						face(i).area = triangle_area(p1, p2, p3);
						triangle_center(p1, p2, p3, face(i).center.data());
						triangle_normal(p1, p2, p3, face(i).n_LR.data(), face(i).n_RL.data());
					}
					else if (cnct.x == XF_FACE::QUADRILATERAL)
					{
						auto na = cnct.n[0], nb = cnct.n[1], nc = cnct.n[2], nd = cnct.n[3];
						auto p1 = node(na).coordinate.data(), p2 = node(nb).coordinate.data(), p3 = node(nc).coordinate.data(), p4 = node(nd).coordinate.data();

						face(i).area = quadrilateral_area(p1, p2, p3, p4);
						quadrilateral_center(p1, p2, p3, p4, face(i).center.data());
						quadrilateral_normal(p1, p2, p3, p4, face(i).n_LR.data(), face(i).n_RL.data());
					}
					else if (cnct.x == XF_FACE::POLYGONAL)
						throw std::runtime_error("Not supported currently!");
					else
						throw std::runtime_error("Internal error!");
				}
			}
		}

		// Parse records of cell
		for (auto curPtr : m_content)
		{
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
					cell_standardization(curCell);

					// Adjacent cells and normal
					curCell.adjCell.resize(curCell.face.size());
					curCell.n.resize(curCell.face.size());
					curCell.S.resize(curCell.face.size());
					for (size_t j = 0; j < curCell.face.size(); ++j)
					{
						const auto f_idx = curCell.face[j];
						const auto &f = face(f_idx);
						const auto c0 = f.leftCell, c1 = f.rightCell;
						if (c0 == i)
						{
							curCell.adjCell[j] = c1;
							curCell.n[j] = f.n_LR;
						}
						else if (c1 == i)
						{
							curCell.adjCell[j] = c0;
							curCell.n[j] = f.n_RL;
						}
						else
							throw std::runtime_error("Internal error.");

						for (int k = 1; k <= dimension(); ++k)
							curCell.S[j](k) = f.area * curCell.n[j](k);
					}

					// Volume and Centroid. 
					// Based on the divergence theorem. See (5.15) and (5.17) of Jiri Blazek's CFD book.
					curCell.volume = 0.0;
					curCell.center.x() = 0.0; curCell.center.y() = 0.0; curCell.center.z() = 0.0;
					for (size_t j = 0; j < curCell.face.size(); ++j)
					{
						const auto cfi = curCell.face.at(j);
						const auto &cf = face(cfi);
						const auto &cf_c = cf.center;
						const auto &cf_S = curCell.S.at(j);
						const auto w = dot_product(cf_c, cf_S);
						curCell.volume += w;
						for (int k = 1; k <= dimension(); ++k)
							curCell.center(k) += w * cf_c(k);
					}
					curCell.volume /= dimension();
					const double cde = (1.0 + dimension()) * curCell.volume;
					for (int k = 1; k <= dimension(); ++k)
						curCell.center(k) /= cde;
				}
			}
		}

		// Parse records of zone
		m_totalZoneNum = 0;
		std::map<int, int> tmp;
		for (auto curPtr : m_content)
		{
			auto curObj = dynamic_cast<XF_RANGE*>(curPtr);
			if (curObj == nullptr)
				continue;

			const auto curZone = curObj->zone();
			if (tmp.find(curZone) == tmp.end())
				tmp[curZone] = m_totalZoneNum++;
			else
				throw std::runtime_error("Duplicated zone detected.");
		}
		m_zone.clear();
		m_zone.resize(numOfZone());
		for (auto curPtr : m_content)
		{
			auto curObj = dynamic_cast<XF_RANGE*>(curPtr);
			if (curObj == nullptr)
				continue;

			const auto curZoneIdx = curObj->zone();
			if (curZoneIdx > numOfZone())
				throw std::runtime_error("Inconsistent zone index detected.");
			else
				zone(curZoneIdx).obj = curObj;
		}
		for (auto e : m_content)
		{
			const auto curObj = dynamic_cast<XF_ZONE*>(e);
			if (curObj == nullptr)
				continue;

			const auto curZoneIdx = curObj->zone();
			zone(curZoneIdx).name = curObj->name();
			zone(curZoneIdx).type = curObj->type();
		}
	}

	void derived2raw()
	{
		// Update size
		m_totalCellNum = m_cell.size();
		m_totalFaceNum = m_face.size();
		m_totalNodeNum = m_node.size();

		// TODO
	}

	void tet_standardization(CELL_ELEM &tet)
	{
		// Check num of total faces
		if (tet.face.size() != 4)
			throw std::runtime_error(R"(Mismatch between cell type ")" + XF_CELL::ELEM_MAPPING_Idx2Str.at(tet.type) + R"(" and num of faces: )" + std::to_string(tet.face.size()));

		// Ensure all faces are triangular
		const auto &f0 = face(tet.face.at(0));
		if (f0.type != XF_FACE::TRIANGULAR)
			throw std::runtime_error("Internal error.");

		const auto &f1 = face(tet.face.at(1));
		if (f1.type != XF_FACE::TRIANGULAR)
			throw std::runtime_error("Internal error.");

		const auto &f2 = face(tet.face.at(2));
		if (f2.type != XF_FACE::TRIANGULAR)
			throw std::runtime_error("Internal error.");

		const auto &f3 = face(tet.face.at(3));
		if (f3.type != XF_FACE::TRIANGULAR)
			throw std::runtime_error("Internal error.");

		// Find all 4 nodes
		size_t n1 = f0.node.at(0);
		if (n1 == 0)
			throw std::runtime_error("Internal error.");

		size_t n2 = f0.node.at(1);
		if (n2 == 0)
			throw std::runtime_error("Internal error.");

		size_t n3 = f0.node.at(2);
		if (n3 == 0)
			throw std::runtime_error("Internal error.");

		size_t n0 = 0;
		for (auto e : f2.node)
			if (e != n1 && e != n2 && e != n3)
			{
				n0 = e;
				break;
			}
		if (n0 == 0)
			throw std::runtime_error("Internal error.");

		// Assign node index
		tet.node.resize(4);
		tet.node.at(0) = n0;
		tet.node.at(1) = n1;
		tet.node.at(2) = n2;
		tet.node.at(3) = n3;
	}

	void pyramid_standardization(CELL_ELEM &pyramid)
	{
		// Check num of total faces
		if (pyramid.face.size() != 5)
			throw std::runtime_error(R"(Mismatch between cell type ")" + XF_CELL::ELEM_MAPPING_Idx2Str.at(pyramid.type) + R"(" and num of faces: )" + std::to_string(pyramid.face.size()));

		// Find the bottom quad and ensure other faces are triangular.
		size_t f0_idx = 0;
		for (auto e : pyramid.face)
		{
			const auto &f = face(e);
			if (f.type == XF_FACE::QUADRILATERAL)
			{
				if (f0_idx == 0)
					f0_idx = e;
				else
					throw std::runtime_error("Internal error.");
			}
			else if (f.type == XF_FACE::TRIANGULAR)
				continue;
			else
				throw std::runtime_error("Internal error.");
		}
		if (f0_idx == 0)
			throw std::runtime_error("Internal error.");

		// Nodes at bottom
		const auto &f0 = face(f0_idx);

		const size_t n0 = f0.node.at(0);
		if (n0 == 0)
			throw std::runtime_error("Internal error.");

		const size_t n1 = f0.node.at(1);
		if (n1 == 0)
			throw std::runtime_error("Internal error.");

		const size_t n2 = f0.node.at(2);
		if (n2 == 0)
			throw std::runtime_error("Internal error.");

		const size_t n3 = f0.node.at(3);
		if (n3 == 0)
			throw std::runtime_error("Internal error.");

		// Find other 4 triangles
		size_t f1_idx = 0;
		for (auto e : pyramid.face)
		{
			if (e == f0_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n0, n3))
			{
				f1_idx = e;
				break;
			}
		}
		if (f1_idx == 0)
			throw std::runtime_error("Internal error.");

		size_t f2_idx = 0;
		for (auto e : pyramid.face)
		{
			if (e == f0_idx || e == f1_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n3, n2))
			{
				f2_idx = e;
				break;
			}
		}
		if (f2_idx == 0)
			throw std::runtime_error("Internal error.");

		size_t f3_idx = 0;
		for (auto e : pyramid.face)
		{
			if (e == f0_idx || e == f1_idx || e == f2_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n2, n1))
			{
				f3_idx = e;
				break;
			}
		}
		if (f3_idx == 0)
			throw std::runtime_error("Internal error.");

		size_t f4_idx = 0;
		for (auto e : pyramid.face)
		{
			if (e != f0_idx && e != f1_idx && e != f2_idx && e != f3_idx)
			{
				f4_idx = e;
				const auto &f = face(e);
				if (!f.node.contains(n1, n0))
					throw std::runtime_error("Internal error.");

				break;
			}
		}
		if (f4_idx == 0)
			throw std::runtime_error("Internal error.");

		// Assign face index
		pyramid.face.at(0) = f0_idx;
		pyramid.face.at(1) = f1_idx;
		pyramid.face.at(2) = f2_idx;
		pyramid.face.at(3) = f3_idx;
		pyramid.face.at(4) = f4_idx;

		// The last node
		size_t n4 = 0;
		const auto &f1 = face(f1_idx);
		for (auto e : f1.node)
			if (e != n0 && e != n3)
			{
				n4 = e;
				break;
			}
		if (n4 == 0)
			throw std::runtime_error("Internal error.");

		// Assign node index
		pyramid.node.resize(5);
		pyramid.node.at(0) = n0;
		pyramid.node.at(1) = n1;
		pyramid.node.at(2) = n2;
		pyramid.node.at(3) = n3;
		pyramid.node.at(4) = n4;
	}

	void prism_standardization(CELL_ELEM &prism)
	{
		// Check num of total faces
		if (prism.face.size() != 5)
			throw std::runtime_error(R"(Mismatch between cell type ")" + XF_CELL::ELEM_MAPPING_Idx2Str.at(prism.type) + R"(" and num of faces: )" + std::to_string(prism.face.size()));

		// Ensure there're only 2 triangle and 3 quad
		size_t f0_idx = 0, f1_idx = 0;
		for (auto e : prism.face)
		{
			const auto &f = face(e);
			if (f.type == XF_FACE::TRIANGULAR)
			{
				if (f0_idx == 0)
					f0_idx = e;
				else if (f1_idx == 0)
					f1_idx = e;
				else
					throw std::runtime_error("There're more than 2 triangular faces in a prism cell.");
			}
			else if (f.type == XF_FACE::QUADRILATERAL)
				continue;
			else
				throw std::runtime_error("Internal error.");
		}
		if (f0_idx == 0 || f1_idx == 0)
			throw std::runtime_error("Missing triangular faces in a prism cell.");

		// 2 triangular faces
		const auto &f0 = face(f0_idx);
		const auto &f1 = face(f1_idx);

		// 3 nodes on the bottom triangular face
		const size_t n0 = f0.node.at(0);
		const size_t n1 = f0.node.at(1);
		const size_t n2 = f0.node.at(2);

		// Find face 4
		size_t f4_idx = 0;
		for (auto e : prism.face)
		{
			const auto &f = face(e);
			if (f.type == XF_FACE::QUADRILATERAL && f.node.contains(n0, n1))
			{
				if (f4_idx == 0)
					f4_idx = e;
				else
					throw std::runtime_error("Inconsistent face composition.");
			}
		}
		if (f4_idx == 0)
			throw std::runtime_error("Missing face 4.");

		const auto &f4 = face(f4_idx);

		// Find face 3
		size_t f3_idx = 0;
		for (auto e : prism.face)
		{
			const auto &f = face(e);
			if (f.type == XF_FACE::QUADRILATERAL && f.node.contains(n1, n2))
			{
				if (f3_idx == 0)
					f3_idx = e;
				else
					throw std::runtime_error("Inconsistent face composition.");
			}
		}
		if (f3_idx == 0 || f3_idx == f4_idx)
			throw std::runtime_error("Missing face 3.");

		const auto &f3 = face(f3_idx);

		// Find face 2
		size_t f2_idx = 0;
		for (auto e : prism.face)
		{
			const auto &f = face(e);
			if (f.type == XF_FACE::QUADRILATERAL)
			{
				if (e != f4_idx && e != f3_idx)
				{
					f2_idx = e;
					break;
				}
			}
		}
		if (f2_idx == 0)
			throw std::runtime_error("Missing face 2.");

		const auto &f2 = face(f2_idx);
		if (!f2.node.contains(n2, n0))
			throw std::runtime_error("Inconsistent face composition.");

		// Assign face index
		prism.face.at(0) = f0_idx;
		prism.face.at(1) = f1_idx;
		prism.face.at(2) = f2_idx;
		prism.face.at(3) = f3_idx;
		prism.face.at(4) = f4_idx;

		// 3 nodes on the top triangular face to be defined
		size_t n3 = 0, n5 = 0, n4 = 0;
		for (auto e : f2.node)
		{
			if (e != n0 && e != n2)
			{
				if (n3 == 0)
					n3 = e;
				else if (n5 == 0)
					n5 = e;
				else
					throw std::runtime_error("Internal error.");
			}
		}
		if (n3 == 0 || n5 == 0 || n3 == n5)
			throw std::runtime_error("Missing nodes on the top");

		// Check if n3 and n5 needs to be swaped
		if (!(f1.node.contains(n3) && f4.node.contains(n3)))
			std::swap(n3, n5);

		// Ensure n5 located on f1, f2 and f3
		if (!(f1.node.contains(n5) && f3.node.contains(n5)))
			throw std::runtime_error("Internal error.");

		// Find n4
		for (auto e : f1.node)
		{
			if (e != n3 && e != n5)
			{
				if (n4 == 0)
					n4 = e;
				else
					throw std::runtime_error("Internal error.");
			}
		}
		if (n4 == 0)
			throw std::runtime_error("Internal error.");

		// Assign node index
		prism.node.resize(6);
		prism.node.at(0) = n0;
		prism.node.at(1) = n1;
		prism.node.at(2) = n2;
		prism.node.at(3) = n3;
		prism.node.at(4) = n4;
		prism.node.at(5) = n5;
	}

	void hex_standardization(CELL_ELEM &hex)
	{
		// Check num of total faces
		if (hex.face.size() != 6)
			throw std::runtime_error(R"(Mismatch between cell type ")" + XF_CELL::ELEM_MAPPING_Idx2Str.at(hex.type) + R"(" and num of faces: )" + std::to_string(hex.face.size()));

		// Ensure all faces are quad
		for (auto e : hex.face)
		{
			if (e == 0)
				throw std::runtime_error("Internal error.");
			const auto &f = face(e);
			if (f.type != XF_FACE::QUADRILATERAL)
				throw std::runtime_error(R"(Inconsistent face type ")" + XF_FACE::MAPPING_Idx2Str.at(f.type) + R"(" in a hex cell.)");
		}

		// Face 4 at bottom
		const size_t f4_idx = hex.face.at(0);
		const auto &f4 = face(f4_idx);
		const size_t n0 = f4.node.at(0);
		const size_t n1 = f4.node.at(1);
		const size_t n2 = f4.node.at(2);
		const size_t n3 = f4.node.at(3);
		if (n0 == 0 || n1 == 0 || n2 == 0 || n3 == 0)
			throw std::runtime_error("Internal error.");

		// Face 0
		size_t f0_idx = 0;
		for (auto e : hex.face)
		{
			if (e == f4_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n3, n0))
			{
				if (f0_idx == 0)
					f0_idx = e;
				else
					throw std::runtime_error("Duplicated face detected.");
			}
		}
		if (f0_idx == 0)
			throw std::runtime_error("Missing face 0");

		const auto &f0 = face(f0_idx);
		if (f0.node.size() != 4)
			throw std::runtime_error("Internal error.");

		// Face 2
		size_t f2_idx = 0;
		for (auto e : hex.face)
		{
			if (e == f4_idx || e == f0_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n0, n1))
			{
				if (f2_idx == 0)
					f2_idx = e;
				else
					throw std::runtime_error("Duplicated face detected.");
			}
		}
		if (f2_idx == 0)
			throw std::runtime_error("Missing face 2");

		const auto &f2 = face(f2_idx);
		if (f2.node.size() != 4)
			throw std::runtime_error("Internal error.");

		// Face 1
		size_t f1_idx = 0;
		for (auto e : hex.face)
		{
			if (e == f4_idx || e == f0_idx || e == f2_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n1, n2))
			{
				if (f1_idx == 0)
					f1_idx = e;
				else
					throw std::runtime_error("Duplicated face detected.");
			}
		}
		if (f1_idx == 0)
			throw std::runtime_error("Missing face 1");

		const auto &f1 = face(f1_idx);
		if (f1.node.size() != 4)
			throw std::runtime_error("Internal error.");

		// Face 3
		size_t f3_idx = 0;
		for (auto e : hex.face)
		{
			if (e == f4_idx || e == f0_idx || e == f2_idx || e == f1_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n2, n3))
			{
				if (f3_idx == 0)
					f3_idx = e;
				else
					throw std::runtime_error("Duplicated face detected.");
			}
		}
		if (f3_idx == 0)
			throw std::runtime_error("Missing face 3");

		const auto &f3 = face(f3_idx);
		if (f3.node.size() != 4)
			throw std::runtime_error("Internal error.");

		// Face 5
		size_t f5_idx = 0;
		for (auto e : hex.face)
		{
			if (e == f4_idx || e == f0_idx || e == f2_idx || e == f1_idx || e == f3_idx)
				continue;

			f5_idx = e;
			break;
		}
		if (f5_idx == 0)
			throw std::runtime_error("Missing face 5");

		const auto &f5 = face(f5_idx);
		if (f5.node.size() != 4)
			throw std::runtime_error("Internal error.");

		// Assign face index
		hex.face.at(0) = f0_idx;
		hex.face.at(1) = f1_idx;
		hex.face.at(2) = f2_idx;
		hex.face.at(3) = f3_idx;
		hex.face.at(4) = f4_idx;
		hex.face.at(5) = f5_idx;

		// 4 nodes at top to be defined
		size_t n4 = 0, n5 = 0, n6 = 0, n7 = 0;

		// Node 4
		for (auto e : f5.node)
		{
			if (f0.node.contains(e) && f2.node.contains(e))
			{
				n4 = e;
				break;
			}
		}
		if (n4 == 0 || n4 == n0)
			throw std::runtime_error("Missing node 4");

		// Node 5
		for (auto e : f2.node)
		{
			if (e == n0 || e == n1 || e == n4)
				continue;

			n5 = e;
			break;
		}
		if (n5 == 0)
			throw std::runtime_error("Missing node 5");

		// Node 6
		for (auto e : f1.node)
		{
			if (e == n1 || e == n2 || e == n5)
				continue;

			n6 = e;
			break;
		}
		if (n6 == 0)
			throw std::runtime_error("Missing node 6");

		// Node 7
		for (auto e : f0.node)
		{
			if (e == n0 || e == n3 || e == n4)
				continue;

			n7 = e;
			break;
		}
		if (n7 == 0)
			throw std::runtime_error("Missing node 7");

		if (!f3.node.contains(n6, n7))
			throw std::runtime_error("Inconsistent node composition.");

		// Assign node index
		hex.node.resize(8);
		hex.node.at(0) = n0;
		hex.node.at(1) = n1;
		hex.node.at(2) = n2;
		hex.node.at(3) = n3;
		hex.node.at(4) = n4;
		hex.node.at(5) = n5;
		hex.node.at(6) = n6;
		hex.node.at(7) = n7;
	}

	void triangle_standardization(CELL_ELEM &tri)
	{
		// Check num of total faces
		if (tri.face.size() != 3)
			throw std::runtime_error("Inconsitent face composition.");

		// Ensure all faces are lines
		for (auto e : tri.face)
		{
			const auto &f = face(e);
			if (e == 0 || f.type != XF_FACE::LINEAR || f.node.size() != 2)
				throw std::runtime_error("Invalid face detected.");
		}

		// Faces
		// Keep the order of faces as it is.
		// Doesn't matter as any order of the 3 faces is valid.
		const auto f0_idx = tri.face.at(0);
		const auto f1_idx = tri.face.at(1);
		const auto f2_idx = tri.face.at(2);

		const auto &f0 = face(f0_idx);
		const auto &f1 = face(f1_idx);
		const auto &f2 = face(f2_idx);

		// Nodes
		const size_t n0 = f0.node.at(0);
		const size_t n1 = f0.node.at(1);
		if (n0 == 0 || n1 == 0 || n0 == n1)
			throw std::runtime_error("Missing node detected.");

		size_t n2 = 0;
		for (auto e : f1.node)
		{
			if (f0.node.contains(e))
				continue;

			n2 = e;
			break;
		}
		if (n2 == 0 || n2 == n0 || n2 == n1 || !f2.node.contains(n2))
			throw std::runtime_error("Node 2 is missing.");

		// Assign node index
		tri.node.resize(3);
		tri.node.at(0) = n0;
		tri.node.at(1) = n1;
		tri.node.at(2) = n2;
	}

	void quad_standardization(CELL_ELEM &quad)
	{
		// Check num of total faces
		if (quad.face.size() != 4)
			throw std::runtime_error("Inconsitent face composition.");

		// Ensure all faces are lines
		for (auto e : quad.face)
		{
			const auto &f = face(e);
			if (e == 0 || f.type != XF_FACE::LINEAR || f.node.size() != 2)
				throw std::runtime_error("Invalid face detected.");
		}

		// Face 0
		const auto f0_idx = quad.face.at(0);
		const auto &f0 = face(f0_idx);

		// Node 0 and 1
		const auto n0 = f0.node.at(0);
		const auto n1 = f0.node.at(1);

		// Face 1
		size_t f1_idx = 0;
		for (auto e : quad.face)
		{
			if (e == f0_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n1))
			{
				f1_idx = e;
				break;
			}
		}
		if (f1_idx == 0)
			throw std::runtime_error("Missing face 1");

		const auto &f1 = face(f1_idx);

		// Node 2
		const size_t n2 = f1.node.at(0) == n1 ? f1.node.at(1) : f1.node.at(0);

		// Face 3
		size_t f3_idx = 0;
		for (auto e : quad.face)
		{
			if (e == f0_idx || e == f1_idx)
				continue;

			const auto &f = face(e);
			if (f.node.contains(n0))
			{
				f3_idx = e;
				break;
			}
		}
		if (f3_idx == 0)
			throw std::runtime_error("Missing face 3");

		const auto &f3 = face(f3_idx);

		// Node 3
		const size_t n3 = f3.node.at(0) == n0 ? f3.node.at(1) : f3.node.at(0);

		// Check face 2
		size_t f2_idx = 0;
		for (auto e : quad.face)
		{
			if (e == f0_idx || e == f1_idx || e == f3_idx)
				continue;

			f2_idx = e;
			break;
		}
		if (f2_idx == 0)
			throw std::runtime_error("Missing face 2");

		const auto &f2 = face(f2_idx);
		if (!f2.node.contains(n2, n3))
			throw std::runtime_error("Inconsistent node and face includance on face 2");

		// Assign face index
		quad.face.at(1) = f1_idx;
		quad.face.at(2) = f2_idx;
		quad.face.at(3) = f3_idx;

		// Assign node index
		quad.node.resize(4);
		quad.node.at(0) = n0;
		quad.node.at(1) = n1;
		quad.node.at(2) = n2;
		quad.node.at(3) = n3;
	}

	void cell_standardization(CELL_ELEM &c)
	{
		switch (c.type)
		{
		case XF_CELL::TETRAHEDRAL:
			tet_standardization(c);
			break;
		case XF_CELL::HEXAHEDRAL:
			hex_standardization(c);
			break;
		case XF_CELL::PYRAMID:
			pyramid_standardization(c);
			break;
		case XF_CELL::WEDGE:
			prism_standardization(c);
			break;
		case XF_CELL::TRIANGULAR:
			triangle_standardization(c);
			break;
		case XF_CELL::QUADRILATERAL:
			quad_standardization(c);
		default:
			throw std::runtime_error("Invalid cell element type: " + std::to_string(c.type) + ", maybe caused by internal error.");
		}
	}
};

#endif
