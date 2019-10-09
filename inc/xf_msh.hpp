#ifndef __XF_MSH_HPP__
#define __XF_MSH_HPP__

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
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <stdexcept>

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

template<typename T>
class XF_Array2D
{
private:
	std::vector<T> m_data;
	size_t m_Nx, m_Ny;

	// Internal 0-based indexing interface.
	size_t idx(size_t i, size_t j) const { return i + m_Nx * j; }

public:
	XF_Array2D(size_t nx = 0, size_t ny = 0) : m_Nx(nx), m_Ny(ny), m_data(nx * ny) {}

	XF_Array2D(size_t nx, size_t ny, const T &val) : m_Nx(nx), m_Ny(ny), m_data(nx * ny, val) {}

	~XF_Array2D() = default;

	void resize(size_t nI, size_t nJ)
	{
		m_Nx = nI;
		m_Ny = nJ;
		m_data.resize(nI * nJ);
	}

	// 0-based indexing
	T &at(size_t i, size_t j) { return m_data[idx(i, j)]; }
	T at(size_t i, size_t j) const { return m_data[idx(i, j)]; }

	// 1-based indexing
	T &operator()(size_t i, size_t j) { return at(i - 1, j - 1); }
	T operator()(size_t i, size_t j) const { return at(i - 1, j - 1); }
};

template<typename T>
class XF_Array3D
{
private:
	std::vector<T> m_data;
	size_t m_Nx, m_Ny, m_Nz;

	// Internal 0-based indexing interface.
	size_t idx(size_t i, size_t j, size_t k) const { return i + m_Nx * (j + m_Ny * k); }

public:
	XF_Array3D(size_t nx = 0, size_t ny = 0, size_t nz = 0) : m_Nx(nx), m_Ny(ny), m_Nz(nz), m_data(nx * ny * nz) {}

	XF_Array3D(size_t nx, size_t ny, size_t nz, const T &val) : m_Nx(nx), m_Ny(ny), m_Nz(nz), m_data(nx * ny * nz, val) {}

	~XF_Array3D() = default;

	void resize(size_t nI, size_t nJ, size_t nK)
	{
		m_Nx = nI;
		m_Ny = nJ;
		m_Nz = nK;
		m_data.resize(nI * nJ * nK);
	}

	// 0-based indexing
	T &at(size_t i, size_t j, size_t k) { return m_data[idx(i, j, k)]; }
	T at(size_t i, size_t j, size_t k) const { return m_data[idx(i, j, k)]; }

	// 1-based indexing
	T &operator()(size_t i, size_t j, size_t k) { return at(i - 1, j - 1, k - 1); }
	T operator()(size_t i, size_t j, size_t k) const { return at(i - 1, j - 1, k - 1); }
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

	XF_BC() = default;

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

class XF_RANGE
{
protected:
	int m_zone;
	int m_first, m_last;

public:
	XF_RANGE(int zone, int first, int last) : m_zone(zone), m_first(first), m_last(last)
	{
		if (m_first > m_last)
			throw std::runtime_error("Invalid node index!");
	}

	virtual ~XF_RANGE() = default;

	int zone() const { return m_zone; }

	int first_index() const { return m_first; }

	int last_index() const { return m_last; }

	int num() const { return (last_index() - first_index() + 1); }
};

class XF_NODE :public XF_SECTION, public XF_RANGE, public XF_DIM
{
private:
	int m_type;
	std::vector<double> m_node;

public:
	enum { VIRTUAL = 0, ANY = 1, BOUNDARY = 2 };

	XF_NODE(int zone, int first, int last, int type, int ND) : XF_SECTION(XF_SECTION::NODE), XF_RANGE(zone, first, last), XF_DIM(ND)
	{
		if (type == 0)
			m_type = XF_NODE::VIRTUAL;
		else if (type == 1)
			m_type = XF_NODE::ANY;
		else if (type == 2)
			m_type = XF_NODE::BOUNDARY;
		else
			throw std::runtime_error("Invalid description of node type!");

		m_node.resize(ND * num());
		std::fill(m_node.begin(), m_node.end(), 0.0);
	}

	~XF_NODE() = default;

	int type() const { return m_type; }

	int ND() const { return dimension(); }

	void get_coordinate(size_t loc_idx, std::vector<double> &dst) const
	{
		size_t stx = STX(loc_idx);
		for (int i = 0; i < m_dim; ++i)
			dst[i] = m_node[stx + i];
	}

	void get_coordinate(size_t loc_idx, double *dst) const
	{
		size_t stx = STX(loc_idx);
		for (int i = 0; i < m_dim; ++i)
			dst[i] = m_node[stx + i];
	}

	void set_coordinate(size_t loc_idx, double x0, double x1, double x2)
	{
		const size_t stx = loc_idx * 3;
		m_node[stx] = x0;
		m_node[stx + 1] = x1;
		m_node[stx + 2] = x2;
	}

	void set_coordinate(size_t loc_idx, double x0, double x1)
	{
		const size_t stx = loc_idx * 2;
		m_node[stx] = x0;
		m_node[stx + 1] = x1;
	}

	void repr(std::ostream &out)
	{
		const int n_dim = ND();
		const int N = num();

		out << "(" << std::dec << identity();
		out << " (" << std::hex << zone() << " " << first_index() << " " << last_index() << " ";
		out << std::dec << type() << " " << n_dim << ")(" << std::endl;

		size_t loc_idx = 0;
		out.precision(12);
		for (int i = 0; i < N; ++i)
		{
			for (int k = 0; k < n_dim; ++k)
				out << " " << m_node[loc_idx + k];
			out << std::endl;
			loc_idx += n_dim;
		}

		out << "))" << std::endl;
	}

	bool is_virtual_node() const { return m_type == XF_NODE::VIRTUAL; }
	bool is_boundary_node() const { return m_type == XF_NODE::BOUNDARY; }
	bool is_internal_node() const { return m_type == XF_NODE::ANY; }

private:
	size_t STX(size_t loc_idx) const { return loc_idx * m_dim; }
};

class XF_CELL :public XF_SECTION, public XF_RANGE
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

	XF_CELL(int zone, int first, int last, int type, int elem_type) : XF_SECTION(XF_SECTION::CELL), XF_RANGE(zone, first, last)
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

	void set(int x, size_t *n, size_t *c)
	{
		this->x = x;
		this->c[0] = c[0];
		this->c[1] = c[1];

		int i;
		for (i = 0; i < x; ++i)
			this->n[i] = n[i];
		while (i < 4)
		{
			this->n[i] = 0;
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

class XF_FACE : public XF_SECTION, public XF_RANGE
{
private:
	int m_bc;
	int m_face;
	std::vector<XF_CONNECTIVITY> m_connectivity;

public:
	enum { MIXED = 0, LINEAR = 2, TRIANGULAR = 3, QUADRILATERAL = 4, POLYGONAL = 5 };

	static const std::map<int, std::string> MAPPING_Idx2Str;

	static const std::map<std::string, int> MAPPING_Str2Idx;

	XF_FACE(int zone, int first, int last, int bc, int face) : XF_SECTION(XF_SECTION::FACE), XF_RANGE(zone, first, last)
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
				for (int i = 0; i < loc_cnect.x; ++i)
					out << " " << loc_cnect.n[i];
				out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
			}
		}
		else
		{
			for (int i = 0; i < N; ++i)
			{
				const auto &loc_cnect = m_connectivity[i];
				for (int i = 0; i < loc_cnect.x; ++i)
					out << " " << loc_cnect.n[i];
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
	XF_ZONE(int zone, const std::string &type, const std::string &name) : XF_SECTION(XF_SECTION::ZONE)
	{
		m_zoneID = zone;
		m_zoneType = type;
		m_zoneName = name;
		m_domainID = 0;
	}

	~XF_ZONE() = default;

	int zone() const { return m_zoneID; }

	const std::string &zone_type() const { return m_zoneType; }

	const std::string &zone_name() const { return m_zoneName; }

	int domain() const { return m_domainID; }

	void repr(std::ostream &out)
	{
		out << std::dec;
		out << "(" << identity() << " (" << m_zoneID << " " << m_zoneType << " " << m_zoneName << ")())" << std::endl;
	}
};

class XF_MSH : public XF_DIM
{
private:
	struct NODE_ELEM
	{
		double coordinate[3];
		bool atBdry;

		double &x() { return coordinate[0]; }
		double &y() { return coordinate[1]; }
		double &z() { return coordinate[2]; }
	};

	struct FACE_ELEM
	{
		int type;
		double center[3];
		double area;
		XF_Array1D<size_t> node;
		size_t leftCell, rightCell;
		bool atBdry;
		double n_LR[3], n_RL[3]; // Surface unit normal
	};

	struct CELL_ELEM
	{
		int type;
		double center[3];
		double volume;
		XF_Array1D<size_t> face;
		XF_Array1D<size_t> node;
		XF_Array1D<size_t> adjCell;
		XF_Array2D<double> unitNormal;
	};

	// Raw
	std::vector<XF_SECTION*> m_content;
	size_t m_totalNodeNum, m_totalCellNum, m_totalFaceNum;

	// Derived
	XF_Array1D<NODE_ELEM> m_node;
	XF_Array1D<FACE_ELEM> m_face;
	XF_Array1D<CELL_ELEM> m_cell;

public:
	XF_MSH() : XF_DIM(3), m_content(0, nullptr)
	{
		m_totalNodeNum = 0;
		m_totalCellNum = 0;
		m_totalFaceNum = 0;
	}

	XF_MSH(const std::string &inp) : XF_DIM(3), m_content(0, nullptr)
	{
		m_totalNodeNum = 0;
		m_totalCellNum = 0;
		m_totalFaceNum = 0;

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

	// 1-based indexing
	NODE_ELEM &node(size_t idx) { return m_node(idx); }
	FACE_ELEM &face(size_t idx) { return m_face(idx); }
	CELL_ELEM &cell(size_t idx) { return m_cell(idx); }

private:
	static void eat(std::istream &in, char c)
	{
		char tmp;

		do {
			in >> tmp;
		} while (tmp != c);
	}

	static void skip_white(std::istream &in)
	{
		char tmp;

		do {
			in >> tmp;
		} while (tmp == ' ' || tmp == '\t' || tmp == '\n');

		if (!in.eof())
			in.unget();
	}

	static void cross_product(double *a, double *b, double *dst)
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
		delta(na, nb, dst);
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
		// Order of nodes follows the right-hand convention.
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

	static double quadrilateral_area(double *n1, double *n2, double *n3, double *n4)
	{
		// 1, 2, 3, 4 are in anti-clockwise direction.
		const double S123 = triangle_area(n1, n2, n3);
		const double S134 = triangle_area(n1, n3, n4);
		return S123 + S134;
	}

	static void quadrilateral_normal(double *n1, double *n2, double *n3, double *n4, double *dst, double *dst_r)
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

		// Parse record of node and face
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
					curObj->get_coordinate(i - cur_first, node(i).coordinate);

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
		}

		// Parse record of cell
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

	void derived2raw()
	{
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


	}
};

#endif
