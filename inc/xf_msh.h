#ifndef __XF_MSH_H__
#define __XF_MSH_H__

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
	T &operator()(size_t i)
	{
		return std::vector<T>::at(i - 1);
	}

	T operator()(size_t i) const
	{
		return std::vector<T>::at(i - 1);
	}
};

template<typename T>
class XF_Array2D
{
private:
	std::vector<T> m_data;
	size_t m_Nx, m_Ny;

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
	T &at(size_t i, size_t j)
	{
		return m_data[idx(i, j)];
	}

	T at(size_t i, size_t j) const
	{
		return m_data[idx(i, j)];
	}

	// 1-based indexing
	T &operator()(size_t i, size_t j)
	{
		return at(i - 1, j - 1);
	}

	T operator()(size_t i, size_t j) const
	{
		return at(i - 1, j - 1);
	}

private:
	// Internal 0-based indexing interface.
	size_t idx(size_t i, size_t j) const
	{
		return i + m_Nx * j;
	}
};

template<typename T>
class XF_Array3D
{
private:
	std::vector<T> m_data;
	size_t m_Nx, m_Ny, m_Nz;

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
	T &at(size_t i, size_t j, size_t k)
	{
		return m_data[idx(i, j, k)];
	}

	T at(size_t i, size_t j, size_t k) const
	{
		return m_data[idx(i, j, k)];
	}

	// 1-based indexing
	T &operator()(size_t i, size_t j, size_t k)
	{
		return at(i - 1, j - 1, k - 1);
	}

	T operator()(size_t i, size_t j, size_t k) const
	{
		return at(i - 1, j - 1, k - 1);
	}

private:
	// Internal 0-based indexing interface.
	size_t idx(size_t i, size_t j, size_t k) const
	{
		return i + m_Nx * (j + m_Ny * k);
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

	~XF_DIM() = default;

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
	XF_RANGE(int zone, int first, int last) :
		m_zone(zone),
		m_first(first),
		m_last(last)
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

	XF_NODE(int zone, int first, int last, int type, int ND) :
		XF_SECTION(XF_SECTION::NODE),
		XF_RANGE(zone, first, last),
		XF_DIM(ND)
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

	bool is_virtual_node() const
	{
		return m_type == XF_NODE::VIRTUAL;
	}

	bool is_boundary_node() const
	{
		return m_type == XF_NODE::BOUNDARY;
	}

	bool is_internal_node() const
	{
		return m_type == XF_NODE::ANY;
	}

private:
	size_t STX(size_t loc_idx) const
	{
		return loc_idx * m_dim;
	}
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

	XF_CELL(int zone, int first, int last, int type, int elem_type) :
		XF_SECTION(XF_SECTION::CELL),
		XF_RANGE(zone, first, last)
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

	int type() const
	{
		return m_type;
	}

	int element_type() const
	{
		return m_elem;
	}

	int elem(size_t loc_idx) const
	{
		return m_mixedElemDesc[loc_idx];
	}

	int &elem(size_t loc_idx)
	{
		return m_mixedElemDesc[loc_idx];
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

	XF_FACE(int zone, int first, int last, int bc, int face) :
		XF_SECTION(XF_SECTION::FACE),
		XF_RANGE(zone, first, last)
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

	const XF_CONNECTIVITY &connectivity(size_t loc_idx) const // 0-based local indexing
	{
		return m_connectivity[loc_idx];
	}

	XF_CONNECTIVITY &connectivity(size_t loc_idx) // 0-based local indexing
	{
		return m_connectivity[loc_idx];
	}

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
	XF_ZONE(int zone, const std::string &type, const std::string &name) :
		XF_SECTION(XF_SECTION::ZONE)
	{
		m_zoneID = zone;
		m_zoneType = type;
		m_zoneName = name;
		m_domainID = 0;
	}

	~XF_ZONE() = default;

	int zone() const
	{
		return m_zoneID;
	}

	const std::string &zone_type() const
	{
		return m_zoneType;
	}

	const std::string &zone_name() const
	{
		return m_zoneName;
	}

	int domain() const
	{
		return m_domainID;
	}

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
	XF_MSH() :
		XF_DIM(3),
		m_content(0, nullptr)
	{
		m_totalNodeNum = 0;
		m_totalCellNum = 0;
		m_totalFaceNum = 0;
	}

	~XF_MSH()
	{
		if (!m_content.empty())
		{
			for (size_t i = 0; i < m_content.size(); ++i)
				delete m_content[i];
		}
	}

	int readFromFile(const std::string &src);

	int writeToFile(const std::string &dst) const;

	size_t numOfNode() const { return m_totalNodeNum; }

	size_t numOfFace() const { return m_totalFaceNum; }

	size_t numOfCell() const { return m_totalCellNum; }

	NODE_ELEM &node(size_t idx) { return m_node(idx); }

	FACE_ELEM &face(size_t idx) { return m_face(idx); }

	CELL_ELEM &cell(size_t idx) { return m_cell(idx); }

private:
	void add_entry(XF_SECTION *e)
	{
		m_content.push_back(e);
	}

	void clear_entry()
	{
		// Release previous contents
		for (auto &ptr : m_content)
			if (ptr)
			{
				delete ptr;
				ptr = nullptr;
			}

		// Clear container
		m_content.clear();
	}

	void raw2derived();

	void derived2raw();
};

#endif
