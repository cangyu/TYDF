#ifndef __XF_MSH_H__
#define __XF_MSH_H__

#include <istream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <string>
#include <cstdint>
#include <algorithm>
#include <cmath>

typedef enum {
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
} XF_BC;

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

class XF_COMMENT : public XF_SECTION
{
private:
	std::string m_info;

public:
	XF_COMMENT(const std::string &info) : XF_SECTION(XF_SECTION::COMMENT), m_info(info) {}

	~XF_COMMENT() = default;

	const std::string &str() const { return m_info; }

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " \"" << str() << "\")" << std::endl;
	}
};

class XF_HEADER :public XF_SECTION
{
private:
	std::string m_msg;

public:
	XF_HEADER(const std::string &msg) : XF_SECTION(XF_SECTION::HEADER), m_msg(msg) {}

	~XF_HEADER() = default;

	const std::string &str() const { return m_msg; }

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " \"" << str() << "\")" << std::endl;
	}
};

class XF_DIMENSION :public XF_SECTION
{
private:
	bool m_is3D;
	int m_dim;

public:
	XF_DIMENSION(int dim) :
		XF_SECTION(XF_SECTION::DIMENSION)
	{
		if (dim == 2)
			m_is3D = false;
		else if (dim == 3)
			m_is3D = true;
		else
			throw("Invalid dimension!");

		m_dim = dim;
	}

	~XF_DIMENSION() = default;

	int ND() const { return m_dim; }

	bool is3D() const { return m_is3D; }

	void repr(std::ostream &out)
	{
		out << "(" << std::dec << identity() << " " << ND() << ")" << std::endl;
	}
};

class XF_MAIN_RECORD :public XF_SECTION
{
protected:
	int m_zone;
	int m_first, m_last;

public:
	XF_MAIN_RECORD(int id, int zone, int first, int last) :
		XF_SECTION(id),
		m_zone(zone),
		m_first(first),
		m_last(last)
	{
		if (m_first > m_last)
			throw("Invalid node index!");
	}

	virtual ~XF_MAIN_RECORD() = default;

	int zone() const { return m_zone; }

	int first_index() const { return m_first; }

	int last_index() const { return m_last; }

	int num() const { return (m_last - m_first + 1); }
};

class XF_NODE :public XF_MAIN_RECORD
{
private:
	int m_type;
	bool m_is3D;
	int m_dim;
	std::vector<double> m_node;

public:
	enum { VIRTUAL = 0, ANY = 1, BOUNDARY = 2 };

	XF_NODE(int zone, int first, int last, int type, int ND) :
		XF_MAIN_RECORD(XF_SECTION::NODE, zone, first, last)
	{
		if (type == 0)
			m_type = XF_NODE::VIRTUAL;
		else if (type == 1)
			m_type = XF_NODE::ANY;
		else if (type == 2)
			m_type = XF_NODE::BOUNDARY;
		else
			throw("Invalid description of node type!");

		if (ND == 2)
			m_is3D = false;
		else if (ND == 3)
			m_is3D = true;
		else
			throw("Invalid specification of node dimension!");
		m_dim = ND;

		m_node.resize(ND * num());
		std::fill(m_node.begin(), m_node.end(), 0.0);
	}

	~XF_NODE() = default;

	int type() const { return m_type; }

	int ND() const { return m_dim; }

	bool is3D() const { return m_is3D; }

	void get_coordinate(size_t loc_idx, std::vector<double> &dst) const
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

class XF_CELL :public XF_MAIN_RECORD
{
private:
	int m_type;
	int m_elem;
	std::vector<int> m_mixedElemDesc; // Only effective when 'm_elem == MIXED'.

public:
	enum { DEAD = 0, FLUID = 1, SOLID = 17 }; // Cell type.

	enum { MIXED = 0, TRIANGULAR = 1, TETRAHEDRAL = 2, QUADRILATERAL = 3, HEXAHEDRAL = 4, PYRAMID = 5, WEDGE = 6, POLYHEDRAL = 7 }; // Cell element type.

	XF_CELL(int zone, int first, int last, int type, int elem_type);

	~XF_CELL() = default;

	int type() const { return m_type; }

	int element_type() const { return m_elem; }

	int elem(size_t loc_idx) const { return m_mixedElemDesc[loc_idx]; }

	int &elem(size_t loc_idx) { return m_mixedElemDesc[loc_idx]; }

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

	static const std::string &cell_name(int x)
	{
		const static std::string sMIXED("MIXED");
		const static std::string sTRIANGULAR("TRIANGULAR");
		const static std::string sTETRAHEDRAL("TETRAHEDRAL");
		const static std::string sQUADRILATERAL("QUADRILATERAL");
		const static std::string sHEXAHEDRAL("HEXAHEDRAL");
		const static std::string sPYRAMID("PYRAMID");
		const static std::string sWEDGE("WEDGE");
		const static std::string sPOLYHEDRAL("POLYHEDRAL");

		switch (x)
		{
		case XF_CELL::MIXED:
			return sMIXED;
		case XF_CELL::TRIANGULAR:
			return sTRIANGULAR;
		case XF_CELL::TETRAHEDRAL:
			return sTETRAHEDRAL;
		case XF_CELL::QUADRILATERAL:
			return sQUADRILATERAL;
		case XF_CELL::HEXAHEDRAL:
			return sHEXAHEDRAL;
		case XF_CELL::PYRAMID:
			return sPYRAMID;
		case XF_CELL::WEDGE:
			return sWEDGE;
		case XF_CELL::POLYHEDRAL:
			return sPOLYHEDRAL;
		default:
			throw("Unsupported CELL!");
		}
	}
};

class XF_CONNECTIVITY
{
public:
	int x; // Length of n
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

class XF_FACE : public XF_MAIN_RECORD
{
private:
	XF_BC m_bc;
	int m_face;
	std::vector<XF_CONNECTIVITY> m_connectivity;

public:
	enum { MIXED = 0, LINEAR = 2, TRIANGULAR = 3, QUADRILATERAL = 4, POLYGONAL = 5 };

	XF_FACE(int zone, int first, int last, int bc, int face);

	~XF_FACE() = default;

	int bc_type() const { return m_bc; }

	int face_type() const { return m_face; }

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
	XF_ZONE(int zone, const std::string &type, const std::string &name) :
		XF_SECTION(XF_SECTION::ZONE),
		m_zoneID(zone),
		m_zoneType(type),
		m_zoneName(name),
		m_domainID(0)
	{
	}

	~XF_ZONE() = default;

	void repr(std::ostream &out)
	{
		out << std::dec;
		out << "(" << identity() << " (" << m_zoneID << " " << m_zoneType << " " << m_zoneName << ")())" << std::endl;
	}
};

class XF_MSH
{
private:
	std::vector<XF_SECTION*> m_content;
	size_t m_totalNodeNum, m_totalCellNum, m_totalFaceNum;
	bool m_is3D;
	int m_dim;

public:
	XF_MSH() :m_content(0, nullptr)
	{
		m_totalNodeNum = 0;
		m_totalCellNum = 0;
		m_totalFaceNum = 0;
		m_is3D = false;
		m_dim = 3;
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

	bool is3D() const { return m_is3D; }

	int dimension() const { return m_dim; }

	size_t numOfNode() const { return m_totalNodeNum; }

	size_t numOfFace() const { return m_totalFaceNum; }

	size_t numOfCell() const { return m_totalCellNum; }

	int computeTopology(
		std::vector<std::vector<double>> &nCoord, // Coordinates of each node.
		std::vector<bool> &nBdryFlag, // If located at boundary of each node.
		std::vector<std::vector<size_t>> &nAdjN, // Adjacent nodes of each node.
		std::vector<std::vector<size_t>> &nDepF, // Dependent faces of each node.
		std::vector<std::vector<size_t>> &nDepC, // Dependent cells of each node.
		std::vector<std::vector<size_t>> &fIncN, // Included nodes of each face, the order of nodes follows right-hand convention.
		std::vector<std::vector<size_t>> &fAdjC, // Adjacent cells of each face, the order of nodes follows right-hand convention.
		std::vector<double> &fArea, // Area of each face.
		std::vector<bool> &fBdryFlag, // If located at boundary of each face.
		std::vector<std::vector<double>> &fCenCoord, // Coordinates of each face centre.
		std::vector<std::vector<double>> &fUNLR, // Unit normal vector of each face, from c0/cl to c1/cr.
		std::vector<std::vector<double>> &fUNRL, // Unit normal vector of each face, from c1/cr to c0/cl.
		std::vector<std::vector<double>> &fNLR, // Normal vector of each face, from c0/cl to c1/cr, magnitude equals to face area.
		std::vector<std::vector<double>> &fNRL, // Normal vector of each face, from c1/cr to c0/cl, magnitude equals to face area.
		std::vector<std::vector<double>> &cCenCoord, // Coordinates of each cell centre.
		std::vector<double> &cVol, // Volumn of each face.
		std::vector<std::vector<size_t>> &cIncN, // Included nodes of each cell.
		std::vector<std::vector<size_t>> &cIncF, // Included faces of each cell.
		std::vector<std::vector<size_t>> &cAdjC, // Adjacent cells of each cell, the order is in accordance with cIncF.
		std::vector<std::vector<std::vector<double>>> &cFUNVec, // Positive unit normal vector of each included face of each cell, the order is in accordance with cIncF.
		std::vector<std::vector<std::vector<double>>> &cFNVec // Positive normal vector of each included face of each cell, the order is in accordance with cIncF, magnitude equals to face area.
	) const
	{
		int ret = 0;

		ret = computeTopology_nodeCoordinates(nCoord);
		if (!ret)
			throw(ret);

		ret = computeTopology_nodeBoundaryFlag(nBdryFlag);
		if (!ret)
			throw(ret);

		ret = computeTopology_nodeAdjacentNodes(nAdjN);
		if (!ret)
			throw(ret);

		ret = computeTopology_nodeDependentFaces(nDepF);
		if (!ret)
			throw(ret);

		ret = computeTopology_nodeDependentCells(nDepC);
		if (!ret)
			throw(ret);

		ret = computeTopology_faceIncludedNodes(fIncN);
		if (!ret)
			throw(ret);

		ret = computeTopology_faceAdjacentCells(fAdjC);
		if (!ret)
			throw(ret);

		ret = computeTopology_faceArea(nCoord, fArea);
		if (!ret)
			throw(ret);

		ret = computeTopology_faceBoundaryFlag(fBdryFlag);
		if (!ret)
			throw(ret);

		ret = computeTopology_faceCenterCoordinates(nCoord, fCenCoord);
		if (!ret)
			throw(ret);

		ret = computeTopology_faceUnitNormalVector(nCoord, fUNLR);
		if (!ret)
			throw(ret);

		const auto NF = numOfFace();
		const auto ND = dimension();
		for (size_t i = 0; i < NF; ++i)
			for (int j = 0; j < ND; ++j)
			{
				fUNRL[i][j] = -fUNLR[i][j];
				fNRL[i][j] = fUNRL[i][j] * fArea[i];
				fNLR[i][j] = fUNLR[i][j] * fArea[i];
			}

		ret = computeTopology_cellIncludedNodes(cIncN);
		if (!ret)
			throw(ret);

		ret = computeTopology_cellIncludedFaces(cIncF);
		if (!ret)
			throw(ret);

		ret = computeTopology_cellAdjacentCells(cIncF, fAdjC, cAdjC);
		if (!ret)
			throw(ret);

		ret = computeTopology_cellFaceNormal(cIncF, fAdjC, fNLR, fNRL, cFNVec);
		if (!ret)
			throw(ret);

		ret = computeTopology_cellFaceUnitNormal(cIncF, fArea, cFNVec, cFUNVec);
		if (!ret)
			throw(ret);

		if (is3D())
			ret = computeTopology_cellVolume3D(cIncF, fCenCoord, cFNVec, cVol);
		else
			ret = computeTopology_cellVolume2D(nCoord, cIncN, cIncF, cVol);
		if (!ret)
			throw(ret);

		ret = computeTopology_cellCentroidCoordinates(nCoord, cIncN, cCenCoord);
		if (!ret)
			throw(ret);

		return ret;
	}

private:
	void add_entry(XF_SECTION *e) { m_content.push_back(e); }

	void clear_entry()
	{
		// Release previous contents
		const size_t N = m_content.size();
		for (size_t i = 0; i < N; ++i)
			if (m_content[i])
				delete m_content[i];

		// Clear container
		m_content.clear();
	}

	int computeTopology_nodeCoordinates(std::vector<std::vector<double>> &dst) const;

	int computeTopology_nodeBoundaryFlag(std::vector<bool> &dst) const;

	int computeTopology_nodeAdjacentNodes(std::vector<std::vector<size_t>> &dst) const;

	int computeTopology_nodeDependentFaces(std::vector<std::vector<size_t>> &dst) const;

	int computeTopology_nodeDependentCells(std::vector<std::vector<size_t>> &dst) const;

	int computeTopology_faceIncludedNodes(std::vector<std::vector<size_t>> &dst) const;

	int computeTopology_faceAdjacentCells(std::vector<std::vector<size_t>> &dst) const;

	int computeTopology_faceArea(
		const std::vector<std::vector<double>> &nCoord,
		std::vector<double> &dst
	) const;

	int computeTopology_faceBoundaryFlag(std::vector<bool> &dst) const;

	int computeTopology_faceCenterCoordinates(
		const std::vector<std::vector<double>> &nCoord,
		std::vector<std::vector<double>> &dst
	) const;

	int computeTopology_faceUnitNormalVector(
		const std::vector<std::vector<double>> &nCoord,
		std::vector<std::vector<double>> &dst
	) const; // From left to right.

	int computeTopology_cellIncludedNodes(std::vector<std::vector<size_t>> &dst) const;

	int computeTopology_cellIncludedFaces(std::vector<std::vector<size_t>> &dst) const;

	int computeTopology_cellAdjacentCells(
		const std::vector<std::vector<size_t>> &cIncF,
		const std::vector<std::vector<size_t>> &fAdjC,
		std::vector<std::vector<size_t>> &dst
	) const;

	int computeTopology_cellFaceNormal(
		const std::vector<std::vector<size_t>> &cIncF,
		const std::vector<std::vector<size_t>> &fAdjC,
		const std::vector<std::vector<double>> &fNLR,
		const std::vector<std::vector<double>> &fNRL,
		std::vector<std::vector<std::vector<double>>> &dst
	) const;

	int computeTopology_cellFaceUnitNormal(
		const std::vector<std::vector<size_t>> &cIncF,
		const std::vector<double> &fArea,
		const std::vector<std::vector<std::vector<double>>> &cFNVec,
		std::vector<std::vector<std::vector<double>>> &dst
	) const;

	int computeTopology_cellVolume2D(
		const std::vector<std::vector<double>> &nCoord,
		const std::vector<std::vector<size_t>> &cIncN,
		const std::vector<std::vector<size_t>> &cIncF,
		std::vector<double> &dst
	) const;

	int computeTopology_cellVolume3D(
		const std::vector<std::vector<size_t>> &cIncF,
		const std::vector<std::vector<double>> &fCenCoord,
		const std::vector<std::vector<std::vector<double>>> &cFNVec,
		std::vector<double> &dst
	)const;

	int computeTopology_cellCentroidCoordinates(const std::vector<std::vector<double>> &nCoord, const std::vector<std::vector<size_t>> &cIncN, std::vector<std::vector<double>> &dst) const;
};

#endif
