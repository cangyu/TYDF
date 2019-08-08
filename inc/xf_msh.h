#ifndef __XF_MSH_H__
#define __XF_MSH_H__

#include <vector>
#include <string>
#include <istream>
#include <ostream>
#include <cstdint>

typedef enum {
	COMMENT = 0,
	HEADER = 1,
	DIMENSION = 2,
	NODE = 10,
	CELL = 12,
	FACE = 13,
	EDGE = 11,
	ZONE = 39
} XF_SECTION;

typedef enum {
	VIRTUAL = 0,
	ANY = 1,
	BOUNDARY = 2
} XF_NODE_TYPE;

typedef enum {
	DEAD = 0,
	FLUID = 1,
	SOLID = 17
} XF_CELL_TYPE;

typedef enum {
	MIXED = 0,
	TRIANGULAR = 1,
	TETRAHEDRAL = 2,
	QUADRILATERAL = 3,
	HEXAHEDRAL = 4,
	PYRAMID = 5,
	WEDGE = 6,
	POLYHEDRAL = 7
} XF_CELL_ELEM_TYPE;

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
} XF_BC_TYPE;

typedef enum {
	F_MIXED = 0,
	F_LINEAR = 2,
	F_TRIANGULAR = 3,
	F_QUADRILATERAL = 4,
	F_POLYGONAL = 5
} XF_FACE_TYPE;

class XF_ENTRY
{
protected:
	XF_SECTION m_identity;

public:
	XF_ENTRY(XF_SECTION id) :
		m_identity(id)
	{
	}

	virtual ~XF_ENTRY() = default;

	virtual void repr(std::ostream &out) = 0;
};

class XF_COMMENT : public XF_ENTRY
{
private:
	std::string m_info;

public:
	XF_COMMENT() :
		XF_ENTRY(XF_SECTION::COMMENT),
		m_info("")
	{
	}

	XF_COMMENT(const std::string &info) :
		XF_ENTRY(XF_SECTION::COMMENT),
		m_info(info)
	{
	}

	~XF_COMMENT() = default;

	const std::string &str() const
	{
		return m_info;
	}

	void repr(std::ostream &out);
};

class XF_HEADER :public XF_ENTRY
{
private:
	std::string m_msg;

public:
	XF_HEADER() :
		XF_ENTRY(XF_SECTION::HEADER),
		m_msg("")
	{
	}

	XF_HEADER(const std::string &msg) :
		XF_ENTRY(XF_SECTION::HEADER),
		m_msg(msg)
	{
	}

	~XF_HEADER() = default;

	const std::string &str() const
	{
		return m_msg;
	}

	void repr(std::ostream &out);
};

class XF_DIMENSION :public XF_ENTRY
{
private:
	bool m_is3D;

public:
	XF_DIMENSION() :
		XF_ENTRY(XF_SECTION::DIMENSION),
		m_is3D(true)
	{
	}

	XF_DIMENSION(int dim);

	~XF_DIMENSION() = default;

	int ND() const
	{
		if (m_is3D)
			return 3;
		else
			return 2;
	}

	void repr(std::ostream &out);
};

class XF_MAIN_RECORD :public XF_ENTRY
{
protected:
	int m_zone;
	int m_first, m_last;

public:
	XF_MAIN_RECORD(XF_SECTION id, int zone, int first, int last) :
		XF_ENTRY(id),
		m_zone(zone),
		m_first(first),
		m_last(last)
	{
		if (m_first > m_last)
			throw("Invalid node index!");
	}

	~XF_MAIN_RECORD() = default;

	int zone() const
	{
		return m_zone;
	}

	int first_index() const
	{
		return m_first;
	}

	int last_index() const
	{
		return m_last;
	}

	int num() const
	{
		return (m_last - m_first + 1);
	}
};

class XF_NODE :public XF_MAIN_RECORD
{
private:
	XF_NODE_TYPE m_type;
	bool m_is3D;
	std::vector<double> m_pnt;

public:
	XF_NODE(int zone, int first, int last, int type, int ND);

	~XF_NODE() = default;

	int type() const
	{
		return m_type;
	}

	int ND() const
	{
		if (m_is3D)
			return 3;
		else
			return 2;
	}

	void record_pnt_coordinate(size_t loc_idx, double x0, double x1, double x2);

	void record_pnt_coordinate(size_t loc_idx, double x0, double x1);

	void repr(std::ostream &out);
};

class XF_CELL :public XF_MAIN_RECORD
{
private:
	XF_CELL_TYPE m_type;
	XF_CELL_ELEM_TYPE m_elem;
	std::vector<XF_CELL_ELEM_TYPE> m_mixedElemDesc;

public:
	XF_CELL(int zone, int first, int last, int type, int elem_type);

	~XF_CELL() = default;

	int type() const
	{
		return m_type;
	}

	int element_type() const
	{
		return m_elem;
	}

	XF_CELL_ELEM_TYPE elem(size_t loc_idx) const
	{
		return m_mixedElemDesc[loc_idx];
	}

	XF_CELL_ELEM_TYPE &elem(size_t loc_idx)
	{
		return m_mixedElemDesc[loc_idx];
	}

	void repr(std::ostream &out);

	static const std::string &cell_name(int x)
	{
		const static std::string MIXED("MIXED");
		const static std::string TRIANGULAR("TRIANGULAR");
		const static std::string TETRAHEDRAL("TETRAHEDRAL");
		const static std::string QUADRILATERAL("QUADRILATERAL");
		const static std::string HEXAHEDRAL("HEXAHEDRAL");
		const static std::string PYRAMID("PYRAMID");
		const static std::string WEDGE("WEDGE");
		const static std::string POLYHEDRAL("POLYHEDRAL");

		switch (x)
		{
		case XF_CELL_ELEM_TYPE::MIXED:
			return MIXED;
		case XF_CELL_ELEM_TYPE::TRIANGULAR:
			return TRIANGULAR;
		case XF_CELL_ELEM_TYPE::TETRAHEDRAL:
			return TETRAHEDRAL;
		case XF_CELL_ELEM_TYPE::QUADRILATERAL:
			return QUADRILATERAL;
		case XF_CELL_ELEM_TYPE::HEXAHEDRAL:
			return HEXAHEDRAL;
		case XF_CELL_ELEM_TYPE::PYRAMID:
			return PYRAMID;
		case XF_CELL_ELEM_TYPE::WEDGE:
			return WEDGE;
		case XF_CELL_ELEM_TYPE::POLYHEDRAL:
			return POLYHEDRAL;
		default:
			throw("Unsupported \"XF_CELL_ELEM_TYPE\" input!");
		}
	}
};

class XF_CONNECTIVITY
{
public:
	uint8_t x;
	size_t n[4];
	size_t c[2];

public:
	XF_CONNECTIVITY() :
		x(1),
		n{ 0, 0, 0, 0 },
		c{ 0, 0 }
	{
	}

	~XF_CONNECTIVITY() = default;

	void set(uint8_t x, size_t *n, size_t *c)
	{
		this->x = x;
		this->c[0] = c[0];
		this->c[1] = c[1];

		uint8_t i;
		for(i=0; i < x; ++i)
			this->n[i] = n[i];
		while(i < 4)
		{
			this->n[i] = 0;
			++i;
		}
	}

	void writeTo(std::ostream &out, bool with_x = false)
	{
		if (with_x)
			out << " " << x;

		for (uint8_t i = 0; i < x; ++i)
			out << " " << n[i];
		out << " " << c[0] << " " << c[1];
	}
};

class XF_FACE : public XF_MAIN_RECORD
{
private:
	XF_BC_TYPE m_bc;
	XF_FACE_TYPE m_face;
	std::vector<XF_CONNECTIVITY> m_connectivity;

public:
	XF_FACE(int zone, int first, int last, int bc, int face);

	~XF_FACE() = default;

	int bc_type() const
	{
		return m_bc;
	}

	int face_type() const
	{
		return m_face;
	}

	const XF_CONNECTIVITY &connectivity(size_t loc_idx) const
	{
		return m_connectivity[loc_idx];
	}

	XF_CONNECTIVITY &connectivity(size_t loc_idx)
	{
		return m_connectivity[loc_idx];
	}

	void repr(std::ostream &out);
};

class XF_ZONE :public XF_ENTRY
{
private:
	int m_zoneID;
	std::string m_zoneType, m_zoneName;
	int m_domainID;

public:
	XF_ZONE(int zone, const std::string &type, const std::string &name) :
		XF_ENTRY(XF_SECTION::ZONE),
		m_zoneID(zone),
		m_zoneType(type),
		m_zoneName(name),
		m_domainID(0)
	{
	}

	~XF_ZONE() = default;

	void repr(std::ostream &out);
};

class XF_MSH
{
private:
	std::vector<XF_ENTRY*> m_content; 
	size_t m_totalNodeNum, m_totalCellNum, m_totalFaceNum;
	bool m_is3D;

public:
	XF_MSH() {}

	~XF_MSH()
	{
		const size_t N = m_content.size();
		for (size_t i = 0; i < N; ++i)
			delete m_content[i];
	}

	int readFromFile(const std::string &src);

	int writeToFile(const std::string &dst);

private:
	void add_entry(XF_ENTRY *e)
	{
		m_content.push_back(e);
	}

	void eat(std::istream &in, char c)
	{
		char tmp;

		do { in >> tmp; } while (tmp != c);
	}
};

#endif
