#ifndef __GT_XF_H__
#define __GT_XF_H__

#include <cstddef>
#include <istream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "common.h"

/**
	If connectivity info concerning each node is required, uncomment the following
	line or define the marco somewhere else in the compiling environment.
 */
 // #define XF_EXTRACT_NODE_CONNECTIVITY

namespace GridTool
{
	namespace XF
	{
		using COMMON::Vector;
		using COMMON::Array1D;
		using COMMON::DIM;

		class SECTION
		{
		private:
			int m_identity;

		public:
			enum {
				COMMENT = 0,
				HEADER = 1,
				DIMENSION = 2,
				NODE = 10,
				CELL = 12,
				FACE = 13,
				EDGE = 11,
				ZONE = 39, ZONE_MESHING = 45
			};

			SECTION() = delete;
			SECTION(int id) : m_identity(id) {}
			SECTION(const SECTION &rhs) = default;
			virtual ~SECTION() = default;

			virtual void repr(std::ostream &out) = 0;

			int identity() const { return m_identity; }
		};

		class BC
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

			static bool isValidIdx(int x);
			static bool isValidStr(const std::string &x);
			static const std::string &idx2str(int x);
			static int str2idx(const std::string &x);
		};

		class STR : public SECTION
		{
		private:
			std::string m_msg;

		public:
			STR() = delete;
			STR(int id, const std::string &msg) : SECTION(id), m_msg(msg) {}
			STR(const STR &rhs) : SECTION(rhs.identity()), m_msg(rhs.m_msg) {}
			virtual ~STR() = default;

			const std::string &str() const { return m_msg; }
			std::string &str() { return m_msg; }

			void repr(std::ostream &out);
		};

		class COMMENT : public STR
		{
		public:
			COMMENT() = delete;
			COMMENT(const std::string &info) : STR(SECTION::COMMENT, info) {}
			COMMENT(const COMMENT &rhs) : STR(SECTION::COMMENT, rhs.str()) {}
			~COMMENT() = default;
		};

		class HEADER : public STR
		{
		public:
			HEADER() = delete;
			HEADER(const std::string &info) : STR(SECTION::HEADER, info) {}
			HEADER(const HEADER &rhs) : STR(SECTION::HEADER, rhs.str()) {}
			~HEADER() = default;
		};

		class DIMENSION :public SECTION, public DIM
		{
		public:
			DIMENSION() = delete;
			DIMENSION(int dim, bool id3d = true) : SECTION(SECTION::DIMENSION), DIM(dim, id3d) {}
			DIMENSION(const DIMENSION &rhs) : SECTION(SECTION::DIMENSION), DIM(rhs.dimension(), rhs.is3D()) {}
			~DIMENSION() = default;

			int ND() const { return dimension(); }

			void repr(std::ostream &out);
		};

		class RANGE : public SECTION
		{
		protected:
			size_t m_zone;
			size_t m_first, m_last;

		public:
			RANGE() = delete;
			RANGE(int id, size_t zone, size_t first, size_t last);
			RANGE(const RANGE &rhs);
			virtual ~RANGE() = default;

			size_t zone() const { return m_zone; }

			size_t first_index() const { return m_first; }

			size_t last_index() const { return m_last; }

			size_t num() const { return(last_index() - first_index() + 1); }
		};

		class NODE : public RANGE, public DIM, public std::vector<Vector>
		{
		private:
			int m_type;

		public:
			enum { VIRTUAL = 0, ANY = 1, BOUNDARY = 2 };

			static bool isValidTypeIdx(int x);
			static bool isValidTypeStr(const std::string &x);
			static const std::string &idx2str(int x);
			static int str2idx(const std::string &x);

			NODE() = delete;
			NODE(size_t zone, size_t first, size_t last, int tp, int ND);
			NODE(const NODE &rhs);
			~NODE() = default;

			bool is_virtual_node() const { return type() == NODE::VIRTUAL; }
			bool is_boundary_node() const { return type() == NODE::BOUNDARY; }
			bool is_internal_node() const { return type() == NODE::ANY; }

			int &type() { return m_type; }
			int type() const { return m_type; }

			int ND() const { return dimension(); }

			void repr(std::ostream &out);
		};

		class CELL : public RANGE, public std::vector<int>
		{
		private:
			int m_type;
			int m_elem;

		public:
			enum {
				DEAD = 0,
				FLUID = 1,
				SOLID = 17
			};

			static bool isValidTypeIdx(int x);
			static bool isValidTypeStr(const std::string &x);
			static const std::string &idx2str_type(int x);
			static int str2idx_type(const std::string &x);

			enum {
				MIXED = 0,
				TRIANGULAR = 1,
				TETRAHEDRAL = 2,
				QUADRILATERAL = 3,
				HEXAHEDRAL = 4,
				PYRAMID = 5,
				WEDGE = 6,
				POLYHEDRAL = 7
			};

			static bool isValidElemIdx(int x);
			static bool isValidElemStr(const std::string &x);
			static const std::string &idx2str_elem(int x);
			static int str2idx_elem(const std::string &x);

			CELL() = delete;
			CELL(size_t zone, size_t first, size_t last, int type, int elem_type);
			CELL(const CELL &rhs);
			~CELL() = default;

			// Type of cells within this section: DEAD cell, FLUID cell or SOLID cell.
			int type() const { return m_type; }
			int &type() { return m_type; }

			// General description of ALL cell elements within this section.
			int element_type() const { return  m_elem; }
			int &element_type() { return m_elem; }

			void repr(std::ostream &out);
		};

		struct CONNECTIVITY
		{
			int x; // Num of nodes.
			size_t n[4]; // At most 4 nodes within a single face, polygon faces are not supported currently.
			size_t c[2]; // Adjacent cells.

			CONNECTIVITY() : x(1), n{ 0, 0, 0, 0 }, c{ 0, 0 } {}
			CONNECTIVITY(const CONNECTIVITY &rhs) = default;
			~CONNECTIVITY() = default;

			size_t cl() const { return c[0]; }
			size_t cr() const { return c[1]; }

			size_t c0() const { return c[0]; }
			size_t c1() const { return c[1]; }

			void set(int x_, size_t *n_, size_t *c_)
			{
				if (x_ > 4)
					throw std::invalid_argument("Too many nodes within a face, polygon face are not supported currently.");
				if (x_ < 1)
					throw std::invalid_argument("Invalid num of nodes within a face.");

				x = x_;
				c[0] = c_[0];
				c[1] = c_[1];

				int i = 0;
				for (; i < x_; ++i)
					n[i] = n_[i];
				while (i < 4)
				{
					n[i] = 0;
					++i;
				}
			}

			// Index of its left-hand side node.
			size_t leftAdj(int loc_idx) const
			{
				if (loc_idx == 0)
					return n[x - 1];
				else
					return n[loc_idx - 1];
			}

			// Index of its right-hand side node.
			size_t rightAdj(int loc_idx) const
			{
				if (loc_idx == x - 1)
					return n[0];
				else
					return n[loc_idx + 1];
			}
		};

		class FACE : public RANGE, public std::vector<CONNECTIVITY>
		{
		private:
			int m_bc;
			int m_face;

		public:
			enum {
				MIXED = 0,
				LINEAR = 2,
				TRIANGULAR = 3,
				QUADRILATERAL = 4,
				POLYGONAL = 5
			};

			static bool isValidIdx(int x);
			static bool isValidStr(const std::string &x);
			static const std::string &idx2str(int x);
			static int str2idx(const std::string &x);

			FACE() = delete;
			FACE(size_t zone, size_t first, size_t last, int bc, int face);
			~FACE() = default;

			// The B.C. of this group of faces.
			int bc_type() const { return m_bc; }
			int &bc_type() { return m_bc; }

			// The shape of this group of faces.
			int face_type() const { return m_face; }
			int &face_type() { return m_face; }

			void repr(std::ostream &out);
		};

		class ZONE :public SECTION
		{
		private:
			size_t m_zoneID;
			std::string m_zoneType, m_zoneName;
			int m_domainID;

		public:
			ZONE() = delete;
			ZONE(int zone, const std::string &type, const std::string &name) : SECTION(SECTION::ZONE), m_zoneID(zone), m_zoneType(type), m_zoneName(name), m_domainID(0) {}
			~ZONE() = default;

			// The index of this zone, may be any non-consecutive positive integer.
			size_t zone() const { return m_zoneID; }

			// The B.C. string literal.
			const std::string &type() const { return m_zoneType; }

			// The name of this zone.
			const std::string &name() const { return m_zoneName; }

			// NOT used.
			int domain() const { return m_domainID; }

			void repr(std::ostream &out);
		};

		class MESH : public DIM
		{
		protected:
			/// Index of node, face, and cell starts from 1 and increase continuously.
			/// But zone is different.
			struct NODE_ELEM
			{
				Vector coordinate;
				bool atBdry;
#ifdef XF_EXTRACT_NODE_CONNECTIVITY
				Array1D<size_t> adjacentNode;
				Array1D<size_t> dependentFace;
				Array1D<size_t> dependentCell;
#endif // XF_EXTRACT_NODE_CONNECTIVITY
			};
			struct FACE_ELEM
			{
				int type;
				Vector center;
				double area;
				Array1D<size_t> includedNode;
				size_t leftCell, rightCell;
				bool atBdry;
				Vector n_LR; // Surface unit normal
				Vector n_RL;
			};
			struct CELL_ELEM
			{
				int type;
				Vector center;
				double volume;
				Array1D<size_t> includedFace;
				Array1D<size_t> includedNode;
				Array1D<size_t> adjacentCell; // The size is equal to that of "includedFace", set to 0 if the adjacent cell is boundary.
				Array1D<Vector> n;
				Array1D<Vector> S;
			};
			struct ZONE_ELEM
			{
				size_t ID; // May not start from 1, and may be given arbitrarily.
				std::string type;
				std::string name;
				RANGE *obj;
			};

		private:
			// Raw
			std::vector<SECTION*> m_content;
			size_t m_totalNodeNum;
			size_t m_totalCellNum;
			size_t m_totalFaceNum;

			// Derived
			Array1D<NODE_ELEM> m_node;
			Array1D<FACE_ELEM> m_face;
			Array1D<CELL_ELEM> m_cell;
			size_t m_totalZoneNum;
			std::map<size_t, size_t> m_zoneMapping;
			Array1D<ZONE_ELEM> m_zone;

		public:
			MESH() : DIM(3), m_totalNodeNum(0), m_totalCellNum(0), m_totalFaceNum(0), m_totalZoneNum(0) {} // 3D by default
			MESH(const std::string &inp, std::ostream &fout);
			MESH(const std::string &f_nmf, const std::string &f_p3d, std::ostream &fout);
			MESH(const MESH &rhs) = delete;
			~MESH() { clear_entry(); }

			/* IO */
			void readFromFile(const std::string &src, std::ostream &fout);
			void writeToFile(const std::string &dst) const;

			/* Num of elements */
			size_t numOfNode() const { return m_totalNodeNum; }
			size_t numOfFace() const { return m_totalFaceNum; }
			size_t numOfCell() const { return m_totalCellNum; }
			size_t numOfZone() const { return m_totalZoneNum; }

			/* 1-based access */
			const NODE_ELEM &node(int id) const { return m_node(id); }
			NODE_ELEM &node(int id) { return m_node(id); }

			const FACE_ELEM &face(int id) const { return m_face(id); }
			FACE_ELEM &face(int id) { return m_face(id); }

			const CELL_ELEM &cell(int id) const { return m_cell(id); }
			CELL_ELEM &cell(int id) { return m_cell(id); }

			/// If "isRealZoneID" is "true", then "id" is the real zone index,
			/// otherwise, "id" is the internal storage index.
			/// Whatever "isRealZoneID" is, "id" is always 1-based for consistency.
			const ZONE_ELEM &zone(int id, bool isRealZoneID = false) const { return  isRealZoneID ? m_zone.at(m_zoneMapping.at(id)) : m_zone(id); }
			ZONE_ELEM &zone(int id, bool isRealZoneID = false) { return  isRealZoneID ? m_zone.at(m_zoneMapping.at(id)) : m_zone(id); }

		private:
			void add_entry(SECTION *e)
			{
				m_content.push_back(e);
			}
			void clear_entry()
			{
				// Release previous contents
				for (auto ptr : m_content)
					if (ptr)
						delete ptr;

				// Clear container
				m_content.clear();
			}

			void raw2derived();
			void derived2raw();

			void cell_standardization(CELL_ELEM &c);
			void tet_standardization(CELL_ELEM &tet);
			void pyramid_standardization(CELL_ELEM &pyramid);
			void prism_standardization(CELL_ELEM &prism);
			void hex_standardization(CELL_ELEM &hex);
			void triangle_standardization(CELL_ELEM &tri);
			void quad_standardization(CELL_ELEM &quad);
		};
	}
}
#endif
