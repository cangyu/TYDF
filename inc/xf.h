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

		class NODE : public RANGE, public DIM
		{
		private:
			int m_type;
			Array1D<Vector> m_node;

		public:
			enum { VIRTUAL = 0, ANY = 1, BOUNDARY = 2 };

			static bool isValidTypeIdx(int x);
			static bool isValidTypeStr(const std::string &x);
			static const std::string &idx2str(int x);
			static int str2idx(const std::string &x);

			NODE() = delete;
			NODE(size_t zone, size_t first, size_t last, int tp, int ND) :
				RANGE(SECTION::NODE, zone, first, last),
				DIM(ND),
				m_type(tp),
				m_node(num())
			{
				if (!isValidTypeIdx(type()))
					throw std::runtime_error("Invalid description of node type!");
			}
			NODE(const NODE &rhs) :
				RANGE(SECTION::NODE, rhs.zone(), rhs.first_index(), rhs.last_index()),
				DIM(rhs.ND(), rhs.is3D()),
				m_type(rhs.type()),
				m_node(rhs.num())
			{
				if (!isValidTypeIdx(type()))
					throw std::runtime_error("Invalid description of node type!");
			}
			~NODE() = default;

			bool is_virtual_node() const { return type() == NODE::VIRTUAL; }
			bool is_boundary_node() const { return type() == NODE::BOUNDARY; }
			bool is_internal_node() const { return type() == NODE::ANY; }

			int &type() { return m_type; }
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

			void repr(std::ostream &out);
		};

		class CELL : public RANGE
		{
		private:
			int m_type;
			int m_elem;
			std::vector<int> m_mixedElemDesc; // Only effective when 'm_elem == MIXED', empty otherwise.

		public:
			enum { DEAD = 0, FLUID = 1, SOLID = 17 };
			enum { MIXED = 0, TRIANGULAR = 1, TETRAHEDRAL = 2, QUADRILATERAL = 3, HEXAHEDRAL = 4, PYRAMID = 5, WEDGE = 6, POLYHEDRAL = 7 };

			static bool isValidTypeIdx(int x);
			static bool isValidTypeStr(const std::string &x);
			static const std::string &idx2str_type(int x);
			static int str2idx_type(const std::string &x);
			static bool isValidElemIdx(int x);
			static bool isValidElemStr(const std::string &x);
			static const std::string &idx2str_elem(int x);
			static int str2idx_elem(const std::string &x);

			CELL() = delete;
			CELL(size_t zone, size_t first, size_t last, int type, int elem_type) :
				RANGE(SECTION::CELL, zone, first, last)
			{
				// Check cell type before assign
				if (!isValidTypeIdx(type))
					throw std::runtime_error("Invalid cell type: " + std::to_string(type));
				else
					m_type = type;

				// Check cell elem before assign
				if (!isValidElemIdx(elem_type))
					throw std::runtime_error("Invalid cell element type: " + std::to_string(elem_type));
				else
					m_elem = elem_type;

				// Special treatment for mixed cell
				if (elem_type == CELL::MIXED)
				{
					m_mixedElemDesc.resize(num());
					std::fill(m_mixedElemDesc.begin(), m_mixedElemDesc.end(), CELL::MIXED);
				}
			}
			~CELL() = default;

			// Type of cells within this section: DEAD cell, FLUOID cell or SOLID cell.
			int type() const { return m_type; }
			int &type() { return m_type; }

			// General description of ALL cell elements within this section.
			int element_type() const { return m_elem; }
			int &element_type() { return m_elem; }

			// Element-Type of each cell within this section.
			// Typically used when "element_type() == MIXED".
			int elem(size_t loc_idx) const
			{
				int et = element_type();
				if (et == CELL::MIXED)
					return m_mixedElemDesc[loc_idx];
				else
					return et;
			}
			int &elem(size_t loc_idx)
			{
				int &et = element_type();
				if (et == CELL::MIXED)
					return m_mixedElemDesc[loc_idx];
				else
					return et;
			}

			void repr(std::ostream &out);
		};

		class CONNECTIVITY
		{
		public:
			int x; // Num of nodes.
			size_t n[4]; // At most 4 nodes within a single face, polygon faces are not supported currently.
			size_t c[2]; // Adjacent cells.

		public:
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

		class FACE : public RANGE
		{
		private:
			int m_bc;
			int m_face;
			std::vector<CONNECTIVITY> m_connectivity;

		public:
			enum { MIXED = 0, LINEAR = 2, TRIANGULAR = 3, QUADRILATERAL = 4, POLYGONAL = 5 };

			static bool isValidIdx(int x);
			static bool isValidStr(const std::string &x);
			static const std::string &idx2str(int x);
			static int str2idx(const std::string &x);

			FACE() = delete;
			FACE(size_t zone, size_t first, size_t last, int bc, int face) :
				RANGE(SECTION::FACE, zone, first, last)
			{
				// Check B.C. before assign
				if (!BC::isValidIdx(bc))
					throw std::runtime_error("Invalid B.C. type: " + std::to_string(bc));
				else
					m_bc = bc;

				// Check face type before assign
				if (!isValidIdx(face))
					throw std::runtime_error("Invalid face type: " + std::to_string(face));
				else if (face == FACE::POLYGONAL)
					throw std::runtime_error("Polygonal face is not supported currently.");
				else
					m_face = face;

				// Resize local storage
				m_connectivity.resize(num());
			}
			~FACE() = default;

			// The B.C. of this group of faces.
			int bc_type() const { return m_bc; }

			// The shape of this group of faces.
			int face_type() const { return m_face; }

			// 0-based local indexing
			const CONNECTIVITY &connectivity(size_t loc_idx) const { return m_connectivity.at(loc_idx); }
			CONNECTIVITY &connectivity(size_t loc_idx) { return m_connectivity.at(loc_idx); }

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
		private:
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
			const ZONE_ELEM &zone(int id, bool isRealZoneID = false) const
			{
				return  isRealZoneID ? m_zone.at(m_zoneMapping.at(id)) : m_zone(id);
			}
			ZONE_ELEM &zone(int id, bool isRealZoneID = false)
			{
				return  isRealZoneID ? m_zone.at(m_zoneMapping.at(id)) : m_zone(id);
			}

		private:
			static double dot_product(const Vector &na, const Vector &nb)
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
