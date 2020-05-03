#ifndef TYDF_XF_H
#define TYDF_XF_H

#include <cstddef>
#include <istream>
#include <ostream>
#include <iostream>
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
#include "common.h"

namespace GridTool::XF
{
    using GridTool::COMMON::Scalar;
    using GridTool::COMMON::Vector;
    using GridTool::COMMON::DIM;
    using GridTool::COMMON::Array1D;
    using GridTool::COMMON::wrong_index;
    using GridTool::COMMON::wrong_string;

    class SECTION
    {
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

    private:
        int m_identity;

    public:
        SECTION() = delete;

        SECTION(int id);

        SECTION(const SECTION &rhs) = default;

        virtual ~SECTION() = default;

        virtual void repr(std::ostream &out) = 0;

        int identity() const;
    };

    class BC
    {
    private:
        struct invalid_bc_idx : public wrong_index
        {
            invalid_bc_idx(int bc) : wrong_index(bc, "is not a valid B.C. index") {}
        };

        struct invalid_bc_str : public wrong_string
        {
            invalid_bc_str(const std::string &bc) : wrong_string(bc, "is not a valid B.C. string") {}
        };

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

        STR(int id, const std::string &msg);

        STR(const STR &rhs);

        virtual ~STR() = default;

        const std::string &str() const;

        std::string &str();

        void repr(std::ostream &out);
    };

    class COMMENT : public STR
    {
    public:
        COMMENT() = delete;

        COMMENT(const std::string &info);

        COMMENT(const COMMENT &rhs);

        ~COMMENT() = default;
    };

    class HEADER : public STR
    {
    public:
        HEADER() = delete;

        HEADER(const std::string &info);

        HEADER(const HEADER &rhs);

        ~HEADER() = default;
    };

    class DIMENSION :public SECTION, public DIM
    {
    public:
        DIMENSION() = delete;

        DIMENSION(int dim, bool id3d = true);

        DIMENSION(const DIMENSION &rhs);

        ~DIMENSION() = default;

        int ND() const;

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

        size_t zone() const;

        size_t first_index() const;

        size_t last_index() const;

        size_t num() const;
    };

    class NODE : public RANGE, public DIM, public std::vector<Vector>
    {
    private:
        struct invalid_node_type_idx : public wrong_index
        {
            invalid_node_type_idx(int x) : wrong_index(x, "is not a valid NODE-TYPE index") {}
        };

        struct invalid_node_type_str : public wrong_string
        {
            invalid_node_type_str(const std::string &s) : wrong_string(s, "is not a valid NODE-TYPE string") {}
        };

    public:
        enum {
            VIRTUAL = 0,
            ANY = 1,
            BOUNDARY = 2
        };

        static bool isValidTypeIdx(int x);

        static bool isValidTypeStr(const std::string &x);

        static const std::string &idx2str(int x);

        static int str2idx(const std::string &x);

    private:
        int m_type;

    public:
        NODE() = delete;

        NODE(size_t zone, size_t first, size_t last, int tp, int ND);

        NODE(const NODE &rhs);

        ~NODE() = default;

        bool is_virtual_node() const;

        bool is_boundary_node() const;

        bool is_internal_node() const;

        int &type();

        int type() const;

        int ND() const;

        void repr(std::ostream &out);
    };

    class CELL : public RANGE, public std::vector<int>
    {
    private:
        struct invalid_cell_type_idx : public wrong_index
        {
            invalid_cell_type_idx(int x) : wrong_index(x, "is not a valid CELL-TYPE index") {}
        };

        struct invalid_cell_type_str : public wrong_string
        {
            invalid_cell_type_str(const std::string &s) : wrong_string(s, "is not a valid CELL-TYPE string") {}
        };

        struct invalid_elem_type_idx : public wrong_index
        {
            invalid_elem_type_idx(int x) : wrong_index(x, "is not a valid CELL-ELEM-TYPE index") {}
        };

        struct invalid_elem_type_str : public wrong_string
        {
            invalid_elem_type_str(const std::string &s) : wrong_string(s, "is not a valid CELL-ELEM-TYPE string") {}
        };

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

    private:
        int m_type;
        int m_elem;

    public:
        CELL() = delete;

        CELL(size_t zone, size_t first, size_t last, int type, int elem_type);

        CELL(const CELL &rhs);

        ~CELL() = default;

        /// Type of cells within this section: DEAD cell, FLUID cell or SOLID cell.
        int type() const;

        int &type();

        /// General description of ALL cell elements within this section.
        int element_type() const;

        int &element_type();

        void repr(std::ostream &out);
    };

    struct CONNECTIVITY
    {
        /// Num of nodes.
        int x;

        /// At most 4 nodes within a single face, 
        /// polygon faces are not supported currently.
        size_t n[4];

        /// Adjacent cells.
        size_t c[2];

        CONNECTIVITY() : x(1), n{ 0, 0, 0, 0 }, c{ 0, 0 } {}
        CONNECTIVITY(const CONNECTIVITY &rhs) = default;
        ~CONNECTIVITY() = default;

        size_t cl() const { return c[0]; }
        size_t cr() const { return c[1]; }

        size_t c0() const { return c[0]; }
        size_t c1() const { return c[1]; }

        void set(int x_, size_t *n_, size_t *c_);

        /// Index of adjacent node.
        size_t leftAdj(int loc_idx) const;
        size_t rightAdj(int loc_idx) const;
    };

    class FACE : public RANGE, public std::vector<CONNECTIVITY>
    {
    private:
        struct invalid_face_type_idx : public wrong_index
        {
            invalid_face_type_idx(int x) : wrong_index(x, "is not a valid FACE-TYPE index") {}
        };

        struct invalid_face_type_str : public wrong_string
        {
            invalid_face_type_str(const std::string &s) : wrong_string(s, "is not a valid FACE-TYPE string") {}
        };

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

    private:
        int m_bc;
        int m_face;

    public:
        FACE() = delete;

        FACE(size_t zone, size_t first, size_t last, int bc, int face);

        FACE(const FACE &rhs);

        ~FACE() = default;

        /// B.C. of this group of faces.
        int bc_type() const;

        int &bc_type();

        /// Shape of this group of faces.
        int face_type() const;

        int &face_type();

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

        ZONE(int zone, const std::string &bc_type, const std::string &name);

        ZONE(const ZONE &rhs) = default;

        ~ZONE() = default;

        /// Index of this zone, may be any 
        /// non-consecutive positive integer.
        size_t zone() const;

        size_t &zone();

        /// B.C. string literal.
        const std::string &type() const;

        std::string &type();

        /// Name of this zone.
        const std::string &name() const;

        std::string &name();

        /// Domain ID, NOT used.
        int domain() const;

        int &domain();

        void repr(std::ostream &out);
    };

    class MESH : public DIM
    {
    protected:
        /// Index of node, face, and cell starts from 1 
        /// and increase continuously. But zone is different.
        struct NODE_ELEM
        {
            Vector coordinate;
            bool atBdry;
            Array1D<size_t> adjacentNode;
            Array1D<size_t> dependentFace;
            Array1D<size_t> dependentCell;
        };
        struct FACE_ELEM
        {
            int type;
            Vector center;
            double area;
            Array1D<size_t> includedNode;
            size_t leftCell, rightCell;
            bool atBdry;

            /// Surface unit normal.
            Vector n_LR, n_RL;
        };
        struct CELL_ELEM
        {
            int type;
            Vector center;
            double volume;
            Array1D<size_t> includedFace;
            Array1D<size_t> includedNode;

            /// The size is equal to that of "includedFace", set to 0 if adjacent cell is boundary.
            Array1D<size_t> adjacentCell;

            Array1D<Vector> n;
            Array1D<Vector> S;
        };
        struct ZONE_ELEM
        {
            /// May not start from 1, and may be given arbitrarily.
            size_t ID;

            std::string type;
            std::string name;
            RANGE *obj;
        };

    private:
        /// Raw
        std::vector<SECTION*> m_content;
        size_t m_totalNodeNum;
        size_t m_totalCellNum;
        size_t m_totalFaceNum;

        /// Derived
        Array1D<NODE_ELEM> m_node;
        Array1D<FACE_ELEM> m_face;
        Array1D<CELL_ELEM> m_cell;
        size_t m_totalZoneNum;
        std::map<size_t, size_t> m_zoneMapping;
        Array1D<ZONE_ELEM> m_zone;

    public:
        MESH();
        MESH(const std::string &inp, std::ostream &fout = std::cout);
        MESH(const std::string &f_nmf, const std::string &f_p3d, std::ostream &fout = std::cout);
        MESH(const MESH &rhs) = delete;
        ~MESH() { clear_entry(); }

        /// IO
        void readFromFile(const std::string &src, std::ostream &fout);
        void writeToFile(const std::string &dst) const;

        /// Num of elements
        size_t numOfNode() const { return m_totalNodeNum; }
        size_t numOfFace() const { return m_totalFaceNum; }
        size_t numOfCell() const { return m_totalCellNum; }
        size_t numOfZone() const { return m_totalZoneNum; }

        /// 1-based access
        const NODE_ELEM &node(size_t id) const { return m_node(id); }
        NODE_ELEM &node(size_t id) { return m_node(id); }

        const FACE_ELEM &face(size_t id) const { return m_face(id); }
        FACE_ELEM &face(size_t id) { return m_face(id); }

        const CELL_ELEM &cell(size_t id) const { return m_cell(id); }
        CELL_ELEM &cell(size_t id) { return m_cell(id); }

        /// If "isRealZoneID" is "true", then "id" is the real zone index,
        /// otherwise, "id" is the internal storage index.
        /// Whatever "isRealZoneID" is, "id" is always 1-based for consistency.
        const ZONE_ELEM &zone(size_t id, bool isRealZoneID = false) const { return  isRealZoneID ? m_zone.at(m_zoneMapping.at(id)) : m_zone(id); }
        ZONE_ELEM &zone(size_t id, bool isRealZoneID = false) { return  isRealZoneID ? m_zone.at(m_zoneMapping.at(id)) : m_zone(id); }

    private:
        void add_entry(SECTION *e);
        void clear_entry();

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
#endif
