#include "../inc/xf.h"

/// Convert a boundary condition string literal to unified form within the scope of this code.
/// Outcome will be composed of LOWER case letters and '-' only!
static void formalize_inplace(std::string &s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    for (auto &e : s)
        if (e == '_')
            e = '-';
}

static std::string formalize(const std::string &s)
{
    std::string ret(s);
    formalize_inplace(ret);
    return ret;
}

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

namespace GridTool::XF
{
    SECTION::SECTION(int id) :
        m_identity(id)
    {
        /// Empty body.
    }

    int SECTION::identity() const
    {
        return m_identity;
    }

    struct BC::invalid_bc_idx : public wrong_index
    {
        invalid_bc_idx(int x) :
            wrong_index(x, "is not a valid B.C. index")
        {
            /// Empty body.
        }
    };

    struct BC::invalid_bc_str : public wrong_string
    {
        invalid_bc_str(const std::string &s) :
            wrong_string(s, "is not a valid B.C. string")
        {
            /// Empty body.
        }
    };

    bool BC::isValidIdx(int x)
    {
        static const std::set<int> candidate_set{
            INTERIOR,
            WALL,
            PRESSURE_INLET,
            PRESSURE_OUTLET,
            SYMMETRY,
            PERIODIC_SHADOW,
            PRESSURE_FAR_FIELD,
            VELOCITY_INLET,
            PERIODIC,
            FAN,
            MASS_FLOW_INLET,
            INTERFACE,
            PARENT,
            OUTFLOW,
            AXIS
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool BC::isValidStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "interior",
            "wall",
            "pressure-inlet", "inlet-vent", "intake-fan",
            "pressure-outlet", "exhaust-fan", "outlet-vent",
            "symmetry",
            "periodic-shadow",
            "pressure-far-field",
            "velocity-inlet",
            "periodic",
            "fan", "porous-jump", "radiator",
            "mass-flow-inlet",
            "interface",
            "parent",
            "outflow",
            "axis"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &BC::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(INTERIOR, "interior"),
            std::pair<int, std::string>(WALL, "wall"),
            std::pair<int, std::string>(PRESSURE_INLET, "pressure-inlet"),
            std::pair<int, std::string>(PRESSURE_OUTLET, "pressure-outlet"),
            std::pair<int, std::string>(SYMMETRY, "symmetry"),
            std::pair<int, std::string>(PERIODIC_SHADOW, "periodic-shadow"),
            std::pair<int, std::string>(PRESSURE_FAR_FIELD, "pressure-far-field"),
            std::pair<int, std::string>(VELOCITY_INLET, "velocity-inlet"),
            std::pair<int, std::string>(PERIODIC, "periodic"),
            std::pair<int, std::string>(FAN, "fan"),
            std::pair<int, std::string>(MASS_FLOW_INLET, "mass-flow-inlet"),
            std::pair<int, std::string>(INTERFACE, "interface"),
            std::pair<int, std::string>(PARENT, "parent"),
            std::pair<int, std::string>(OUTFLOW, "outflow"),
            std::pair<int, std::string>(AXIS, "axis")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_bc_idx(x);
        else
            return it->second;
    }

    int BC::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("interior", INTERIOR),
            std::pair<std::string, int>("wall", WALL),
            std::pair<std::string, int>("pressure-inlet", PRESSURE_INLET),
            std::pair<std::string, int>("inlet-vent", INLET_VENT),
            std::pair<std::string, int>("intake-fan", INTAKE_FAN),
            std::pair<std::string, int>("pressure-outlet", PRESSURE_OUTLET),
            std::pair<std::string, int>("exhaust-fan", EXHAUST_FAN),
            std::pair<std::string, int>("outlet-vent", OUTLET_VENT),
            std::pair<std::string, int>("symmetry", SYMMETRY),
            std::pair<std::string, int>("periodic-shadow", PERIODIC_SHADOW),
            std::pair<std::string, int>("pressure-far-field", PRESSURE_FAR_FIELD),
            std::pair<std::string, int>("velocity-inlet", VELOCITY_INLET),
            std::pair<std::string, int>("periodic", PERIODIC),
            std::pair<std::string, int>("fan", FAN),
            std::pair<std::string, int>("porous-jump", POROUS_JUMP),
            std::pair<std::string, int>("radiator", RADIATOR),
            std::pair<std::string, int>("mass-flow-inlet", MASS_FLOW_INLET),
            std::pair<std::string, int>("interface", INTERFACE),
            std::pair<std::string, int>("parent", PARENT),
            std::pair<std::string, int>("outflow", OUTFLOW),
            std::pair<std::string, int>("axis", AXIS)
        };

        auto it = mapping_set.find(formalize(x));
        if (it == mapping_set.end())
            throw invalid_bc_str(x);
        else
            return it->second;
    }

    STR::STR(int id, const std::string &msg) :
        SECTION(id),
        m_msg(msg)
    {
        /// Empty body.
    }

    STR::STR(const STR &rhs) :
        SECTION(rhs.identity()),
        m_msg(rhs.m_msg)
    {
        /// Empty body.
    }

    const std::string &STR::str() const
    {
        return m_msg;
    }

    std::string &STR::str()
    {
        return m_msg;
    }

    void STR::repr(std::ostream &out)
    {
        out << "(" << std::dec << identity() << " \"" << str() << "\")" << std::endl;
    }

    COMMENT::COMMENT(const std::string &info) :
        STR(SECTION::COMMENT, info)
    {
        /// Empty body.
    }

    COMMENT::COMMENT(const COMMENT &rhs) :
        STR(SECTION::COMMENT, rhs.str())
    {
        /// Empty body.
    }

    HEADER::HEADER(const std::string &info) :
        STR(SECTION::HEADER, info)
    {
        /// Empty body.
    }

    HEADER::HEADER(const HEADER &rhs) :
        STR(SECTION::HEADER, rhs.str())
    {
        /// Empty body.
    }

    DIMENSION::DIMENSION(int dim, bool id3d) :
        SECTION(SECTION::DIMENSION),
        DIM(dim, id3d)
    {
        /// Empty body.
    }

    DIMENSION::DIMENSION(const DIMENSION &rhs) :
        SECTION(SECTION::DIMENSION),
        DIM(rhs.dimension(), rhs.is3D())
    {
        /// Empty body.
    }

    int DIMENSION::ND() const
    {
        return dimension();
    }

    void DIMENSION::repr(std::ostream &out)
    {
        out << "(" << std::dec << identity() << " " << ND() << ")" << std::endl;
    }

    RANGE::RANGE(int id, size_t zone, size_t first, size_t last) :
        SECTION(id),
        m_zone(zone),
        m_first(first),
        m_last(last)
    {
        if (first > last)
            throw std::invalid_argument("Invalid range in constructor.");
    }

    RANGE::RANGE(const RANGE &rhs) :
        SECTION(rhs.identity()),
        m_zone(rhs.zone()),
        m_first(rhs.first_index()),
        m_last(rhs.last_index())
    {
        if (first_index() > last_index())
            throw std::invalid_argument("Invalid range in copy-constructor.");
    }

    size_t RANGE::zone() const
    {
        return m_zone;
    }

    size_t RANGE::first_index() const
    {
        return m_first;
    }

    size_t RANGE::last_index() const
    {
        return m_last;
    }

    size_t RANGE::num() const
    {
        const size_t ret = last_index() - first_index() + 1;
        return ret;
    }

    struct NODE::invalid_node_type_idx : public wrong_index
    {
        invalid_node_type_idx(int x) :
            wrong_index(x, "is not a valid NODE-TYPE index")
        {
            /// Empty body.
        }
    };

    struct NODE::invalid_node_type_str : public wrong_string
    {
        invalid_node_type_str(const std::string &s) :
            wrong_string(s, "is not a valid NODE-TYPE string")
        {
            /// Empty body.
        }
    };

    bool NODE::isValidTypeIdx(int x)
    {
        return x == VIRTUAL || x == ANY || x == BOUNDARY;
    }

    bool NODE::isValidTypeStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "virtual",
            "any",
            "boundary"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &NODE::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(VIRTUAL, "virtual"),
            std::pair<int, std::string>(ANY, "any"),
            std::pair<int, std::string>(BOUNDARY, "boundary")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_node_type_idx(x);
        else
            return it->second;
    }

    int NODE::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("virtual", VIRTUAL),
            std::pair<std::string, int>("any", ANY),
            std::pair<std::string, int>("boundary", BOUNDARY),
        };

        auto it = mapping_set.find(formalize(x));
        if (it == mapping_set.end())
            throw invalid_node_type_str(x);
        else
            return it->second;
    }

    NODE::NODE(size_t zone, size_t first, size_t last, int tp, int ND) :
        RANGE(SECTION::NODE, zone, first, last),
        DIM(ND, ND == 3),
        std::vector<Vector>(num()),
        m_type(tp)
    {
        if (!isValidTypeIdx(type()))
            throw std::invalid_argument("Invalid description of node type in constructor.");
    }

    NODE::NODE(const NODE &rhs) :
        RANGE(SECTION::NODE, rhs.zone(), rhs.first_index(), rhs.last_index()),
        DIM(rhs.ND(), rhs.is3D()),
        std::vector<Vector>(rhs.begin(), rhs.end()),
        m_type(rhs.type())
    {
        if (!isValidTypeIdx(type()))
            throw std::invalid_argument("Invalid description of node type in copy-constructor.");
    }

    bool NODE::is_virtual_node() const
    {
        return type() == NODE::VIRTUAL;
    }

    bool NODE::is_boundary_node() const
    {
        return type() == NODE::BOUNDARY;
    }

    bool NODE::is_internal_node() const
    {
        return type() == NODE::ANY;
    }

    int &NODE::type()
    {
        return m_type;
    }

    int NODE::type() const
    {
        return m_type;
    }

    int NODE::ND() const
    {
        return dimension();
    }

    void NODE::repr(std::ostream &out)
    {
        out << "(" << std::dec << identity();
        out << " (" << std::hex << zone() << " " << first_index() << " " << last_index() << " ";
        out << std::dec << type() << " " << ND() << ")(" << std::endl;

        out.precision(12);
        const size_t N = num();
        for (size_t i = 0; i < N; ++i)
        {
            const auto &node = at(i);
            for (int k = 0; k < m_dim; ++k)
                out << " " << node.at(k);
            out << std::endl;
        }
        out << "))" << std::endl;
    }

    struct CELL::invalid_cell_type_idx : public wrong_index
    {
        invalid_cell_type_idx(int x) :
            wrong_index(x, "is not a valid CELL-TYPE index")
        {
            /// Empty body.
        }
    };

    struct CELL::invalid_cell_type_str : public wrong_string
    {
        invalid_cell_type_str(const std::string &s) :
            wrong_string(s, "is not a valid CELL-TYPE string")
        {
            /// Empty body.
        }
    };

    struct CELL::invalid_elem_type_idx : public wrong_index
    {
        invalid_elem_type_idx(int x) :
            wrong_index(x, "is not a valid CELL-ELEM-TYPE index")
        {
            /// Empty body.
        }
    };

    struct CELL::invalid_elem_type_str : public wrong_string
    {
        invalid_elem_type_str(const std::string &s) :
            wrong_string(s, "is not a valid CELL-ELEM-TYPE string")
        {
            /// Empty body.
        }
    };

    bool CELL::isValidTypeIdx(int x)
    {
        static const std::set<int> candidate_set{
            DEAD,
            FLUID,
            SOLID
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool CELL::isValidTypeStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "dead",
            "fluid",
            "solid"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &CELL::idx2str_type(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(FLUID, "fluid"),
            std::pair<int, std::string>(SOLID, "solid"),
            std::pair<int, std::string>(DEAD, "dead")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_cell_type_idx(x);
        else
            return it->second;
    }

    int CELL::str2idx_type(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("fluid", FLUID),
            std::pair<std::string, int>("solid", SOLID),
            std::pair<std::string, int>("dead", DEAD),
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_cell_type_str(x);
        else
            return it->second;
    }

    bool CELL::isValidElemIdx(int x)
    {
        static const std::set<int> candidate_set{
            MIXED,
            TRIANGULAR,
            TETRAHEDRAL,
            QUADRILATERAL,
            HEXAHEDRAL,
            PYRAMID,
            WEDGE,
            POLYHEDRAL
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool CELL::isValidElemStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "mixed",
            "triangular",
            "tetrahedral",
            "quadrilateral",
            "hexahedral",
            "pyramid",
            "wedge",
            "prism",
            "polyhedral"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &CELL::idx2str_elem(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(MIXED, "mixed"),
            std::pair<int, std::string>(TRIANGULAR, "triangular"),
            std::pair<int, std::string>(TETRAHEDRAL, "tetrahedral"),
            std::pair<int, std::string>(QUADRILATERAL, "quadrilateral"),
            std::pair<int, std::string>(HEXAHEDRAL, "hexahedral"),
            std::pair<int, std::string>(PYRAMID, "pyramid"),
            std::pair<int, std::string>(WEDGE, "wedge"),
            std::pair<int, std::string>(POLYHEDRAL, "polyhedral")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_elem_type_idx(x);
        else
            return it->second;
    }

    int CELL::str2idx_elem(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("mixed", MIXED),
            std::pair<std::string, int>("triangular", TRIANGULAR),
            std::pair<std::string, int>("tri", TRIANGULAR),
            std::pair<std::string, int>("tetrahedral", TETRAHEDRAL),
            std::pair<std::string, int>("tet", TETRAHEDRAL),
            std::pair<std::string, int>("quadrilateral", QUADRILATERAL),
            std::pair<std::string, int>("quad", QUADRILATERAL),
            std::pair<std::string, int>("hexahedral", HEXAHEDRAL),
            std::pair<std::string, int>("hex", HEXAHEDRAL),
            std::pair<std::string, int>("pyramid", PYRAMID),
            std::pair<std::string, int>("wedge", WEDGE),
            std::pair<std::string, int>("prism", WEDGE),
            std::pair<std::string, int>("polyhedral", POLYHEDRAL)
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_elem_type_str(x);
        else
            return it->second;
    }

    CELL::CELL(size_t zone, size_t first, size_t last, int type, int elem_type) :
        RANGE(SECTION::CELL, zone, first, last),
        std::vector<int>(num(), elem_type),
        m_type(type),
        m_elem(elem_type)
    {
        if (!isValidTypeIdx(type))
            throw invalid_cell_type_idx(type);

        if (!isValidElemIdx(elem_type))
            throw invalid_elem_type_idx(elem_type);
    }

    CELL::CELL(const CELL &rhs) :
        RANGE(SECTION::CELL, rhs.zone(), rhs.first_index(), rhs.last_index()),
        std::vector<int>(rhs.begin(), rhs.end()),
        m_type(rhs.type()),
        m_elem(rhs.element_type())
    {
        if (num() != rhs.num())
            throw std::runtime_error("Default copy operation is inconsistent.");
    }

    int CELL::type() const
    {
        return m_type;
    }

    int &CELL::type()
    {
        return m_type;
    }

    int CELL::element_type() const
    {
        return  m_elem;
    }

    int &CELL::element_type()
    {
        return m_elem;
    }

    void CELL::repr(std::ostream &out)
    {
        static const size_t NumPerLine = 40;

        out << "(" << std::dec << identity() << " (";
        out << std::hex;
        out << zone() << " " << first_index() << " " << last_index() << " ";
        out << m_type << " " << m_elem << ")";

        if (m_elem != CELL::MIXED)
            out << ")" << std::endl;
        else
        {
            out << "(";
            const size_t N = num();
            for (size_t i = 0; i < N; ++i)
            {
                if (i % NumPerLine == 0)
                    out << std::endl;
                out << " " << at(i);
            }
            out << std::endl << "))" << std::endl;
        }
    }

    CONNECTIVITY::CONNECTIVITY() : x(1), n{ 0, 0, 0, 0 }, c{ 0, 0 } {}

    size_t CONNECTIVITY::cl() const
    {
        return c1();
    }

    size_t CONNECTIVITY::cr() const
    {
        return c0();
    }

    size_t CONNECTIVITY::c0() const
    {
        return c[0];
    }

    size_t CONNECTIVITY::c1() const
    {
        return c[1];
    }

    void CONNECTIVITY::set(int x_, const size_t *n_, const size_t *c_)
    {
        if (x_ > 4 || x_ < 1)
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

    size_t CONNECTIVITY::leftAdj(int loc_idx) const
    {
        if (loc_idx == 0)
            return n[x - 1];
        else
            return n[loc_idx - 1];
    }

    size_t CONNECTIVITY::rightAdj(int loc_idx) const
    {
        if (loc_idx == x - 1)
            return n[0];
        else
            return n[loc_idx + 1];
    }

    struct FACE::polygon_not_supported : public std::invalid_argument
    {
        polygon_not_supported() :
            std::invalid_argument("Polygonal faces are NOT supported currently!")
        {
            /// Empty body.
        }
    };

    struct FACE::invalid_face_type_idx : public wrong_index
    {
        invalid_face_type_idx(int x) :
            wrong_index(x, "is not a valid FACE-TYPE index")
        {
            /// Empty body.
        }
    };

    struct FACE::invalid_face_type_str : public wrong_string
    {
        invalid_face_type_str(const std::string &s) :
            wrong_string(s, "is not a valid FACE-TYPE string")
        {
            /// Empty body.
        }
    };

    bool FACE::isValidIdx(int x)
    {
        static const std::set<int> candidate_set{
            MIXED,
            LINEAR,
            TRIANGULAR,
            QUADRILATERAL,
            POLYGONAL
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    bool FACE::isValidStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "mixed",
            "linear",
            "triangular",
            "quadrilateral",
            "polygonal"
        };

        return candidate_set.find(x) != candidate_set.end();
    }

    const std::string &FACE::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(MIXED, "mixed"),
            std::pair<int, std::string>(LINEAR, "linear"),
            std::pair<int, std::string>(TRIANGULAR, "triangular"),
            std::pair<int, std::string>(QUADRILATERAL, "quadrilateral"),
            std::pair<int, std::string>(POLYGONAL, "polygonal")
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_face_type_idx(x);
        else
            return it->second;
    }

    int FACE::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("mixed", MIXED),
            std::pair<std::string, int>("linear", LINEAR),
            std::pair<std::string, int>("line", LINEAR),
            std::pair<std::string, int>("triangular", TRIANGULAR),
            std::pair<std::string, int>("tri", TRIANGULAR),
            std::pair<std::string, int>("quadrilateral", QUADRILATERAL),
            std::pair<std::string, int>("quad", QUADRILATERAL),
            std::pair<std::string, int>("polygonal", POLYGONAL)
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_face_type_str(x);
        else
            return it->second;
    }

    FACE::FACE(size_t zone, size_t first, size_t last, int bc, int face) :
        RANGE(SECTION::FACE, zone, first, last),
        std::vector<CONNECTIVITY>(num()),
        m_bc(bc),
        m_face(face)
    {
        if (!BC::isValidIdx(bc))
            throw BC::invalid_bc_idx(bc);

        if (!isValidIdx(face))
            throw invalid_face_type_idx(face);
        if (face == POLYGONAL)
            throw FACE::polygon_not_supported();
    }

    FACE::FACE(const FACE &rhs) :
        RANGE(SECTION::FACE, rhs.zone(), rhs.first_index(), rhs.last_index()),
        std::vector<CONNECTIVITY>(rhs.begin(), rhs.end()),
        m_bc(rhs.bc_type()),
        m_face(rhs.face_type())
    {
        if (!BC::isValidIdx(bc_type()))
            throw std::runtime_error("Invalid B.C. not detected in previous construction.");

        if (!isValidIdx(face_type()))
            throw std::runtime_error("Invalid FACE-TYPE not detected in previous construction.");
        if (face_type() == POLYGONAL)
            throw FACE::polygon_not_supported();
    }

    int FACE::bc_type() const
    {
        return m_bc;
    }

    int &FACE::bc_type()
    {
        return m_bc;
    }

    int FACE::face_type() const
    {
        return m_face;
    }

    int &FACE::face_type()
    {
        return m_face;
    }

    void FACE::repr(std::ostream &out)
    {
        out << "(" << std::dec << identity() << " (";
        out << std::hex;
        out << zone() << " " << first_index() << " " << last_index() << " ";
        out << bc_type() << " " << face_type() << ")(" << std::endl;

        const size_t N = num();
        if (m_face == MIXED)
        {
            for (size_t i = 0; i < N; ++i)
            {
                const auto &loc_cnect = at(i);
                out << " " << loc_cnect.x;
                for (int j = 0; j < loc_cnect.x; ++j)
                    out << " " << loc_cnect.n[j];
                out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
            }
        }
        else
        {
            for (size_t i = 0; i < N; ++i)
            {
                const auto &loc_cnect = at(i);
                for (int j = 0; j < loc_cnect.x; ++j)
                    out << " " << loc_cnect.n[j];
                out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
            }
        }

        out << "))" << std::endl;
    }

    struct ZONE::invalid_zone_type_idx : public wrong_index
    {
        invalid_zone_type_idx(int x) :
            wrong_index(x, "is not a valid ZONE-TYPE index")
        {
            /// Empty body.
        }
    };

    struct ZONE::invalid_zone_type_str : public wrong_string
    {
        invalid_zone_type_str(const std::string &s) :
            wrong_string(s, "is not a valid specification of ZONE-TYPE")
        {
            /// Empty body.
        }
    };

    bool ZONE::isValidIdx(int x)
    {
        const bool ret = ZONE::DEGASSING <= x && x <= ZONE::WRAPPER;
        return ret;
    }

    bool ZONE::isValidStr(const std::string &x)
    {
        static const std::set<std::string> candidate_set{
            "degassing",
            "exhaust-fan",
            "fan",
            "fluid",
            "geometry",
            "inlet-vent",
            "intake-fan",
            "interface",
            "interior",
            "internal",
            "mass-flow-inlet",
            "outflow",
            "outlet-vent",
            "parent-face",
            "porous-jump",
            "pressure-far-field",
            "pressure-inlet",
            "pressure-outlet",
            "radiator",
            "solid",
            "symmetry",
            "velocity-inlet",
            "wall",
            "wrapper"
        };

        return candidate_set.find(formalize(x)) != candidate_set.end();
    }

    const std::string &ZONE::idx2str(int x)
    {
        static const std::map<int, std::string> mapping_set{
            std::pair<int, std::string>(DEGASSING, "degassing"),
            std::pair<int, std::string>(EXHAUST_FAN, "exhaust-fan"),
            std::pair<int, std::string>(FAN, "fan"),
            std::pair<int, std::string>(FLUID, "fluid"),
            std::pair<int, std::string>(GEOMETRY, "geometry"),
            std::pair<int, std::string>(INLET_VENT, "inlet_vent"),
            std::pair<int, std::string>(INTAKE_FAN, "intake_fan"),
            std::pair<int, std::string>(INTERFACE, "interface"),
            std::pair<int, std::string>(INTERIOR, "interior"),
            std::pair<int, std::string>(INTERNAL, "internal"),
            std::pair<int, std::string>(MASS_FLOW_INLET, "mass_flow_inlet"),
            std::pair<int, std::string>(OUTFLOW, "outflow"),
            std::pair<int, std::string>(OUTLET_VENT, "outlet_vent"),
            std::pair<int, std::string>(PARENT_FACE, "parent_face"),
            std::pair<int, std::string>(POROUS_JUMP, "porous_jump"),
            std::pair<int, std::string>(PRESSURE_PAR_FIELD, "pressure_par_field"),
            std::pair<int, std::string>(PRESSURE_INLET, "pressure_inlet"),
            std::pair<int, std::string>(PRESSURE_OUTLET, "pressure_outlet"),
            std::pair<int, std::string>(RADIATOR, "radiator"),
            std::pair<int, std::string>(SOLID, "solid"),
            std::pair<int, std::string>(SYMMETRY, "symmetry"),
            std::pair<int, std::string>(VELOCITY_INLET, "velocity_inlet"),
            std::pair<int, std::string>(WALL, "wall"),
            std::pair<int, std::string>(WRAPPER, "wrapper"),
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_zone_type_idx(x);
        else
            return it->second;
    }

    int ZONE::str2idx(const std::string &x)
    {
        static const std::map<std::string, int> mapping_set{
            std::pair<std::string, int>("degassing", DEGASSING),
            std::pair<std::string, int>("exhaust_fan", EXHAUST_FAN),
            std::pair<std::string, int>("fan", FAN),
            std::pair<std::string, int>("fluid", FLUID),
            std::pair<std::string, int>("geometry", GEOMETRY),
            std::pair<std::string, int>("inlet_vent", INLET_VENT),
            std::pair<std::string, int>("intake_fan", INTAKE_FAN),
            std::pair<std::string, int>("interface", INTERFACE),
            std::pair<std::string, int>("interior", INTERIOR),
            std::pair<std::string, int>("internal", INTERNAL),
            std::pair<std::string, int>("mass_flow_inlet", MASS_FLOW_INLET),
            std::pair<std::string, int>("outflow", OUTFLOW),
            std::pair<std::string, int>("outlet_vent", OUTLET_VENT),
            std::pair<std::string, int>("parent_face", PARENT_FACE),
            std::pair<std::string, int>("porous_jump", POROUS_JUMP),
            std::pair<std::string, int>("pressure_par_field", PRESSURE_PAR_FIELD),
            std::pair<std::string, int>("pressure_inlet", PRESSURE_INLET),
            std::pair<std::string, int>("pressure_outlet", PRESSURE_OUTLET),
            std::pair<std::string, int>("radiator", RADIATOR),
            std::pair<std::string, int>("solid", SOLID),
            std::pair<std::string, int>("symmetry", SYMMETRY),
            std::pair<std::string, int>("velocity_inlet", VELOCITY_INLET),
            std::pair<std::string, int>("wall", WALL),
            std::pair<std::string, int>("wrapper", WRAPPER)
        };

        auto it = mapping_set.find(x);
        if (it == mapping_set.end())
            throw invalid_zone_type_str(x);
        else
            return it->second;
    }

    ZONE::ZONE(int zone, const std::string &zt, const std::string &name, int id) :
        SECTION(SECTION::ZONE),
        m_zoneID(zone),
        m_zoneType(formalize(zt)),
        m_zoneName(name),
        m_domainID(id)
    {
        if (!isValidStr(zt))
            throw invalid_zone_type_str(zt);
    }

    size_t ZONE::zone() const
    {
        return m_zoneID;
    }

    size_t &ZONE::zone()
    {
        return m_zoneID;
    }

    const std::string &ZONE::type() const
    {
        return m_zoneType;
    }

    std::string &ZONE::type()
    {
        return m_zoneType;
    }

    const std::string &ZONE::name() const
    {
        return m_zoneName;
    }

    std::string &ZONE::name()
    {
        return m_zoneName;
    }

    int ZONE::domain() const
    {
        return m_domainID;
    }

    int &ZONE::domain()
    {
        return m_domainID;
    }

    void ZONE::repr(std::ostream &out)
    {
        out << std::dec << "(" << identity() << " (" << zone() << " " << type() << " " << name() << ")())" << std::endl;
    }

    struct MESH::internal_error : public std::runtime_error
    {
        internal_error(int err) :
            std::runtime_error("Internal error occurred with error code: " + std::to_string(err) + ".")
        {
            /// Empty body.
        }

        internal_error(const std::string &msg) :
            std::runtime_error("Internal error occurred with error message: \"" + msg + "\".")
        {
            /// Empty body.
        }

        internal_error(int err, const std::string &msg) :
            std::runtime_error("Internal error occurred with error code: " + std::to_string(err) + " and error message: \"" + msg + "\".")
        {
            /// Empty body.
        }
    };

    MESH::MESH() :
        DIM(3), /// 3D by default.
        m_totalNodeNum(0),
        m_totalCellNum(0),
        m_totalFaceNum(0),
        m_totalZoneNum(0)
    {
        /// Empty body.
    }

    MESH::MESH(const std::string &inp, std::ostream &fout) :
        DIM(3), /// 3D by default, may be modified when input mesh is loaded.
        m_totalNodeNum(0),
        m_totalCellNum(0),
        m_totalFaceNum(0),
        m_totalZoneNum(0)
    {
        readFromFile(inp, fout);
    }

    MESH::~MESH()
    {
        clear_entry();
    }

    size_t MESH::numOfNode() const
    {
        return m_totalNodeNum;
    }

    size_t MESH::numOfFace() const
    {
        return m_totalFaceNum;
    }

    size_t MESH::numOfCell() const
    {
        return m_totalCellNum;
    }

    size_t MESH::numOfZone() const
    {
        return m_totalZoneNum;
    }

    const MESH::NODE_ELEM &MESH::node(size_t id) const
    {
        return m_node((int)id);
    }

    MESH::NODE_ELEM &MESH::node(size_t id)
    {
        return m_node((int)id);
    }

    const MESH::FACE_ELEM &MESH::face(size_t id) const
    {
        return m_face((int)id);
    }

    MESH::FACE_ELEM &MESH::face(size_t id)
    {
        return m_face((int)id);
    }

    const MESH::CELL_ELEM &MESH::cell(size_t id) const
    {
        return m_cell((int)id);
    }

    MESH::CELL_ELEM &MESH::cell(size_t id)
    {
        return m_cell((int)id);
    }

    const MESH::ZONE_ELEM &MESH::zone(size_t id, bool isRealZoneID) const
    {
        if (isRealZoneID)
        {
            const auto real_idx = m_zoneMapping.at(id);
            return m_zone.at(real_idx);
        }
        else
            return m_zone((int)id);
    }

    MESH::ZONE_ELEM &MESH::zone(size_t id, bool isRealZoneID)
    {
        if (isRealZoneID)
        {
            const auto real_idx = m_zoneMapping.at(id);
            return m_zone.at(real_idx);
        }
        else
            return m_zone((int)id);
    }

    void MESH::raw2derived()
    {
        /************************* Allocate storage ***************************/
        m_node.resize(numOfNode());
        m_face.resize(numOfFace());
        m_cell.resize(numOfCell());

        /************************ Set initial values **************************/
        for (auto &e : m_node)
        {
            e.coordinate.z() = 0.0;
            e.atBdry = false;
            e.adjacentNode.clear();
            e.dependentFace.clear();
            e.dependentCell.clear();
        }

        for (auto &e : m_face)
        {
            e.center.z() = 0.0;
            e.includedNode.clear();
            e.n_LR.z() = 0.0;
            e.n_RL.z() = 0.0;
        }

        for (auto &e : m_cell)
        {
            e.center.z() = 0.0;
            e.includedFace.clear();
            e.includedNode.clear();
            e.adjacentCell.clear();
            e.n.clear();
            e.S.clear();
        }

        /************************* Parse node and face ************************/
        /// Basic records
        for (auto curPtr : m_content)
        {
            if (curPtr->identity() == SECTION::NODE)
            {
                auto curObj = dynamic_cast<NODE*>(curPtr);
                if (curObj == nullptr)
                    throw internal_error(-1);

                /// Node type within this zone
                const bool flag = curObj->is_boundary_node();

                /// 1-based global node index
                const size_t cur_first = curObj->first_index();
                const size_t cur_last = curObj->last_index();
                for (size_t i = cur_first; i <= cur_last; ++i)
                {
                    /// Node Coordinates
                    node(i).coordinate = curObj->at(i - cur_first);

                    /// Node on boundary or not
                    node(i).atBdry = flag;
                }
            }

            if (curPtr->identity() == SECTION::FACE)
            {
                auto curObj = dynamic_cast<FACE*>(curPtr);
                if (curObj == nullptr)
                    throw internal_error(-2);

                /// 1-based global face index
                const size_t cur_first = curObj->first_index();
                const size_t cur_last = curObj->last_index();

                /// Face type of this zone
                const int ft = curObj->face_type();

                for (size_t i = cur_first; i <= cur_last; ++i)
                {
                    const auto &cnct = curObj->at(i - cur_first);

                    /// Check consistency of face-type
                    if (ft == 0)
                        face(i).type = cnct.x;
                    else if (cnct.x != ft)
                        throw internal_error("local face shape is inconsistent with global specification");
                    else
                        face(i).type = ft;

                    /// Nodes within this face.
                    /// 1-based node index are stored.
                    /// Right-hand convention is preserved.
                    face(i).includedNode.assign(cnct.n, cnct.n + cnct.x);

                    /// Adjacent cells.
                    /// 1-based cell index are stored, 0 stands for boundary.
                    /// Right-hand convention is preserved.
                    const size_t lc = cnct.cl(), rc = cnct.cr();
                    face(i).leftCell = lc;
                    face(i).rightCell = rc;
                    if (lc != 0)
                        cell(lc).includedFace.push_back(i);
                    if (rc != 0)
                        cell(rc).includedFace.push_back(i);

                    /// Face on boundary or not
                    face(i).atBdry = (cnct.c0() == 0 || cnct.c1() == 0);

                    /// Face area, center and unit normal vectors
                    if (cnct.x == FACE::LINEAR)
                    {
                        const size_t na = cnct.n[0], nb = cnct.n[1];
                        const auto &p1 = node(na).coordinate;
                        const auto &p2 = node(nb).coordinate;

                        face(i).area = GridTool::COMMON::line_length(p1, p2);
                        GridTool::COMMON::line_center(p1, p2, face(i).center);
                        GridTool::COMMON::line_normal(p1, p2, face(i).n_LR, face(i).n_RL);
                    }
                    else if (cnct.x == FACE::TRIANGULAR)
                    {
                        const size_t na = cnct.n[0], nb = cnct.n[1], nc = cnct.n[2];
                        const auto &p1 = node(na).coordinate;
                        const auto &p2 = node(nb).coordinate;
                        const auto &p3 = node(nc).coordinate;

                        face(i).area = GridTool::COMMON::triangle_area(p1, p2, p3);
                        GridTool::COMMON::triangle_center(p1, p2, p3, face(i).center);
                        GridTool::COMMON::triangle_normal(p1, p2, p3, face(i).n_LR, face(i).n_RL);
                    }
                    else if (cnct.x == FACE::QUADRILATERAL)
                    {
                        const size_t na = cnct.n[0], nb = cnct.n[1], nc = cnct.n[2], nd = cnct.n[3];
                        const auto &p1 = node(na).coordinate;
                        const auto &p2 = node(nb).coordinate;
                        const auto &p3 = node(nc).coordinate;
                        const auto &p4 = node(nd).coordinate;

                        face(i).area = GridTool::COMMON::quadrilateral_area(p1, p2, p3, p4);
                        GridTool::COMMON::quadrilateral_center(p1, p2, p3, p4, face(i).center);
                        GridTool::COMMON::quadrilateral_normal(p1, p2, p3, p4, face(i).n_LR, face(i).n_RL);
                    }
                    else if (cnct.x == FACE::POLYGONAL)
                        throw FACE::polygon_not_supported();
                    else
                        throw internal_error("face shape not recognized");
                }
            }
        }

        /// Adjacent nodes, dependent faces, and dependent cells of each node.
        /// Step1: Count all occurance
        for (auto curPtr : m_content)
        {
            if (curPtr->identity() == SECTION::FACE)
            {
                auto curObj = dynamic_cast<FACE*>(curPtr);
                if (curObj == nullptr)
                    throw internal_error(-3);

                /// 1-based index
                const size_t cur_first = curObj->first_index();
                const size_t cur_last = curObj->last_index();

                for (size_t i = cur_first; i <= cur_last; ++i)
                {
                    const auto &cnct = curObj->at(i - cur_first);
                    const size_t loc_leftCell = cnct.cl();
                    const size_t loc_rightCell = cnct.cr();

                    for (int j = 0; j < cnct.x; ++j)
                    {
                        const size_t idx_ = cnct.n[j];
                        auto &curNode = node(idx_);

                        const size_t loc_leftNode = cnct.leftAdj(j);
                        const size_t loc_rightNode = cnct.rightAdj(j);

                        /// Adjacent nodes
                        curNode.adjacentNode.push_back(loc_leftNode);
                        if (cnct.x > 2)
                            curNode.adjacentNode.push_back(loc_rightNode);

                        /// Dependent faces
                        curNode.dependentFace.push_back(i);

                        /// Dependent cells
                        if (loc_leftCell != 0)
                            curNode.dependentCell.push_back(loc_leftCell);
                        if (loc_rightCell != 0)
                            curNode.dependentCell.push_back(loc_rightCell);
                    }
                }
            }
        }

        /// Step2: Remove duplication
        for (size_t i = 1; i <= numOfNode(); ++i)
        {
            auto &curNode = node(i);

            const std::set<size_t> st1(curNode.adjacentNode.begin(), curNode.adjacentNode.end());
            curNode.adjacentNode.assign(st1.begin(), st1.end());

            const std::set<size_t> st2(curNode.dependentCell.begin(), curNode.dependentCell.end());
            curNode.dependentCell.assign(st2.begin(), st2.end());
        }

        /*********************** Parse records of cell ************************/
        for (auto curPtr : m_content)
        {
            if (curPtr->identity() == SECTION::CELL)
            {
                auto curObj = dynamic_cast<CELL*>(curPtr);
                if (curObj == nullptr)
                    throw internal_error(-4);

                /// 1-based global face index
                const size_t cur_first = curObj->first_index();
                const size_t cur_last = curObj->last_index();

                for (size_t i = cur_first; i <= cur_last; ++i)
                {
                    auto &curCell = cell(i);

                    /// Element type of cells in this zone
                    curCell.type = curObj->at(i - cur_first);

                    /// Organize order of included nodes and faces
                    cell_standardization(curCell);

                    /// Adjacent cells and normal
                    /// Allocate storage.
                    curCell.adjacentCell.resize(curCell.includedFace.size());
                    curCell.n.resize(curCell.includedFace.size());
                    curCell.S.resize(curCell.includedFace.size());

                    for (size_t j = 0; j < curCell.includedFace.size(); ++j)
                    {
                        const size_t f_idx = curCell.includedFace[j];
                        const auto &f = face(f_idx);
                        const auto c0 = f.leftCell, c1 = f.rightCell;
                        if (c0 == i)
                        {
                            curCell.adjacentCell[j] = c1;
                            curCell.n[j] = f.n_LR;
                        }
                        else if (c1 == i)
                        {
                            curCell.adjacentCell[j] = c0;
                            curCell.n[j] = f.n_RL;
                        }
                        else
                            throw internal_error(-5);

                        for (int k = 1; k <= dimension(); ++k)
                            curCell.S[j](k) = f.area * curCell.n[j](k);
                    }

                    /// Volume and Centroid.
                    /// Based on the divergence theorem. See (5.15) and (5.17) of Jiri Blazek's CFD book.
                    curCell.volume = 0.0;
                    curCell.center.x() = 0.0; curCell.center.y() = 0.0; curCell.center.z() = 0.0;
                    for (size_t j = 0; j < curCell.includedFace.size(); ++j)
                    {
                        const auto cfi = curCell.includedFace.at(j);
                        const auto &cf = face(cfi);
                        const auto &cf_c = cf.center;
                        const auto &cf_S = curCell.S.at(j);
                        const auto w = cf_c.dot(cf_S);
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

        /*********************** Parse records of zone ************************/
        m_totalZoneNum = 0;
        m_zoneMapping.clear();
        for (auto curPtr : m_content) // Determine the total num of zones.
        {
            auto curObj = dynamic_cast<RANGE*>(curPtr);
            if (curObj == nullptr)
                continue;

            const auto curZone = curObj->zone();
            if (m_zoneMapping.find(curZone) == m_zoneMapping.end())
                m_zoneMapping[curZone] = m_totalZoneNum++;
            else
                throw internal_error("Duplicated zone detected.");
        }
        m_zone.clear();
        m_zone.resize(numOfZone());
        for (auto curPtr : m_content) // Assign zone contents
        {
            auto curObj = dynamic_cast<RANGE*>(curPtr);
            if (curObj == nullptr)
                continue;

            const auto curZoneIdx = curObj->zone();
            auto &curZone = m_zone.at(m_zoneMapping.at(curZoneIdx));
            curZone.ID = curZoneIdx;
            curZone.obj = curObj;
        }
        for (auto e : m_content)
        {
            const auto curObj = dynamic_cast<ZONE*>(e);
            if (curObj == nullptr)
                continue;

            const auto curZoneIdx = curObj->zone();
            auto &curZone = m_zone.at(m_zoneMapping.at(curZoneIdx));
            curZone.name = curObj->name();
            curZone.type = curObj->type();
        }
    }

    void MESH::add_entry(SECTION *e)
    {
        m_content.push_back(e);
    }

    void MESH::clear_entry()
    {
        /// Release previous contents.
        for (auto ptr : m_content)
            delete ptr; /// NO side-effect when deleting a nullptr.

        /// Clear container.
        m_content.clear();
    }

    void MESH::readFromFile(const std::string &src, std::ostream &fout)
    {
        // Open grid file
        std::ifstream fin(src);
        if (fin.fail())
            throw std::runtime_error("Failed to open input grid file: \"" + src + "\".");

        // Clear existing records if any.
        clear_entry();

        // Read contents
        while (!fin.eof())
        {
            skip_white(fin);
            eat(fin, '(');
            int ti = -1;
            fin >> std::dec >> ti;
            if (ti == SECTION::COMMENT)
            {
                eat(fin, '\"');
                std::string ts;
                char tc;
                while ((tc = fin.get()) != '\"')
                    ts.push_back(tc);
                eat(fin, ')');
                add_entry(new COMMENT(ts));
                skip_white(fin);
            }
            else if (ti == SECTION::HEADER)
            {
                eat(fin, '\"');
                std::string ts;
                char tc;
                while ((tc = fin.get()) != '\"')
                    ts.push_back(tc);
                eat(fin, ')');
                add_entry(new HEADER(ts));
                skip_white(fin);
            }
            else if (ti == SECTION::DIMENSION)
            {
                int nd = 0;
                fin >> std::dec >> nd;
                eat(fin, ')');
                add_entry(new DIMENSION(nd));
                skip_white(fin);
                m_dim = nd;
                m_is3D = (nd == 3);
            }
            else if (ti == SECTION::NODE)
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
                    fout << "Total number of nodes: " << m_totalNodeNum << std::endl;
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
                    auto e = new NODE(zone, first, last, tp, nd);
                    eat(fin, ')');
                    eat(fin, '(');
                    fout << "Reading " << e->num() << " nodes in zone " << zone << " (from " << first << " to " << last << "), whose type is \"" << NODE::idx2str(tp) << "\"  ... ";

                    if (nd != dimension())
                        throw std::runtime_error("Inconsistent with previous DIMENSION declaration!");

                    if (nd == 3)
                    {
                        for (int i = first; i <= last; ++i)
                        {
                            auto &ce = e->at(i - first);
                            fin >> ce.x() >> ce.y() >> ce.z();
                        }
                    }
                    else
                    {
                        for (int i = first; i <= last; ++i)
                        {
                            auto &ce = e->at(i - first);
                            fin >> ce.x() >> ce.y();
                        }
                    }
                    eat(fin, ')');
                    eat(fin, ')');
                    fout << "Done!" << std::endl;
                    add_entry(e);
                }
                skip_white(fin);
            }
            else if (ti == SECTION::CELL)
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
                    fout << "Total number of cells: " << m_totalCellNum << std::endl;
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
                    auto e = new CELL(zone, first, last, tp, elem);
                    eat(fin, ')');

                    if (elem == 0)
                    {
                        fout << "Reading " << e->num() << " mixed cells in zone " << zone << " (from " << first << " to " << last << ") ... ";
                        eat(fin, '(');
                        for (int i = first; i <= last; ++i)
                        {
                            fin >> elem;
                            if (CELL::isValidElemIdx(elem))
                                e->at(i - first) = elem;
                            else
                                throw std::runtime_error("Invalid CELL-ELEM-TYPE: \"" + std::to_string(elem) + "\"");
                        }
                        eat(fin, ')');
                        fout << "Done!" << std::endl;
                    }
                    else
                        fout << e->num() << " " << CELL::idx2str_elem(elem) << " in zone " << zone << " (from " << first << " to " << last << ")" << std::endl;

                    eat(fin, ')');
                    add_entry(e);
                }
                skip_white(fin);
            }
            else if (ti == SECTION::FACE)
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
                    fout << "Total number of faces: " << m_totalFaceNum << std::endl;
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
                    auto e = new FACE(zone, first, last, bc, face);
                    eat(fin, ')');
                    eat(fin, '(');
                    fout << "Reading " << e->num() << " " << FACE::idx2str(face) << " faces in zone " << zone << " (from " << first << " to " << last << "), whose B.C. is \"" << BC::idx2str(bc) << "\" ... ";

                    size_t tmp_n[4];
                    size_t tmp_c[2];
                    if (face == FACE::MIXED)
                    {
                        int x = -1;
                        for (size_t i = first; i <= last; ++i)
                        {
                            // Read connectivity record
                            fin >> x;
                            if (x <= 1 || x >= 5)
                                throw std::invalid_argument("Invalid node num in the mixed face.");
                            for (int j = 0; j < x; ++j)
                                fin >> tmp_n[j];
                            fin >> tmp_c[0] >> tmp_c[1];

                            // Store current connectivity info
                            e->at(i - first).set(x, tmp_n, tmp_c);
                        }
                    }
                    else
                    {
                        for (size_t i = first; i <= last; ++i)
                        {
                            // Read connectivity record
                            for (int j = 0; j < face; ++j)
                                fin >> tmp_n[j];
                            fin >> tmp_c[0] >> tmp_c[1];

                            // Store current connectivity info
                            e->at(i - first).set(face, tmp_n, tmp_c);
                        }
                    }
                    eat(fin, ')');
                    eat(fin, ')');
                    fout << "Done!" << std::endl;
                    add_entry(e);
                }
                skip_white(fin);
            }
            else if (ti == SECTION::ZONE || ti == SECTION::ZONE_MESHING)
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
                auto e = new ZONE(zone, ztp, zname);
                add_entry(e);
                skip_white(fin);
                fout << "ZONE " << e->zone() << ", named " << R"(")" << e->name() << R"(", )" << "is " << R"(")" << e->type() << R"(")" << std::endl;
            }
            else
                throw std::runtime_error("Unsupported section index: " + std::to_string(ti));
        }

        // Close grid file
        fin.close();

        // Re-orginize grid connectivities in a much easier way,
        // and compute some derived quantities.
        fout << "Converting into high-level representation ... ";
        raw2derived();
        fout << "Done!" << std::endl;
    }

    void MESH::writeToFile(const std::string &dst) const
    {
        if (numOfCell() == 0)
            throw std::runtime_error("Invalid num of cells.");
        if (numOfFace() == 0)
            throw std::runtime_error("Invalid num of faces.");
        if (numOfNode() == 0)
            throw std::runtime_error("Invalid num of nodes.");
        if (m_content.empty())
            throw std::runtime_error("Invalid num of contents.");

        /// Open grid file
        std::ofstream fout(dst);
        if (fout.fail())
            throw std::runtime_error("Failed to open output grid file: " + dst);

        /// Write until dimension declaration
        size_t i = 0;
        while (true)
        {
            m_content[i]->repr(fout);
            bool flag = dynamic_cast<DIMENSION*>(m_content[i]) != nullptr;
            ++i;
            if (flag)
                break;
        }

        /// Declaration of NODE, FACE, CELL
        fout << "(" << std::dec << SECTION::NODE << " (";
        fout << std::hex << 0 << " " << 1 << " " << m_totalNodeNum << " ";
        fout << std::dec << 0 << " " << (m_is3D ? 3 : 2) << "))" << std::endl;
        fout << "(" << std::dec << SECTION::CELL << " (";
        fout << std::hex << 0 << " " << 1 << " " << m_totalCellNum << " ";
        fout << std::dec << 0 << " " << 0 << "))" << std::endl;
        fout << "(" << std::dec << SECTION::FACE << " (";
        fout << std::hex << 0 << " " << 1 << " " << m_totalFaceNum << " ";
        fout << std::dec << 0 << " " << 0 << "))" << std::endl;

        /// Contents
        for (; i < m_content.size(); ++i)
            m_content[i]->repr(fout);

        /// Close grid file
        fout.close();
    }

    void MESH::cell_standardization(CELL_ELEM &c)
    {
        switch (c.type)
        {
        case CELL::TETRAHEDRAL:
            tet_standardization(c);
            break;
        case CELL::HEXAHEDRAL:
            hex_standardization(c);
            break;
        case CELL::PYRAMID:
            pyramid_standardization(c);
            break;
        case CELL::WEDGE:
            prism_standardization(c);
            break;
        case CELL::TRIANGULAR:
            triangle_standardization(c);
            break;
        case CELL::QUADRILATERAL:
            quad_standardization(c);
            break;
        default:
            throw CELL::invalid_cell_type_idx(c.type); /// May caused by internal error.
        }
    }

    void MESH::tet_standardization(CELL_ELEM &tet)
    {
        // Check num of total faces
        if (tet.includedFace.size() != 4)
            throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::idx2str_elem(tet.type) + R"(" and num of faces: )" + std::to_string(tet.includedFace.size()));

        // Ensure all faces are triangular
        const auto &f0 = face(tet.includedFace.at(0));
        if (f0.type != FACE::TRIANGULAR)
            throw std::runtime_error("Internal error.");

        const auto &f1 = face(tet.includedFace.at(1));
        if (f1.type != FACE::TRIANGULAR)
            throw std::runtime_error("Internal error.");

        const auto &f2 = face(tet.includedFace.at(2));
        if (f2.type != FACE::TRIANGULAR)
            throw std::runtime_error("Internal error.");

        const auto &f3 = face(tet.includedFace.at(3));
        if (f3.type != FACE::TRIANGULAR)
            throw std::runtime_error("Internal error.");

        // Find all 4 nodes
        size_t n1 = f0.includedNode.at(0);
        if (n1 == 0)
            throw std::runtime_error("Internal error.");

        size_t n2 = f0.includedNode.at(1);
        if (n2 == 0)
            throw std::runtime_error("Internal error.");

        size_t n3 = f0.includedNode.at(2);
        if (n3 == 0)
            throw std::runtime_error("Internal error.");

        size_t n0 = 0;
        for (auto e : f2.includedNode)
            if (e != n1 && e != n2 && e != n3)
            {
                n0 = e;
                break;
            }
        if (n0 == 0)
            throw std::runtime_error("Internal error.");

        // Assign node index
        tet.includedNode.resize(4);
        tet.includedNode.at(0) = n0;
        tet.includedNode.at(1) = n1;
        tet.includedNode.at(2) = n2;
        tet.includedNode.at(3) = n3;
    }

    void MESH::pyramid_standardization(CELL_ELEM &pyramid)
    {
        // Check num of total faces
        if (pyramid.includedFace.size() != 5)
            throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::idx2str_elem(pyramid.type) + R"(" and num of faces: )" + std::to_string(pyramid.includedFace.size()));

        // Find the bottom quad and ensure other faces are triangular.
        size_t f0_idx = 0;
        for (auto e : pyramid.includedFace)
        {
            const auto &f = face(e);
            if (f.type == FACE::QUADRILATERAL)
            {
                if (f0_idx == 0)
                    f0_idx = e;
                else
                    throw std::runtime_error("Internal error.");
            }
            else if (f.type == FACE::TRIANGULAR)
                continue;
            else
                throw std::runtime_error("Internal error.");
        }
        if (f0_idx == 0)
            throw std::runtime_error("Internal error.");

        // Nodes at bottom
        const auto &f0 = face(f0_idx);

        const size_t n0 = f0.includedNode.at(0);
        if (n0 == 0)
            throw std::runtime_error("Internal error.");

        const size_t n1 = f0.includedNode.at(1);
        if (n1 == 0)
            throw std::runtime_error("Internal error.");

        const size_t n2 = f0.includedNode.at(2);
        if (n2 == 0)
            throw std::runtime_error("Internal error.");

        const size_t n3 = f0.includedNode.at(3);
        if (n3 == 0)
            throw std::runtime_error("Internal error.");

        // Find other 4 triangles
        size_t f1_idx = 0;
        for (auto e : pyramid.includedFace)
        {
            if (e == f0_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n0, n3))
            {
                f1_idx = e;
                break;
            }
        }
        if (f1_idx == 0)
            throw std::runtime_error("Internal error.");

        size_t f2_idx = 0;
        for (auto e : pyramid.includedFace)
        {
            if (e == f0_idx || e == f1_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n3, n2))
            {
                f2_idx = e;
                break;
            }
        }
        if (f2_idx == 0)
            throw std::runtime_error("Internal error.");

        size_t f3_idx = 0;
        for (auto e : pyramid.includedFace)
        {
            if (e == f0_idx || e == f1_idx || e == f2_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n2, n1))
            {
                f3_idx = e;
                break;
            }
        }
        if (f3_idx == 0)
            throw std::runtime_error("Internal error.");

        size_t f4_idx = 0;
        for (auto e : pyramid.includedFace)
        {
            if (e != f0_idx && e != f1_idx && e != f2_idx && e != f3_idx)
            {
                f4_idx = e;
                const auto &f = face(e);
                if (!f.includedNode.contains(n1, n0))
                    throw std::runtime_error("Internal error.");

                break;
            }
        }
        if (f4_idx == 0)
            throw std::runtime_error("Internal error.");

        // Assign face index
        pyramid.includedFace.at(0) = f0_idx;
        pyramid.includedFace.at(1) = f1_idx;
        pyramid.includedFace.at(2) = f2_idx;
        pyramid.includedFace.at(3) = f3_idx;
        pyramid.includedFace.at(4) = f4_idx;

        // The last node
        size_t n4 = 0;
        const auto &f1 = face(f1_idx);
        for (auto e : f1.includedNode)
            if (e != n0 && e != n3)
            {
                n4 = e;
                break;
            }
        if (n4 == 0)
            throw std::runtime_error("Internal error.");

        // Assign node index
        pyramid.includedNode.resize(5);
        pyramid.includedNode.at(0) = n0;
        pyramid.includedNode.at(1) = n1;
        pyramid.includedNode.at(2) = n2;
        pyramid.includedNode.at(3) = n3;
        pyramid.includedNode.at(4) = n4;
    }

    void MESH::prism_standardization(CELL_ELEM &prism)
    {
        // Check num of total faces
        if (prism.includedFace.size() != 5)
            throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::idx2str_elem(prism.type) + R"(" and num of faces: )" + std::to_string(prism.includedFace.size()));

        // Ensure there're only 2 triangle and 3 quad
        size_t f0_idx = 0, f1_idx = 0;
        for (auto e : prism.includedFace)
        {
            const auto &f = face(e);
            if (f.type == FACE::TRIANGULAR)
            {
                if (f0_idx == 0)
                    f0_idx = e;
                else if (f1_idx == 0)
                    f1_idx = e;
                else
                    throw std::runtime_error("There're more than 2 triangular faces in a prism cell.");
            }
            else if (f.type == FACE::QUADRILATERAL)
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
        const size_t n0 = f0.includedNode.at(0);
        const size_t n1 = f0.includedNode.at(1);
        const size_t n2 = f0.includedNode.at(2);

        // Find face 4
        size_t f4_idx = 0;
        for (auto e : prism.includedFace)
        {
            const auto &f = face(e);
            if (f.type == FACE::QUADRILATERAL && f.includedNode.contains(n0, n1))
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
        for (auto e : prism.includedFace)
        {
            const auto &f = face(e);
            if (f.type == FACE::QUADRILATERAL && f.includedNode.contains(n1, n2))
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
        for (auto e : prism.includedFace)
        {
            const auto &f = face(e);
            if (f.type == FACE::QUADRILATERAL)
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
        if (!f2.includedNode.contains(n2, n0))
            throw std::runtime_error("Inconsistent face composition.");

        // Assign face index
        prism.includedFace.at(0) = f0_idx;
        prism.includedFace.at(1) = f1_idx;
        prism.includedFace.at(2) = f2_idx;
        prism.includedFace.at(3) = f3_idx;
        prism.includedFace.at(4) = f4_idx;

        // 3 nodes on the top triangular face to be defined
        size_t n3 = 0, n5 = 0, n4 = 0;
        for (auto e : f2.includedNode)
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
        if (!(f1.includedNode.contains(n3) && f4.includedNode.contains(n3)))
            std::swap(n3, n5);

        // Ensure n5 located on f1, f2 and f3
        if (!(f1.includedNode.contains(n5) && f3.includedNode.contains(n5)))
            throw std::runtime_error("Internal error.");

        // Find n4
        for (auto e : f1.includedNode)
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
        prism.includedNode.resize(6);
        prism.includedNode.at(0) = n0;
        prism.includedNode.at(1) = n1;
        prism.includedNode.at(2) = n2;
        prism.includedNode.at(3) = n3;
        prism.includedNode.at(4) = n4;
        prism.includedNode.at(5) = n5;
    }

    void MESH::hex_standardization(CELL_ELEM &hex)
    {
        // Check num of total faces
        if (hex.includedFace.size() != 6)
            throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::idx2str_elem(hex.type) + R"(" and num of faces: )" + std::to_string(hex.includedFace.size()));

        // Ensure all faces are quad
        for (auto e : hex.includedFace)
        {
            if (e == 0)
                throw std::runtime_error("Internal error.");
            const auto &f = face(e);
            if (f.type != FACE::QUADRILATERAL)
                throw std::runtime_error(R"(Inconsistent face type ")" + FACE::idx2str(f.type) + R"(" in a hex cell.)");
        }

        // Face 4 at bottom
        const size_t f4_idx = hex.includedFace.at(0);
        const auto &f4 = face(f4_idx);
        const size_t n0 = f4.includedNode.at(0);
        const size_t n1 = f4.includedNode.at(1);
        const size_t n2 = f4.includedNode.at(2);
        const size_t n3 = f4.includedNode.at(3);
        if (n0 == 0 || n1 == 0 || n2 == 0 || n3 == 0)
            throw std::runtime_error("Internal error.");

        // Face 0
        size_t f0_idx = 0;
        for (auto e : hex.includedFace)
        {
            if (e == f4_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n3, n0))
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
        if (f0.includedNode.size() != 4)
            throw std::runtime_error("Internal error.");

        // Face 2
        size_t f2_idx = 0;
        for (auto e : hex.includedFace)
        {
            if (e == f4_idx || e == f0_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n0, n1))
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
        if (f2.includedNode.size() != 4)
            throw std::runtime_error("Internal error.");

        // Face 1
        size_t f1_idx = 0;
        for (auto e : hex.includedFace)
        {
            if (e == f4_idx || e == f0_idx || e == f2_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n1, n2))
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
        if (f1.includedNode.size() != 4)
            throw std::runtime_error("Internal error.");

        // Face 3
        size_t f3_idx = 0;
        for (auto e : hex.includedFace)
        {
            if (e == f4_idx || e == f0_idx || e == f2_idx || e == f1_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n2, n3))
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
        if (f3.includedNode.size() != 4)
            throw std::runtime_error("Internal error.");

        // Face 5
        size_t f5_idx = 0;
        for (auto e : hex.includedFace)
        {
            if (e == f4_idx || e == f0_idx || e == f2_idx || e == f1_idx || e == f3_idx)
                continue;

            f5_idx = e;
            break;
        }
        if (f5_idx == 0)
            throw std::runtime_error("Missing face 5");

        const auto &f5 = face(f5_idx);
        if (f5.includedNode.size() != 4)
            throw std::runtime_error("Internal error.");

        // Assign face index
        hex.includedFace.at(0) = f0_idx;
        hex.includedFace.at(1) = f1_idx;
        hex.includedFace.at(2) = f2_idx;
        hex.includedFace.at(3) = f3_idx;
        hex.includedFace.at(4) = f4_idx;
        hex.includedFace.at(5) = f5_idx;

        // 4 nodes at top to be defined
        size_t n4 = 0, n5 = 0, n6 = 0, n7 = 0;

        // Node 4
        for (auto e : f5.includedNode)
        {
            if (f0.includedNode.contains(e) && f2.includedNode.contains(e))
            {
                n4 = e;
                break;
            }
        }
        if (n4 == 0 || n4 == n0)
            throw std::runtime_error("Missing node 4");

        // Node 5
        for (auto e : f2.includedNode)
        {
            if (e == n0 || e == n1 || e == n4)
                continue;

            n5 = e;
            break;
        }
        if (n5 == 0)
            throw std::runtime_error("Missing node 5");

        // Node 6
        for (auto e : f1.includedNode)
        {
            if (e == n1 || e == n2 || e == n5)
                continue;

            n6 = e;
            break;
        }
        if (n6 == 0)
            throw std::runtime_error("Missing node 6");

        // Node 7
        for (auto e : f0.includedNode)
        {
            if (e == n0 || e == n3 || e == n4)
                continue;

            n7 = e;
            break;
        }
        if (n7 == 0)
            throw std::runtime_error("Missing node 7");

        if (!f3.includedNode.contains(n6, n7))
            throw std::runtime_error("Inconsistent node composition.");

        // Assign node index
        hex.includedNode.resize(8);
        hex.includedNode.at(0) = n0;
        hex.includedNode.at(1) = n1;
        hex.includedNode.at(2) = n2;
        hex.includedNode.at(3) = n3;
        hex.includedNode.at(4) = n4;
        hex.includedNode.at(5) = n5;
        hex.includedNode.at(6) = n6;
        hex.includedNode.at(7) = n7;
    }

    void MESH::triangle_standardization(CELL_ELEM &tri)
    {
        // Check num of total faces
        if (tri.includedFace.size() != 3)
            throw std::runtime_error("Inconsitent face composition.");

        // Ensure all faces are lines
        for (auto e : tri.includedFace)
        {
            const auto &f = face(e);
            if (e == 0 || f.type != FACE::LINEAR || f.includedNode.size() != 2)
                throw std::runtime_error("Invalid face detected.");
        }

        // Faces
        // Keep the order of faces as it is.
        // Doesn't matter as any order of the 3 faces is valid.
        const auto f0_idx = tri.includedFace.at(0);
        const auto f1_idx = tri.includedFace.at(1);
        const auto f2_idx = tri.includedFace.at(2);

        const auto &f0 = face(f0_idx);
        const auto &f1 = face(f1_idx);
        const auto &f2 = face(f2_idx);

        // Nodes
        const size_t n0 = f0.includedNode.at(0);
        const size_t n1 = f0.includedNode.at(1);
        if (n0 == 0 || n1 == 0 || n0 == n1)
            throw std::runtime_error("Missing node detected.");

        size_t n2 = 0;
        for (auto e : f1.includedNode)
        {
            if (f0.includedNode.contains(e))
                continue;

            n2 = e;
            break;
        }
        if (n2 == 0 || n2 == n0 || n2 == n1 || !f2.includedNode.contains(n2))
            throw std::runtime_error("Node 2 is missing.");

        // Assign node index
        tri.includedNode.resize(3);
        tri.includedNode.at(0) = n0;
        tri.includedNode.at(1) = n1;
        tri.includedNode.at(2) = n2;
    }

    void MESH::quad_standardization(CELL_ELEM &quad)
    {
        // Check num of total faces
        if (quad.includedFace.size() != 4)
            throw std::runtime_error("Inconsitent face composition.");

        // Ensure all faces are lines
        for (auto e : quad.includedFace)
        {
            const auto &f = face(e);
            if (e == 0 || f.type != FACE::LINEAR || f.includedNode.size() != 2)
                throw std::runtime_error("Invalid face detected.");
        }

        // Face 0
        const auto f0_idx = quad.includedFace.at(0);
        const auto &f0 = face(f0_idx);

        // Node 0 and 1
        const auto n0 = f0.includedNode.at(0);
        const auto n1 = f0.includedNode.at(1);

        // Face 1
        size_t f1_idx = 0;
        for (auto e : quad.includedFace)
        {
            if (e == f0_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n1))
            {
                f1_idx = e;
                break;
            }
        }
        if (f1_idx == 0)
            throw std::runtime_error("Missing face 1");

        const auto &f1 = face(f1_idx);

        // Node 2
        const size_t n2 = f1.includedNode.at(0) == n1 ? f1.includedNode.at(1) : f1.includedNode.at(0);

        // Face 3
        size_t f3_idx = 0;
        for (auto e : quad.includedFace)
        {
            if (e == f0_idx || e == f1_idx)
                continue;

            const auto &f = face(e);
            if (f.includedNode.contains(n0))
            {
                f3_idx = e;
                break;
            }
        }
        if (f3_idx == 0)
            throw std::runtime_error("Missing face 3");

        const auto &f3 = face(f3_idx);

        // Node 3
        const size_t n3 = f3.includedNode.at(0) == n0 ? f3.includedNode.at(1) : f3.includedNode.at(0);

        // Check face 2
        size_t f2_idx = 0;
        for (auto e : quad.includedFace)
        {
            if (e == f0_idx || e == f1_idx || e == f3_idx)
                continue;

            f2_idx = e;
            break;
        }
        if (f2_idx == 0)
            throw std::runtime_error("Missing face 2");

        const auto &f2 = face(f2_idx);
        if (!f2.includedNode.contains(n2, n3))
            throw std::runtime_error("Inconsistent node and face includance on face 2");

        // Assign face index
        quad.includedFace.at(1) = f1_idx;
        quad.includedFace.at(2) = f2_idx;
        quad.includedFace.at(3) = f3_idx;

        // Assign node index
        quad.includedNode.resize(4);
        quad.includedNode.at(0) = n0;
        quad.includedNode.at(1) = n1;
        quad.includedNode.at(2) = n2;
        quad.includedNode.at(3) = n3;
    }
}
