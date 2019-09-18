#ifndef __NMF_H__
#define __NMF_H__

#include <vector>
#include <cstddef>
#include <string>
#include <map>
#include <utility>
#include <stdexcept>

class NMF_Range
{
private:
	size_t m_blk; // Block index, 1-based.
	size_t m_face; // Face index, ranges from 1 to 6.
	size_t m_s1; // Primary direction starting index, 1-based.
	size_t m_e1; // Primary direction ending index, 1-based.
	size_t m_s2; // Secondary direction starting index, 1-based.
	size_t m_e2; // Secondary direction ending index, 1-based.

public:
	NMF_Range()
	{
		m_blk = 0;
		m_face = 0;
		m_s1 = 0;
		m_e1 = 0;
		m_s2 = 0;
		m_e2 = 0;
	}

	NMF_Range(size_t *src)
	{
		m_blk = src[0];
		m_face = src[1];
		m_s1 = src[2];
		m_e1 = src[3];
		m_s2 = src[4];
		m_e2 = src[5];
	}

	NMF_Range(size_t b, size_t f, size_t s1, size_t e1, size_t s2, size_t e2)
	{
		m_blk = b;
		m_face = f;
		m_s1 = s1;
		m_e1 = e1;
		m_s2 = s2;
		m_e2 = e2;
	}

	~NMF_Range() = default;

	size_t B() const
	{
		return m_blk;
	}

	size_t F() const
	{
		return m_face;
	}

	size_t S1() const
	{
		return m_s1;
	}

	size_t E1() const
	{
		return m_e1;
	}

	size_t S2() const
	{
		return m_s2;
	}

	size_t E2() const
	{
		return m_e2;
	}

	// Check if given index is within this range.
	bool constains(size_t pri, size_t sec) const
	{
		const bool t1 = (m_s1 <= pri) && (pri <= m_e1);
		const bool t2 = (m_s2 <= sec) && (sec <= m_e2);
		return t1 && t2;
	}

	// Nodes in primary direction.
	size_t pri_node_num() const
	{
		return m_e1 - m_s1 + 1;
	}

	// Nodes in secondary direction.
	size_t sec_node_num() const
	{
		return m_e2 - m_s2 + 1;
	}

	// Total nodes on this interface.
	size_t node_num() const
	{
		return pri_node_num() * sec_node_num();
	}

	// Total faces/edges on this interface.
	size_t face_num() const
	{
		const size_t n_pri = (pri_node_num() - 1) * sec_node_num();
		const size_t n_sec = (sec_node_num() - 1) * pri_node_num();
		return n_pri + n_sec;
	}

	// Total quad cells on this interface.
	size_t cell_num() const
	{
		return (pri_node_num() - 1) * (sec_node_num() - 1);
	}
};

typedef NMF_Range *pNMF_Range;

typedef enum {
	COLLAPSED = 1,
	ONE_TO_ONE = 2,
	PATCHED = 3,
	POLE_DIR1 = 4,
	POLE_DIR2 = 5,
	SYM_X = 6,
	SYM_Y = 7,
	SYM_Z = 8,
	UNPROCESSED = 9,
	WALL = 10,
	SYM = 11
}NMF_BC;

const std::string &NMF_BC_STR(NMF_BC x)
{
	static const std::map<NMF_BC, std::string> BC_NAME
	{
		std::pair<NMF_BC, std::string>(NMF_BC::COLLAPSED, "COLLAPSED"),
		std::pair<NMF_BC, std::string>(NMF_BC::ONE_TO_ONE, "ONE_TO_ONE"),
		std::pair<NMF_BC, std::string>(NMF_BC::PATCHED, "PATCHED"),
		std::pair<NMF_BC, std::string>(NMF_BC::POLE_DIR1, "POLE_DIR1"),
		std::pair<NMF_BC, std::string>(NMF_BC::POLE_DIR2, "POLE_DIR2"),
		std::pair<NMF_BC, std::string>(NMF_BC::SYM_X, "SYM_X"),
		std::pair<NMF_BC, std::string>(NMF_BC::SYM_Y, "SYM_Y"),
		std::pair<NMF_BC, std::string>(NMF_BC::SYM_Z, "SYM_Z"),
		std::pair<NMF_BC, std::string>(NMF_BC::UNPROCESSED, "UNPROCESSED"),
		std::pair<NMF_BC, std::string>(NMF_BC::WALL, "WALL"),
		std::pair<NMF_BC, std::string>(NMF_BC::SYM, "SYM")
	};

	auto it = BC_NAME.find(x);
	if (it == BC_NAME.end())
		throw std::runtime_error("Unsupported B.C. enumeration: " + std::to_string(x));
	else
		return it->second;
}

NMF_BC NMF_BC_ENUM(const std::string &x)
{
	static const std::map<std::string, NMF_BC> BC_ENUM
	{
		// COLLAPSED
		std::pair<std::string, NMF_BC>("COLLAPSED", NMF_BC::COLLAPSED),
		std::pair<std::string, NMF_BC>("Collapsed", NMF_BC::COLLAPSED),
		std::pair<std::string, NMF_BC>("collapsed", NMF_BC::COLLAPSED),
		// ONE_TO_ONE
		std::pair<std::string, NMF_BC>("ONE_TO_ONE", NMF_BC::ONE_TO_ONE),
		std::pair<std::string, NMF_BC>("One_To_One", NMF_BC::ONE_TO_ONE),
		std::pair<std::string, NMF_BC>("One_to_One", NMF_BC::ONE_TO_ONE),
		std::pair<std::string, NMF_BC>("one_to_one", NMF_BC::ONE_TO_ONE),
		// PATCHED
		std::pair<std::string, NMF_BC>("PATCHED", NMF_BC::PATCHED),
		std::pair<std::string, NMF_BC>("Patched", NMF_BC::PATCHED),
		std::pair<std::string, NMF_BC>("patched", NMF_BC::PATCHED),
		// POLE_DIR1
		std::pair<std::string, NMF_BC>("POLE_DIR1", NMF_BC::POLE_DIR1),
		std::pair<std::string, NMF_BC>("Pole_Dir1", NMF_BC::POLE_DIR1),
		std::pair<std::string, NMF_BC>("pole_dir1", NMF_BC::POLE_DIR1),
		// POLE_DIR2
		std::pair<std::string, NMF_BC>("POLE_DIR2", NMF_BC::POLE_DIR2),
		std::pair<std::string, NMF_BC>("Pole_Dir2", NMF_BC::POLE_DIR2),
		std::pair<std::string, NMF_BC>("pole_dir2", NMF_BC::POLE_DIR2),
		// SYM_X
		std::pair<std::string, NMF_BC>("SYM_X", NMF_BC::SYM_X),
		std::pair<std::string, NMF_BC>("Sym_X", NMF_BC::SYM_X),
		std::pair<std::string, NMF_BC>("Sym_x", NMF_BC::SYM_X),
		std::pair<std::string, NMF_BC>("sym_X", NMF_BC::SYM_X),
		std::pair<std::string, NMF_BC>("sym_x", NMF_BC::SYM_X),
		// SYM_Y
		std::pair<std::string, NMF_BC>("SYM_Y", NMF_BC::SYM_Y),
		std::pair<std::string, NMF_BC>("Sym_Y", NMF_BC::SYM_Y),
		std::pair<std::string, NMF_BC>("Sym_y", NMF_BC::SYM_Y),
		std::pair<std::string, NMF_BC>("sym_Y", NMF_BC::SYM_Y),
		std::pair<std::string, NMF_BC>("sym_y", NMF_BC::SYM_Y),
		// SYM_Z
		std::pair<std::string, NMF_BC>("SYM_Z", NMF_BC::SYM_Z),
		std::pair<std::string, NMF_BC>("Sym_Z", NMF_BC::SYM_Z),
		std::pair<std::string, NMF_BC>("Sym_z", NMF_BC::SYM_Z),
		std::pair<std::string, NMF_BC>("sym_Z", NMF_BC::SYM_Z),
		std::pair<std::string, NMF_BC>("sym_z", NMF_BC::SYM_Z),
		// UNPROCESSED
		std::pair<std::string, NMF_BC>("UNPROCESSED", NMF_BC::UNPROCESSED),
		std::pair<std::string, NMF_BC>("Unprocessed", NMF_BC::UNPROCESSED),
		std::pair<std::string, NMF_BC>("unprocessed", NMF_BC::UNPROCESSED),
		// WALL
		std::pair<std::string, NMF_BC>("WALL", NMF_BC::WALL),
		std::pair<std::string, NMF_BC>("Wall", NMF_BC::WALL),
		std::pair<std::string, NMF_BC>("wall", NMF_BC::WALL),
		// SYM
		std::pair<std::string, NMF_BC>("SYM", NMF_BC::SYM),
		std::pair<std::string, NMF_BC>("Sym", NMF_BC::SYM),
		std::pair<std::string, NMF_BC>("sym", NMF_BC::SYM),
		std::pair<std::string, NMF_BC>("SYMMETRY", NMF_BC::SYM),
		std::pair<std::string, NMF_BC>("Symmetry", NMF_BC::SYM),
		std::pair<std::string, NMF_BC>("symmetry", NMF_BC::SYM)
	};

	auto it = BC_ENUM.find(x);
	if (it == BC_ENUM.end())
		throw std::runtime_error("Unsupported B.C. name: " + x);
	else
		return it->second;
}

class NMF_Entry
{
private:
	NMF_BC m_bc;
	pNMF_Range m_rg1, m_rg2;
	bool m_swap;

public:
	NMF_Entry(const std::string &t, size_t *s) :
		m_rg1(new NMF_Range(s))
	{
		m_rg2 = nullptr;
		m_bc = NMF_BC_ENUM(t);
		m_swap = false;
	}

	NMF_Entry(const std::string &t, uint32_t *s1, uint32_t *s2, bool f) :
		m_rg1(new NMF_Range(s1)),
		m_rg2(new NMF_Range(s2))
	{
		m_bc = NMF_BC_ENUM(t);
		m_swap = f;
	}

	~NMF_Entry()
	{
		if (m_rg1)
			delete m_rg1;

		if (m_rg2)
			delete m_rg2;
	}

	NMF_BC Type() const
	{
		return m_bc;
	}

	size_t B1() const
	{
		return m_rg1->B();
	}

	size_t F1() const
	{
		return m_rg1->F();
	}

	size_t B2() const
	{
		if (m_rg2)
			return m_rg2->B();
		else
			return 0;
	}

	size_t F2() const
	{
		if (m_rg2)
			return m_rg2->F();
		else
			return 0;
	}

	pNMF_Range Range1()
	{
		return m_rg1;
	}

	pNMF_Range Range2()
	{
		return m_rg2;
	}

	bool Swap() const
	{
		return m_swap;
	}

	int contains(size_t bs, size_t fs, size_t lpri, size_t lsec)
	{
		// Left
		if (B1() == bs && F1() == fs && Range1()->constains(lpri, lsec))
			return 1;
		else if (B2() == bs && F2() == fs && Range2()->constains(lpri, lsec))
			return 2;
		else
			return 0;
	}

	vector<uint32_t> opposite_logical_desc(uint32_t bs, uint32_t fs, uint32_t lpri, uint32_t lsec, uint8_t side)
	{
		if (!m_rg2)
			throw "Internal Error!";

		vector<uint32_t> ret(4, 0);
		uint32_t pri_off, sec_off;

		//Calc Offset in each direction
		if (side == 1)
		{
			ret[0] = m_rg2->m_blk;
			ret[1] = m_rg2->m_face;
			pri_off = lpri - m_rg1->m_s1;
			sec_off = lsec - m_rg1->m_s2;
		}
		else if (side == 2)
		{
			ret[0] = m_rg1->m_blk;
			ret[1] = m_rg1->m_face;
			pri_off = lpri - m_rg2->m_s1;
			sec_off = lsec - m_rg2->m_s2;
		}
		else
			throw "Invalid side indication.";

		//Swap on necessary, using bit manipulation
		if (swap)
		{
			pri_off ^= sec_off;
			sec_off ^= pri_off;
			pri_off ^= sec_off;
		}

		//Mark opposite position
		if (side == 1)
		{
			ret[2] = m_rg2->m_s1 + pri_off;
			ret[3] = m_rg2->m_s2 + sec_off;
		}
		else
		{
			ret[2] = m_rg1->m_s1 + pri_off;
			ret[3] = m_rg1->m_s2 + sec_off;
		}

		return ret;
	}
};

#endif
