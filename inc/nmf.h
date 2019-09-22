#ifndef __NMF_H__
#define __NMF_H__

#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cctype>
#include <algorithm>
#include <cstddef>
#include <string>
#include <map>
#include <utility>
#include <stdexcept>

class NMF_Block
{
private:
	template<typename T>
	class Array1D : public std::vector<T>
	{
	public:
		Array1D(size_t n) : std::vector<T>(n) {}

		Array1D(size_t n, const T &val) : std::vector<T>(n, val) {}

		// 1-based indexing
		T &operator()(size_t i) { return std::vector<T>::at(i - 1); }
	};

	class HEX_CELL
	{
	private:
		size_t m_cell;
		size_t m_node[8];
		size_t m_face[6];

	public:
		HEX_CELL() : m_cell(0), m_node{ 0 }, m_face{ 0 } {}

		size_t CellSeq() const { return m_cell; }

		size_t &CellSeq() { return m_cell; }

		size_t NodeSeq(int n) const { return m_node[n - 1]; } // 1-based indexing

		size_t &NodeSeq(int n) { return m_node[n - 1]; } // 1-based indexing

		size_t FaceSeq(int n) const { return m_face[n - 1]; } // 1-based indexing

		size_t &FaceSeq(int n) { return m_face[n - 1]; } // 1-based indexing
	};

public:
	NMF_Block(size_t nI, size_t nJ) :
		m_hex((nI - 1)*(nJ - 1))
	{
		m_nI = nI;
		m_nJ = nJ;
		m_nK = 1;
	}

	NMF_Block(size_t nI, size_t nJ, size_t nK) :
		m_hex((nI - 1)*(nJ - 1)*(nK - 1))
	{
		m_nI = nI;
		m_nJ = nJ;
		m_nK = nK;
	}

	size_t IDIM() const { return m_nI; }

	size_t &IDIM() { return m_nI; }

	size_t JDIM() const { return m_nJ; }

	size_t &JDIM() { return m_nJ; }

	size_t KDIM() const { return m_nK; }

	size_t &KDIM() { return m_nK; }

	int dimension() const { return KDIM() == 1 ? 2 : 3; }

	bool is3D() const { return KDIM() > 1; }

	size_t node_num() const
	{
		return IDIM() * JDIM() * KDIM();
	}

	size_t cell_num() const
	{
	    if(is3D())
		    return (IDIM() - 1) * (JDIM() - 1) *(KDIM() - 1);
	    else
            return (IDIM() - 1) * (JDIM() - 1);
	}

	HEX_CELL &cell(size_t i, size_t j)
	{
		// Convert 1-based index to 0-based
		const size_t i0 = i - 1;
		const size_t j0 = j - 1;

		// Access
		return m_hex.at(i0 + (m_nI - 1) * j0);
	}

	HEX_CELL &cell(size_t i, size_t j, size_t k)
	{
		// Convert 1-based index to 0-based
		const size_t i0 = i - 1;
		const size_t j0 = j - 1;
		const size_t k0 = k - 1;

		// Access
		return m_hex.at(i0 + (m_nI - 1) * (j0 + (m_nJ - 1) * k0));
	}

private:
	size_t m_nI, m_nJ, m_nK;
	Array1D<HEX_CELL> m_hex;
};

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

	size_t B() const { return m_blk; }

	size_t &B() { return m_blk; }

	size_t F() const { return m_face; }

	size_t &F() { return m_face; }

	size_t S1() const { return m_s1; }

	size_t &S1() { return m_s1; }

	size_t E1() const { return m_e1; }

	size_t &E1() { return m_e1; }

	size_t S2() const { return m_s2; }

	size_t &S2() { return m_s2; }

	size_t E2() const { return m_e2; }

	size_t &E2() { return m_e2; }

	// Check if given index is within this range.
	bool constains(size_t pri, size_t sec) const
	{
		const bool t1 = (m_s1 <= pri) && (pri <= m_e1);
		const bool t2 = (m_s2 <= sec) && (sec <= m_e2);
		return t1 && t2;
	}

	// Nodes in primary direction.
	size_t pri_node_num() const { return m_e1 - m_s1 + 1; }

	// Nodes in secondary direction.
	size_t sec_node_num() const { return m_e2 - m_s2 + 1; }

	// Total nodes on this interface.
	size_t node_num() const { return pri_node_num() * sec_node_num(); }

	// Total faces/edges on this interface.
	size_t face_num() const
	{
		const size_t n_pri = (pri_node_num() - 1) * sec_node_num();
		const size_t n_sec = (sec_node_num() - 1) * pri_node_num();
		return n_pri + n_sec;
	}

	// Total quad cells on this interface.
	size_t cell_num() const { return (pri_node_num() - 1) * (sec_node_num() - 1); }
};

class NMF_BC
{
public:
	enum {
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
		SYM = 11,
		INFLOW = 12,
		OUTFLOW = 13
	};

	static const std::map<int, std::string> MAPPING_Idx2Str;

	static const std::map<std::string, int> MAPPING_Str2Idx;
};

class NMF_Entry
{
private:
	int m_bc;
	NMF_Range m_rg1, m_rg2;
	bool m_swap;

public:
	NMF_Entry(const std::string &t, size_t *s) :
		m_rg1(s),
		m_rg2(),
		m_swap(false)
	{
		// Check before assign
		auto it = NMF_BC::MAPPING_Str2Idx.find(t);
		if (it == NMF_BC::MAPPING_Str2Idx.end())
			throw std::runtime_error("Unsupported B.C. name: " + t);
		else
			m_bc = it->second;
	}

	NMF_Entry(const std::string &t, size_t *s1, size_t *s2, bool f) :
		m_rg1(s1),
		m_rg2(s2),
		m_swap(f)
	{
		// Check before assign
		auto it = NMF_BC::MAPPING_Str2Idx.find(t);
		if (it == NMF_BC::MAPPING_Str2Idx.end())
			throw std::runtime_error("Unsupported B.C. name: " + t);
		else
			m_bc = it->second;
	}

	int Type() const { return m_bc; }

	int &Type() { return m_bc; }

	size_t B1() const { return m_rg1.B(); }

	size_t &B1() { return m_rg1.B(); }

	size_t F1() const { return m_rg1.F(); }

	size_t &F1() { return m_rg1.F(); }

	size_t B2() const { return m_rg2.B(); }

	size_t &B2() { return m_rg2.B(); }

	size_t F2() const { return m_rg2.F(); }

	size_t &F2() { return m_rg2.F(); }

	NMF_Range &Range1() { return m_rg1; }

	NMF_Range &Range2() { return m_rg2; }

	bool Swap() const { return m_swap; }

	bool &Swap() { return m_swap; }

	int contains(size_t bs, size_t fs, size_t lpri, size_t lsec)
	{
		// Left
		if (B1() == bs && F1() == fs && Range1().constains(lpri, lsec))
			return 1;
		else if (B2() == bs && F2() == fs && Range2().constains(lpri, lsec))
			return 2;
		else
			return 0;
	}
};

class NMF
{
public:
	size_t nBlk() const { return m_blk.size(); }

	int readFromFile(const std::string &path);

	int writeToFile(const std::string &path);

	size_t nHex() const
	{
		size_t ret = 0;
		for (const auto & e : m_blk)
			ret += e.cell_num();
		return ret;
	}

	int compute_topology();

private:
	std::vector<NMF_Block> m_blk;
	std::vector<NMF_Entry> m_entry;
};

#endif
