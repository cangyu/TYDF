#ifndef __NMF_H__
#define __NMF_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <cstdint>
#include <unordered_map>
#include <cctype>
#include <algorithm>
#include <cstring>
#include <list>
#include <cstddef>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <stdexcept>

class NMF_Block
{
private:
	size_t m_nI;
	size_t m_nJ;
	size_t m_nK;

public:
	NMF_Block()
	{
		m_nI = 0;
		m_nJ = 0;
		m_nK = 0;
	}

	~NMF_Block() = default;

	size_t IDIM() const { return m_nI; }

	size_t &IDIM() { return m_nI; }

	size_t JDIM() const { return m_nJ; }

	size_t &JDIM() { return m_nJ; }

	size_t KDIM() const { return m_nK; }

	size_t &KDIM() { return m_nK; }
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

	~NMF_Range() = default;

	size_t B() const { return m_blk; }

	size_t F() const { return m_face; }

	size_t S1() const { return m_s1; }

	size_t E1() const { return m_e1; }

	size_t S2() const { return m_s2; }

	size_t E2() const { return m_e2; }

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

	NMF_BC() = default;

	~NMF_BC() = default;
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

	~NMF_Entry() = default;

	int Type() const { return m_bc; }

	size_t B1() const { return m_rg1.B(); }

	size_t F1() const { return m_rg1.F(); }

	size_t B2() const { return m_rg2.B(); }

	size_t F2() const { return m_rg2.F(); }

	NMF_Range &Range1() { return m_rg1; }

	NMF_Range &Range2() { return m_rg2; }

	bool Swap() const { return m_swap; }

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
private:
	// Raw content
	std::vector<NMF_Block> m_blk;
	std::vector<NMF_Entry> m_entry;

	// Topology info
	// TODO

public:
	NMF() = default;

	~NMF() = default;

	size_t nBlk() const { return m_blk.size(); }

	int readFromFile(const std::string &path);

	int writeToFile(const std::string &path);

	int compute_topology();
};

#endif
