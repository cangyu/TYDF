#ifndef __NMF_HPP__
#define __NMF_HPP__

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
#include <set>
#include <utility>
#include <stdexcept>

namespace NMF
{

template<typename T>
class Array1D : public std::vector<T>
{
public:
	Array1D(size_t n) : std::vector<T>(n) {}
	Array1D(size_t n, const T &val) : std::vector<T>(n, val) {}

	// 1-based indexing
	T &operator()(size_t i) { return std::vector<T>::at(i - 1); }
	const T &operator()(size_t i) const { return std::vector<T>::at(i - 1); }
};

class Block
{
private:
	class HEX_CELL
	{
	private:
		size_t m_cell; // 1-based sequence
		size_t m_node[8]; // 1-based sequence
		size_t m_face[6]; // 1-based sequence

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
	Block(size_t nI, size_t nJ, size_t nK) : m_hex((nI - 1)*(nJ - 1)*(nK - 1))
	{
		m_nI = nI;
		m_nJ = nJ;
		m_nK = nK;
	}

	size_t IDIM() const { return m_nI; }
	size_t JDIM() const { return m_nJ; }
	size_t KDIM() const { return m_nK; }
	
	size_t &IDIM() { return m_nI; }
	size_t &JDIM() { return m_nJ; }
	size_t &KDIM() { return m_nK; }

	size_t node_num() const { return IDIM() * JDIM() * KDIM(); }
	size_t face_num() const
	{
		size_t ret = 0;
		ret += (IDIM() - 1) * JDIM() * KDIM();
		ret += IDIM() * (JDIM() - 1) * KDIM();
		ret += IDIM() * JDIM() * (KDIM() - 1);
		return ret;
	}
	size_t cell_num() const { return (IDIM() - 1) * (JDIM() - 1) * (KDIM() - 1); }

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

class Range
{
private:
	size_t m_blk; // Block index, 1-based.
	size_t m_face; // Face index, ranges from 1 to 6.
	size_t m_s1; // Primary direction starting index, 1-based.
	size_t m_e1; // Primary direction ending index, 1-based.
	size_t m_s2; // Secondary direction starting index, 1-based.
	size_t m_e2; // Secondary direction ending index, 1-based.

public:
	Range()
	{
		m_blk = 0;
		m_face = 0;
		m_s1 = 0;
		m_e1 = 0;
		m_s2 = 0;
		m_e2 = 0;
	}

	Range(size_t *src)
	{
		m_blk = src[0];
		m_face = src[1];
		m_s1 = src[2];
		m_e1 = src[3];
		m_s2 = src[4];
		m_e2 = src[5];
	}

	Range(size_t b, size_t f, size_t s1, size_t e1, size_t s2, size_t e2)
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

	// Total edges on this interface.
	size_t edge_num() const
	{
		const size_t n_pri = (pri_node_num() - 1) * sec_node_num();
		const size_t n_sec = (sec_node_num() - 1) * pri_node_num();
		return n_pri + n_sec;
	}

	// Total quad cells on this interface.
	size_t face_num() const { return (pri_node_num() - 1) * (sec_node_num() - 1); }
};

class BC
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

	static bool isValidBCIdx(int x)
	{
		static const std::set<int> candidate_set{ 
			COLLAPSED, ONE_TO_ONE, PATCHED, 
			POLE_DIR1, POLE_DIR2, 
			SYM_X, SYM_Y, SYM_Z, 
			UNPROCESSED, WALL, SYM, 
			INFLOW, OUTFLOW 
		};

		return candidate_set.find(x) != candidate_set.end();
	}

	static bool isValidBCStr(const std::string &x)
	{
		static const std::set<std::string> candidate_set{ 
			"COLLAPSED", "ONE_TO_ONE", "PATCHED", 
			"POLE_DIR1", "POLE_DIR2", 
			"SYM_X", "SYM_Y", "SYM_Z",
			 "UNPROCESSED", "WALL", "SYM", 
			 "INFLOW", "OUTFLOW" 
		};

		return candidate_set.find(x) != candidate_set.end();
	}

	static const std::string &idx2str(int x)
	{
		static const std::map<int, std::string> mapping_set{
			std::pair<int, std::string>(COLLAPSED, "COLLAPSED"),
			std::pair<int, std::string>(ONE_TO_ONE, "ONE_TO_ONE"),
			std::pair<int, std::string>(PATCHED, "PATCHED"),
			std::pair<int, std::string>(POLE_DIR1, "POLE_DIR1"),
			std::pair<int, std::string>(POLE_DIR2, "POLE_DIR2"),
			std::pair<int, std::string>(SYM_X, "SYM_X"),
			std::pair<int, std::string>(SYM_Y, "SYM_Y"),
			std::pair<int, std::string>(SYM_Z, "SYM_Z"),
			std::pair<int, std::string>(UNPROCESSED, "UNPROCESSED"),
			std::pair<int, std::string>(WALL, "WALL"),
			std::pair<int, std::string>(SYM, "SYM"),
			std::pair<int, std::string>(INFLOW, "INFLOW"),
			std::pair<int, std::string>(OUTFLOW, "OUTFLOW")
		};

		auto it = mapping_set.find(x);
		if(it == mapping_set.end())
			throw std::runtime_error("Not found!");
		else
			return it->second;
	}

	static int str2idx(const std::string &x)
	{
		static const std::map<std::string, int> mapping_set{
			std::pair<std::string, int>("COLLAPSED", COLLAPSED),
			std::pair<std::string, int>("ONE_TO_ONE", ONE_TO_ONE),
			std::pair<std::string, int>("PATCHED", PATCHED),
			std::pair<std::string, int>("POLE_DIR1", POLE_DIR1),
			std::pair<std::string, int>("POLE_DIR2", POLE_DIR2),
			std::pair<std::string, int>("SYM_X", SYM_X),
			std::pair<std::string, int>("SYM_Y", SYM_Y),
			std::pair<std::string, int>("SYM_Z", SYM_Z),
			std::pair<std::string, int>("UNPROCESSED", UNPROCESSED),
			std::pair<std::string, int>("WALL", WALL),
			std::pair<std::string, int>("SYM", SYM),
			std::pair<std::string, int>("SYMMETRY", SYM),
			std::pair<std::string, int>("INFLOW", INFLOW),
			std::pair<std::string, int>("OUTFLOW", OUTFLOW)
		};

		auto it = mapping_set.find(x);
		if(it == mapping_set.end())
			throw std::runtime_error("Not found!");
		else
			return it->second;
	}

	BC() = delete;

	~BC() = default;
};

class Entry
{
private:
	int m_bc;
	Range m_rg1, m_rg2;
	bool m_swap;

public:
	Entry(const std::string &t, size_t *s) :
		m_rg1(s),
		m_rg2(),
		m_swap(false)
	{
		// Check before assign
		auto it = BC::MAPPING_Str2Idx.find(t);
		if (it == BC::MAPPING_Str2Idx.end())
			throw std::runtime_error("Unsupported B.C. name: " + t);
		else
			m_bc = it->second;
	}

	Entry(const std::string &t, size_t *s1, size_t *s2, bool f) :
		m_rg1(s1),
		m_rg2(s2),
		m_swap(f)
	{
		// Check before assign
		auto it = BC::MAPPING_Str2Idx.find(t);
		if (it == BC::MAPPING_Str2Idx.end())
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

	Range &Range1() { return m_rg1; }

	Range &Range2() { return m_rg2; }

	bool Swap() const { return m_swap; }

	bool &Swap() { return m_swap; }

	size_t node_num() const { return m_rg1.node_num(); }

	size_t face_num() const { return m_rg1.face_num(); }

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

class MAPPING
{
public:
	MAPPING() = default;

	~MAPPING() = default;

	size_t nBlk() const { return m_blk.size(); }

	int readFromFile(const std::string &path)
	{
		std::string s;
		std::stringstream ss;

		//Load file
		std::ifstream mfp(path);
		if (mfp.fail())
			throw std::runtime_error("Can not open target input file: " + path);

		//Skip header
		do {
			std::getline(mfp, s, '\n');
		} while (s.find('#') != std::string::npos);

		//Read block dimension info
		size_t NumOfBlk;
		ss << s;
		ss >> NumOfBlk;
		if (NumOfBlk == 0)
			throw std::runtime_error("Invalid num of blocks: " + std::to_string(NumOfBlk));

		m_blk.clear();
		for (size_t i = 1; i <= NumOfBlk; i++)
		{
			size_t idx, i_max, j_max, k_max;
			std::getline(mfp, s, '\n');
			ss.clear();
			ss << s;
			ss >> idx >> i_max >> j_max >> k_max;

			if (idx != i)
				throw std::runtime_error("Invalid order of block: " + std::to_string(idx));
			if (i_max == 0)
				throw std::runtime_error("Invalid I dimension: " + std::to_string(i_max));
			if (j_max == 0)
				throw std::runtime_error("Invalid J dimension: " + std::to_string(j_max));
			if (k_max == 0)
				throw std::runtime_error("Invalid K dimension: " + std::to_string(k_max));

			m_blk.emplace_back(i_max, j_max, k_max);
		}

		//Skip separators
		do {
			std::getline(mfp, s, '\n');
		} while (s.find('#') != std::string::npos);

		//Read bc mappings
		m_entry.clear();
		do {
			transform(s.begin(), s.end(), s.begin(), ::toupper);
			ss.clear();
			ss << s;
			std::string bc_str;
			ss >> bc_str;
			size_t connectivity[2][6] = { 0 };
			if (BC::MAPPING_Str2Idx.at(bc_str) == BC::ONE_TO_ONE)
			{
				for (int i = 0; i < 2; i++)
					for (int j = 0; j < 6; j++)
						ss >> connectivity[i][j];

				std::string swp;
				ss >> swp;

				m_entry.emplace_back(bc_str, connectivity[0], connectivity[1], swp == "TRUE");
			}
			else
			{
				for (int i = 0; i < 6; i++)
					ss >> connectivity[0][i];

				m_entry.emplace_back(bc_str, connectivity[0]);
			}
		} while (std::getline(mfp, s, '\n'));

		// Finalize
		mfp.close();
		return 0;
	}

	int writeToFile(const std::string &path)
	{
		// Open target file
		std::ofstream f_out(path);
		if (f_out.fail())
			throw std::runtime_error("Can not open target output file: " + path);

		// Header
		f_out << "# ======================== Neutral Map File generated by the Grid-Glue software ==============================" << std::endl;
		f_out << "# ============================================================================================================" << std::endl;
		f_out << "# Block#    IDIM    JDIM    KDIM" << std::endl;
		f_out << "# ------------------------------------------------------------------------------------------------------------" << std::endl;

		// Block info
		const size_t NumOfBlk = nBlk();
		f_out << std::setw(8) << std::right << NumOfBlk << std::endl;
		for (size_t i = 0; i < NumOfBlk; i++)
		{
			f_out << std::setw(8) << std::right << i + 1;
			f_out << std::setw(8) << std::right << m_blk[i].IDIM();
			f_out << std::setw(8) << std::right << m_blk[i].JDIM();
			f_out << std::setw(8) << std::right << m_blk[i].KDIM();
			f_out << std::endl;
		}

		// Interface and boundary info
		f_out << "# ============================================================================================================" << std::endl;
		f_out << "# Type           B1    F1       S1    E1       S2    E2       B2    F2       S1    E1       S2    E2      Swap" << std::endl;
		f_out << "# ------------------------------------------------------------------------------------------------------------" << std::endl;
		for (auto & e : m_entry)
		{
			f_out << std::setw(13) << std::left << BC::MAPPING_Idx2Str.at(e.Type());
			f_out << std::setw(6) << std::right << e.B1();
			f_out << std::setw(6) << std::right << e.F1();
			f_out << std::setw(9) << std::right << e.Range1().S1();
			f_out << std::setw(6) << std::right << e.Range1().E1();
			f_out << std::setw(9) << std::right << e.Range1().S2();
			f_out << std::setw(6) << std::right << e.Range1().E2();
			if (e.Type() == BC::ONE_TO_ONE)
			{
				f_out << std::setw(9) << std::right << e.B2();
				f_out << std::setw(6) << std::right << e.F2();
				f_out << std::setw(9) << std::right << e.Range2().S1();
				f_out << std::setw(6) << std::right << e.Range2().E1();
				f_out << std::setw(9) << std::right << e.Range2().S2();
				f_out << std::setw(6) << std::right << e.Range2().E2();
				f_out << std::setw(10) << std::right << (e.Swap() ? "TRUE" : "FALSE");
			}
			f_out << std::endl;
		}

		// Finalize
		f_out.close();
		return 0;
	}
	size_t nHex() const
	{
		size_t ret = 0;
		for (const auto & blk : m_blk)
			ret += blk.cell_num();
		return ret;
	}

	size_t nQuad() const
	{
		size_t ret = 0;
		for (const auto & blk : m_blk)
			ret += blk.cell_num();

		// Sustract duplicated interface
		for (const auto & e : m_entry)
			if (e.Type() == BC::ONE_TO_ONE)
				ret -= e.face_num();

		return ret;
	}

	int compute_topology()
	{
		// Indexing of cells
		size_t cnt = 0;
		for (auto &blk : m_blk)
		{
			const size_t cI = blk.IDIM();
			const size_t cJ = blk.JDIM();
			const size_t cK = blk.KDIM();

			for (size_t k = 1; k < cK; ++k)
				for (size_t j = 1; j < cJ; ++j)
					for (size_t i = 1; i < cI; ++i)
						blk.cell(i, j, k).CellSeq() = ++cnt;
		}

		const size_t totalCellNum = nHex();
		assert(cnt == totalCellNum);

		// Indexing of faces
		cnt = 0;
		for (auto &blk : m_blk)
		{
			const size_t cI = blk.IDIM();
			const size_t cJ = blk.JDIM();
			const size_t cK = blk.KDIM();

			

		}
		const size_t totalFaceNum = nQuad();
		//assert(cnt == totalFaceNum);

		// Indexing of nodes
		cnt = 0;
		// TODO

		return 0;
	}

private:
	std::vector<Block> m_blk;
	std::vector<Entry> m_entry;
};

}

#endif
