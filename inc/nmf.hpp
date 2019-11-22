#ifndef __NMF_HPP__
#define __NMF_HPP__

#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cctype>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <utility>
#include <stdexcept>
#include <regex>

namespace NMF
{
	template<typename T>
	class Array1D : public std::vector<T>
	{
	public:
		Array1D() : std::vector<T>() {}
		Array1D(size_t n) : std::vector<T>(n) {}
		Array1D(size_t n, const T &val) : std::vector<T>(n, val) {}
		Array1D(const Array1D &rhs) = default;
		~Array1D() = default;

		// 1-based indexing
		T &operator()(long long i)
		{
			const long long N = this->size();

			if (1 <= i && i <= N)
				return this->at(i - 1);
			else if (-N <= i && i <= -1)
				return this->at(N + i);
			else
				throw std::invalid_argument("\"" + std::to_string(i) + "\" is not a valid 1-based index.");
		}
		const T &operator()(long long i) const
		{
			const long long N = this->size();

			if (1 <= i && i <= N)
				return this->at(i - 1);
			else if (-N <= i && i <= -1)
				return this->at(N + i);
			else
				throw std::invalid_argument("\"" + std::to_string(i) + "\" is not a valid 1-based index.");
		}
	};

	class BC
	{
	public:
		static void str_formalize(std::string &s)
		{
			std::transform(s.begin(), s.end(), s.begin(), ::toupper);
			for (auto &e : s)
				if (e == '-')
					e = '_';
		}

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
				COLLAPSED,
				ONE_TO_ONE,
				PATCHED,
				POLE_DIR1,
				POLE_DIR2,
				SYM_X,
				SYM_Y,
				SYM_Z,
				UNPROCESSED,
				WALL,
				SYM,
				INFLOW,
				OUTFLOW
			};

			return candidate_set.find(x) != candidate_set.end();
		}

		static bool isValidBCStr(const std::string &x)
		{
			static const std::set<std::string> candidate_set{
				"COLLAPSED",
				"ONE_TO_ONE",
				"PATCHED",
				"POLE_DIR1", "POLE_DIR2",
				"SYM_X", "SYM_Y", "SYM_Z",
				"UNPROCESSED",
				"WALL",
				"SYM",
				"INFLOW",
				"OUTFLOW"
			};

			std::string x_(x);
			str_formalize(x_);
			return candidate_set.find(x_) != candidate_set.end();
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
			if (it == mapping_set.end())
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

			std::string x_(x);
			str_formalize(x_);

			auto it = mapping_set.find(x_);
			if (it == mapping_set.end())
				throw std::runtime_error("Not found!");
			else
				return it->second;
		}

		BC() = delete;
		BC(const BC &rhs) = delete;
		~BC() = default;
	};

	class CELL
	{
	protected:
		size_t m_cell; // 1-based sequence
		Array1D<size_t> m_node; // 1-based sequence
		Array1D<size_t> m_face; // 1-based sequence

	public:
		CELL() = delete;
		CELL(bool is3D) : m_cell(0), m_node((is3D ? 8 : 4), 0), m_face((is3D ? 6 : 4), 0) {}
		CELL(const CELL &rhs) : m_cell(rhs.m_cell), m_node(rhs.m_node), m_face(rhs.m_face) {}
		virtual ~CELL() = default;

		size_t CellSeq() const { return m_cell; }
		size_t &CellSeq() { return m_cell; }

		// 1-based indexing of node
		size_t NodeSeq(int n) const { return m_node.at(n - 1); }
		size_t &NodeSeq(int n) { return m_node.at(n - 1); }

		// 1-based indexing of face
		size_t FaceSeq(int n) const { return m_face.at(n - 1); }
		size_t &FaceSeq(int n) { return m_face.at(n - 1); }
	};

	class QUAD_CELL : public CELL
	{
	public:
		QUAD_CELL() : CELL(false) {}
		QUAD_CELL(const QUAD_CELL &rhs) = default;
		~QUAD_CELL() = default;
	};

	class HEX_CELL : public CELL
	{
	public:
		HEX_CELL() : CELL(true) {}
		HEX_CELL(const HEX_CELL &rhs) = default;
		~HEX_CELL() = default;
	};

	class BLOCK
	{
	protected:
		size_t idx; // 1-based global index
		size_t m_nI, m_nJ, m_nK; // Dimensions

	public:
		BLOCK() = delete;
		BLOCK(int nI, int nJ) :
			idx(0),
			m_nI(nI),
			m_nJ(nJ),
			m_nK(1)
		{
			if (nI < 2 || nJ < 2)
				throw std::invalid_argument("Invalid dimension.");
		}
		BLOCK(int nI, int nJ, int nK) :
			idx(0),
			m_nI(nI),
			m_nJ(nJ),
			m_nK(nK)
		{
			if (nI < 2 || nJ < 2 || nK < 2)
				throw std::invalid_argument("Invalid dimension.");
		}
		BLOCK(const BLOCK &rhs) = default;
		virtual ~BLOCK() = default;

		size_t index() const { return idx; }
		size_t &index() { return idx; }

		size_t IDIM() const { return m_nI; }
		size_t JDIM() const { return m_nJ; }
		size_t KDIM() const { return m_nK; }

		bool is3D() const
		{
			if (m_nK == 1)
				return false;
			else
				return true;
		}
		short dimension() const
		{
			return is3D() ? 3 : 2;
		}

		// Here the terms 'node', 'face' and 'cell' are 
		// consistent with that in ANSYS Fluent convention.
		size_t node_num() const
		{
			if (is3D())
				return IDIM() * JDIM() * KDIM();
			else
				return IDIM() * JDIM();
		}
		size_t face_num() const
		{
			size_t ret = 0;
			if (is3D())
			{
				ret += (IDIM() - 1) * (JDIM() - 1) * KDIM();
				ret += IDIM() * (JDIM() - 1) * (KDIM() - 1);
				ret += (IDIM() - 1) *JDIM() * (KDIM() - 1);
			}
			else
			{
				ret += (IDIM() - 1) * JDIM();
				ret += IDIM() * (JDIM() - 1);
			}
			return ret;
		}
		size_t cell_num() const
		{
			if (is3D())
				return (IDIM() - 1) * (JDIM() - 1) * (KDIM() - 1);
			else
				return (IDIM() - 1) * (JDIM() - 1);
		}

		size_t block_internal_node_num() const
		{
			if (is3D())
				return (IDIM() - 2) * (JDIM() - 2) * (KDIM() - 2);
			else
				return (IDIM() - 2) * (JDIM() - 2);
		}
		size_t surface_internal_node_num(int idx) const
		{
			size_t ret = 0;
			if (is3D())
			{
				switch (idx)
				{
				case 1:
				case 2:
					ret = (IDIM() - 2) * (JDIM() - 2);
					break;
				case 3:
				case 4:
					ret = (JDIM() - 2) * (KDIM() - 2);
					break;
				case 5:
				case 6:
					ret = (KDIM() - 2) * (IDIM() - 2);
					break;
				default:
					throw std::invalid_argument("\"" + std::to_string(idx) + "\" is not a valid surface index of a 3D block.");
				}
			}
			else
			{
				switch (idx)
				{
				case 1:
				case 2:
					ret = (JDIM() - 2);
					break;
				case 3:
				case 4:
					ret = (IDIM() - 2);
					break;
				default:
					throw std::invalid_argument("\"" + std::to_string(idx) + "\" is not a valid surface index of a 2D block.");
				}
			}
			return ret;
		}
	};

	class Block2D : public BLOCK
	{
	public:
		static const short NumOfVertex = 4;
		static const short NumOfFrame = 4;

		struct FRAME
		{
			short local_index = 0; // Ranges from 1 to 'NumOfFrame', set to 0 when uninitialized.
			int global_index = 0; // 1-based global index, set to 0 when uninitialized.
			Block2D *dependentBlock = nullptr;
			Block2D *neighbourBlock = nullptr;
		};
		struct VERTEX
		{
			short local_index = 0; // Ranges from 1 to 'NumOfVertex', set to 0 when uninitialized.
			int global_index = 0; // 1-based global index, set to 0 when uninitialized.
			Block2D *dependentBlock = nullptr;
			std::array<FRAME*, 2> dependentFrame{ nullptr, nullptr };
		};

	private:
		Array1D<QUAD_CELL> m_cell;
		Array1D<FRAME> m_frame;
		Array1D<VERTEX> m_vertex;

	public:
		Block2D() = delete;
		Block2D(int nI, int nJ) :
			BLOCK(nI, nJ),
			m_cell(cell_num()),
			m_vertex(NumOfVertex),
			m_frame(NumOfFrame)
		{
			for (size_t i = 0; i < m_frame.size(); ++i)
				m_frame[i].local_index = i + 1;

			for (size_t i = 0; i < m_vertex.size(); ++i)
				m_vertex[i].local_index = i + 1;
		}
		Block2D(const Block2D &rhs) = default;
		~Block2D() = default;

		// Access internal cell through 1-based index.
		// Indexing convention:
		//      "i" ranges from 1 to IDIM()-1;
		//      "j" ranges from 1 to JDIM()-1;
		// When the IJ-axis follows the right-hand convention, (i, j) represents
		// the left-most, bottom-most node of the selected cell.
		QUAD_CELL &cell(size_t i, size_t j)
		{
			const size_t i0 = i - 1, j0 = j - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (m_nI - 1) * j0;
			return m_cell.at(idx);
		}
		const QUAD_CELL &cell(size_t i, size_t j) const
		{
			const size_t i0 = i - 1, j0 = j - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (m_nI - 1) * j0;
			return m_cell.at(idx);
		}

		// Access frame through 1-based index.
		// The indexing convention follows OpenFOAM specification.
		FRAME &frame(int n)
		{
			if (1 <= n && n <= NumOfFrame)
				return m_frame.at(n - 1);
			else if (-NumOfFrame <= n && n <= -1)
				return m_frame.at(NumOfFrame + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid frame index for a 2D block.");
		}
		const FRAME &frame(int n) const
		{
			if (1 <= n && n <= NumOfFrame)
				return m_frame.at(n - 1);
			else if (-NumOfFrame <= n && n <= -1)
				return m_frame.at(NumOfFrame + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid frame index for a 2D block.");
		}

		// Access vertex through 1-based index.
		// The indexing convention follows OpenFOAM specification.
		VERTEX &vertex(int n)
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid vertex index for a 2D block.");
		}
		const VERTEX &vertex(int n) const
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid vertex index for a 2D block.");
		}
	};

	class Block3D : public BLOCK
	{
	public:
		static const short NumOfVertex = 8;
		static const short NumOfFrame = 12;
		static const short NumOfSurf = 6;

		struct SURF;
		struct FRAME;
		struct VERTEX
		{
			short local_index = 0; // Ranges from 1 to 8, set to 0 when uninitialized.
			int global_index = 0; // 1-based global index, set to 0 when uninitialized.
			Block3D *dependentBlock = nullptr;
			std::array<SURF*, 3> dependentSurf = { nullptr, nullptr, nullptr };
			std::array<FRAME*, 3> dependentFrame = { nullptr, nullptr, nullptr };
		};
		struct FRAME
		{
			short local_index = 0; // Ranges from 1 to 12, set to 0 when uninitialized.
			int global_index = 0; // 1-based global index, set to 0 when uninitialized.
			Block3D *dependentBlock = nullptr;
			std::array<SURF*, 2> dependentSurf{ nullptr, nullptr };
			std::array<VERTEX*, 2> includedVertex{ nullptr, nullptr };
		};
		struct SURF
		{
			short local_index = 0; // Ranges from 1 to 6, set to 0 when uninitialized.
			int global_index = 0; // 1-based global index, set to 0 when uninitialized.
			std::string name = "";
			Block3D *dependentBlock = nullptr;
			SURF *neighbourSurf = nullptr;
			std::array<FRAME*, 4> includedFrame{ nullptr, nullptr, nullptr, nullptr };
			std::array<FRAME*, 4> counterpartFrame{ nullptr, nullptr, nullptr, nullptr };
			std::array<VERTEX*, 4> includedVertex{ nullptr, nullptr, nullptr, nullptr };
			std::array<VERTEX*, 4> counterpartVertex{ nullptr, nullptr, nullptr, nullptr };
		};

	private:
		Array1D<HEX_CELL> m_cell;
		Array1D<VERTEX> m_vertex;
		Array1D<FRAME> m_frame;
		Array1D<SURF> m_surf;

	public:
		Block3D() = delete;
		Block3D(int nI, int nJ, int nK) :
			BLOCK(nI, nJ, nK),
			m_cell(cell_num()),
			m_vertex(NumOfVertex),
			m_frame(NumOfFrame),
			m_surf(NumOfSurf)
		{
			for (short i = 0; i < NumOfVertex; ++i)
			{
				auto &v = m_vertex[i];
				v.local_index = i + 1;
				v.dependentBlock = this;
			}
			for (short i = 0; i < NumOfFrame; ++i)
			{
				auto &e = m_frame[i];
				e.local_index = i + 1;
				e.dependentBlock = this;
			}
			for (short i = 0; i < NumOfSurf; ++i)
			{
				auto &s = m_surf[i];
				s.local_index = i + 1;
				s.dependentBlock = this;
			}
			establish_connections();
		}
		Block3D(const Block3D &rhs) = default;
		~Block3D() = default;

		// Access internal cell through 1-based index.
		// Indexing convention:
		//      "i" ranges from 1 to IDIM()-1;
		//      "j" ranges from 1 to JDIM()-1;
		//      "k" ranges from 1 to KDIM()-1;
		// When the IJK-axis follows the right-hand convention, (i, j, k) represents
		// the left-most, bottom-most and back-most node of the selected cell.
		HEX_CELL &cell(size_t i, size_t j, size_t k)
		{
			const size_t i0 = i - 1, j0 = j - 1, k0 = k - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (m_nI - 1) * (j0 + (m_nJ - 1) * k0);
			return m_cell.at(idx);
		}
		const HEX_CELL &cell(size_t i, size_t j, size_t k) const
		{
			const size_t i0 = i - 1, j0 = j - 1, k0 = k - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (m_nI - 1) * (j0 + (m_nJ - 1) * k0);
			return m_cell.at(idx);
		}

		// Access vertex through 1-based index.
		// The indexing convention follows OpenFOAM specification.
		VERTEX &vertex(int n)
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid vertex index for a 3D block.");
		}
		const VERTEX &vertex(int n) const
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid vertex index for a 3D block.");
		}

		// Access the frame edges through 1-based index.
		// The indexing convention follows OpenFOAM specification.
		FRAME &frame(int n)
		{
			if (1 <= n && n <= NumOfFrame)
				return m_frame.at(n - 1);
			else if (-NumOfFrame <= n && n <= -1)
				return m_frame.at(NumOfFrame + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid frame index for a 3D block.");
		}
		const FRAME &frame(int n) const
		{
			if (1 <= n && n <= NumOfFrame)
				return m_frame.at(n - 1);
			else if (-NumOfFrame <= n && n <= -1)
				return m_frame.at(NumOfFrame + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid frame index for a 3D block.");
		}

		// Access surrounding surface through 1-based index.
		// The index follows NMF convection.
		SURF &surf(int n)
		{
			if (1 <= n && n <= NumOfSurf)
				return m_surf.at(n - 1);
			else if (-NumOfSurf <= n && n <= -1)
				return m_surf.at(NumOfSurf + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid 1-based surface index for a block.");
		}
		const SURF &surf(int n) const
		{
			if (1 <= n && n <= NumOfSurf)
				return m_surf.at(n - 1);
			else if (-NumOfSurf <= n && n <= -1)
				return m_surf.at(NumOfSurf + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid 1-based surface index for a block.");
		}

		size_t frame_node_num(int idx) const
		{
			switch (idx - 1)
			{
			case 4:
			case 5:
			case 6:
			case 7:
				return IDIM();
			case 8:
			case 9:
			case 10:
			case 11:
				return JDIM();
			case 0:
			case 1:
			case 2:
			case 3:
				return KDIM();
			default:
				throw std::invalid_argument("Invalid frame index.");
			}
		}
		size_t frame_internal_node_num(int idx) const
		{
			const auto n = frame_node_num(idx);
			if (n >= 2)
				return n - 2;
			else
				throw std::runtime_error("Internal error: too less nodes not detected when constructing.");
		}

	private:
		void establish_connections()
		{
			// Connection between frame edges and surrounding surfaces.
			// The order matters a lot!
			m_surf[0].includedFrame = { &m_frame[4], &m_frame[11], &m_frame[7], &m_frame[8] };
			m_surf[1].includedFrame = { &m_frame[5], &m_frame[10], &m_frame[6], &m_frame[9] };
			m_surf[2].includedFrame = { &m_frame[8], &m_frame[3], &m_frame[9], &m_frame[0] };
			m_surf[3].includedFrame = { &m_frame[11], &m_frame[2], &m_frame[10], &m_frame[1] };
			m_surf[4].includedFrame = { &m_frame[0], &m_frame[5], &m_frame[1], &m_frame[4] };
			m_surf[5].includedFrame = { &m_frame[3], &m_frame[6], &m_frame[2], &m_frame[7] };

			m_frame[0].dependentSurf = { &m_surf[2], &m_surf[4] };
			m_frame[1].dependentSurf = { &m_surf[4], &m_surf[3] };
			m_frame[2].dependentSurf = { &m_surf[3], &m_surf[5] };
			m_frame[3].dependentSurf = { &m_surf[5], &m_surf[2] };
			m_frame[4].dependentSurf = { &m_surf[0], &m_surf[4] };
			m_frame[5].dependentSurf = { &m_surf[4], &m_surf[1] };
			m_frame[6].dependentSurf = { &m_surf[1], &m_surf[5] };
			m_frame[7].dependentSurf = { &m_surf[5], &m_surf[0] };
			m_frame[8].dependentSurf = { &m_surf[0], &m_surf[2] };
			m_frame[9].dependentSurf = { &m_surf[2], &m_surf[1] };
			m_frame[10].dependentSurf = { &m_surf[1], &m_surf[3] };
			m_frame[11].dependentSurf = { &m_surf[3], &m_surf[0] };

			// Connection between surrounding surfaces and block vertexes.
			// The order matters a lot!
			m_surf[0].includedVertex = { &m_vertex[0], &m_vertex[3], &m_vertex[7], &m_vertex[4] };
			m_surf[1].includedVertex = { &m_vertex[1], &m_vertex[2], &m_vertex[6], &m_vertex[5] };
			m_surf[2].includedVertex = { &m_vertex[0], &m_vertex[4], &m_vertex[5], &m_vertex[1] };
			m_surf[3].includedVertex = { &m_vertex[3], &m_vertex[7], &m_vertex[6], &m_vertex[2] };
			m_surf[4].includedVertex = { &m_vertex[0], &m_vertex[1], &m_vertex[2], &m_vertex[3] };
			m_surf[5].includedVertex = { &m_vertex[4], &m_vertex[5], &m_vertex[6], &m_vertex[7] };

			m_vertex[0].dependentSurf = { &m_surf[4], &m_surf[0], &m_surf[2] };
			m_vertex[1].dependentSurf = { &m_surf[4], &m_surf[2], &m_surf[1] };
			m_vertex[2].dependentSurf = { &m_surf[4], &m_surf[1], &m_surf[3] };
			m_vertex[3].dependentSurf = { &m_surf[4], &m_surf[3], &m_surf[0] };
			m_vertex[4].dependentSurf = { &m_surf[5], &m_surf[0], &m_surf[2] };
			m_vertex[5].dependentSurf = { &m_surf[5], &m_surf[2], &m_surf[1] };
			m_vertex[6].dependentSurf = { &m_surf[5], &m_surf[1], &m_surf[3] };
			m_vertex[7].dependentSurf = { &m_surf[5], &m_surf[3], &m_surf[0] };

			// Connection between frame edges and block vertexes.
			m_frame[0].includedVertex = { &m_vertex[0], &m_vertex[1] };
			m_frame[1].includedVertex = { &m_vertex[3], &m_vertex[2] };
			m_frame[2].includedVertex = { &m_vertex[7], &m_vertex[6] };
			m_frame[3].includedVertex = { &m_vertex[4], &m_vertex[5] };
			m_frame[4].includedVertex = { &m_vertex[0], &m_vertex[3] };
			m_frame[5].includedVertex = { &m_vertex[1], &m_vertex[2] };
			m_frame[6].includedVertex = { &m_vertex[5], &m_vertex[6] };
			m_frame[7].includedVertex = { &m_vertex[4], &m_vertex[7] };
			m_frame[8].includedVertex = { &m_vertex[0], &m_vertex[4] };
			m_frame[9].includedVertex = { &m_vertex[1], &m_vertex[5] };
			m_frame[10].includedVertex = { &m_vertex[2], &m_vertex[6] };
			m_frame[11].includedVertex = { &m_vertex[3], &m_vertex[7] };

			m_vertex[0].dependentFrame = { &m_frame[8], &m_frame[4], &m_frame[0] };
			m_vertex[1].dependentFrame = { &m_frame[9], &m_frame[0], &m_frame[5] };
			m_vertex[2].dependentFrame = { &m_frame[10], &m_frame[5], &m_frame[1] };
			m_vertex[3].dependentFrame = { &m_frame[11], &m_frame[1], &m_frame[4] };
			m_vertex[4].dependentFrame = { &m_frame[8], &m_frame[7], &m_frame[3] };
			m_vertex[5].dependentFrame = { &m_frame[9], &m_frame[3], &m_frame[6] };
			m_vertex[6].dependentFrame = { &m_frame[10], &m_frame[6], &m_frame[2] };
			m_vertex[7].dependentFrame = { &m_frame[11], &m_frame[2], &m_frame[7] };
		}
	};

	class Mapping2D
	{
	public:
		Mapping2D() = default;
		Mapping2D(const Mapping2D &rhs) = default;
		~Mapping2D() = default;
	};

	class Mapping3D
	{
	protected:
		class ENTRY
		{
		protected:
			class RANGE
			{
			private:
				size_t m_blk; // Block index, 1-based.
				size_t m_face; // Face index, ranges from 1 to 6.
				size_t m_s1; // Primary direction starting index, 1-based.
				size_t m_e1; // Primary direction ending index, 1-based.
				size_t m_s2; // Secondary direction starting index, 1-based.
				size_t m_e2; // Secondary direction ending index, 1-based.

			public:
				RANGE() : m_blk(0), m_face(0), m_s1(0), m_e1(0), m_s2(0), m_e2(0) {}
				RANGE(size_t *src) :
					m_blk(src[0]),
					m_face(src[1]),
					m_s1(src[2]),
					m_e1(src[3]),
					m_s2(src[4]),
					m_e2(src[5])
				{
					check_param();
				}
				RANGE(size_t b, size_t f, size_t s1, size_t e1, size_t s2, size_t e2) :
					m_blk(b),
					m_face(f),
					m_s1(s1),
					m_e1(e1),
					m_s2(s2),
					m_e2(e2)
				{
					check_param();
				}
				RANGE(const RANGE &rhs) = default;
				~RANGE() = default;

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
				size_t pri_node_num() const
				{
					return m_e1 - m_s1 + 1;
				}

				// Nodes in secondary direction.
				size_t sec_node_num() const
				{
					return m_e2 - m_s2 + 1;
				}

				// True - Asscending;
				// False - Descending.
				bool pri_trend() const
				{
					return E1() > S1();
				}

				// True - Asscending;
				// False - Descending.
				bool sec_trend() const
				{
					return E2() > S2();
				}

				// Total nodes on this interface.
				size_t node_num() const
				{
					return pri_node_num() * sec_node_num();
				}

				// Total edges on this interface.
				size_t edge_num() const
				{
					const size_t n_pri = (pri_node_num() - 1) * sec_node_num();
					const size_t n_sec = (sec_node_num() - 1) * pri_node_num();
					return n_pri + n_sec;
				}

				// Total quad cells on this interface.
				size_t face_num() const
				{
					return (pri_node_num() - 1) * (sec_node_num() - 1);
				}

			private:
				void check_param()
				{
					if (B() == 0)
						throw std::invalid_argument("Invalid block index, must be positive.");
					if (F() == 0)
						throw std::invalid_argument("Invalid face index, must be positive.");
					if (S1() == E1())
						throw std::invalid_argument("Starting index and ending index in primary direction must differ,");
					if (S2() == E2())
						throw std::invalid_argument("Starting index and ending index in secondary direction must differ,");
				}
			};

		private:
			int m_bc;
			RANGE m_rg1;

		public:
			ENTRY() : m_bc(-1) {}
			ENTRY(const std::string &t, size_t B, size_t F, size_t S1, size_t E1, size_t S2, size_t E2) :
				m_rg1(B, F, S1, E1, S2, E2)
			{
				if (BC::isValidBCStr(t))
					m_bc = BC::str2idx(t);
				else
					throw std::runtime_error("Unsupported B.C. name: \"" + t + "\"");
			}
			ENTRY(const std::string &t, size_t *s) :
				m_rg1(s)
			{
				if (BC::isValidBCStr(t))
					m_bc = BC::str2idx(t);
				else
					throw std::runtime_error("Unsupported B.C. name: \"" + t + "\"");
			}
			ENTRY(const ENTRY &rhs) = default;
			virtual ~ENTRY() = default;

			int Type() const { return m_bc; }
			int &Type() { return m_bc; }

			RANGE &Range1() { return m_rg1; }
			const RANGE &Range1() const { return m_rg1; }

			virtual int contains(size_t bs, size_t fs, size_t lpri, size_t lsec) const = 0;
		};
		class SingleSideEntry : public ENTRY
		{
		public:
			SingleSideEntry() = default;
			SingleSideEntry(const std::string &t, size_t B, size_t F, size_t S1, size_t E1, size_t S2, size_t E2) : ENTRY(t, B, F, S1, E1, S2, E2) {}
			SingleSideEntry(const std::string &t, size_t *s) : ENTRY(t, s) {}
			SingleSideEntry(const SingleSideEntry &rhs) = default;
			~SingleSideEntry() = default;

			int contains(size_t bs, size_t fs, size_t lpri, size_t lsec) const
			{
				const auto &rg = this->Range1();
				return (rg.B() == bs && rg.F() == fs && rg.constains(lpri, lsec)) ? 1 : 0;
			}
		};
		class DoubleSideEntry : public ENTRY
		{
		private:
			RANGE m_rg2;
			bool m_swap;

		public:
			DoubleSideEntry() : m_swap(false) {}
			DoubleSideEntry(const std::string &t, size_t *s1, size_t *s2, bool f) : ENTRY(t, s1), m_rg2(s2), m_swap(f) {}
			DoubleSideEntry(const DoubleSideEntry &rhs) = default;
			~DoubleSideEntry() = default;

			RANGE &Range2() { return m_rg2; }
			const RANGE &Range2() const { return m_rg2; }

			bool Swap() const { return m_swap; }
			bool &Swap() { return m_swap; }

			int contains(size_t bs, size_t fs, size_t lpri, size_t lsec) const
			{
				const auto &rg1 = this->Range1();
				const auto &rg2 = this->Range2();

				if (rg1.B() == bs && rg1.F() == fs && rg1.constains(lpri, lsec))
					return 1;
				else if (rg2.B() == bs && rg2.F() == fs && rg2.constains(lpri, lsec))
					return 2;
				else
					return 0;
			}
		};

	private:
		Array1D<Block3D*> m_blk;
		Array1D<ENTRY*> m_entry;
		Array1D<Array1D<Block3D::VERTEX*>> m_vertex;
		Array1D<Array1D<Block3D::FRAME*>> m_frame;
		Array1D<Array1D<Block3D::SURF*>> m_surf;

	public:
		Mapping3D() = default;
		Mapping3D(const std::string &inp)
		{
			readFromFile(inp);
			compute_topology();
		}
		Mapping3D(const Mapping3D &rhs) :
			m_blk(rhs.nBlock(), nullptr),
			m_entry(rhs.m_entry.size(), nullptr)
		{
			// Copy blocks
			for (size_t i = 0; i < m_blk.size(); ++i)
				m_blk[i] = new Block3D(*rhs.m_blk[i]);

			// Copy entries
			for (size_t i = 0; i < m_entry.size(); ++i)
			{
				auto ptr1 = dynamic_cast<SingleSideEntry*>(rhs.m_entry[i]);
				auto ptr2 = dynamic_cast<DoubleSideEntry*>(rhs.m_entry[i]);

				if (ptr2)
					m_entry[i] = new DoubleSideEntry(*ptr2);
				else
					m_entry[i] = new SingleSideEntry(*ptr1);
			}
		}
		~Mapping3D()
		{
			release_all();
		}

		int readFromFile(const std::string &path)
		{
			std::string s;
			std::stringstream ss;

			//Open file
			std::ifstream mfp(path);
			if (mfp.fail())
				throw std::runtime_error("Can not open target input file: \"" + path + "\".");

			//Skip header
			do {
				std::getline(mfp, s);
			} while (isBlankLine(s) || checkStarting(s, '#'));

			//Read block nums
			int NumOfBlk = -1;
			static const std::regex pattern1(R"(\s*(\d+)\s*)");
			std::smatch res1;
			if (std::regex_match(s, res1, pattern1))
			{
				NumOfBlk = std::stoi(res1[1].str());
				if (NumOfBlk > 0)
				{
					// NOT release all existing resources until 
					// it is ensured that this input file is valid.
					release_all();

					// Re-Allocate storage for new recordings
					m_blk.resize(NumOfBlk, nullptr);
				}
				else
					throw std::runtime_error("Invalid num of blocks: \"" + res1[1].str() + "\".");
			}
			else
				throw std::runtime_error("Failed to match the single line, where only the num of blocks is specified.");

			// Read dimension info of each block
			static const std::regex pattern2(R"(\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*)");
			for (int i = 0; i < NumOfBlk; i++)
			{
				std::getline(mfp, s);
				std::smatch res2;
				if (std::regex_match(s, res2, pattern2))
				{
					const int idx = std::stoi(res2[1].str());
					if (idx < 1 || idx > NumOfBlk)
						throw std::runtime_error("Invalid order of block: " + std::to_string(idx));

					const int i_max = std::stoi(res2[2].str());
					if (i_max < 1)
						throw std::runtime_error("Invalid I dimension: " + std::to_string(i_max));

					const int j_max = std::stoi(res2[3].str());
					if (j_max < 1)
						throw std::runtime_error("Invalid J dimension: " + std::to_string(j_max));

					const int k_max = std::stoi(res2[4].str());
					if (k_max < 1)
						throw std::runtime_error("Invalid K dimension: " + std::to_string(k_max));

					auto &e = m_blk(idx);
					e = new Block3D(i_max, j_max, k_max);
					e->index() = idx;
				}
				else
					throw std::runtime_error("Failed to match 4 integers.");
			}

			//Skip separators
			do {
				std::getline(mfp, s);
			} while (isBlankLine(s) || checkStarting(s, '#'));

			//Read connections
			do {
				BC::str_formalize(s);
				ss.clear();
				ss << s;
				std::string bc_str;
				ss >> bc_str;
				size_t connectivity[2][6] = { 0 };
				if (BC::str2idx(bc_str) == BC::ONE_TO_ONE)
				{
					for (int i = 0; i < 2; i++)
						for (int j = 0; j < 6; j++)
							ss >> connectivity[i][j];

					std::string swp;
					ss >> swp;

					auto ptr = new DoubleSideEntry(bc_str, connectivity[0], connectivity[1], swp == "TRUE");
					m_entry.push_back(ptr);
				}
				else
				{
					for (int i = 0; i < 6; i++)
						ss >> connectivity[0][i];

					auto ptr = new SingleSideEntry(bc_str, connectivity[0]);
					m_entry.push_back(ptr);
				}
			} while (std::getline(mfp, s));

			// Close input file
			mfp.close();

			// Finalize
			return 0;
		}

		void compute_topology()
		{
			connecting();

			const int nfm = coloring_frame();
			if (nfm < Block3D::NumOfFrame)
				throw std::runtime_error("Internal error occured when counting frames.");

			m_frame.resize(nfm);
			for (auto &e : m_frame)
				e.clear();
			for (auto b : m_blk)
				for (int i = 1; i <= Block3D::NumOfFrame; ++i)
				{
					auto &e = b->frame(i);
					m_frame(e.global_index).push_back(&e);
				}

			const int nvt = coloring_vertex();
			if (nvt < Block3D::NumOfVertex)
				throw std::runtime_error("Internal error occured when counting vertexes.");

			m_vertex.resize(nvt);
			for (auto &e : m_vertex)
				e.clear();
			for (auto b : m_blk)
				for (int i = 1; i <= Block3D::NumOfVertex; ++i)
				{
					auto &v = b->vertex(i);
					m_vertex(v.global_index).push_back(&v);
				}
		}

		void summary(std::ostream &out = std::cout)
		{
			std::cout << "=================================== SUMMARY ===================================" << std::endl;
			std::cout << "Total num of blocks: " << nBlock() << std::endl;
			size_t nSa = 0, nSi = 0, nSb = 0;
			nSurface(nSa, nSi, nSb);
			std::cout << "Total num of surfaces: " << nSa << ", among which " << nSi << " are internal, " << nSb << " located at boundary" << std::endl;
			std::cout << "Total num of frames: " << nFrame() << std::endl;
			std::cout << "Total num of vertexs: " << nVertex() << std::endl;
			std::cout << "Total num of HEX cells: " << nCell() << std::endl;
			std::cout << "Total num of QUAD faces: " << nFace() << " (duplication removed)" << std::endl;
			std::cout << "Total num of nodes: " << nNode() << " (duplication removed)" << std::endl;
			std::cout << "-------------------------------------------------------------------------------" << std::endl;
			static const std::string sep("    ");
			for (auto &b : m_blk)
			{
				std::cout << "Block" << b->index() << ":" << std::endl;
				std::cout << sep << "I=" << b->IDIM() << ", J=" << b->JDIM() << ", K=" << b->KDIM() << std::endl;
				std::cout << sep << "Num of HEX cells: " << b->cell_num() << std::endl;
				std::cout << sep << "Num of QUAD faces: " << b->face_num() << std::endl;
				std::cout << sep << "Num of nodes: " << b->node_num() << std::endl;
				std::cout << sep << "Local Frame Index:  ";
				for (int i = 1; i <= Block3D::NumOfFrame; ++i)
					std::cout << std::setw(4) << std::right << b->frame(i).local_index;
				std::cout << std::endl;
				std::cout << sep << "Global Frame Index: ";
				for (int i = 1; i <= Block3D::NumOfFrame; ++i)
					std::cout << std::setw(4) << std::right << b->frame(i).global_index;
				std::cout << std::endl;
			}
			std::cout << "-------------------------------------------------------------------------------" << std::endl;
			for (int i = 1; i <= nFrame(); ++i)
			{
				std::cout << "Frame" << i << ": " << std::endl;
				std::cout << sep << "Num of nodes: " << -1 << std::endl;
				std::cout << sep << "Occurance: ";
				for (auto e : m_frame(i))
					std::cout << "(" << e->dependentBlock->index() << ", " << e->local_index << ") ";
				std::cout << std::endl;

			}
			std::cout << "===================================== END =====================================" << std::endl;
		}

		int numbering()
		{
			// Indexing of cells
			size_t cnt = 0;
			for (auto &blk : m_blk)
			{
				const size_t cI = blk->IDIM();
				const size_t cJ = blk->JDIM();
				const size_t cK = blk->KDIM();

				for (size_t k = 1; k < cK; ++k)
					for (size_t j = 1; j < cJ; ++j)
						for (size_t i = 1; i < cI; ++i)
							blk->cell(i, j, k).CellSeq() = ++cnt;
			}

			const size_t totalCellNum = nCell();
			if (cnt != totalCellNum)
				throw std::length_error("Inconsistent num of cells detected.");

			// Indexing of faces
			cnt = 0;
			for (auto &blk : m_blk)
			{
				const size_t cI = blk->IDIM();
				const size_t cJ = blk->JDIM();
				const size_t cK = blk->KDIM();

				// TODO
			}
			const size_t totalFaceNum = nFace();
			if (cnt != totalFaceNum)
				throw std::length_error("Inconsistent num of cells detected.");

			// Indexing of nodes
			cnt = 0;
			// TODO

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
			const size_t NumOfBlk = nBlock();
			f_out << std::setw(8) << std::right << NumOfBlk << std::endl;
			for (size_t i = 0; i < NumOfBlk; i++)
			{
				f_out << std::setw(8) << std::right << i + 1;
				f_out << std::setw(8) << std::right << m_blk[i]->IDIM();
				f_out << std::setw(8) << std::right << m_blk[i]->JDIM();
				f_out << std::setw(8) << std::right << m_blk[i]->KDIM();
				f_out << std::endl;
			}

			// Interface and boundary info
			f_out << "# ============================================================================================================" << std::endl;
			f_out << "# Type           B1    F1       S1    E1       S2    E2       B2    F2       S1    E1       S2    E2      Swap" << std::endl;
			f_out << "# ------------------------------------------------------------------------------------------------------------" << std::endl;
			for (auto & e : m_entry)
			{
				f_out << std::setw(13) << std::left << BC::idx2str(e->Type());
				f_out << std::setw(6) << std::right << e->Range1().B();
				f_out << std::setw(6) << std::right << e->Range1().F();
				f_out << std::setw(9) << std::right << e->Range1().S1();
				f_out << std::setw(6) << std::right << e->Range1().E1();
				f_out << std::setw(9) << std::right << e->Range1().S2();
				f_out << std::setw(6) << std::right << e->Range1().E2();
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = dynamic_cast<DoubleSideEntry*>(e);
					f_out << std::setw(9) << std::right << p->Range2().B();
					f_out << std::setw(6) << std::right << p->Range2().F();
					f_out << std::setw(9) << std::right << p->Range2().S1();
					f_out << std::setw(6) << std::right << p->Range2().E1();
					f_out << std::setw(9) << std::right << p->Range2().S2();
					f_out << std::setw(6) << std::right << p->Range2().E2();
					f_out << std::setw(10) << std::right << (p->Swap() ? "TRUE" : "FALSE");
				}
				f_out << std::endl;
			}

			// Close output file
			f_out.close();

			// Finalize
			return 0;
		}

		size_t nBlock() const
		{
			return m_blk.size();
		}

		void nSurface(size_t &_all, size_t &_inter, size_t &_bdry) const
		{
			_all = Block3D::NumOfSurf * nBlock();
			_inter = 0;
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					_all -= 1;
					_inter += 1;
				}
			}
			_bdry = _all - _inter;
		}

		size_t nFrame() const
		{
			return m_frame.size();
		}

		size_t nVertex() const
		{
			return m_vertex.size();
		}

		size_t nCell() const
		{
			size_t ret = 0;
			for (const auto & blk : m_blk)
				ret += blk->cell_num();
			return ret;
		}

		size_t nFace() const
		{
			size_t ret = 0;
			for (const auto &blk : m_blk)
				ret += blk->face_num();

			// Substract duplicated interface
			for (const auto &e : m_entry)
				if (e->Type() == BC::ONE_TO_ONE)
					ret -= e->Range1().face_num();

			return ret;
		}

		size_t nNode() const
		{
			size_t ret = nVertex();
			for (auto b : m_blk)
			{
				ret += b->block_internal_node_num();
				for (int i = 1; i <= Block3D::NumOfSurf; ++i)
					ret += b->surface_internal_node_num(i);
			}
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					const auto b = e->Range1().B();
					const auto f = e->Range1().F();
					ret -= m_blk(b)->surface_internal_node_num(f);
				}
			}
			for (const auto &f : m_frame)
			{
				auto ff = f[0];
				int b = ff->dependentBlock->index();
				int f = ff->local_index;
				ret += m_blk(b)->frame_internal_node_num(f);
			}
			return ret;
		}

		Block3D &block(int n)
		{
			if (1 <= n && n <= nBlock())
				return *m_blk[n - 1];
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid 1-based index.");
		}
		const Block3D &block(int n) const
		{
			if (1 <= n && n <= nBlock())
				return *m_blk[n - 1];
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid 1-based index.");
		}

		void add_block(size_t _nI, size_t _nJ, size_t _nK)
		{
			auto b = new Block3D(_nI, _nJ, _nK);
			b->index() = nBlock() + 1;
			m_blk.push_back(b);
		}
		void add_block(const Block3D &x)
		{
			auto b = new Block3D(x);
			b->index() = nBlock() + 1;
			m_blk.push_back(b);
		}

		void add_entry(const std::string &_bc, size_t b, size_t f, size_t s1, size_t e1, size_t s2, size_t e2)
		{
			auto e = new SingleSideEntry(_bc, b, f, s1, e1, s2, e2);
			m_entry.push_back(e);
		}
		void add_entry(const SingleSideEntry &x)
		{
			auto e = new SingleSideEntry(x);
			m_entry.push_back(e);
		}
		void add_entry(const DoubleSideEntry &x)
		{
			auto e = new DoubleSideEntry(x);
			m_entry.push_back(e);
		}

	private:
		void release_all()
		{
			// Release memory used for blocks.
			for (auto e : m_blk)
				if (e)
					delete e;

			// Release memory used for entries.
			for (auto e : m_entry)
				if (e)
					delete e;

			m_blk.clear();
			m_entry.clear();
		}

		static bool isWhite(char c)
		{
			return c == '\n' || c == ' ' || c == '\t';
		}

		static bool isBlankLine(const std::string &s)
		{
			for (const auto &e : s)
				if (!isWhite(e))
					return false;
			return true;
		}

		static bool checkStarting(const std::string &s, char c)
		{
			for (const auto &e : s)
			{
				if (isWhite(e))
					continue;
				else
					return e == c;
			}
			return false;
		}

		void connecting()
		{
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = dynamic_cast<DoubleSideEntry*>(e);
					auto B1 = m_blk(p->Range1().B());
					auto B2 = m_blk(p->Range2().B());
					auto F1 = &B1->surf(p->Range1().F());
					auto F2 = &B2->surf(p->Range2().F());

					// Surface connectivity
					F1->neighbourSurf = F2;
					F2->neighbourSurf = F1;

					// Counterpart concerning frames on each surface.
					// There're 4 possible mapping cases.
					if (p->Swap())
					{
						// When the primary directions of F1 and F2 are not aligned,
						// the primary direction of F1 goes parallel with the secondary
						// direction of F2, and the secondary direction of F1 goes
						// parallel with the primary direction of F2. However, under
						// the right-hand convention, there're 2 further possibilities:
						// one is the primary direction of F1 and the secondary direction
						// of F2 not only go parallel, but also run in the same direction,
						// in this case, the remaining pair MUST runs in different direction;
						// the other is the primary direction of F1 and the secondary
						// direction of F2 go parallel, but run in different direction,
						// in this case, the remaining pair MUST run in same direction.
						// The exact case is determined by checking the trend of
						// corresponding ranges.

						if (p->Range1().pri_trend() == p->Range2().sec_trend()) // case 1
						{
							F1->counterpartFrame[0] = F2->includedFrame[1];
							F1->counterpartFrame[1] = F2->includedFrame[2];
							F1->counterpartFrame[2] = F2->includedFrame[3];
							F1->counterpartFrame[3] = F2->includedFrame[0];

							F2->counterpartFrame[0] = F1->includedFrame[3];
							F2->counterpartFrame[1] = F1->includedFrame[0];
							F2->counterpartFrame[2] = F1->includedFrame[1];
							F2->counterpartFrame[3] = F1->includedFrame[2];

							F1->counterpartVertex[0] = F2->includedVertex[1];
							F1->counterpartVertex[1] = F2->includedVertex[2];
							F1->counterpartVertex[2] = F2->includedVertex[3];
							F1->counterpartVertex[3] = F2->includedVertex[0];

							F2->counterpartVertex[0] = F1->includedVertex[3];
							F2->counterpartVertex[1] = F1->includedVertex[0];
							F2->counterpartVertex[2] = F1->includedVertex[1];
							F2->counterpartVertex[3] = F1->includedVertex[2];
						}
						else // case 2
						{
							F1->counterpartFrame[0] = F2->includedFrame[3];
							F1->counterpartFrame[1] = F2->includedFrame[0];
							F1->counterpartFrame[2] = F2->includedFrame[1];
							F1->counterpartFrame[3] = F2->includedFrame[2];

							F2->counterpartFrame[0] = F1->includedFrame[1];
							F2->counterpartFrame[1] = F1->includedFrame[2];
							F2->counterpartFrame[2] = F1->includedFrame[3];
							F2->counterpartFrame[3] = F1->includedFrame[0];

							F1->counterpartVertex[0] = F2->includedVertex[3];
							F1->counterpartVertex[1] = F2->includedVertex[0];
							F1->counterpartVertex[2] = F2->includedVertex[1];
							F1->counterpartVertex[3] = F2->includedVertex[2];

							F2->counterpartVertex[0] = F1->includedVertex[1];
							F2->counterpartVertex[1] = F1->includedVertex[2];
							F2->counterpartVertex[2] = F1->includedVertex[3];
							F2->counterpartVertex[3] = F1->includedVertex[0];
						}
					}
					else
					{
						// Even the primary directions of F1 and F2 are aligned, 
						// they may be in opposite directions. This is further deteced 
						// by the compare the trend from S1 to E1 in F1 with that in F2. 
						// If these 2 trends are the same, the 2 primary directions are 
						// not only parallel but also in the same directions, otherwise, 
						// they are only parallel, but runs in different directions.

						if (p->Range1().pri_trend() != p->Range2().pri_trend()) // parallel, but opposite.
						{
							F1->counterpartFrame[0] = F2->includedFrame[2];
							F1->counterpartFrame[1] = F2->includedFrame[3];
							F1->counterpartFrame[2] = F2->includedFrame[0];
							F1->counterpartFrame[3] = F2->includedFrame[1];

							F2->counterpartFrame[0] = F1->includedFrame[2];
							F2->counterpartFrame[1] = F1->includedFrame[3];
							F2->counterpartFrame[2] = F1->includedFrame[0];
							F2->counterpartFrame[3] = F1->includedFrame[1];

							F1->counterpartVertex[0] = F2->includedVertex[2];
							F1->counterpartVertex[1] = F2->includedVertex[3];
							F1->counterpartVertex[2] = F2->includedVertex[0];
							F1->counterpartVertex[3] = F2->includedVertex[1];

							F2->counterpartVertex[0] = F1->includedVertex[2];
							F2->counterpartVertex[1] = F1->includedVertex[3];
							F2->counterpartVertex[2] = F1->includedVertex[0];
							F2->counterpartVertex[3] = F1->includedVertex[1];
						}
						else // parallel, and same direction.
						{
							F1->counterpartFrame[0] = F2->includedFrame[0];
							F1->counterpartFrame[1] = F2->includedFrame[1];
							F1->counterpartFrame[2] = F2->includedFrame[2];
							F1->counterpartFrame[3] = F2->includedFrame[3];

							F2->counterpartFrame[0] = F1->includedFrame[0];
							F2->counterpartFrame[1] = F1->includedFrame[1];
							F2->counterpartFrame[2] = F1->includedFrame[2];
							F2->counterpartFrame[3] = F1->includedFrame[3];

							F1->counterpartVertex[0] = F2->includedVertex[0];
							F1->counterpartVertex[1] = F2->includedVertex[1];
							F1->counterpartVertex[2] = F2->includedVertex[2];
							F1->counterpartVertex[3] = F2->includedVertex[3];

							F2->counterpartVertex[0] = F1->includedVertex[0];
							F2->counterpartVertex[1] = F1->includedVertex[1];
							F2->counterpartVertex[2] = F1->includedVertex[2];
							F2->counterpartVertex[3] = F1->includedVertex[3];
						}
					}
				}
			}
		}

		int coloring_frame()
		{
			int global_cnt = 0;
			for (auto b : m_blk)
			{
				for (size_t j = 1; j <= Block3D::NumOfFrame; ++j)
				{
					auto e = &b->frame(j);
					if (e->global_index != 0)
						continue;

					++global_cnt;
					std::queue<Block3D::FRAME*> q;
					q.push(e);

					// BFS
					while (!q.empty())
					{
						auto ce = q.front();
						q.pop();
						ce->global_index = global_cnt;

						for (auto sf : ce->dependentSurf)
						{
							if (!sf)
								throw std::runtime_error("Dependent surface of a frame should NOT be empty.");

							if (sf->neighbourSurf)
							{
								Block3D::FRAME *t = nullptr;
								for (int ii = 0; ii < 4; ++ii)
									if (sf->includedFrame[ii] == ce)
									{
										t = sf->counterpartFrame[ii];
										break;
									}
								if (!t)
									throw std::runtime_error("Internal error.");

								if (t->global_index == 0)
									q.push(t);
							}
						}
					}
				}
			}
			return global_cnt; // The total num of block frames.
		}

		int coloring_vertex()
		{
			int global_cnt = 0;
			for (auto b : m_blk)
			{
				for (size_t j = 1; j <= Block3D::NumOfVertex; ++j)
				{
					auto v = &b->vertex(j);
					if (v->global_index != 0)
						continue;

					++global_cnt;
					std::stack<Block3D::VERTEX*> s;
					s.push(v);

					// DFS
					while (!s.empty())
					{
						auto cv = s.top();
						s.pop();
						cv->global_index = global_cnt;

						for (auto sf : cv->dependentSurf)
						{
							if (!sf)
								throw std::runtime_error("Dependent surface of a vertex should NOT be empty.");

							if (sf->neighbourSurf)
							{
								Block3D::VERTEX *t = nullptr;
								for (int ii = 0; ii < 4; ++ii)
									if (sf->includedVertex[ii] == cv)
									{
										t = sf->counterpartVertex[ii];
										break;
									}
								if (!t)
									throw std::runtime_error("Counterpart vertex should exist.");

								if (t->global_index == 0)
									s.push(t);
							}
						}
					}
				}
			}
			return global_cnt;
		}
	};
}

#endif
