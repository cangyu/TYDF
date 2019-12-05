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
#include <exception>
#include <stdexcept>
#include <regex>
#include "common.hpp"

namespace NMF
{
	using COMMON::wrong_index;
	using COMMON::Array1D;

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
		size_t m_cell; // 1-based cell index.

	public:
		CELL(size_t idx = 0) : m_cell(idx) {}
		CELL(const CELL &rhs) = default;
		virtual ~CELL() = default;

		size_t CellSeq() const { return m_cell; }
		size_t &CellSeq() { return m_cell; }

		// 1-based indexing of node
		virtual size_t NodeSeq(size_t n) const = 0;
		virtual size_t &NodeSeq(size_t n) = 0;

		// 1-based indexing of face
		virtual size_t FaceSeq(size_t n) const = 0;
		virtual size_t &FaceSeq(size_t n) = 0;
	};

	class QUAD_CELL : public CELL
	{
	private:
		std::array<size_t, 4> m_node; // 1-based node sequence.
		std::array<size_t, 4> m_face; // 1-based face sequence.

	public:
		QUAD_CELL(size_t idx = 0) : CELL(idx), m_node{ 0, 0, 0, 0 }, m_face{ 0, 0, 0, 0 } {}
		QUAD_CELL(const QUAD_CELL &rhs) = default;
		~QUAD_CELL() = default;

		// 1-based indexing of node
		size_t NodeSeq(size_t n) const { return m_node.at(n - 1); }
		size_t &NodeSeq(size_t n) { return m_node.at(n - 1); }

		// 1-based indexing of face
		size_t FaceSeq(size_t n) const { return m_face.at(n - 1); }
		size_t &FaceSeq(size_t n) { return m_face.at(n - 1); }
	};

	class HEX_CELL : public CELL
	{
	private:
		std::array<size_t, 8> m_node; // 1-based node sequence.
		std::array<size_t, 6> m_face; // 1-based face sequence.

	public:
		HEX_CELL(size_t idx = 0) : CELL(idx), m_node{ 0, 0, 0, 0, 0, 0, 0, 0 }, m_face{ 0, 0, 0, 0, 0, 0 } {}
		HEX_CELL(const HEX_CELL &rhs) = default;
		~HEX_CELL() = default;

		// 1-based indexing of node
		size_t NodeSeq(size_t n) const { return m_node.at(n - 1); }
		size_t &NodeSeq(size_t n) { return m_node.at(n - 1); }

		// 1-based indexing of face
		size_t FaceSeq(size_t n) const { return m_face.at(n - 1); }
		size_t &FaceSeq(size_t n) { return m_face.at(n - 1); }
	};

	class BLOCK
	{
	protected:
		size_t m_idx; // 1-based global index
		std::string m_name;
		std::array<size_t, 3> m_dim; // Dimensions

	public:
		BLOCK() = delete;
		BLOCK(size_t nI, size_t nJ) :
			m_idx(0),
			m_name(""),
			m_dim{ nI, nJ, 1 }
		{
			if (nI < 2 || nJ < 2)
				throw std::invalid_argument("Invalid dimension.");
		}
		BLOCK(size_t nI, size_t nJ, size_t nK) :
			m_idx(0),
			m_name(""),
			m_dim{ nI, nJ, nK }
		{
			if (nI < 2 || nJ < 2 || nK < 2)
				throw std::invalid_argument("Invalid dimension.");
		}
		BLOCK(const BLOCK &rhs) : m_idx(rhs.index()), m_name(rhs.name()), m_dim{ rhs.IDIM(), rhs.JDIM(), rhs.KDIM() } {}
		virtual ~BLOCK() = default;

		size_t index() const { return m_idx; }
		size_t &index() { return m_idx; }

		const std::string &name() const { return m_name; }
		std::string &name() { return m_name; }

		size_t IDIM() const { return m_dim[0]; }
		size_t JDIM() const { return m_dim[1]; }
		size_t KDIM() const { return m_dim[2]; }

		bool is3D() const { return (KDIM() != 1); }
		short dimension() const { return is3D() ? 3 : 2; }

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

		virtual size_t surface_internal_node_num(short s_idx) const = 0;
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
			std::string name = "";
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
		Array1D<QUAD_CELL*> m_cell;
		Array1D<FRAME> m_frame;
		Array1D<VERTEX> m_vertex;

	public:
		Block2D() = delete;
		Block2D(int nI, int nJ) :
			BLOCK(nI, nJ),
			m_cell(cell_num(), nullptr),
			m_vertex(NumOfVertex),
			m_frame(NumOfFrame)
		{
			for (short i = 0; i < NumOfFrame; ++i)
				m_frame[i].local_index = i + 1;

			for (short i = 0; i < NumOfVertex; ++i)
				m_vertex[i].local_index = i + 1;
		}
		Block2D(const Block2D &rhs) :
			BLOCK(rhs.IDIM(), rhs.JDIM()),
			m_cell(rhs.m_cell.size(), nullptr),
			m_frame(rhs.m_frame),
			m_vertex(rhs.m_vertex)
		{
			for (size_t i = 0; i < m_cell.size(); ++i)
			{
				auto p = rhs.m_cell[i];
				if (p)
					m_cell[i] = new QUAD_CELL(*p);
			}
		}
		~Block2D()
		{
			release_cell_storage();
		}

		void release_cell_storage()
		{
			for (auto &c : m_cell)
				if (c != nullptr)
				{
					delete c;
					c = nullptr;
				}
		}
		void allocate_cell_storage()
		{
			for (auto &c : m_cell)
			{
				c = new QUAD_CELL();
				if (!c)
					throw std::bad_alloc();
			}
		}

		// Access internal cell through 1-based index.
		// Indexing convention:
		//      "i" ranges from 1 to IDIM()-1;
		//      "j" ranges from 1 to JDIM()-1;
		// When the IJ-axis follows the right-hand convention, (i, j) represents
		// the left-most, bottom-most node of the selected cell.
		QUAD_CELL &cell(size_t i, size_t j)
		{
			const size_t i0 = i - 1, j0 = j - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (IDIM() - 1) * j0;
			return *m_cell.at(idx);
		}
		const QUAD_CELL &cell(size_t i, size_t j) const
		{
			const size_t i0 = i - 1, j0 = j - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (IDIM() - 1) * j0;
			return *m_cell.at(idx);
		}

		// Access frame through 1-based index.
		// The indexing convention follows OpenFOAM specification.
		FRAME &frame(short n)
		{
			if (1 <= n && n <= NumOfFrame)
				return m_frame.at(n - 1);
			else if (-NumOfFrame <= n && n <= -1)
				return m_frame.at(NumOfFrame + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid frame index for a 2D block.");
		}
		const FRAME &frame(short n) const
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
		VERTEX &vertex(short n)
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid vertex index for a 2D block.");
		}
		const VERTEX &vertex(short n) const
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid vertex index for a 2D block.");
		}

		size_t surface_internal_node_num(short s_idx) const
		{
			switch (s_idx)
			{
			case 1:
			case 2:
				return (JDIM() - 2);
			case 3:
			case 4:
				return (IDIM() - 2);
			default:
				throw std::invalid_argument("\"" + std::to_string(s_idx) + "\" is not a valid surface index of a 2D block.");
			}
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
			std::array<bool, 4> counterpartFrameIsOpposite{ false, false, false, false };
			std::array<VERTEX*, 4> includedVertex{ nullptr, nullptr, nullptr, nullptr };
			std::array<VERTEX*, 4> counterpartVertex{ nullptr, nullptr, nullptr, nullptr };
		};

	protected:
		class not_a_surface : public wrong_index
		{
		public:
			not_a_surface(short x) :
				wrong_index(x)
			{
				m_msg += "is not a valid surface index of a 3D block.";
			}
		};
		class not_a_frame : public wrong_index
		{
		public:
			not_a_frame(short x) :
				wrong_index(x)
			{
				m_msg += "is not a valid frame index of a 3D block.";
			}
		};
		class not_a_vertex : public wrong_index
		{
		public:
			not_a_vertex(short x) :
				wrong_index(x)
			{
				m_msg += "is not a valid vertex index of a 3D block.";
			}
		};

	private:
		Array1D<HEX_CELL*> m_cell;
		Array1D<VERTEX> m_vertex;
		Array1D<FRAME> m_frame;
		Array1D<SURF> m_surf;

	public:
		Block3D() = delete;
		Block3D(int nI, int nJ, int nK) :
			BLOCK(nI, nJ, nK),
			m_cell(cell_num(), nullptr),
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
		Block3D(const Block3D &rhs) :
			BLOCK(rhs.IDIM(), rhs.JDIM(), rhs.KDIM()),
			m_cell(cell_num(), nullptr),
			m_vertex(rhs.m_vertex),
			m_frame(rhs.m_frame),
			m_surf(rhs.m_surf)
		{
			for (size_t i = 0; i < m_cell.size(); ++i)
			{
				auto p = rhs.m_cell[i];
				if (p)
					m_cell[i] = new HEX_CELL(*p);
			}
		}
		~Block3D()
		{
			release_cell_storage();
		}

		void release_cell_storage()
		{
			for (auto &c : m_cell)
				if (c != nullptr)
				{
					delete c;
					c = nullptr;
				}
		}
		void allocate_cell_storage()
		{
			for (auto &c : m_cell)
			{
				c = new HEX_CELL();
				if (!c)
					throw std::bad_alloc();
			}
		}

		/// Access internal cell through 1-based index.
		/// Indexing convention:
		///      "i" ranges from 1 to IDIM()-1;
		///      "j" ranges from 1 to JDIM()-1;
		///      "k" ranges from 1 to KDIM()-1;
		/// When the IJK-axis follows the right-hand convention, (i, j, k) represents
		/// the left-most, bottom-most and back-most node of the selected cell.
		HEX_CELL &cell(size_t i, size_t j, size_t k)
		{
			const size_t i0 = i - 1, j0 = j - 1, k0 = k - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (IDIM() - 1) * (j0 + (JDIM() - 1) * k0);
			return *m_cell.at(idx);
		}
		const HEX_CELL &cell(size_t i, size_t j, size_t k) const
		{
			const size_t i0 = i - 1, j0 = j - 1, k0 = k - 1; // Convert 1-based index to 0-based
			const size_t idx = i0 + (IDIM() - 1) * (j0 + (JDIM() - 1) * k0);
			return *m_cell.at(idx);
		}

		/// Access vertex through 1-based index.
		/// The indexing convention follows OpenFOAM specification.
		VERTEX &vertex(short n)
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw not_a_vertex(n);
		}
		const VERTEX &vertex(short n) const
		{
			if (1 <= n && n <= NumOfVertex)
				return m_vertex.at(n - 1);
			else if (-NumOfVertex <= n && n <= -1)
				return m_vertex.at(NumOfVertex + n);
			else
				throw not_a_vertex(n);
		}

		/// Access the frame edges through 1-based index.
		/// The indexing convention follows OpenFOAM specification.
		FRAME &frame(short n)
		{
			if (1 <= n && n <= NumOfFrame)
				return m_frame.at(n - 1);
			else if (-NumOfFrame <= n && n <= -1)
				return m_frame.at(NumOfFrame + n);
			else
				throw not_a_frame(n);
		}
		const FRAME &frame(short n) const
		{
			if (1 <= n && n <= NumOfFrame)
				return m_frame.at(n - 1);
			else if (-NumOfFrame <= n && n <= -1)
				return m_frame.at(NumOfFrame + n);
			else
				throw not_a_frame(n);
		}

		/// Access surrounding surface through 1-based index.
		/// The index follows NMF convection.
		SURF &surf(short n)
		{
			if (1 <= n && n <= NumOfSurf)
				return m_surf.at(n - 1);
			else if (-NumOfSurf <= n && n <= -1)
				return m_surf.at(NumOfSurf + n);
			else
				throw not_a_surface(n);
		}
		const SURF &surf(short n) const
		{
			if (1 <= n && n <= NumOfSurf)
				return m_surf.at(n - 1);
			else if (-NumOfSurf <= n && n <= -1)
				return m_surf.at(NumOfSurf + n);
			else
				throw not_a_surface(n);
		}

		size_t surface_internal_node_num(short s) const
		{
			switch (s)
			{
			case 1:
			case 2:
				return (IDIM() - 2) * (JDIM() - 2);
			case 3:
			case 4:
				return (JDIM() - 2) * (KDIM() - 2);
			case 5:
			case 6:
				return (KDIM() - 2) * (IDIM() - 2);
			default:
				throw not_a_surface(s);
			}
		}

		size_t frame_node_num(short idx) const
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
				throw not_a_frame(idx);
			}
		}

		size_t frame_internal_node_num(short idx) const
		{
			const auto n = frame_node_num(idx);
			if (n >= 2)
				return n - 2;
			else
				throw std::runtime_error("Internal error: too less nodes not detected when constructing.");
		}

		size_t surface_node_num(short idx) const
		{
			switch (idx - 1)
			{
			case 0:
			case 1:
				return IDIM() * JDIM();
			case 2:
			case 3:
				return JDIM() * KDIM();
			case 4:
			case 5:
				return KDIM() * IDIM();
			default:
				throw not_a_surface(idx);
			}
		}

		size_t surface_face_num(short idx) const
		{
			switch (idx - 1)
			{
			case 0:
			case 1:
				return (IDIM() - 1) * (JDIM() - 1);
			case 2:
			case 3:
				return (JDIM() - 1) * (KDIM() - 1);
			case 4:
			case 5:
				return (KDIM() - 1) * (IDIM() - 1);
			default:
				throw not_a_surface(idx);
			}
		}

		size_t shell_face_num() const
		{
			size_t ret = 0;
			for (short i = 1; i <= NumOfSurf; ++i)
				ret += surface_face_num(i);
			return ret;
		}

		size_t &surface_face_index(short f, size_t pri, size_t sec)
		{
			HEX_CELL *p = nullptr;
			switch (f)
			{
			case 1:
				p = &cell(pri, sec, 1);
				break;
			case 2:
				p = &cell(pri, sec, KDIM() - 1);
				break;
			case 3:
				p = &cell(1, pri, sec);
				break;
			case 4:
				p = &cell(IDIM() - 1, pri, sec);
				break;
			case 5:
				p = &cell(sec, 1, pri);
				break;
			case 6:
				p = &cell(sec, JDIM() - 1, pri);
				break;
			default:
				throw not_a_surface(f);
			}
			return p->FaceSeq(f);
		}

		size_t &vertex_node_index(short v)
		{
			HEX_CELL *p = nullptr;
			switch (v)
			{
			case 1:
				p = &cell(1, 1, 1);
				break;
			case 2:
				p = &cell(1, 1, KDIM() - 1);
				break;
			case 3:
				p = &cell(IDIM() - 1, 1, KDIM() - 1);
				break;
			case 4:
				p = &cell(IDIM() - 1, 1, 1);
				break;
			case 5:
				p = &cell(1, JDIM() - 1, 1);
				break;
			case 6:
				p = &cell(1, JDIM() - 1, KDIM() - 1);
				break;
			case 7:
				p = &cell(IDIM() - 1, JDIM() - 1, KDIM() - 1);
				break;
			case 8:
				p = &cell(IDIM() - 1, JDIM() - 1, 1);
				break;
			default:
				throw not_a_vertex(v);
			}
			return p->NodeSeq(v);
		}

		void interior_node_occurance(size_t i, size_t j, size_t k, std::vector<size_t*> &oc)
		{
			oc[0] = &cell(i - 1, j - 1, k - 1).NodeSeq(7);
			oc[1] = &cell(i, j - 1, k - 1).NodeSeq(6);
			oc[2] = &cell(i - 1, j, k - 1).NodeSeq(3);
			oc[3] = &cell(i, j, k - 1).NodeSeq(2);
			oc[4] = &cell(i - 1, j - 1, k).NodeSeq(8);
			oc[5] = &cell(i, j - 1, k).NodeSeq(5);
			oc[6] = &cell(i - 1, j, k).NodeSeq(4);
			oc[7] = &cell(i, j, k).NodeSeq(1);
		}

		void surface_node_coordinate(short f, size_t pri_seq, size_t sec_seq, size_t &i, size_t &j, size_t &k)
		{
			switch (f)
			{
			case 1:
				k = 1;
				i = pri_seq;
				j = sec_seq;
				break;
			case 2:
				k = KDIM();
				i = pri_seq;
				j = sec_seq;
				break;
			case 3:
				i = 1;
				j = pri_seq;
				k = sec_seq;
				break;
			case 4:
				i = IDIM();
				j = pri_seq;
				k = sec_seq;
				break;
			case 5:
				j = 1;
				k = pri_seq;
				i = sec_seq;
				break;
			case 6:
				j = JDIM();
				k = pri_seq;
				i = sec_seq;
				break;
			default:
				throw not_a_surface(f);
			}
		}

		void surface_internal_node_occurance(short f, size_t pri, size_t sec, std::vector<size_t*> &oc)
		{
			size_t i = 0, j = 0, k = 0;
			surface_node_coordinate(f, pri, sec, i, j, k);
			switch (f)
			{
			case 1:
				oc[0] = &cell(i, j, k).NodeSeq(1);
				oc[1] = &cell(i - 1, j, k).NodeSeq(4);
				oc[2] = &cell(i, j - 1, k).NodeSeq(5);
				oc[3] = &cell(i - 1, j - 1, k).NodeSeq(8);
				break;
			case 2:
				oc[0] = &cell(i, j, k - 1).NodeSeq(2);
				oc[1] = &cell(i - 1, j, k - 1).NodeSeq(3);
				oc[2] = &cell(i, j - 1, k - 1).NodeSeq(6);
				oc[3] = &cell(i - 1, j - 1, k - 1).NodeSeq(7);
				break;
			case 3:
				oc[0] = &cell(i, j, k).NodeSeq(1);
				oc[1] = &cell(i, j - 1, k).NodeSeq(5);
				oc[2] = &cell(i, j, k - 1).NodeSeq(2);
				oc[3] = &cell(i, j - 1, k - 1).NodeSeq(6);
				break;
			case 4:
				oc[0] = &cell(i - 1, j, k).NodeSeq(4);
				oc[1] = &cell(i - 1, j - 1, k).NodeSeq(8);
				oc[2] = &cell(i - 1, j, k - 1).NodeSeq(3);
				oc[3] = &cell(i - 1, j - 1, k - 1).NodeSeq(7);
				break;
			case 5:
				oc[0] = &cell(i, j, k).NodeSeq(1);
				oc[1] = &cell(i - 1, j, k).NodeSeq(4);
				oc[2] = &cell(i, j, k - 1).NodeSeq(2);
				oc[3] = &cell(i - 1, j, k - 1).NodeSeq(3);
				break;
			case 6:
				oc[0] = &cell(i, j - 1, k).NodeSeq(5);
				oc[1] = &cell(i - 1, j - 1, k).NodeSeq(8);
				oc[2] = &cell(i, j - 1, k - 1).NodeSeq(6);
				oc[3] = &cell(i - 1, j - 1, k - 1).NodeSeq(7);
				break;
			default:
				throw not_a_surface(f);
			}
		}

		void frame_internal_node_occurace(short f, size_t idx, std::vector<size_t*> &oc)
		{
			switch (f - 1)
			{
			case 0:
				oc[0] = &cell(1, 1, idx).NodeSeq(1);
				oc[1] = &cell(1, 1, idx - 1).NodeSeq(2);
				break;
			case 1:
				oc[0] = &cell(IDIM() - 1, 1, idx).NodeSeq(4);
				oc[1] = &cell(IDIM() - 1, 1, idx - 1).NodeSeq(3);
				break;
			case 2:
				oc[0] = &cell(IDIM() - 1, JDIM() - 1, idx).NodeSeq(8);
				oc[1] = &cell(IDIM() - 1, JDIM() - 1, idx - 1).NodeSeq(7);
				break;
			case 3:
				oc[0] = &cell(1, JDIM() - 1, idx).NodeSeq(5);
				oc[1] = &cell(1, JDIM() - 1, idx - 1).NodeSeq(6);
				break;
			case 4:
				oc[0] = &cell(idx, 1, 1).NodeSeq(1);
				oc[1] = &cell(idx - 1, 1, 1).NodeSeq(4);
				break;
			case 5:
				oc[0] = &cell(idx, 1, KDIM() - 1).NodeSeq(2);
				oc[1] = &cell(idx - 1, 1, KDIM() - 1).NodeSeq(3);
				break;
			case 6:
				oc[0] = &cell(idx, JDIM() - 1, KDIM() - 1).NodeSeq(6);
				oc[1] = &cell(idx - 1, JDIM() - 1, KDIM() - 1).NodeSeq(7);
				break;
			case 7:
				oc[0] = &cell(idx, JDIM() - 1, 1).NodeSeq(5);
				oc[1] = &cell(idx - 1, JDIM() - 1, 1).NodeSeq(8);
				break;
			case 8:
				oc[0] = &cell(1, idx, 1).NodeSeq(1);
				oc[1] = &cell(1, idx - 1, 1).NodeSeq(5);
				break;
			case 9:
				oc[0] = &cell(1, idx, KDIM() - 1).NodeSeq(2);
				oc[1] = &cell(1, idx - 1, KDIM() - 1).NodeSeq(6);
				break;
			case 10:
				oc[0] = &cell(IDIM() - 1, idx, KDIM() - 1).NodeSeq(3);
				oc[1] = &cell(IDIM() - 1, idx - 1, KDIM() - 1).NodeSeq(7);
				break;
			case 11:
				oc[0] = &cell(IDIM() - 1, idx, 1).NodeSeq(4);
				oc[1] = &cell(IDIM() - 1, idx - 1, 1).NodeSeq(8);
				break;
			default:
				throw not_a_frame(f);
			}
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
				short m_face; // Face index, ranges from 1 to 6.
				size_t m_s1; // Primary direction starting index, 1-based.
				size_t m_e1; // Primary direction ending index, 1-based.
				size_t m_s2; // Secondary direction starting index, 1-based.
				size_t m_e2; // Secondary direction ending index, 1-based.

			public:
				RANGE() : m_blk(0), m_face(0), m_s1(0), m_e1(0), m_s2(0), m_e2(0) {}
				RANGE(size_t b, short f, size_t s1, size_t e1, size_t s2, size_t e2) :
					m_blk(b), m_face(f), m_s1(s1), m_e1(e1), m_s2(s2), m_e2(e2)
				{
					check_param();
				}
				RANGE(const RANGE &rhs) = default;
				~RANGE() = default;

				size_t B() const { return m_blk; }
				size_t &B() { return m_blk; }

				short F() const { return m_face; }
				short &F() { return m_face; }

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
					size_t ret = 1;

					if (pri_trend())
						ret += (E1() - S1());
					else
						ret += (S1() - E1());

					return ret;
				}

				// Nodes in secondary direction.
				size_t sec_node_num() const
				{
					size_t ret = 1;

					if (sec_trend())
						ret += (E2() - S2());
					else
						ret += (S2() - E2());

					return ret;
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
					if (F() <= 0)
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
			ENTRY() : m_bc(-1), m_rg1() {}
			ENTRY(const std::string &t, size_t B, short F, size_t S1, size_t E1, size_t S2, size_t E2) :
				m_rg1(B, F, S1, E1, S2, E2)
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

			virtual int contains(size_t bs, short fs, size_t lpri, size_t lsec) const = 0;
		};
		class SingleSideEntry : public ENTRY
		{
		public:
			SingleSideEntry() = default;
			SingleSideEntry(const std::string &t, size_t B, short F, size_t S1, size_t E1, size_t S2, size_t E2) : ENTRY(t, B, F, S1, E1, S2, E2) {}
			SingleSideEntry(const SingleSideEntry &rhs) = default;
			~SingleSideEntry() = default;

			int contains(size_t bs, short fs, size_t lpri, size_t lsec) const
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
			DoubleSideEntry() : m_rg2(), m_swap(false) {}
			DoubleSideEntry(const std::string &t, size_t B1, short F1, size_t S11, size_t E11, size_t S12, size_t E12, size_t B2, short F2, size_t S21, size_t E21, size_t S22, size_t E22, bool f) :
				ENTRY(t, B1, F1, S11, E11, S12, E12),
				m_rg2(B2, F2, S21, E21, S22, E22),
				m_swap(f)
			{
				if (!check_dim_consistency())
					throw std::invalid_argument("Inconsistent dimensions bwtween 2 surfaces.");
			}
			DoubleSideEntry(const DoubleSideEntry &rhs) = default;
			~DoubleSideEntry() = default;

			RANGE &Range2() { return m_rg2; }
			const RANGE &Range2() const { return m_rg2; }

			bool Swap() const { return m_swap; }
			bool &Swap() { return m_swap; }

			int contains(size_t bs, short fs, size_t lpri, size_t lsec) const
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

		private:
			bool check_dim_consistency()
			{
				if (Swap())
				{
					const bool cond1 = Range1().pri_node_num() == Range2().sec_node_num();
					const bool cond2 = Range1().sec_node_num() == Range2().pri_node_num();
					return cond1 && cond2;
				}
				else
				{
					const bool cond1 = Range1().pri_node_num() == Range2().pri_node_num();
					const bool cond2 = Range1().sec_node_num() == Range2().sec_node_num();
					return cond1 && cond2;
				}
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

			// TODO
		}
		~Mapping3D()
		{
			release_all();
		}

		void readFromFile(const std::string &path)
		{
			std::string s;
			std::stringstream ss;

			// Open file
			std::ifstream mfp(path);
			if (mfp.fail())
				throw std::runtime_error("Can not open target input file: \"" + path + "\".");

			// Skip header
			do {
				std::getline(mfp, s);
			} while (isBlankLine(s) || checkStarting(s, '#'));

			// Read block nums
			static const std::regex pattern1(R"(\s*(\d+)\s*)");
			std::smatch res1;
			if (std::regex_match(s, res1, pattern1))
			{
				int NumOfBlk = std::stoi(res1[1].str());
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
			const auto NumOfBlk = nBlock();
			for (size_t i = 0; i < NumOfBlk; i++)
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

			// Skip separators
			while (std::getline(mfp, s))
			{
				if (!isBlankLine(s) && !checkStarting(s, '#'))
					break;
			}

			// Read connections
			if (!mfp.eof())
			{
				do {
					BC::str_formalize(s);
					ss.clear();
					ss << s;
					std::string bc_str;
					ss >> bc_str;
					if (BC::str2idx(bc_str) == BC::ONE_TO_ONE)
					{
						size_t cB[2];
						short cF[2];
						size_t cS1[2], cE1[2], cS2[2], cE2[2];
						std::string swp;
						for (int i = 0; i < 2; i++)
							ss >> cB[i] >> cF[i] >> cS1[i] >> cE1[i] >> cS2[i] >> cE2[i];
						ss >> swp;
						add_entry(bc_str, cB[0], cF[0], cS1[0], cE1[0], cS2[0], cE2[0], cB[1], cF[1], cS1[1], cE1[1], cS2[1], cE2[1], swp == "TRUE");
					}
					else
					{
						size_t cB;
						short cF;
						size_t cS1, cE1, cS2, cE2;
						ss >> cB >> cF >> cS1 >> cE1 >> cS2 >> cE2;
						add_entry(bc_str, cB, cF, cS1, cE1, cS2, cE2);
					}
				} while (std::getline(mfp, s));
			}

			// Close input file
			mfp.close();
		}

		void compute_topology()
		{
			connecting();

			const int nsf = coloring_surface();
			if (nsf < Block3D::NumOfSurf)
				throw std::runtime_error("Internal error occured when counting surfaces.");
			m_surf.resize(nsf);
			for (auto &e : m_surf)
				e.clear();
			for (auto b : m_blk)
				for (short i = 1; i <= Block3D::NumOfSurf; ++i)
				{
					auto &s = b->surf(i);
					m_surf(s.global_index).push_back(&s);
				}

			const int nfm = coloring_frame();
			if (nfm < Block3D::NumOfFrame)
				throw std::runtime_error("Internal error occured when counting frames.");
			m_frame.resize(nfm);
			for (auto &e : m_frame)
				e.clear();
			for (auto b : m_blk)
				for (short i = 1; i <= Block3D::NumOfFrame; ++i)
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
				for (short i = 1; i <= Block3D::NumOfVertex; ++i)
				{
					auto &v = b->vertex(i);
					m_vertex(v.global_index).push_back(&v);
				}
		}

		void summary(std::ostream &out = std::cout)
		{
			out << "======================================== SUMMARY =======================================" << std::endl;
			out << "Num of blocks: " << nBlock() << std::endl;
			size_t nSa = 0, nSi = 0, nSb = 0;
			nSurface(nSa, nSi, nSb);
			out << "Num of surfaces: " << nSa << ", among which " << nSi << " are internal, " << nSb << " at boundary" << std::endl;
			out << "Num of frames: " << nFrame() << std::endl;
			out << "Num of vertexs: " << nVertex() << std::endl;
			out << "Num of HEX cells: " << nCell() << std::endl;
			size_t nFa = 0, nFi = 0, nFb = 0;
			nFace(nFa, nFi, nFb);
			out << "Num of QUAD faces: " << nFa << ", among which " << nFi << " are internal, " << nFb << " at boundary" << std::endl;
			out << "Num of nodes: " << nNode() << " (duplication removed)" << std::endl;
			out << "------------------------------------------Blocks----------------------------------------";
			for (auto b : m_blk)
			{
				out << "\nIndex: " << b->index() << std::endl;
				out << "I=" << b->IDIM() << ", J=" << b->JDIM() << ", K=" << b->KDIM() << std::endl;
				out << "Num of HEX cells: " << b->cell_num() << std::endl;
				out << "Num of QUAD faces: " << b->face_num() << std::endl;
				out << "Num of nodes: " << b->node_num() << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Local Vertex Index:" << std::resetiosflags(std::ios::left);
				for (int i = 1; i <= Block3D::NumOfVertex; ++i)
					out << std::setw(5) << b->vertex(i).local_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Global Vertex Index:" << std::resetiosflags(std::ios::left);
				for (int i = 1; i <= Block3D::NumOfVertex; ++i)
					out << std::setw(5) << b->vertex(i).global_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Local Frame Index:" << std::resetiosflags(std::ios::left);
				for (int i = 1; i <= Block3D::NumOfFrame; ++i)
					out << std::setw(5) << b->frame(i).local_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Global Frame Index:" << std::resetiosflags(std::ios::left);
				for (int i = 1; i <= Block3D::NumOfFrame; ++i)
					out << std::setw(5) << b->frame(i).global_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Local Surface Index:" << std::resetiosflags(std::ios::left);
				for (int i = 1; i <= Block3D::NumOfSurf; ++i)
					out << std::setw(5) << b->surf(i).local_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Global Surface Index:" << std::resetiosflags(std::ios::left);
				for (int i = 1; i <= Block3D::NumOfSurf; ++i)
					out << std::setw(5) << b->surf(i).global_index;
				out << std::endl;
			}
			out << "-----------------------------------------Surfaces---------------------------------------";
			for (int i = 1; i <= nSa; ++i)
			{
				out << "\nGlobal Index: " << i << std::endl;
				auto sf_rep = m_surf(i)[0];
				out << "Num of faces: " << sf_rep->dependentBlock->surface_face_num(sf_rep->local_index) << std::endl;
				out << "Num of nodes: " << sf_rep->dependentBlock->surface_node_num(sf_rep->local_index) << std::endl;
				out << "Occurance:";
				for (auto e : m_surf(i))
					out << " (" << e->dependentBlock->index() << ", " << e->local_index << ")";
				out << std::endl;
			}
			out << "------------------------------------------Frames----------------------------------------";
			for (int i = 1; i <= nFrame(); ++i)
			{
				out << "\nGlobal Index: " << i << std::endl;
				auto f_rep = m_frame(i)[0];
				out << "Num of nodes: " << f_rep->dependentBlock->frame_node_num(f_rep->local_index) << std::endl;
				out << "Occurance:";
				for (auto e : m_frame(i))
					out << " (" << e->dependentBlock->index() << ", " << e->local_index << ")";
				out << std::endl;
			}
			out << "-----------------------------------------Vertexes---------------------------------------";
			for (int i = 1; i <= nVertex(); ++i)
			{
				out << "\nGlobal Index: " << i << std::endl;
				out << "Occurance:";
				for (auto e : m_vertex(i))
					out << " (" << e->dependentBlock->index() << ", " << e->local_index << ")";
				out << std::endl;
			}
			out << "========================================== END =========================================" << std::endl;
		}

		void numbering()
		{
			for (auto b : m_blk)
				b->allocate_cell_storage();

			/* Index of cells */
			numbering_cell();

			/* Index of faces */
			numbering_face();

			/* Index of nodes */
			numbering_node();
		}

		void writeToFile(const std::string &path)
		{
			// Open target file
			std::ofstream f_out(path);
			if (f_out.fail())
				throw std::runtime_error("Can not open target output file: \"" + path + "\".");

			// Header
			f_out << "# ============================= Neutral Map File generated by the Grid-Glue software =========================" << std::endl;
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
					auto p = static_cast<DoubleSideEntry*>(e);
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

		void nFace(size_t &_all, size_t &_inner, size_t &_bdry) const
		{
			// Counting all faces.
			_all = 0;
			for (auto b : m_blk)
				_all += b->face_num();

			// Boundary
			_bdry = 0;
			for (auto b : m_blk)
				_bdry += b->shell_face_num();

			// Faces at internal interfaces.
			size_t ii = 0;
			for (const auto &e : m_entry)
				if (e->Type() == BC::ONE_TO_ONE)
					ii += e->Range1().face_num();

			// Substract duplicated
			_all -= ii;
			_bdry -= ii * 2;

			// Inner
			_inner = _all - _bdry;
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
			for (const auto &rec : m_frame)
			{
				auto ff = rec[0];
				auto b = ff->dependentBlock->index();
				auto f = ff->local_index;
				ret += m_blk(b)->frame_internal_node_num(f);
			}
			return ret;
		}

		Block3D &block(size_t n)
		{
			if (1 <= n && n <= nBlock())
				return *m_blk[n - 1];
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid 1-based index.");
		}
		const Block3D &block(size_t n) const
		{
			if (1 <= n && n <= nBlock())
				return *m_blk[n - 1];
			else
				throw std::invalid_argument("\"" + std::to_string(n) + "\" is not a valid 1-based index.");
		}

		void add_block(int _nI, int _nJ, int _nK)
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

		void add_entry(const std::string &_bc, size_t b, short f, size_t s1, size_t e1, size_t s2, size_t e2)
		{
			auto e = new SingleSideEntry(_bc, b, f, s1, e1, s2, e2);
			m_entry.push_back(e);
		}
		void add_entry(const SingleSideEntry &x)
		{
			auto e = new SingleSideEntry(x);
			m_entry.push_back(e);
		}
		void add_entry(const std::string &_bc, size_t b1, short f1, size_t s11, size_t e11, size_t s12, size_t e12, size_t b2, short f2, size_t s21, size_t e21, size_t s22, size_t e22, bool swp)
		{
			auto e = new DoubleSideEntry(_bc, b1, f1, s11, e11, s12, e12, b2, f2, s21, e21, s22, e22, swp);
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

		static void distribute_index(size_t s, size_t e, std::vector<size_t> &dst)
		{
			if (s > e)
			{
				const size_t n = s - e + 1;
				dst.resize(n);
				size_t val = s;

				for (size_t i = 0; i < n; ++i)
					dst[i] = val--;
			}
			else
			{
				const size_t n = e - s + 1;
				dst.resize(n);
				size_t val = s;
				for (size_t i = 0; i < n; ++i)
					dst[i] = val++;
			}
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

							F1->counterpartFrameIsOpposite[0] = false;
							F1->counterpartFrameIsOpposite[1] = true;
							F1->counterpartFrameIsOpposite[2] = false;
							F1->counterpartFrameIsOpposite[3] = true;

							F2->counterpartFrame[0] = F1->includedFrame[3];
							F2->counterpartFrame[1] = F1->includedFrame[0];
							F2->counterpartFrame[2] = F1->includedFrame[1];
							F2->counterpartFrame[3] = F1->includedFrame[2];

							F2->counterpartFrameIsOpposite[0] = true;
							F2->counterpartFrameIsOpposite[1] = false;
							F2->counterpartFrameIsOpposite[2] = true;
							F2->counterpartFrameIsOpposite[3] = false;

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

							F1->counterpartFrameIsOpposite[0] = true;
							F1->counterpartFrameIsOpposite[1] = false;
							F1->counterpartFrameIsOpposite[2] = true;
							F1->counterpartFrameIsOpposite[3] = false;

							F2->counterpartFrame[0] = F1->includedFrame[1];
							F2->counterpartFrame[1] = F1->includedFrame[2];
							F2->counterpartFrame[2] = F1->includedFrame[3];
							F2->counterpartFrame[3] = F1->includedFrame[0];

							F2->counterpartFrameIsOpposite[0] = false;
							F2->counterpartFrameIsOpposite[1] = true;
							F2->counterpartFrameIsOpposite[2] = false;
							F2->counterpartFrameIsOpposite[3] = true;

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

						if (p->Range1().pri_trend() != p->Range2().pri_trend()) // Parallel, but goes in opposite direction.
						{
							F1->counterpartFrame[0] = F2->includedFrame[2];
							F1->counterpartFrame[1] = F2->includedFrame[3];
							F1->counterpartFrame[2] = F2->includedFrame[0];
							F1->counterpartFrame[3] = F2->includedFrame[1];

							F1->counterpartFrameIsOpposite[0] = true;
							F1->counterpartFrameIsOpposite[1] = true;
							F1->counterpartFrameIsOpposite[2] = true;
							F1->counterpartFrameIsOpposite[3] = true;

							F2->counterpartFrame[0] = F1->includedFrame[2];
							F2->counterpartFrame[1] = F1->includedFrame[3];
							F2->counterpartFrame[2] = F1->includedFrame[0];
							F2->counterpartFrame[3] = F1->includedFrame[1];

							F2->counterpartFrameIsOpposite[0] = true;
							F2->counterpartFrameIsOpposite[1] = true;
							F2->counterpartFrameIsOpposite[2] = true;
							F2->counterpartFrameIsOpposite[3] = true;

							F1->counterpartVertex[0] = F2->includedVertex[2];
							F1->counterpartVertex[1] = F2->includedVertex[3];
							F1->counterpartVertex[2] = F2->includedVertex[0];
							F1->counterpartVertex[3] = F2->includedVertex[1];

							F2->counterpartVertex[0] = F1->includedVertex[2];
							F2->counterpartVertex[1] = F1->includedVertex[3];
							F2->counterpartVertex[2] = F1->includedVertex[0];
							F2->counterpartVertex[3] = F1->includedVertex[1];
						}
						else // Parallel, and goes in the same direction.
						{
							F1->counterpartFrame[0] = F2->includedFrame[0];
							F1->counterpartFrame[1] = F2->includedFrame[1];
							F1->counterpartFrame[2] = F2->includedFrame[2];
							F1->counterpartFrame[3] = F2->includedFrame[3];

							F1->counterpartFrameIsOpposite[0] = false;
							F1->counterpartFrameIsOpposite[1] = false;
							F1->counterpartFrameIsOpposite[2] = false;
							F1->counterpartFrameIsOpposite[3] = false;

							F2->counterpartFrame[0] = F1->includedFrame[0];
							F2->counterpartFrame[1] = F1->includedFrame[1];
							F2->counterpartFrame[2] = F1->includedFrame[2];
							F2->counterpartFrame[3] = F1->includedFrame[3];

							F2->counterpartFrameIsOpposite[0] = false;
							F2->counterpartFrameIsOpposite[1] = false;
							F2->counterpartFrameIsOpposite[2] = false;
							F2->counterpartFrameIsOpposite[3] = false;

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

		int coloring_surface()
		{
			int global_cnt = 0;
			for (auto b : m_blk)
			{
				for (short j = 1; j <= Block3D::NumOfSurf; ++j)
				{
					auto s = &b->surf(j);
					if (s->global_index != 0)
						continue;

					s->global_index = ++global_cnt;
					if (s->neighbourSurf)
						s->neighbourSurf->global_index = s->global_index;
				}
			}
			return global_cnt;
		}

		int coloring_frame()
		{
			int global_cnt = 0;
			for (auto b : m_blk)
			{
				for (short j = 1; j <= Block3D::NumOfFrame; ++j)
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
				for (short j = 1; j <= Block3D::NumOfVertex; ++j)
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

		void numbering_cell()
		{
			const auto totalCellNum = nCell();

			size_t cnt = 0;
			for (auto b : m_blk)
				for (size_t k = 1; k < b->KDIM(); ++k)
					for (size_t j = 1; j < b->JDIM(); ++j)
						for (size_t i = 1; i < b->IDIM(); ++i)
							b->cell(i, j, k).CellSeq() = ++cnt;

			if (cnt != totalCellNum)
				throw std::length_error("Inconsistent num of cells.");
		}

		void numbering_face()
		{
			size_t totalFaceNum = 0, innerFaceNum = 0, bdryFaceNum = 0;
			nFace(totalFaceNum, innerFaceNum, bdryFaceNum);

			size_t cnt = 0;
			for (auto b : m_blk)
			{
				/* Internal faces */
				// K - direction
				for (size_t k = 1; k <= b->KDIM() - 2; ++k)
					for (size_t j = 1; j <= b->JDIM() - 1; ++j)
						for (size_t i = 1; i <= b->IDIM() - 1; ++i)
							b->cell(i, j, k + 1).FaceSeq(1) = b->cell(i, j, k).FaceSeq(2) = ++cnt;

				// I - direction
				for (size_t i = 1; i <= b->IDIM() - 2; ++i)
					for (size_t k = 1; k <= b->KDIM() - 1; ++k)
						for (size_t j = 1; j <= b->JDIM() - 1; ++j)
							b->cell(i + 1, j, k).FaceSeq(3) = b->cell(i, j, k).FaceSeq(4) = ++cnt;

				// J - direction
				for (size_t j = 1; j <= b->JDIM() - 2; ++j)
					for (size_t i = 1; i <= b->IDIM() - 1; ++i)
						for (size_t k = 1; k <= b->KDIM() - 1; ++k)
							b->cell(i, j + 1, k).FaceSeq(5) = b->cell(i, j, k).FaceSeq(6) = ++cnt;

				/* External faces */
				// Single-Sided
				for (short f = 1; f <= Block3D::NumOfSurf; ++f)
				{
					auto &sf = b->surf(f);
					if (!sf.neighbourSurf)
					{
						if (sf.local_index == 1)
						{
							for (size_t j = 1; j <= b->JDIM() - 1; ++j)
								for (size_t i = 1; i <= b->IDIM() - 1; ++i)
									b->cell(i, j, 1).FaceSeq(1) = ++cnt;
						}
						else if (sf.local_index == 2)
						{
							for (size_t j = 1; j <= b->JDIM() - 1; ++j)
								for (size_t i = 1; i <= b->IDIM() - 1; ++i)
									b->cell(i, j, b->KDIM() - 1).FaceSeq(2) = ++cnt;
						}
						else if (sf.local_index == 3)
						{
							for (size_t k = 1; k <= b->KDIM() - 1; ++k)
								for (size_t j = 1; j <= b->JDIM() - 1; ++j)
									b->cell(1, j, k).FaceSeq(3) = ++cnt;
						}
						else if (sf.local_index == 4)
						{
							for (size_t k = 1; k <= b->KDIM() - 1; ++k)
								for (size_t j = 1; j <= b->JDIM() - 1; ++j)
									b->cell(b->IDIM() - 1, j, k).FaceSeq(4) = ++cnt;
						}
						else if (sf.local_index == 5)
						{
							for (size_t i = 1; i <= b->IDIM() - 1; ++i)
								for (size_t k = 1; k <= b->KDIM() - 1; ++k)
									b->cell(i, 1, k).FaceSeq(5) = ++cnt;
						}
						else if (sf.local_index == 6)
						{
							for (size_t i = 1; i <= b->IDIM() - 1; ++i)
								for (size_t k = 1; k <= b->KDIM() - 1; ++k)
									b->cell(i, b->JDIM() - 1, k).FaceSeq(6) = ++cnt;
						}
						else
							throw std::invalid_argument("Internal error: Wrong local index of block surface.");
					}
				}
			}

			// Double-Sided
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = static_cast<DoubleSideEntry*>(e);
					const auto &rg1 = p->Range1();
					const auto &rg2 = p->Range2();
					auto b1 = &block(rg1.B());
					const auto f1 = rg1.F();
					auto b2 = &block(rg2.B());
					const auto f2 = rg2.F();

					std::vector<size_t> b1_dim_pri, b1_dim_sec, b2_dim_pri, b2_dim_sec;
					distribute_index(rg1.S1(), rg1.E1(), b1_dim_pri);
					distribute_index(rg1.S2(), rg1.E2(), b1_dim_sec);
					distribute_index(rg2.S1(), rg2.E1(), b2_dim_pri);
					distribute_index(rg2.S2(), rg2.E2(), b2_dim_sec);
					if (p->Swap())
						std::swap(b2_dim_pri, b2_dim_sec);

					if (b1_dim_pri.size() != b2_dim_pri.size() || b1_dim_sec.size() != b2_dim_sec.size())
						throw std::runtime_error("Inconsistent num of nodes.");

					const auto n1 = rg1.pri_node_num();
					const auto n2 = rg1.sec_node_num();
					if (p->Swap())
					{
						for (size_t l1 = 1; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 1; l2 <= n2 - 1; ++l2)
							{
								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];
								b1->surface_face_index(f1, b1i1, b1i2) = b2->surface_face_index(f2, b2i2, b2i1) = ++cnt;
							}
					}
					else
					{
						for (size_t l1 = 1; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 1; l2 <= n2 - 1; ++l2)
							{
								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];
								b1->surface_face_index(f1, b1i1, b1i2) = b2->surface_face_index(f2, b2i1, b2i2) = ++cnt;
							}
					}
				}
			}

			if (cnt != totalFaceNum)
				throw std::length_error("Inconsistent num of faces detected.");
		}

		void numbering_node()
		{
			const auto totalNodeNum = nNode();

			size_t cnt = 0;

			// Block interior
			for (auto b : m_blk)
			{
				std::vector<size_t*> boc(8, nullptr);
				for (size_t k = 2; k <= b->KDIM() - 1; ++k)
					for (size_t j = 2; j <= b->JDIM() - 1; ++j)
						for (size_t i = 2; i <= b->IDIM() - 1; ++i)
						{
							++cnt;
							b->interior_node_occurance(i, j, k, boc);
							for (auto r : boc)
								*r = cnt;
						}
			}

			// Vertex
			for (const auto &e : m_vertex)
			{
				++cnt;
				for (auto r : e)
					r->dependentBlock->vertex_node_index(r->local_index) = cnt;
			}

			// Interior of double-sided surface
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = static_cast<DoubleSideEntry*>(e);
					const auto &rg1 = p->Range1();
					const auto &rg2 = p->Range2();
					auto b1 = &block(rg1.B());
					auto b2 = &block(rg2.B());
					const auto f1 = rg1.F();
					const auto f2 = rg2.F();
					const auto n1 = rg1.pri_node_num();
					const auto n2 = rg1.sec_node_num();

					std::vector<size_t> b1_dim_pri, b1_dim_sec, b2_dim_pri, b2_dim_sec;
					distribute_index(rg1.S1(), rg1.E1(), b1_dim_pri);
					distribute_index(rg1.S2(), rg1.E2(), b1_dim_sec);
					distribute_index(rg2.S1(), rg2.E1(), b2_dim_pri);
					distribute_index(rg2.S2(), rg2.E2(), b2_dim_sec);
					if (p->Swap())
						std::swap(b2_dim_pri, b2_dim_sec);

					if (b1_dim_pri.size() != b2_dim_pri.size() || b1_dim_sec.size() != b2_dim_sec.size())
						throw std::runtime_error("Inconsistent num of nodes.");

					std::vector<size_t*> sioc(4, nullptr);
					if (p->Swap())
					{
						for (size_t l1 = 2; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 2; l2 <= n2 - 1; ++l2)
							{
								++cnt;

								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];

								b1->surface_internal_node_occurance(f1, b1i1, b1i2, sioc);
								for (auto r : sioc)
									*r = cnt;

								b2->surface_internal_node_occurance(f2, b2i2, b2i1, sioc);
								for (auto r : sioc)
									*r = cnt;
							}
					}
					else
					{
						for (size_t l1 = 2; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 2; l2 <= n2 - 1; ++l2)
							{
								++cnt;

								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];

								b1->surface_internal_node_occurance(f1, b1i1, b1i2, sioc);
								for (auto r : sioc)
									*r = cnt;

								b2->surface_internal_node_occurance(f2, b2i1, b2i2, sioc);
								for (auto r : sioc)
									*r = cnt;
							}
					}
				}
			}

			// Interior of single-sided surface
			for (auto b : m_blk)
			{
				std::vector<size_t*> sioc(4, nullptr);
				for (short i = 1; i <= Block3D::NumOfSurf; ++i)
				{
					auto &f = b->surf(i);
					if (!f.neighbourSurf)
					{
						size_t pri_end = 0, sec_end = 0;
						switch (i)
						{
						case 1:
						case 2:
							pri_end = b->IDIM() - 1;
							sec_end = b->JDIM() - 1;
							break;
						case 3:
						case 4:
							pri_end = b->JDIM() - 1;
							sec_end = b->KDIM() - 1;
							break;
						case 5:
						case 6:
							pri_end = b->KDIM() - 1;
							sec_end = b->IDIM() - 1;
							break;
						default:
							break;
						}
						for (size_t sec = 2; sec <= sec_end; ++sec)
							for (size_t pri = 2; pri <= pri_end; ++pri)
							{
								++cnt;
								b->surface_internal_node_occurance(i, pri, sec, sioc);
								for (auto r : sioc)
									*r = cnt;
							}
					}
				}
			}

			// Frame
			for (const auto &e : m_frame)
			{
				std::vector<size_t*> fnoc(2, nullptr);

				// Process the first frame in this group
				auto r = e[0];
				auto b = r->dependentBlock;
				const auto itn = b->frame_internal_node_num(r->local_index);
				for (size_t lidx = 0; lidx < itn; ++lidx)
				{
					b->frame_internal_node_occurace(r->local_index, lidx + 2, fnoc);
					*fnoc[0] = *fnoc[1] = ++cnt;
				}

				// Process remaining frames
				for (size_t i = 1; i < e.size(); ++i)
				{
					r = e[i];
					b = r->dependentBlock;

					// TODO
				}
			}

			if (cnt != totalNodeNum)
				throw std::length_error("Inconsistent num of nodes detected.");
		}
	};
}

#endif
