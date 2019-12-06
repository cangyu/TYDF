#ifndef __GT_NMF_H__
#define __GT_NMF_H__

#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cctype>
#include <algorithm>
#include <cstddef>
#include <vector>
#include <utility>
#include <stdexcept>
#include <regex>
#include "common.h"

#define BC_ENUM { UNPROCESSED, ONE_TO_ONE, SYM, WALL, INFLOW, OUTFLOW, FAR }
#define BC_STR { "UNPROCESSED", "ONE_TO_ONE", "SYM", "WALL", "INFLOW", "OUTFLOW", "FAR" }

namespace GridTool
{
	namespace NMF
	{
		using COMMON::wrong_index;
		using COMMON::Array1D;

		class BC
		{
		public:
			enum BC_ENUM;

			static bool isValidBCIdx(int x);
			static bool isValidBCStr(const std::string &x);

			static const std::string &idx2str(int x);
			static int str2idx(const std::string &x);

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
			BLOCK(size_t nI, size_t nJ) : m_idx(0), m_name(""), m_dim{ nI, nJ, 1 }
			{
				if (nI < 2 || nJ < 2)
					throw std::invalid_argument("Invalid dimension.");
			}
			BLOCK(size_t nI, size_t nJ, size_t nK) : m_idx(0), m_name(""), m_dim{ nI, nJ, nK }
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
			virtual size_t node_num() const = 0;
			virtual size_t face_num() const = 0;
			virtual size_t cell_num() const = 0;
			virtual size_t block_internal_node_num() const = 0;
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

			size_t node_num() const
			{
				return IDIM() * JDIM();
			}

			size_t face_num() const
			{
				size_t ret = 0;
				ret += (IDIM() - 1) * JDIM();
				ret += IDIM() * (JDIM() - 1);
				return ret;
			}

			size_t cell_num() const
			{
				return (IDIM() - 1) * (JDIM() - 1);
			}

			size_t block_internal_node_num() const
			{
				return (IDIM() - 2) * (JDIM() - 2);
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

			size_t node_num() const
			{
				return IDIM() * JDIM() * KDIM();
			}

			size_t face_num() const
			{
				size_t ret = 0;
				ret += (IDIM() - 1) * (JDIM() - 1) * KDIM();
				ret += IDIM() * (JDIM() - 1) * (KDIM() - 1);
				ret += (IDIM() - 1) *JDIM() * (KDIM() - 1);
				return ret;
			}

			size_t cell_num() const
			{
				return (IDIM() - 1) * (JDIM() - 1) * (KDIM() - 1);
			}

			size_t block_internal_node_num() const
			{
				return (IDIM() - 2) * (JDIM() - 2) * (KDIM() - 2);
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
			class not_a_block : public wrong_index
			{
			public:
				not_a_block(size_t x) : wrong_index(x)
				{
					m_msg += "is not a valid 1-based block index.";
				}
			};

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

			void readFromFile(const std::string &path);

			void compute_topology();

			void summary(std::ostream &out);

			void numbering();

			void writeToFile(const std::string &path);

			size_t nBlock() const { return m_blk.size(); }

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

			size_t nFrame() const { return m_frame.size(); }

			size_t nVertex() const { return m_vertex.size(); }

			size_t nCell() const
			{
				size_t ret = 0;
				for (auto b : m_blk)
					ret += b->cell_num();
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
					throw not_a_block(n);
			}
			const Block3D &block(size_t n) const
			{
				if (1 <= n && n <= nBlock())
					return *m_blk[n - 1];
				else
					throw not_a_block(n);
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

			void connecting();

			int coloring_surface();
			int coloring_frame();
			int coloring_vertex();

			void numbering_cell();
			void numbering_face();
			void numbering_node();
		};
	}
}

#endif
