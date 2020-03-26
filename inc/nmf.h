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
			BLOCK(size_t nI, size_t nJ);
			BLOCK(size_t nI, size_t nJ, size_t nK);
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

			size_t node_num() const { return IDIM() * JDIM(); }

			size_t face_num() const { return (IDIM() - 1) * JDIM() + IDIM() * (JDIM() - 1); }

			size_t cell_num() const { return (IDIM() - 1) * (JDIM() - 1); }

			size_t block_internal_node_num() const { return (IDIM() - 2) * (JDIM() - 2); }

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
			struct not_a_surface : public wrong_index
			{
				not_a_surface(short x) : wrong_index(x, "is not a valid surface index of a 3D block") {}
			};
			struct not_a_frame : public wrong_index
			{
				not_a_frame(short x) : wrong_index(x, "is not a valid frame index of a 3D block") {}
			};
			struct not_a_vertex : public wrong_index
			{
				not_a_vertex(short x) : wrong_index(x, "is not a valid vertex index of a 3D block") {}
			};

		private:
			Array1D<HEX_CELL*> m_cell;
			Array1D<VERTEX> m_vertex;
			Array1D<FRAME> m_frame;
			Array1D<SURF> m_surf;

		public:
			Block3D() = delete;
			Block3D(int nI, int nJ, int nK);
			Block3D(const Block3D &rhs);
			~Block3D() { release_cell_storage(); }

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

			size_t node_num() const { return IDIM() * JDIM() * KDIM(); }

			size_t face_num() const
			{
				size_t ret = 0;
				ret += (IDIM() - 1) * (JDIM() - 1) * KDIM();
				ret += IDIM() * (JDIM() - 1) * (KDIM() - 1);
				ret += (IDIM() - 1) *JDIM() * (KDIM() - 1);
				return ret;
			}

			size_t cell_num() const { return (IDIM() - 1) * (JDIM() - 1) * (KDIM() - 1); }

			size_t block_internal_node_num() const { return (IDIM() - 2) * (JDIM() - 2) * (KDIM() - 2); }

			size_t surface_internal_node_num(short s) const;

			size_t frame_node_num(short idx) const;

			size_t frame_internal_node_num(short idx) const
			{
				const auto n = frame_node_num(idx);
				if (n >= 2)
					return n - 2;
				else
					throw std::runtime_error("Internal error: too less nodes not detected when constructing.");
			}

			size_t surface_node_num(short idx) const;

			size_t surface_face_num(short idx) const;

			size_t shell_face_num() const
			{
				size_t ret = 0;
				for (short i = 1; i <= NumOfSurf; ++i)
					ret += surface_face_num(i);
				return ret;
			}

			size_t &surface_face_index(short f, size_t pri, size_t sec);

			size_t &vertex_node_index(short v);

			void interior_node_occurance(size_t i, size_t j, size_t k, std::vector<size_t*> &oc);

			void surface_node_coordinate(short f, size_t pri_seq, size_t sec_seq, size_t &i, size_t &j, size_t &k);

			void surface_internal_node_occurance(short f, size_t pri, size_t sec, std::vector<size_t*> &oc);

			void frame_internal_node_occurace(short f, size_t idx, std::vector<size_t*> &oc);

			size_t node_index(size_t i, size_t j, size_t k);

		private:
            void setup_dependence();

            void establish_connections();

			short isVertexNode(size_t i, size_t j, size_t k) const;

			void isFrameInternalNode(size_t i, size_t j, size_t k, short &f, size_t &idx) const;

			void isSurfaceInternalNode(size_t i, size_t j, size_t k, short &s, size_t &pri, size_t &sec) const;
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
			struct not_a_block : public wrong_index
			{
				not_a_block(size_t x) : wrong_index(x, "is not a valid 1-based block index") {}
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
					bool pri_trend() const { return E1() > S1(); }

					// True - Asscending;
					// False - Descending.
					bool sec_trend() const { return E2() > S2(); }

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
				ENTRY(const std::string &t, size_t B, short F, size_t S1, size_t E1, size_t S2, size_t E2);
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
				DoubleSideEntry(const std::string &t, size_t B1, short F1, size_t S11, size_t E11, size_t S12, size_t E12, size_t B2, short F2, size_t S21, size_t E21, size_t S22, size_t E22, bool f);
				DoubleSideEntry(const DoubleSideEntry &rhs) = default;
				~DoubleSideEntry() = default;

				RANGE &Range2() { return m_rg2; }
				const RANGE &Range2() const { return m_rg2; }

				bool Swap() const { return m_swap; }
				bool &Swap() { return m_swap; }

				int contains(size_t bs, short fs, size_t lpri, size_t lsec) const;

			private:
				bool check_dim_consistency() const;
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

				compute_topology();
			}

			~Mapping3D() { release_all(); }

			void readFromFile(const std::string &path);

			void compute_topology();

			void summary(std::ostream &out);

			void numbering();

			void writeToFile(const std::string &path);

			size_t nBlock() const { return m_blk.size(); }

			void nSurface(size_t &_all, size_t &_inter, size_t &_bdry) const;

			size_t nFrame() const { return m_frame.size(); }

			size_t nVertex() const { return m_vertex.size(); }

			size_t nCell() const
			{
				size_t ret = 0;
				for (auto b : m_blk)
					ret += b->cell_num();
				return ret;
			}

			void nFace(size_t &_all, size_t &_inner, size_t &_bdry) const;

			size_t nNode() const;

			// 1-based indexing
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
			void release_all();

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
