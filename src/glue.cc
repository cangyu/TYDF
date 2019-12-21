#include "../inc/nmf.h"
#include "../inc/plot3d.h"
#include "../inc/xf.h"

static const size_t major = 2;
static const size_t minor = 0;
static const size_t patch = 0;

static std::string version()
{
	return std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(patch);
}

namespace GridTool
{
	namespace XF
	{
		MESH::MESH(const std::string &f_nmf, const std::string &f_p3d, std::ostream &fout) :
			DIM(3),
			m_totalNodeNum(0),
			m_totalCellNum(0),
			m_totalFaceNum(0),
			m_totalZoneNum(0)
		{
			// Load mapping file
			auto nmf = new NMF::Mapping3D(f_nmf);
			nmf->compute_topology();
			nmf->numbering();

			// Load grid file
			auto p3d = new PLOT3D::GRID(f_p3d);

			// Check consistency
			const size_t NBLK = nmf->nBlock();
			if (NBLK != p3d->numOfBlock())
				throw std::invalid_argument("Inconsistent num of blocks between NMF and PLOT3D.");
			for (size_t n = 1; n <= NBLK; ++n)
			{
				const auto &b = nmf->block(n);
				auto g = p3d->block(n - 1);
				if (b.IDIM() != g->nI())
					throw std::invalid_argument("Inconsistent num of nodes in I dimension of Block " + std::to_string(n) + ".");
				if (b.JDIM() != g->nJ())
					throw std::invalid_argument("Inconsistent num of nodes in J dimension of Block " + std::to_string(n) + ".");
				if (b.KDIM() != g->nK())
					throw std::invalid_argument("Inconsistent num of nodes in K dimension of Block " + std::to_string(n) + ".");
			}

			// Allocate storage
			m_totalNodeNum = nmf->nNode();
			m_totalCellNum = nmf->nCell();
			size_t innerFaceNum = 0, bdryFaceNum = 0;
			nmf->nFace(m_totalFaceNum, innerFaceNum, bdryFaceNum);
			m_node.resize(numOfNode());
			m_face.resize(numOfFace());
			m_cell.resize(numOfCell());

			// Copy node info
			std::vector<bool> visited(m_node.size(), false);
			for (size_t n = 1; n <= NBLK; ++n)
			{
				auto &b = nmf->block(n);
				auto &g = *p3d->block(n - 1);

				const size_t nI = b.IDIM();
				const size_t nJ = b.JDIM();
				const size_t nK = b.KDIM();

				for (size_t k = 1; k <= nK; ++k)
					for (size_t j = 1; j <= nJ; ++j)
						for (size_t i = 1; i <= nI; ++i)
						{
							const auto idx = b.node_index(i, j, k); // Global 1-based index, already assigned.
							if (!visited[idx - 1])
							{
								node(idx).coordinate = g(i, j, k);
								visited[idx - 1] = true;
							}
						}
			}

			// Copy cell info
			for (size_t n = 1; n <= NBLK; ++n)
			{
				auto &b = nmf->block(n);

				const size_t nI = b.IDIM();
				const size_t nJ = b.JDIM();
				const size_t nK = b.KDIM();

				for (size_t k = 1; k < nK; ++k)
					for (size_t j = 1; j < nJ; ++j)
						for (size_t i = 1; i < nI; ++i)
						{
							const auto &nc = b.cell(i, j, k);
							auto &fc = cell(nc.CellSeq());

							fc.type = CELL::HEXAHEDRAL;

							fc.includedFace.resize(NMF::Block3D::NumOfSurf);
							for (short r = 1; r <= NMF::Block3D::NumOfSurf; ++r)
								fc.includedFace(r) = nc.FaceSeq(r);

							fc.includedNode.resize(NMF::Block3D::NumOfVertex);
							for (short r = 1; r <= NMF::Block3D::NumOfVertex; ++r)
								fc.includedNode(r) = nc.NodeSeq(r);
						}
			}

			// Copy face info
			visited.resize(m_face.size(), false);
			std::fill(visited.begin(), visited.end(), false);
			for (size_t n = 1; n <= NBLK; ++n)
			{
				auto &b = nmf->block(n);

				const size_t nI = b.IDIM();
				const size_t nJ = b.JDIM();
				const size_t nK = b.KDIM();

				// TODO
			}

			// Convert to primary form
			derived2raw();

			// Finalize
			delete nmf;
			delete p3d;
		}

		void MESH::derived2raw()
		{
			clear_entry();

			add_entry(new HEADER("Grid-Glue V" + version()));
			add_entry(new DIMENSION(3));

			auto pnt = new NODE(1, 1, numOfNode(), NODE::ANY, 3);

			// TODO
		}
	}
}
