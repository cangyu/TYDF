#include "../inc/nmf.h"
#include "../inc/plot3d.h"
#include "../inc/xf.h"

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
			// Load mapping file.
			auto nmf = new NMF::Mapping3D(f_nmf);
			nmf->compute_topology();
			nmf->numbering();

			// Load grid file.
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

			// Allocate storage.
			m_totalNodeNum = nmf->nNode();
			m_totalCellNum = nmf->nCell();
			size_t innerFaceNum = 0, bdryFaceNum = 0;
			nmf->nFace(m_totalFaceNum, innerFaceNum, bdryFaceNum);
			m_node.resize(numOfNode());
			m_face.resize(numOfFace());
			m_cell.resize(numOfCell());

			// Copy coordinates.
			for (size_t n = 1; n <= NBLK; ++n)
			{
				const auto &b = nmf->block(n);
				auto g = p3d->block(n - 1);

				const size_t nI = b.IDIM();
				const size_t nJ = b.JDIM();
				const size_t nK = b.KDIM();

				for (size_t k = 1; k < nK; ++k)
					for (size_t j = 1; j < nJ; ++j)
						for (size_t i = 1; i < nI; ++i)
						{
							const auto &c = b.cell(i, j, k);
							for (short l = 1; l <= NMF::Block3D::NumOfVertex; ++l)
							{
								auto ci = c.NodeSeq(l);
							}
						}
			}

			// Copy connections.


			// Convert to primary form.
			derived2raw();

			// Finalize.
			delete nmf;
			delete p3d;
		}
	}
}
