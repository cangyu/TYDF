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
			if (nmf->nBlock() != p3d->numOfBlock())
				throw std::invalid_argument("Inconsistent num of blocks.");

			// Allocate storage.
			m_totalNodeNum = nmf->nNode();
			m_totalCellNum = nmf->nCell();
			size_t innerFaceNum = 0, bdryFaceNum = 0;
			nmf->nFace(m_totalFaceNum, innerFaceNum, bdryFaceNum);
			m_node.resize(numOfNode());
			m_face.resize(numOfFace());
			m_cell.resize(numOfCell());

			// Copy coordinates.

			// Convert to primary form.
			derived2raw();

			// Finalize.
			delete nmf;
			delete p3d;
		}
	}
}
