#include "../inc/glue.h"

void glue_3d(const std::string &f_nmf, const std::string &f_p3d, std::string &f_msh)
{
	auto nmf = new NMF::Mapping3D(f_nmf);
	nmf->compute_topology();
	nmf->numbering();
	auto p3d = new PLOT3D::GRID(f_p3d);
	auto msh = new XF::MESH();
	glue_3d(*nmf, *p3d, *msh);
	msh->writeToFile(f_msh);

	delete nmf;
	delete p3d;
	delete msh;
}

void glue_3d(const NMF::Mapping3D &src1, const PLOT3D::GRID &src2, XF::MESH &dst)
{
	if (src1.nBlock() != src2.numOfBlock())
		throw std::invalid_argument("Inconsistent num of blocks.");


	// TODO
}
