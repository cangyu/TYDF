#ifndef __GLUE_HPP__
#define __GLUE_HPP__

#include "nmf.hpp"
#include "plot3d.hpp"
#include "xf.hpp"

void glue(const std::string &f_nmf, const std::string &f_p3d, std::string &f_msh)
{
	auto nmf = new NMF::Mapping3D(f_nmf);
	nmf->compute_topology();
	nmf->numbering();
	auto p3d = new PLOT3D::GRID(f_p3d);
	auto msh = new XF::MESH();
	glue(*nmf, *p3d, *msh);
	msh->writeToFile(f_msh);

	delete nmf;
	delete p3d;
	delete msh;
}

void glue(const NMF::Mapping3D &src1, const PLOT3D::GRID &src2, XF::MESH &dst)
{
	if (src1.nBlock() != src2.numOfBlock())
		throw std::invalid_argument("Inconsistent num of blocks.");

}

#endif
