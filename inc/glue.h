#ifndef __GLUE_HPP__
#define __GLUE_HPP__

#include "nmf.hpp"
#include "plot3d.hpp"
#include "xf.hpp"

void glue_3d(const std::string &f_nmf, const std::string &f_p3d, std::string &f_msh);

void glue_3d(const NMF::Mapping3D &src1, const PLOT3D::GRID &src2, XF::MESH &dst);

#endif
