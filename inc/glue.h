#ifndef __GLUE_H__
#define __GLUE_H__

#include <string>
#include "nmf.h"
#include "plot3d.h"
#include "xf_msh.h"

int glue(const std::string &f_nmf, const std::string &f_p3d, std::string &f_msh);
int glue(const NMF &src1, const PLOT3D &src2, XF_MSH &dst);

#endif
