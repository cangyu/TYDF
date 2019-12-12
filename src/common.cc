#include <stdexcept>
#include "../inc/common.h"

namespace GridTool
{
	namespace COMMON
	{
		double relaxation(double a, double b, double x)
		{
			return (1.0 - x) * a + x * b;
		}

		DIM::DIM(int dim, bool is3d) :
			m_is3D(is3d)
		{
			if (dim == 2 || dim == 3)
				m_dim = dim;
			else
				throw wrong_dimension(dim);

			if (dim == 3 && !is3d)
				throw std::invalid_argument("Inconsistent dimensions.");
		}
	}
}