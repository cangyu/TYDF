#ifndef __GT_SPACING_H__
#define __GT_SPACING_H__

#include <vector>
#include <cstddef>

/*
	All nodes distribute through [0, 1] by default unless specified.
*/

namespace GridTool
{
	namespace SPACING
	{
		std::vector<double> uniform(size_t n);
		std::vector<double> uniform(double start, double end, size_t n);
		std::vector<double> chebshev(double start, double end, size_t n);
		std::vector<double> chebshev_multi(const std::vector<double> &seg, const std::vector<size_t> &n);
		std::vector<double> single_exponential(size_t n, double a);
		std::vector<double> double_exponential(size_t n, double a1, double a2, double a3);
		std::vector<double> hyperbolic_tangent(size_t n, double b);
		std::vector<double> hyperbolic_sine(size_t n, double c);
		std::vector<double> linear_expand(const std::vector<double> &seq, double start, double end);
	}
}

#endif
