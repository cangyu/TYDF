#include <cmath>
#include <algorithm>
#include "../inc/spacing.h"

namespace GridTool
{
	namespace SPACING
	{
		std::vector<double> uniform(size_t n)
		{
			std::vector<double> ret(n, 0.0);

			return ret;
		}

		std::vector<double> uniform(double start, double end, size_t n)
		{
			std::vector<double> ret(n, start);

			return ret;
		}

		std::vector<double> expansion(const std::vector<double> &seq, double start, double end)
		{
			std::vector<double> ret;

			return ret;
		}

		std::vector<double> chebshev(double start, double end, size_t n)
		{
			std::vector<double> ret;

			return ret;
		}

		std::vector<double> chebshev_multi(const std::vector<double> &seg, const std::vector<size_t> &n)
		{
			std::vector<double> ret;

			return ret;
		}

		std::vector<double> single_exponential(size_t n, double a)
		{
			std::vector<double> ret;

			return ret;
		}

		std::vector<double> double_exponential(size_t n, double a1, double a2, double a3)
		{
			std::vector<double> ret;

			return ret;
		}

		std::vector<double> hyperbolic_tangent(size_t n, double b)
		{
			std::vector<double> ret;

			return ret;
		}

		std::vector<double> hyperbolic_sine(size_t n, double c)
		{
			std::vector<double> ret;

			return ret;
		}
	}
}
