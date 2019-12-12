#include <cmath>
#include <algorithm>
#include "../inc/common.h"
#include "../inc/spacing.h"

namespace GridTool
{
	namespace SPACING
	{
		using COMMON::relaxation;

		int uniform(size_t n, std::vector<double> &dst)
		{
			if (n == 0)
				return -1;
			if (n == 1)
				return -2;

			dst.resize(n);
			for (size_t i = 0; i < n; ++i)
			{
				const double ratio = 1.0 * i / (n - 1);
				dst[i] = relaxation(0, 1, ratio);
			}
			return 0;
		}

		int uniform(double start, double end, size_t n, std::vector<double> &dst)
		{
			if (n == 0)
				return -1;
			if (n == 1)
				return -2;

			dst.resize(n);
			for (size_t i = 0; i < n; ++i)
			{
				const double ratio = 1.0 * i / (n - 1);
				dst[i] = relaxation(start, end, ratio);
			}
			return 0;
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
