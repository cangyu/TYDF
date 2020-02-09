#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include "../inc/common.h"
#include "../inc/spacing.h"

struct invalid_node_num : public std::invalid_argument
{
    explicit invalid_node_num(int n) : std::invalid_argument("At least 2 nodes are required, got \"" + std::to_string(n) + "\".") {}
};

namespace GridTool
{
    namespace SPACING
    {
        using COMMON::relaxation;

        /**
         * Uniform distribution through [0,1].
         * @param n Num of nodes.
         * @param dst Target distribution.
         */
        void uniform(int n, DIST_ARR &dst)
        {
            if (n < 2)
                throw invalid_node_num(n);

            if (dst.size() != n)
                dst.resize(n);

            for (int i = 0; i < n; ++i)
                dst[i] = 1.0 * i / (n - 1);
        }

        /**
         * Uniform distribution through [a,b].
         * @param a The starting position.
         * @param b The ending position.
         * @param n Num of nodes.
         * @param dst Target distribution.
         */
        void uniform(double a, double b, int n, DIST_ARR &dst)
        {
            uniform(n, dst);

            for (auto &e : dst)
                e = relaxation(a, b, e);
        }

        /**
         * Chebshev distribution through [a,b].
         * @param a The starting position.
         * @param b The ending position.
         * @param n Num of nodes.
         * @param dst Target distribution.
         */
        void chebshev(double a, double b, int n, DIST_ARR &dst)
        {
            static const double pi = std::acos(-1.0);

            uniform(pi, 0.0, n, dst);

            for (auto &e : dst)
            {
                const auto ratio = 0.5 * (1.0 + std::cos(e));
                e = relaxation(a, b, ratio);
            }
        }

        /**
         * Multi-segment Chebshev distribution.
         * @param seg Splitting values.
         * @param num Num of nodes within each segment.
         * @param dst Target distribution.
         */
        void chebshev(const std::vector<double> &seg, const std::vector<int> &num, DIST_ARR &dst)
        {
            if (seg.size() < 2)
                throw std::invalid_argument("Too less splitting values.");
            if (seg.size() != num.size())
                throw std::invalid_argument("Inconsistent splitting intervals.");

            int n = 1;
            for (auto e : num)
            {
                if (e < 2)
                    throw invalid_node_num(e);

                n += (e - 1);
            }

            if (dst.size() != n)
                dst.resize(n);

            dst[0] = seg[0];
            size_t pos = 1;

            for (size_t i = 0; i < num.size(); ++i)
            {
                const auto cn = num[i];
                std::vector<double> ccd(cn, 0.0);
                chebshev(seg[i], seg[i + 1], cn, ccd);
                for (int j = 1; j < cn; ++j)
                    dst[pos++] = ccd[j];
            }

            if (pos != n)
                throw std::runtime_error("Unexpected terminating position.");
        }

        /**
         * Single exponential distribution of n nodes in [0, 1].
         * @param n Number of nodes.
         * @param a The control parameter.
         *          For a>0, nodes aggregate towards the starting position;
         *          For a=0, nodes are uniform, no aggregation;
         *          For a<0, nodes aggregate towards the ending position.
         * @param dst Target distribution.
         */
        void single_exponential(int n, double a, DIST_ARR &dst)
        {
            uniform(n, dst);

            if (std::abs(a) < 1e-12)
                throw std::invalid_argument("\"a\" shouldn't be 0.");
            else
            {
                const double t = std::exp(a) - 1.0;
                for (auto &e : dst)
                    e = (std::exp(a * e) - 1.0) / t;
            }
        }

        /**
         * Double exponential distribution of n nodes in [0, 1].
         * The curve will go through (a3, a1).
         * @param n Number of nodes.
         * @param a1 Horizontal coordinate of the control point.
         * @param a2 The control parameter.
         *           For a2>0, nodes aggregate towards boundary;
         *			 For a2=0, nodes are uniform, no aggregation;
         *           For a2<0, nodes aggregate towards center.
         * @param a3 Vertical coordinate of the control point.
         * @param dst Target distribution.
         */
        void double_exponential(int n, double a1, double a2, double a3, DIST_ARR &dst)
        {
            uniform(n, dst);

            if (std::abs(a2) < 1e-12)
                throw std::invalid_argument("\"a2\" shouldn't be 0.");

            const double p = (1 - a1) * a3 * (std::exp(a2) - 1) / ((1 - a3) * a1 * a2 * std::exp(a2));
            if (p <= 0.0)
                throw std::invalid_argument("\"p\" should be positive.");

            // Solve "a4" using Newton-Raphson iteration.
            static const size_t MAX_ITER = 50;
            static const double PRECISION = -6.0;
            double a4 = 5 * std::log(p);
            bool ok = false;
            for (size_t i = 0; i < MAX_ITER; ++i)
            {
                const double f = std::exp(a4) - 1.0 - p * a4;
                if (std::log10(std::abs(f)) < PRECISION)
                {
                    ok = true;
                    break;
                }

                const double df = std::exp(a4) - p;
                a4 -= f / df;
            }
            if (!ok)
                throw std::runtime_error("Newton-Raphson iteration failed to converge.");

            const double ea21 = std::exp(a2) - 1.0;
            const double ea41 = std::exp(a4) - 1.0;
            for (auto &e : dst)
            {
                if (e <= a3)
                    e = a1 * (std::exp(a2 / a3 * e) - 1.0) / ea21;
                else
                    e = a1 + (1.0 - a1) * (std::exp(a4 / (1.0 - a3) * (e - a3)) - 1.0) / ea41;
            }
        }

        /**
         *  Hyperbolic tangent distribution of n nodes in [0, 1].
         * @param n Number of nodes.
         * @param b The control parameter.
         * @param dst Target distribution.
         */
        void hyperbolic_tangent(int n, double b, DIST_ARR &dst)
        {
            uniform(n, dst);

            const double tb = std::tanh(b);
            for (auto &e : dst)
                e = 1.0 + std::tanh(b * (e - 1.0)) / tb;
        }

        /**
         *  Hyperbolic sine distribution of n nodes in [0, 1].
         * @param n Number of nodes.
         * @param c The control parameter.
         * @param dst Target distribution.
         */
        void hyperbolic_sine(int n, double c, DIST_ARR &dst)
        {
            uniform(n, dst);

            const double sc = std::sinh(c);
            for (auto &e : dst)
                e = 1.0 + std::sinh(c * (e - 1.0)) / sc;
        }
    }
}
