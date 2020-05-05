#ifndef TYDF_SPACING_H
#define TYDF_SPACING_H

#include <vector>

typedef std::vector<double> DIST_ARR;

namespace GridTool::SPACING
{
    /// All nodes distribute through [0, 1] by default unless specified.
    void uniform(int n, DIST_ARR &dst);
    void uniform(double a, double b, int n, DIST_ARR &dst);
    void chebshev(double a, double b, int n, DIST_ARR &dst);
    void chebshev(const std::vector<double> &seg, const std::vector<int> &num, DIST_ARR &dst);
    void single_exponential(int n, double a, DIST_ARR &dst);
    void double_exponential(int n, double a1, double a2, double a3, DIST_ARR &dst);
    void hyperbolic_tangent(int n, double b, DIST_ARR &dst);
    void hyperbolic_sine(int n, double c, DIST_ARR &dst);
}

#endif
