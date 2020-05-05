#ifndef TYDF_SPACING_H
#define TYDF_SPACING_H

#include <vector>

typedef std::vector<double> DIST_ARR;

namespace GridTool::SPACING
{
    /// All nodes distribute through [0, 1] by default unless specified.

    /**
     * Uniform distribution through [0,1].
     * @param n Num of nodes.
     * @param dst Target distribution.
     */
    void uniform(int n, DIST_ARR &dst);

    /**
     * Uniform distribution through [a,b].
     * @param a The starting position.
     * @param b The ending position.
     * @param n Num of nodes.
     * @param dst Target distribution.
     */
    void uniform(double a, double b, int n, DIST_ARR &dst);

    /**
     * Chebshev distribution through [a,b].
     * @param a The starting position.
     * @param b The ending position.
     * @param n Num of nodes.
     * @param dst Target distribution.
     */
    void chebshev(double a, double b, int n, DIST_ARR &dst);

    /**
     * Multi-segment Chebshev distribution.
     * @param seg Splitting values.
     * @param num Num of nodes within each segment.
     * @param dst Target distribution.
     */
    void chebshev(const std::vector<double> &seg, const std::vector<int> &num, DIST_ARR &dst);

    /**
     * Single exponential distribution of n nodes in [0, 1].
     * @param n Number of nodes.
     * @param a The control parameter.
     *          For a>0, nodes aggregate towards the starting position;
     *          For a=0, nodes are uniform, no aggregation;
     *          For a<0, nodes aggregate towards the ending position.
     * @param dst Target distribution.
     */
    void single_exponential(int n, double a, DIST_ARR &dst);

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
    void double_exponential(int n, double a1, double a2, double a3, DIST_ARR &dst);

    /**
     * Hyperbolic tangent distribution of n nodes in [0, 1].
     * @param n Number of nodes.
     * @param b The control parameter.
     * @param dst Target distribution.
     */
    void hyperbolic_tangent(int n, double b, DIST_ARR &dst);

    /**
     * Hyperbolic sine distribution of n nodes in [0, 1].
     * @param n Number of nodes.
     * @param c The control parameter.
     * @param dst Target distribution.
     */
    void hyperbolic_sine(int n, double c, DIST_ARR &dst);
}

#endif
