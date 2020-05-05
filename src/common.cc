#include "../inc/common.h"

namespace GridTool::COMMON
{
    Scalar relaxation(Scalar a, Scalar b, Scalar x)
    {
        return (1.0 - x) * a + x * b;
    }

    void delta(const Vector &na, const Vector &nb, Vector &dst)
    {
        dst = nb;
        dst -= na;
    }

    Scalar line_length(const Vector &na, const Vector &nb)
    {
        Vector tmp;
        delta(na, nb, tmp);
        const Scalar L = tmp.norm();

        return L;
    }

    void line_center(const Vector &na, const Vector &nb, Vector &dst)
    {
        dst = nb;
        dst += na;
        dst /= 2.0;
    }

    void line_normal(const Vector &na, const Vector &nb, Vector &dst_LR, Vector &dst_RL)
    {
        delta(na, nb, dst_RL);

        /// Rotate 90 deg in clockwise direction
        std::swap(dst_RL.x(), dst_RL.y());
        dst_RL.y() *= -1.0;

        dst_RL.normalize();

        dst_LR = dst_RL;
        dst_LR *= -1.0;
    }

    Scalar triangle_area(const Vector &na, const Vector &nb, const Vector &nc)
    {
        const Scalar c = line_length(na, nb);
        const Scalar a = line_length(nb, nc);
        const Scalar b = line_length(nc, na);
        const Scalar p = 0.5*(a + b + c);
        const Scalar S = std::sqrt(p*(p - a)*(p - b)*(p - c)); /// Heron's formula

        return S;
    }

    void triangle_center(const Vector &na, const Vector &nb, const Vector &nc, Vector &dst)
    {
        dst = nc;
        dst += nb;
        dst += na;
        dst /= 3.0;
    }

    void triangle_normal(const Vector &na, const Vector &nb, const Vector &nc, Vector &dst_LR, Vector &dst_RL)
    {
        Vector rab, rac;
        delta(na, nb, rab);
        delta(na, nc, rac);
        dst_LR = rab.cross(rac); /// Take cross product to find normal direction
        dst_LR.normalize();
        dst_RL = dst_LR;
        dst_RL *= -1.0;
    }

    Scalar quadrilateral_area(const Vector &n1, const Vector &n2, const Vector &n3, const Vector &n4)
    {
        const Scalar S123 = triangle_area(n1, n2, n3);
        const Scalar S134 = triangle_area(n1, n3, n4);
        const Scalar S = S123 + S134;

        return S;
    }

    void quadrilateral_center(const Vector &n1, const Vector &n2, const Vector &n3, const Vector &n4, Vector &dst)
    {
        const Scalar S123 = triangle_area(n1, n2, n3);
        const Scalar S134 = triangle_area(n1, n3, n4);

        Vector rc123, rc134;
        triangle_center(n1, n2, n3, rc123);
        triangle_center(n1, n3, n4, rc134);

        const Scalar alpha = S123 / (S123 + S134);
        const Scalar beta = 1.0 - alpha;

        rc123 *= alpha;
        rc134 *= beta;

        dst = rc123;
        dst += rc134;
    }

    void quadrilateral_normal(const Vector &n1, const Vector &n2, const Vector &n3, const Vector &n4, Vector &dst_LR, Vector &dst_RL)
    {
        Vector ra, rb;
        delta(n2, n4, ra);
        delta(n1, n3, rb);
        dst_RL = ra.cross(rb);
        dst_RL.normalize();
        dst_LR = dst_RL;
        dst_LR *= -1.0;
    }
}
