#include "../inc/common.h"

namespace GridTool::COMMON
{
    Scalar relaxation(Scalar a, Scalar b, Scalar x)
    {
        return (1.0 - x) * a + x * b;
    }

    struct DIM::wrong_dimension : public wrong_index
    {
        wrong_dimension(int dim) :
            wrong_index(dim, "is not a valid dimension")
        {
            /// Empty body.
        }
    };

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

    bool DIM::is3D() const
    {
        return m_is3D;
    }

    int DIM::dimension() const
    {
        return m_dim;
    }

    struct Vector::not_vector_component : public wrong_index
    {
        not_vector_component(short x) :
            wrong_index(x, "is not a valid index of certain vector component")
        {
            /// Empty body.
        }
    };

    Vector::Vector() :
        std::array<Scalar, 3>{0.0, 0.0, 0.0}
    {
        /// Empty body.
    }

    Vector::Vector(Scalar val) :
        std::array<Scalar, 3>{val, val, val}
    {
        /// Empty body.
    }

    Vector::Vector(Scalar v1, Scalar v2, Scalar v3) :
        std::array<Scalar, 3>{v1, v2, v3}
    {
        /// Empty body.
    }

    Vector::Vector(const Vector &obj) :
        std::array<Scalar, 3>{obj.x(), obj.y(), obj.z()}
    {
        /// Empty body.
    }

    const Scalar &Vector::operator()(short idx) const
    {
        switch (idx)
        {
        case 1:
            return x();
        case 2:
            return y();
        case 3:
            return z();
        case -3:
            return x();
        case -2:
            return y();
        case -1:
            return z();
        default:
            throw not_vector_component(idx);
        }
    }

    Scalar &Vector::operator()(short idx)
    {
        switch (idx)
        {
        case 1:
            return x();
        case 2:
            return y();
        case 3:
            return z();
        case -3:
            return x();
        case -2:
            return y();
        case -1:
            return z();
        default:
            throw not_vector_component(idx);
        }
    }

    const Scalar &Vector::x() const
    {
        return at(0);
    }

    const Scalar &Vector::y() const
    {
        return at(1);
    }

    const Scalar &Vector::z() const
    {
        return at(2);
    }

    Scalar &Vector::x()
    {
        return at(0);
    }

    Scalar &Vector::y()
    {
        return at(1);
    }

    Scalar &Vector::z()
    {
        return at(2);
    }

    Vector &Vector::operator=(const Vector &rhs)
    {
        x() = rhs.x();
        y() = rhs.y();
        z() = rhs.z();
        return *this;
    }

    Vector &Vector::operator+=(const Vector &rhs)
    {
        x() += rhs.x();
        y() += rhs.y();
        z() += rhs.z();
        return *this;
    }

    Vector &Vector::operator-=(const Vector &rhs)
    {
        x() -= rhs.x();
        y() -= rhs.y();
        z() -= rhs.z();
        return *this;
    }

    Vector &Vector::operator*=(Scalar a)
    {
        x() *= a;
        y() *= a;
        z() *= a;
        return *this;
    }

    Vector &Vector::operator/=(Scalar a)
    {
        x() /= a;
        y() /= a;
        z() /= a;
        return *this;
    }

    Scalar Vector::dot(const Vector &b) const
    {
        Scalar ret = 0.0;
        ret += x() * b.x();
        ret += y() * b.y();
        ret += z() * b.z();
        return ret;
    }

    Vector Vector::cross(const Vector &b) const
    {
        Vector ret;
        ret.x() = y() * b.z() - z() * b.y();
        ret.y() = z() * b.x() - x() * b.z();
        ret.z() = x() * b.y() - y() * b.x();
        return ret;
    }

    Scalar Vector::norm() const
    {
        Scalar ret = 0.0;
        ret += std::pow(x(), 2);
        ret += std::pow(y(), 2);
        ret += std::pow(z(), 2);
        return std::sqrt(ret);
    }

    void Vector::normalize()
    {
        const Scalar L = norm();
        this->operator/=(L);
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
