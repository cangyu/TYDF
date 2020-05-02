#include "../inc/common.h"

namespace GridTool::COMMON
{
    Scalar relaxation(Scalar a, Scalar b, Scalar x)
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

    bool DIM::is3D() const
    {
        return m_is3D;
    }

    int DIM::dimension() const
    {
        return m_dim;
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
        case -1:
            return z();
        case -2:
            return y();
        case -3:
            return x();
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
        case -1:
            return z();
        case -2:
            return y();
        case -3:
            return x();
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
}
