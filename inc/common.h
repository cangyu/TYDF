#ifndef __GT_COMMON_H__
#define __GT_COMMON_H__

#include <cmath>
#include <array>
#include <string>
#include <exception>

namespace GridTool
{
	namespace COMMON
	{
		class wrong_index : public std::exception
		{
		protected:
			std::string m_msg;

		public:
			wrong_index(long long idx) : std::exception(), m_msg("\"" + std::to_string(idx) + "\" ") {}
			virtual ~wrong_index() = default;
		};

		typedef double Scalar;

		class Vector : public std::array<Scalar, 3>
		{
		protected:
			class not_vector_component : public wrong_index
			{
			public:
				not_vector_component(short x) : wrong_index(x)
				{
					m_msg += "is not a valid index of certain vector component.";
				}
			};

		public:
			Vector() : std::array<Scalar, 3>{0.0, 0.0, 0.0} {}
			Vector(Scalar val) : std::array<Scalar, 3>{val, val, val} {}
			Vector(Scalar v1, Scalar v2, Scalar v3) : std::array<Scalar, 3>{v1, v2, v3} {}
			Vector(const Vector &obj) : std::array<Scalar, 3>{obj.x(), obj.y(), obj.z()} {}
			~Vector() = default;

			// 1-based indexing
			const Scalar &operator()(short idx) const
			{
				switch (idx)
				{
				case 1:
					return x();
				case 2:
					return y();
				case 3:
					return z();
				default:
					throw not_vector_component(idx);
				}
			}
			Scalar &operator()(short idx)
			{
				switch (idx)
				{
				case 1:
					return x();
				case 2:
					return y();
				case 3:
					return z();
				default:
					throw not_vector_component(idx);
				}
			}

			// Access through component
			const Scalar &x() const { return at(0); }
			Scalar &x() { return at(0); }

			const Scalar &y() const { return at(1); }
			Scalar &y() { return at(1); }

			const Scalar &z() const { return at(2); }
			Scalar &z() { return at(2); }

			// Operators
			Vector &operator=(const Vector &rhs)
			{
				x() = rhs.x();
				y() = rhs.y();
				z() = rhs.z();
				return *this;
			}
			Vector &operator+=(const Vector &rhs)
			{
				x() += rhs.x();
				y() += rhs.y();
				z() += rhs.z();
				return *this;
			}
			Vector &operator-=(const Vector &rhs)
			{
				x() -= rhs.x();
				y() -= rhs.y();
				z() -= rhs.z();
				return *this;
			}
			Vector &operator*=(Scalar a)
			{
				x() *= a;
				y() *= a;
				z() *= a;
				return *this;
			}
			Vector &operator/=(Scalar a)
			{
				x() /= a;
				y() /= a;
				z() /= a;
				return *this;
			}

			// Mathematical operations
			Scalar dot(const Vector &b) const
			{
				Scalar ret = 0.0;
				ret += x() * b.x();
				ret += y() * b.y();
				ret += z() * b.z();
				return ret;
			}

			Vector cross(const Vector &b) const
			{
				Vector ret;
				ret.x() = y() * b.z() - z() * b.y();
				ret.y() = z() * b.x() - x() * b.z();
				ret.z() = x() * b.y() - y() * b.x();
				return ret;
			}

			Scalar norm() const
			{
				Scalar ret = 0.0;
				ret += std::pow(x(), 2);
				ret += std::pow(y(), 2);
				ret += std::pow(z(), 2);
				return std::sqrt(ret);
			}

			void normalize()
			{
				const Scalar L = norm();
				this->operator/=(L);
			}
		};

		template <typename T>
		class Array1D : public std::vector<T>
		{
		protected:
			class not_1_based : public wrong_index
			{
			public:
				not_1_based(size_t x) : wrong_index(x)
				{
					m_msg += "is not a valid 1-based index.";
				}
			};

		public:
			Array1D(size_t n = 0) : std::vector<T>(n) {}
			Array1D(size_t n, const T &val) : std::vector<T>(n, val) {}
			Array1D(const Array1D &obj) = default;
			~Array1D() = default;

			// 1-based indexing
			T &operator()(size_t i)
			{
				if (1 <= i && i <= this->size())
					return this->at(i - 1);
				else
					throw not_1_based(i);
			}
			const T &operator()(size_t i) const
			{
				if (1 <= i && i <= this->size())
					return this->at(i - 1);
				else
					throw not_1_based(i);
			}

			// Check includances
			bool contains(const T &x) const
			{
				const size_t N = this->size();
				for (size_t i = 0; i < N; ++i)
					if (x == this->at(i))
						return true;

				return false;
			}
			bool contains(const T &a, const T &b) const
			{
				bool flag_a = false, flag_b = false;
				const size_t N = this->size();
				for (size_t i = 0; i < N; ++i)
				{
					const T &x = this->at(i);
					if (!flag_a && a == x)
						flag_a = true;
					if (!flag_b && b == x)
						flag_b = true;

					if (flag_a && flag_b)
						return true;
				}
				return false;
			}
		};
	}
}

#endif