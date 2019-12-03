#ifndef __COMMON_HPP__
#define __COMMON_HPP__

namespace COMMON
{
	typedef double Scalar;

	class Vector : public std::array<Scalar, 3>
	{
	public:
		Vector() : std::array<Scalar, 3>{ 0.0, 0.0, 0.0 } {}
		Vector(Scalar val) : std::array<Scalar, 3>{ val, val, val } {}
		Vector(Scalar v1, Scalar v2, Scalar v3) : std::array<Scalar, 3>{ v1, v2, v3 } {}
		Vector(const Vector &obj) : std::array<Scalar, 3>{ obj.x(), obj.y(), obj.z() } {}
		~Vector() = default;

		// Access through component
		Scalar x() const { return at(0); }
		Scalar y() const { return at(1); }
		Scalar z() const { return at(2); }

		Scalar &x() { return at(0); }
		Scalar &y() { return at(1); }
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

		// 1-based indexing
		Scalar operator()(int idx) const
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
				throw std::invalid_argument("\"" + std::to_string(idx) + "\" is not a valid 1-based index.");
			}
		}
		Scalar &operator()(int idx)
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
				throw std::invalid_argument("\"" + std::to_string(idx) + "\" is not a valid 1-based index.");
			}
		}
	};

	class wrong_index : public std::exception
	{
	protected:
		std::string m_msg;

	public:
		wrong_index(long long idx) : std::exception(), m_msg("\"" + std::to_string(idx) + "\" ") {}
		virtual ~wrong_index() = default;

		char const* what() const
		{
			return m_msg.c_str();
		}
	};

	template<typename T>
	class Array1D : public std::vector<T>
	{
	protected:
		class not_1_based : public wrong_index
		{
		public:
			not_1_based(size_t x) :
				wrong_index(x)
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

#endif
