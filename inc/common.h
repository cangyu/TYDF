#ifndef __GT_COMMON_H__
#define __GT_COMMON_H__

#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <exception>
#include <stdexcept>

namespace GridTool
{
	namespace COMMON
	{
		double relaxation(double a, double b, double x);

		struct wrong_index : public std::logic_error
		{
			wrong_index(long long idx, const std::string &reason) : std::logic_error("\"" + std::to_string(idx) + "\" " + reason + ".") {}
			virtual ~wrong_index() = default;
		};

		class DIM
		{
		private:
			struct wrong_dimension : public wrong_index
			{
				wrong_dimension(int dim) : wrong_index(dim, "is not a valid dimension") {}
			};

		protected:
			bool m_is3D;
			int m_dim;

		public:
			DIM() = delete;
			DIM(int dim, bool is3d = true);
			DIM(const DIM &rhs) = default;
			virtual ~DIM() = default;

			bool is3D() const { return m_is3D; }
			int dimension() const { return m_dim; }
		};

		typedef double Scalar;

		class Vector : public std::array<Scalar, 3>
		{
		protected:
			struct not_vector_component : public wrong_index
			{
				not_vector_component(short x) : wrong_index(x, "is not a valid index of certain vector component") {}
			};

		public:
			Vector() : std::array<Scalar, 3>{0.0, 0.0, 0.0} {}
			Vector(Scalar val) : std::array<Scalar, 3>{val, val, val} {}
			Vector(Scalar v1, Scalar v2, Scalar v3) : std::array<Scalar, 3>{v1, v2, v3} {}
			Vector(const Vector &obj) : std::array<Scalar, 3>{obj.x(), obj.y(), obj.z()} {}
			~Vector() = default;

			// 1-based indexing
			const Scalar &operator()(short idx) const;
			Scalar &operator()(short idx);

			// Access through component
			const Scalar &x() const { return at(0); }
			const Scalar &y() const { return at(1); }
			const Scalar &z() const { return at(2); }
			Scalar &x() { return at(0); }
			Scalar &y() { return at(1); }
			Scalar &z() { return at(2); }

			// Operators
			Vector &operator=(const Vector &rhs);
			Vector &operator+=(const Vector &rhs);
			Vector &operator-=(const Vector &rhs);
			Vector &operator*=(Scalar a);
			Vector &operator/=(Scalar a);

			// Mathematical operations
			Scalar dot(const Vector &b) const;
			Vector cross(const Vector &b) const;
			Scalar norm() const;
			void normalize();
		};

		template <typename T>
		class Array1D : public std::vector<T>
		{
		public:
			Array1D(size_t n = 0) : std::vector<T>(n) {}
			Array1D(size_t n, const T &val) : std::vector<T>(n, val) {}
			Array1D(const Array1D &obj) = default;
			~Array1D() = default;

			// 1-based indexing
			T &operator()(size_t i) { return this->at(i - 1); }
			const T &operator()(size_t i) const { return this->at(i - 1); }

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

		template<typename T>
		class ArrayND
		{
		protected:
			size_t m_Nx, m_Ny, m_Nz;
			std::vector<T> m_data;
			size_t m_NXY;

		public:
			ArrayND() = delete;
			ArrayND(size_t nx, const T &val) :
				m_Nx(nx),
				m_Ny(1),
				m_Nz(1),
				m_data(nx, val),
				m_NXY(nx)
			{
				if (nI() == 0)
					throw std::invalid_argument("0 in I-dim.");
			}
			ArrayND(size_t nx, size_t ny, const T &val) :
				m_Nx(nx),
				m_Ny(ny),
				m_Nz(1),
				m_data(nx*ny, val),
				m_NXY(nx*ny)
			{
				if (nI() == 0)
					throw std::invalid_argument("0 in I-dim.");
				if (nJ() == 0)
					throw std::invalid_argument("0 in J-dim.");
			}
			ArrayND(size_t nx, size_t ny, size_t nz, const T &val) :
				m_Nx(nx),
				m_Ny(ny),
				m_Nz(nz),
				m_data(nx*ny*nz, val),
				m_NXY(nx*ny)
			{
				if (nI() == 0)
					throw std::invalid_argument("0 in I-dim.");
				if (nJ() == 0)
					throw std::invalid_argument("0 in J-dim.");
				if (nK() == 0)
					throw std::invalid_argument("0 in K-dim.");
			}
			ArrayND(const ArrayND &rhs) :
				m_Nx(rhs.m_Nx),
				m_Ny(rhs.m_Ny),
				m_Nz(rhs.m_Nz),
				m_data(rhs.m_data.begin(), rhs.m_data.end()),
				m_NXY(rhs.m_NXY)
			{
				if (nI() == 0)
					throw std::runtime_error("0 in I-dim.");
				if (nJ() == 0)
					throw std::runtime_error("0 in J-dim.");
				if (nK() == 0)
					throw std::runtime_error("0 in K-dim.");
			}
			virtual ~ArrayND() = default;

			size_t nI() const { return m_Nx; }
			size_t nJ() const { return m_Ny; }
			size_t nK() const { return m_Nz; }

		private:
			// Calculate 0-based internal index
			size_t idx(size_t i, size_t j) const { return i + m_Nx * j; }
			size_t idx(size_t i, size_t j, size_t k) const { return i + m_Nx * j + m_NXY * k; }

		public:
			/********************************* 2D *********************************/
			// 0-based indexing
			T& at(size_t i, size_t j) { return m_data[idx(i, j)]; }
			T at(size_t i, size_t j) const { return m_data[idx(i, j)]; }

			// 1-based indexing
			T& operator()(size_t i, size_t j) { return at(i - 1, j - 1); }
			T operator()(size_t i, size_t j) const { return at(i - 1, j - 1); }

			/********************************* 3D *********************************/
			// 0-based indexing
			T& at(size_t i, size_t j, size_t k) { return m_data[idx(i, j, k)]; }
			T at(size_t i, size_t j, size_t k) const { return m_data[idx(i, j, k)]; }

			// 1-based indexing
			T& operator()(size_t i, size_t j, size_t k) { return at(i - 1, j - 1, k - 1); }
			T operator()(size_t i, size_t j, size_t k) const { return at(i - 1, j - 1, k - 1); }
		};
	}
}

#endif
