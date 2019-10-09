#ifndef __PLOT3D_HPP__
#define __PLOT3D_HPP__

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

class PLOT3D_BLK
{
private:
	template<typename T>
	class Array
	{
	private:
		size_t m_Nx, m_Ny, m_Nz;
		std::vector<T> m_data;
		size_t m_NXY;

	public:
		Array(size_t nx, size_t ny, const T &val) :
			m_Nx(nx),
			m_Ny(ny),
			m_Nz(1),
			m_data(nx*ny, val)
		{
			if (nx == 0 || ny == 0)
			{
				std::string msg("Invalid size: nx=" + std::to_string(nx) + ", ny=" + std::to_string(ny));
				throw std::runtime_error(msg);
			}

			m_NXY = nx * ny;
		}

		Array(size_t nx, size_t ny, size_t nz, const T &val) :
			m_Nx(nx),
			m_Ny(ny),
			m_Nz(nz),
			m_data(nx*ny*nz, val)
		{
			if (nx == 0 || ny == 0 || nz == 0)
			{
				std::string msg("Invalid size: nx=" + std::to_string(nx) + ", ny=" + std::to_string(ny) + ", nz=" + std::to_string(nz));
				throw std::runtime_error(msg);
			}

			m_NXY = nx * ny;
		}

		~Array() = default;

		size_t nI() const
		{
			return m_Nx;
		}

		size_t nJ() const
		{
			return m_Ny;
		}

		size_t nK() const
		{
			return m_Nz;
		}

	private:
		size_t idx(size_t i, size_t j) const
		{
			return i + m_Nx * j; // 0-based
		}

		size_t idx(size_t i, size_t j, size_t k) const
		{
			return i + m_Nx * j + m_NXY * k; // 0-based
		}

	public:
		/********************************* 2D *********************************/
		// 0-based indexing
		T& at(size_t i, size_t j)
		{
			return m_data[idx(i, j)];
		}

		T at(size_t i, size_t j) const
		{
			return m_data[idx(i, j)];
		}

		// 1-based indexing
		T& operator()(size_t i, size_t j)
		{
			return at(i - 1, j - 1);
		}

		T operator()(size_t i, size_t j) const
		{
			return at(i - 1, j - 1);
		}

		/********************************* 3D *********************************/
		// 0-based indexing
		T& at(size_t i, size_t j, size_t k)
		{
			return m_data[idx(i, j, k)];
		}

		T at(size_t i, size_t j, size_t k) const
		{
			return m_data[idx(i, j, k)];
		}

		// 1-based indexing
		T& operator()(size_t i, size_t j, size_t k)
		{
			return at(i - 1, j - 1, k - 1);
		}

		T operator()(size_t i, size_t j, size_t k) const
		{
			return at(i - 1, j - 1, k - 1);
		}
	};

public:
	typedef Array<double> *pCoordBlk;

	PLOT3D_BLK(size_t nI, size_t nJ) :
		m_nI(nI),
		m_nJ(nJ),
		m_nK(1),
		m_3d(false),
		m_dim(2)
	{
		if (nI == 0 || nJ == 0)
			throw std::runtime_error("Invalid size.");

		m_nIJK = m_nIJ = m_nI * m_nJ;

		m_x = new Array<double>(nI, nJ, 0.0);
		m_y = new Array<double>(nI, nJ, 0.0);
		m_z = nullptr;
	}

	PLOT3D_BLK(size_t nI, size_t nJ, size_t nK) :
		m_nI(nI),
		m_nJ(nJ),
		m_nK(nK),
		m_3d(true),
		m_dim(3)
	{
		if (nI == 0 || nJ == 0 || nK == 0)
			throw std::runtime_error("Invalid size.");

		m_nIJ = m_nI * m_nJ;
		m_nIJK = m_nI * m_nJ * m_nK;

		m_x = new Array<double>(nI, nJ, nK, 0.0);
		m_y = new Array<double>(nI, nJ, nK, 0.0);
		m_z = new Array<double>(nI, nJ, nK, 0.0);
	}

	~PLOT3D_BLK()
	{
		if (m_x)
			delete m_x;
		if (m_y)
			delete m_y;
		if (m_z)
			delete m_z;
	}

	bool is3D() const
	{
		return m_3d;
	}

	int dimension() const
	{
		return m_dim;
	}

	size_t node_num() const
	{
		return m_nIJK;
	}

	size_t cell_num() const
	{
		size_t ret = 0;

		if (m_3d)
			ret = (m_nI - 1)*(m_nJ - 1)*(m_nK - 1);
		else
			ret = (m_nI - 1)*(m_nJ - 1);

		return ret;
	}

	size_t face_num() const
	{
		size_t ret = 0;

		if (m_3d)
		{
			ret += (m_nI - 1) * (m_nJ - 1) * m_nK;
			ret += (m_nJ - 1) * (m_nK - 1) * m_nI;
			ret += (m_nK - 1) * (m_nI - 1) * m_nJ;
		}
		else
		{
			ret += (m_nI - 1) * m_nJ;
			ret += (m_nJ - 1) * m_nI;
		}

		return ret;
	}

	size_t boundary_face_num() const
	{
		size_t ret = 0;

		if (m_3d)
		{
			ret += (m_nI - 1) * (m_nJ - 1);
			ret += (m_nJ - 1) * (m_nK - 1);
			ret += (m_nK - 1) * (m_nI - 1);
			ret *= 2;
		}
		else
		{
			ret += (m_nI - 1);
			ret += (m_nJ - 1);
			ret *= 2;
		}

		return ret;
	}

	size_t internal_face_num() const
	{
		return face_num() - boundary_face_num();
	}

	size_t nI() const
	{
		return m_nI;
	}

	size_t nJ() const
	{
		return m_nJ;
	}

	size_t nK() const
	{
		return m_nK;
	}

	pCoordBlk x()
	{
		return m_x;
	}

	pCoordBlk y()
	{
		return m_y;
	}

	pCoordBlk z()
	{
		return m_z;
	}

private:
	size_t m_nI, m_nJ, m_nK;
	bool m_3d;
	int m_dim;
	pCoordBlk m_x, m_y, m_z;
	size_t m_nIJ, m_nIJK;
};

typedef PLOT3D_BLK *pPLOT3D_BLK;

class PLOT3D
{
private:
	size_t m_nBLK;
	bool m_3d;
	int m_dim;
	std::vector<pPLOT3D_BLK> m_blk;

public:
	PLOT3D()
	{
		m_nBLK = 0;
		m_3d = true;
		m_dim = 3;
		m_blk.resize(0, nullptr);
	}

	~PLOT3D()
	{
		for (size_t i = 0; i < m_blk.size(); ++i)
			if (m_blk[i])
				delete m_blk[i];
	}

	bool is3D() const
	{
		return m_3d;
	}

	int dimension() const
	{
		return m_dim;
	}

	size_t nBLK() const
	{
		return m_nBLK;
	}

	// Access block through 0-based indexing
	pPLOT3D_BLK at(size_t idx)
	{
		return m_blk[idx];
	}

	// Access block through 1-based indexing
	pPLOT3D_BLK operator()(size_t idx)
	{
		return at(idx - 1);
	}

	int readFromFile(const std::string &src)
	{
		std::string s;
		std::stringstream ss;

		// Open input grid file
		std::ifstream fin(src);
		if (!fin)
			return -1;

		// Read block num
		std::getline(fin, s);
		ss << s;
		ss >> m_nBLK;

		// Read dimensions of each block
		for (size_t n = 0; n < m_nBLK; ++n)
		{
			size_t IMAX, JMAX, KMAX;

			std::getline(fin, s);
			ss.clear();
			ss << s;
			ss >> IMAX >> JMAX;
			if (!(ss >> KMAX))
				m_blk.push_back(new PLOT3D_BLK(IMAX, JMAX));
			else
				m_blk.push_back(new PLOT3D_BLK(IMAX, JMAX, KMAX));
		}

		// Read coordinates of each block
		for (size_t n = 0; n < m_nBLK; ++n)
		{
			auto blk = m_blk[n];
			const size_t NX = blk->nI(), NY = blk->nJ(), NZ = blk->nK();

			if (blk->is3D())
			{
				auto x = blk->x();
				for (size_t k = 0; k < NZ; ++k)
					for (size_t j = 0; j < NY; ++j)
						for (size_t i = 0; i < NX; ++i)
							fin >> x->at(i, j, k);

				auto y = blk->y();
				for (size_t k = 0; k < NZ; ++k)
					for (size_t j = 0; j < NY; ++j)
						for (size_t i = 0; i < NX; ++i)
							fin >> y->at(i, j, k);

				auto z = blk->z();
				for (size_t k = 0; k < NZ; ++k)
					for (size_t j = 0; j < NY; ++j)
						for (size_t i = 0; i < NX; ++i)
							fin >> z->at(i, j, k);
			}
			else
			{
				auto x = blk->x();
				for (size_t j = 0; j < NY; ++j)
					for (size_t i = 0; i < NX; ++i)
						fin >> x->at(i, j);

				auto y = blk->y();
				for (size_t j = 0; j < NY; ++j)
					for (size_t i = 0; i < NX; ++i)
						fin >> y->at(i, j);
			}
		}

		// Check dimension consistency
		m_3d = m_blk[0]->is3D();
		m_dim = m_blk[0]->dimension();
		for (size_t n = 1; n < m_nBLK; ++n)
		{
			auto blk = m_blk[n];
			if (blk->is3D() != m_3d || blk->dimension() != m_dim)
				return -2;
		}

		// Close file
		fin.close();
		return 0;
	}

	int writeToFile(const std::string &dst) const
	{
		// Open output file
		std::ofstream fout(dst);
		if (!fout)
			return -1;

		// Write num of blocks
		fout << "\t" << m_nBLK << std::endl;

		// Write dimensions of each block
		for (size_t n = 0; n < m_nBLK; ++n)
		{
			auto blk = m_blk[n];
			fout << "\t" << blk->nI();
			fout << "\t" << blk->nJ();
			if (blk->is3D())
				fout << "\t" << blk->nK();
			fout << std::endl;
		}

		// Write coordinates of each block
		for (size_t n = 0; n < m_nBLK; ++n)
		{
			auto blk = m_blk[n];
			const size_t NX = blk->nI(), NY = blk->nJ(), NZ = blk->nK();

			if (blk->is3D())
			{
				auto x = blk->x();
				for (size_t k = 0; k < NZ; ++k)
					for (size_t j = 0; j < NY; ++j)
					{
						for (size_t i = 0; i < NX; ++i)
							fout << "\t" << x->at(i, j, k);
						fout << std::endl;
					}

				auto y = blk->y();
				for (size_t k = 0; k < NZ; ++k)
					for (size_t j = 0; j < NY; ++j)
					{
						for (size_t i = 0; i < NX; ++i)
							fout << "\t" << y->at(i, j, k);
						fout << std::endl;
					}

				auto z = blk->z();
				for (size_t k = 0; k < NZ; ++k)
					for (size_t j = 0; j < NY; ++j)
					{
						for (size_t i = 0; i < NX; ++i)
							fout << "\t" << z->at(i, j, k);
						fout << std::endl;
					}
			}
			else
			{
				auto x = blk->x();
				for (size_t j = 0; j < NY; ++j)
				{
					for (size_t i = 0; i < NX; ++i)
						fout << "\t" << x->at(i, j);
					fout << std::endl;
				}

				auto y = blk->y();
				for (size_t j = 0; j < NY; ++j)
				{
					for (size_t i = 0; i < NX; ++i)
						fout << "\t" << y->at(i, j);
					fout << std::endl;
				}
			}
		}

		// Close file
		fout.close();
		return 0;
	}
};

#endif
