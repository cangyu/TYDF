#include <fstream>
#include <sstream>
#include "plot3d.h"

int PLOT3D::readFromFile(const std::string &src)
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

int PLOT3D::writeToFile(const std::string &dst) const
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
