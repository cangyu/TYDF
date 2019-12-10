#include "../inc/plot3d.h"

namespace GridTool
{
	namespace PLOT3D
	{
		void GRID::readFromFile(const std::string &src)
		{
			std::string s;
			std::stringstream ss;

			// Open input grid file
			std::ifstream fin(src);
			if (!fin)
				throw std::runtime_error("Failed to read the input grid.");

			// Read block num
			std::getline(fin, s);
			ss << s;
			int blk_num = 0;
			ss >> blk_num;
			if (blk_num <= 0)
				throw std::invalid_argument("Invalid num of blocks.");

			// Read dimensions of each block
			for (int n = 0; n < blk_num; ++n)
			{
				size_t IMAX = 0, JMAX = 0, KMAX = 0;
				std::getline(fin, s);
				ss.clear();
				ss << s;
				ss >> IMAX >> JMAX;
				if (!(ss >> KMAX))
					m_blk.push_back(new BLK(IMAX, JMAX, false));
				else if (KMAX == 1)
					m_blk.push_back(new BLK(IMAX, JMAX, true));
				else
					m_blk.push_back(new BLK(IMAX, JMAX, KMAX));
			}

			// Read coordinates of each block
			for (auto blk : m_blk)
			{
				const size_t NX = blk->nI(), NY = blk->nJ(), NZ = blk->nK();
				if (blk->dimension() == 3)
				{
					auto &x = blk->x();
					for (size_t k = 0; k < NZ; ++k)
						for (size_t j = 0; j < NY; ++j)
							for (size_t i = 0; i < NX; ++i)
								fin >> x.at(i, j, k);

					auto &y = blk->y();
					for (size_t k = 0; k < NZ; ++k)
						for (size_t j = 0; j < NY; ++j)
							for (size_t i = 0; i < NX; ++i)
								fin >> y.at(i, j, k);

					auto &z = blk->z();
					for (size_t k = 0; k < NZ; ++k)
						for (size_t j = 0; j < NY; ++j)
							for (size_t i = 0; i < NX; ++i)
								fin >> z.at(i, j, k);
				}
				else
				{
					auto &x = blk->x();
					for (size_t j = 0; j < NY; ++j)
						for (size_t i = 0; i < NX; ++i)
							fin >> x.at(i, j);

					auto &y = blk->y();
					for (size_t j = 0; j < NY; ++j)
						for (size_t i = 0; i < NX; ++i)
							fin >> y.at(i, j);

					if (blk->is3D())
					{
						auto &z = blk->z();
						for (size_t j = 0; j < NY; ++j)
							for (size_t i = 0; i < NX; ++i)
								fin >> z.at(i, j);
					}
				}
			}

			// Close file
			fin.close();

			// Check dimension consistency
			m_is3D = m_blk[0]->is3D();
			m_dim = m_blk[0]->dimension();
			for (size_t n = 1; n < m_blk.size(); ++n)
			{
				auto blk = m_blk[n];
				if (blk->is3D() != m_is3D || blk->dimension() != m_dim)
					throw std::runtime_error("Inconsistent dimension of blocks");
			}
		}

		void GRID::writeToFile(const std::string &dst) const
		{
			// Open output file
			std::ofstream fout(dst);
			if (!fout)
				throw std::runtime_error("Failed to open the target output grid file.");

			// Write num of blocks
			fout << "\t" << numOfBlock() << std::endl;

			// Write dimensions of each block
			for (auto blk : m_blk)
			{
				fout << "\t" << blk->nI();
				fout << "\t" << blk->nJ();
				if (blk->is3D())
					fout << "\t" << blk->nK();
				fout << std::endl;
			}

			// Write coordinates of each block
			for (auto blk : m_blk)
			{
				const size_t NX = blk->nI(), NY = blk->nJ(), NZ = blk->nK();
				if (blk->dimension() == 3)
				{
					const auto &x = blk->x();
					for (size_t k = 0; k < NZ; ++k)
						for (size_t j = 0; j < NY; ++j)
						{
							for (size_t i = 0; i < NX; ++i)
								fout << "\t" << x.at(i, j, k);
							fout << std::endl;
						}

					const auto &y = blk->y();
					for (size_t k = 0; k < NZ; ++k)
						for (size_t j = 0; j < NY; ++j)
						{
							for (size_t i = 0; i < NX; ++i)
								fout << "\t" << y.at(i, j, k);
							fout << std::endl;
						}

					const auto &z = blk->z();
					for (size_t k = 0; k < NZ; ++k)
						for (size_t j = 0; j < NY; ++j)
						{
							for (size_t i = 0; i < NX; ++i)
								fout << "\t" << z.at(i, j, k);
							fout << std::endl;
						}
				}
				else
				{
					const auto &x = blk->x();
					for (size_t j = 0; j < NY; ++j)
					{
						for (size_t i = 0; i < NX; ++i)
							fout << "\t" << x.at(i, j);
						fout << std::endl;
					}

					const auto &y = blk->y();
					for (size_t j = 0; j < NY; ++j)
					{
						for (size_t i = 0; i < NX; ++i)
							fout << "\t" << y.at(i, j);
						fout << std::endl;
					}

					if (blk->is3D())
					{
						const auto &z = blk->z();
						for (size_t j = 0; j < NY; ++j)
						{
							for (size_t i = 0; i < NX; ++i)
								fout << "\t" << z.at(i, j);
							fout << std::endl;
						}
					}
				}
			}

			// Close file
			fout.close();
		}
	}
}