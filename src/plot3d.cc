#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "../inc/plot3d.h"

static void formatted_writer(std::ostream &fout, double val, size_t &counter)
{
    static const size_t NumPerLine = 6;
    static const std::string sep("\t");

    fout << sep << val;
    if (++counter == NumPerLine)
    {
        fout << std::endl;
        counter = 0;
    }
}

struct invalid_dimension_size : public std::invalid_argument
{
    invalid_dimension_size(char dim, size_t n) : std::invalid_argument("Invalid " + std::to_string(dim) + " dimension: " + std::to_string(n) + ".") {}
};

namespace GridTool
{
    namespace PLOT3D
    {
        BLK::BLK(size_t nI, size_t nJ, bool is3D) : DIM(2, is3D), ArrayND<Vector>(nI, nJ, { 0.0, 0.0, 0.0 })
        {
            if (nI == 0)
                throw invalid_dimension_size('I', nI);
            if (nJ == 0)
                throw invalid_dimension_size('J', nJ);
        }

        BLK::BLK(size_t nI, size_t nJ, size_t nK) : DIM(3), ArrayND<Vector>(nI, nJ, nK, { 0.0, 0.0, 0.0 })
        {
            if (nI == 0)
                throw invalid_dimension_size('I', nI);
            if (nJ == 0)
                throw invalid_dimension_size('J', nJ);
            if (nK == 0)
                throw invalid_dimension_size('K', nK);
        }

        size_t BLK::node_num() const
        {
            if (is3D())
                return nI()*nJ()*nK();
            else
                return nI()*nJ();
        }

        size_t BLK::cell_num() const
        {
            if (is3D())
                return (nI() - 1)*(nJ() - 1)*(nK() - 1);
            else
                return (nI() - 1)*(nJ() - 1);
        }

        size_t BLK::face_num() const
        {
            size_t ret = 0;
            if (is3D())
            {
                ret += (nI() - 1) * (nJ() - 1) * nK();
                ret += (nJ() - 1) * (nK() - 1) * nI();
                ret += (nK() - 1) * (nI() - 1) * nJ();
            }
            else
            {
                ret += (nI() - 1) * nJ();
                ret += (nJ() - 1) * nI();
            }
            return ret;
        }

        size_t BLK::boundary_face_num() const
        {
            size_t ret = 0;
            if (is3D())
            {
                ret += (nI() - 1) * (nJ() - 1);
                ret += (nJ() - 1) * (nK() - 1);
                ret += (nK() - 1) * (nI() - 1);
            }
            else
            {
                ret += (nI() - 1);
                ret += (nJ() - 1);
            }
            ret *= 2;
            return ret;
        }

        size_t BLK::internal_face_num() const
        {
            return face_num() - boundary_face_num();
        }

        GRID::GRID() : DIM(3), m_blk(0) {}

        GRID::GRID(const std::string &fn) : DIM(3) { readFromFile(fn); }

        GRID::GRID(const GRID &rhs) :
            DIM(rhs.dimension(), rhs.is3D()),
            m_blk(rhs.m_blk.size(), nullptr)
        {
            for (size_t i = 0; i < numOfBlock(); ++i)
                m_blk[i] = new BLK(*rhs.m_blk[i]);
        }

        GRID::~GRID() { release_all(); }

        void GRID::readFromFile(const std::string &src)
        {
            std::string s;
            std::stringstream ss;

            // Open input grid file.
            std::ifstream fin(src);
            if (!fin)
                throw std::runtime_error("Failed to read the input grid.");

            // Read block num.
            std::getline(fin, s);
            ss << s;
            int blk_num = 0;
            ss >> blk_num;
            if (blk_num <= 0)
                throw std::invalid_argument("Invalid num of blocks.");

            // Drop previous contents.
            release_all();

            // Read dimensions of each block,and allocate new storage.
            m_blk.resize(blk_num, nullptr);
            for (int n = 0; n < blk_num; ++n)
            {
                int IMAX = 0, JMAX = 0, KMAX = 0;
                std::getline(fin, s);
                ss.clear();
                ss << s;
                ss >> IMAX >> JMAX;
                if (IMAX <= 0)
                    throw std::invalid_argument("Invalid I dimension of Block " + std::to_string(n + 1) + ".");
                if (JMAX <= 0)
                    throw std::invalid_argument("Invalid J dimension of Block " + std::to_string(n + 1) + ".");

                if (!(ss >> KMAX))
                    m_blk[n] = new BLK((size_t)IMAX, (size_t)JMAX, false);
                else if (KMAX == 1)
                    m_blk[n] = new BLK((size_t)IMAX, (size_t)JMAX, true);
                else
                {
                    if (JMAX <= 0)
                        throw std::invalid_argument("Invalid K dimension of Block " + std::to_string(n + 1) + ".");
                    else
                        m_blk[n] = new BLK((size_t)IMAX, (size_t)JMAX, (size_t)KMAX);
                }
            }

            // Read coordinates of each block.
            for (auto b : m_blk)
            {
                const size_t NX = b->nI(), NY = b->nJ(), NZ = b->nK();
                if (b->dimension() == 3)
                {
                    for (size_t k = 0; k < NZ; ++k)
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                fin >> b->at(i, j, k).x();

                    for (size_t k = 0; k < NZ; ++k)
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                fin >> b->at(i, j, k).y();

                    for (size_t k = 0; k < NZ; ++k)
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                fin >> b->at(i, j, k).z();
                }
                else
                {
                    for (size_t j = 0; j < NY; ++j)
                        for (size_t i = 0; i < NX; ++i)
                            fin >> b->at(i, j).x();

                    for (size_t j = 0; j < NY; ++j)
                        for (size_t i = 0; i < NX; ++i)
                            fin >> b->at(i, j).y();

                    if (b->is3D())
                    {
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                fin >> b->at(i, j).z();
                    }
                }
            }

            // Close file.
            fin.close();

            // Update grid global DIM attributes, and check dimension consistency.
            m_is3D = m_blk[0]->is3D();
            m_dim = m_blk[0]->dimension();
            for (size_t n = 1; n < m_blk.size(); ++n)
            {
                auto blk = m_blk[n];
                if (blk->is3D() != m_is3D || blk->dimension() != m_dim)
                    throw std::runtime_error("Inconsistent DIM properties of Block " + std::to_string(n + 1) + ".");
            }
        }

        void GRID::writeToFile(const std::string &dst) const
        {
            // Open output file.
            std::ofstream fout(dst);
            if (!fout)
                throw std::runtime_error("Failed to open the target output grid file.");

            // Write num of blocks.
            fout << "\t" << numOfBlock() << std::endl;

            // Write dimensions of each block.
            for (auto b : m_blk)
            {
                fout << "\t" << b->nI();
                fout << "\t" << b->nJ();
                if (b->is3D())
                    fout << "\t" << b->nK();
                fout << std::endl;
            }

            // Write coordinates of each block
            size_t cnt = 0;
            for (auto b : m_blk)
            {
                const size_t NX = b->nI(), NY = b->nJ(), NZ = b->nK();
                if (b->dimension() == 3)
                {
                    for (size_t k = 0; k < NZ; ++k)
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                formatted_writer(fout, b->at(i, j, k).x(), cnt);

                    for (size_t k = 0; k < NZ; ++k)
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                formatted_writer(fout, b->at(i, j, k).y(), cnt);

                    for (size_t k = 0; k < NZ; ++k)
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                formatted_writer(fout, b->at(i, j, k).z(), cnt);
                }
                else
                {
                    for (size_t j = 0; j < NY; ++j)
                        for (size_t i = 0; i < NX; ++i)
                            formatted_writer(fout, b->at(i, j).x(), cnt);

                    for (size_t j = 0; j < NY; ++j)
                        for (size_t i = 0; i < NX; ++i)
                            formatted_writer(fout, b->at(i, j).y(), cnt);

                    if (b->is3D())
                    {
                        for (size_t j = 0; j < NY; ++j)
                            for (size_t i = 0; i < NX; ++i)
                                formatted_writer(fout, b->at(i, j).z(), cnt);
                    }
                }
            }

            // Close file.
            fout.close();
        }

        void GRID::release_all()
        {
            for (auto e : m_blk)
                if (e)
                    delete e;

            m_blk.clear();
        }
    }
}
