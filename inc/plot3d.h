#ifndef __GT_PLOT3D_H__
#define __GT_PLOT3D_H__

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include "common.h"

namespace GridTool
{
	namespace PLOT3D
	{
		using COMMON::Vector;
		using COMMON::DIM;
		using COMMON::ArrayND;

		class BLK : public DIM, public ArrayND<Vector>
		{
		public:
			BLK(size_t nI, size_t nJ, bool is3D) :
				DIM(2, is3D),
				ArrayND<Vector>(nI, nJ, { 0.0, 0.0, 0.0 })
			{
				if (nI == 0)
					throw std::invalid_argument("Invalid I dimension.");
				if (nJ == 0)
					throw std::invalid_argument("Invalid J dimension.");
			}
			BLK(size_t nI, size_t nJ, size_t nK) :
				DIM(3),
				ArrayND<Vector>(nI, nJ, nK, { 0.0, 0.0, 0.0 })
			{
				if (nI == 0)
					throw std::invalid_argument("Invalid I dimension.");
				if (nJ == 0)
					throw std::invalid_argument("Invalid J dimension.");
				if (nK == 0)
					throw std::invalid_argument("Invalid K dimension.");
			}
			BLK(const BLK &rhs) = default;
			~BLK() = default;

			size_t node_num() const;
			size_t cell_num() const;
			size_t face_num() const;

			size_t boundary_face_num() const;
			size_t internal_face_num() const;
		};

		class GRID : public DIM
		{
		private:
			std::vector<BLK *> m_blk;

		public:
			GRID() : DIM(3), m_blk(0) {}
			GRID(const std::string &fn) : DIM(3)
			{
				readFromFile(fn);
			}
			GRID(const GRID &rhs) :
				DIM(rhs.dimension(), rhs.is3D()),
				m_blk(rhs.m_blk.size(), nullptr)
			{
				for (size_t i = 0; i < numOfBlock(); ++i)
					m_blk[i] = new BLK(*rhs.m_blk[i]);
			}
			~GRID()
			{
				release_all();
			}

			size_t numOfBlock() const { return m_blk.size(); }

			// IO
			void readFromFile(const std::string &src);
			void writeToFile(const std::string &dst) const;

			// 0-based indexing
			BLK *block(size_t loc_idx) { return m_blk[loc_idx]; }

		private:
			void release_all()
			{
				for (auto e : m_blk)
					if (e)
						delete e;
				m_blk.clear();
			}
		};
	}
}
#endif
