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
		using COMMON::ArrayND;
		using COMMON::DIM;

		class BLK : public DIM
		{
		private:
			size_t m_nI, m_nJ, m_nK;
			ArrayND<double> *m_x, *m_y, *m_z;

		public:
			BLK(size_t nI, size_t nJ, bool is3D) :
				DIM(2, is3D),
				m_nI(nI),
				m_nJ(nJ),
				m_nK(1)
			{
				if (nI == 0 || nJ == 0)
					throw std::runtime_error("Invalid size.");

				m_x = new ArrayND<double>(nI, nJ, 0.0);
				m_y = new ArrayND<double>(nI, nJ, 0.0);
				m_z = is3D ? new ArrayND<double>(nI, nJ, 0.0) : nullptr; // May be shell mesh
			}

			BLK(size_t nI, size_t nJ, size_t nK) :
				DIM(3, true),
				m_nI(nI),
				m_nJ(nJ),
				m_nK(nK)
			{
				if (nI == 0 || nJ == 0 || nK == 0)
					throw std::runtime_error("Invalid size.");

				m_x = new ArrayND<double>(nI, nJ, nK, 0.0);
				m_y = new ArrayND<double>(nI, nJ, nK, 0.0);
				m_z = new ArrayND<double>(nI, nJ, nK, 0.0);
			}

			~BLK()
			{
				if (m_x)
					delete m_x;
				if (m_y)
					delete m_y;
				if (m_z)
					delete m_z;
			}

			size_t nI() const { return m_nI; }
			size_t nJ() const { return m_nJ; }
			size_t nK() const { return m_nK; }

			size_t node_num() const
			{
				if (is3D())
					return nI()*nJ()*nK();
				else
					return nI()*nJ();
			}

			size_t cell_num() const
			{
				if (is3D())
					return (nI() - 1)*(nJ() - 1)*(nK() - 1);
				else
					return (nI() - 1)*(nJ() - 1);
			}

			size_t face_num() const
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

			size_t boundary_face_num() const
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

			size_t internal_face_num() const { return face_num() - boundary_face_num(); }

			const ArrayND<double> &x() const { return *m_x; }
			const ArrayND<double> &y() const { return *m_y; }
			const ArrayND<double> &z() const { return *m_z; }

			ArrayND<double> &x() { return *m_x; }
			ArrayND<double> &y() { return *m_y; }
			ArrayND<double> &z() { return *m_z; }
		};

		class GRID : public DIM
		{
		private:
			std::vector<BLK *> m_blk;

			void release_all()
			{
				for (auto e : m_blk)
					if (e)
						delete e;
				m_blk.clear();
			}

		public:
			GRID() : DIM(3, true), m_blk(0) {}

			GRID(const std::string &fn) :
				DIM(3, true),
				m_blk(0)
			{
				readFromFile(fn);
			}

			~GRID()
			{
				release_all();
			}

			size_t numOfBlock() const { return m_blk.size(); }

			void readFromFile(const std::string &src);

			void writeToFile(const std::string &dst) const;

			// 0-based indexing
			BLK *block(size_t loc_idx) { return m_blk[loc_idx]; }
		};
	}
}
#endif
