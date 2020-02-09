#ifndef __GT_PLOT3D_H__
#define __GT_PLOT3D_H__

#include <vector>
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
            BLK(size_t nI, size_t nJ, bool is3D);
            BLK(size_t nI, size_t nJ, size_t nK);
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
            GRID();
            GRID(const std::string &fn);
            GRID(const GRID &rhs);
            ~GRID();

            size_t numOfBlock() const { return m_blk.size(); }

            // IO
            void readFromFile(const std::string &src);
            void writeToFile(const std::string &dst) const;

            // 0-based indexing
            BLK *block(size_t loc_idx) { return m_blk[loc_idx]; }

        private:
            void release_all();
        };
    }
}
#endif
