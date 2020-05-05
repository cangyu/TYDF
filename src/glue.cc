#include "../inc/nmf.h"
#include "../inc/plot3d.h"
#include "../inc/xf.h"

static std::string version_str()
{
    static const size_t major = 2;
    static const size_t minor = 0;
    static const size_t patch = 0;

    return "V" + std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(patch);
}

namespace GridTool::XF
{
    MESH::MESH(const std::string &f_nmf, const std::string &f_p3d, std::ostream &fout) :
        DIM(3),
        m_totalNodeNum(0),
        m_totalCellNum(0),
        m_totalFaceNum(0),
        m_totalZoneNum(0)
    {
        /// Load mapping file.
        auto nmf = new NMF::Mapping3D(f_nmf);
        nmf->compute_topology();
        nmf->numbering();

        /// Load grid file.
        auto p3d = new PLOT3D::GRID(f_p3d);

        /// Check consistency.
        const size_t NBLK = nmf->nBlock();
        if (NBLK != p3d->numOfBlock())
            throw std::invalid_argument("Inconsistent num of blocks between NMF and PLOT3D.");
        for (size_t n = 1; n <= NBLK; ++n)
        {
            const auto &b = nmf->block(n);
            auto g = p3d->block(n - 1);
            if (b.IDIM() != g->nI())
                throw std::invalid_argument("Inconsistent num of nodes in I dimension of Block " + std::to_string(n) + ".");
            if (b.JDIM() != g->nJ())
                throw std::invalid_argument("Inconsistent num of nodes in J dimension of Block " + std::to_string(n) + ".");
            if (b.KDIM() != g->nK())
                throw std::invalid_argument("Inconsistent num of nodes in K dimension of Block " + std::to_string(n) + ".");
        }

        /// Allocate storage.
        m_totalNodeNum = nmf->nNode();
        m_totalCellNum = nmf->nCell();
        size_t innerFaceNum = 0, bdryFaceNum = 0;
        nmf->nFace(m_totalFaceNum, innerFaceNum, bdryFaceNum);
        m_node.resize(numOfNode());
        m_face.resize(numOfFace());
        m_cell.resize(numOfCell());

        /// Copy node info.
        std::vector<bool> visited(m_node.size(), false);
        for (size_t n = 1; n <= NBLK; ++n)
        {
            auto &b = nmf->block(n);
            auto &g = *p3d->block(n - 1);

            const size_t nI = b.IDIM();
            const size_t nJ = b.JDIM();
            const size_t nK = b.KDIM();

            for (size_t k = 1; k <= nK; ++k)
                for (size_t j = 1; j <= nJ; ++j)
                    for (size_t i = 1; i <= nI; ++i)
                    {
                        /// Global 1-based index, already assigned.
                        const auto idx = b.node_index(i, j, k);

                        if (!visited[idx - 1])
                        {
                            node(idx).coordinate = g(i, j, k);
                            visited[idx - 1] = true;
                        }
                    }
        }

        /// Copy cell info.
        for (size_t n = 1; n <= NBLK; ++n)
        {
            auto &b = nmf->block(n);

            const size_t nI = b.IDIM();
            const size_t nJ = b.JDIM();
            const size_t nK = b.KDIM();

            for (size_t k = 1; k < nK; ++k)
                for (size_t j = 1; j < nJ; ++j)
                    for (size_t i = 1; i < nI; ++i)
                    {
                        const auto &nc = b.cell(i, j, k);
                        const auto idx = nc.CellSeq();
                        auto &fc = cell(idx);

                        fc.type = CELL::HEXAHEDRAL;

                        fc.includedFace.resize(NMF::Block3D::NumOfSurf);
                        for (short r = 1; r <= NMF::Block3D::NumOfSurf; ++r)
                            fc.includedFace(r) = nc.FaceSeq(r);

                        fc.includedNode.resize(NMF::Block3D::NumOfVertex);
                        for (short r = 1; r <= NMF::Block3D::NumOfVertex; ++r)
                            fc.includedNode(r) = nc.NodeSeq(r);
                    }
        }

        /// Copy face info.
        visited.resize(m_face.size(), false);
        std::fill(visited.begin(), visited.end(), false);
        for (size_t n = 1; n <= NBLK; ++n)
        {
            auto &b = nmf->block(n);

            const size_t nI = b.IDIM();
            const size_t nJ = b.JDIM();
            const size_t nK = b.KDIM();

            /// Internal I direction
            for (size_t k = 1; k < nK; ++k)
                for (size_t j = 1; j < nJ; ++j)
                    for (size_t i = 2; i < nI; ++i)
                    {
                        const auto &curCell = b.cell(i, j, k);
                        const auto &adjCell = b.cell(i - 1, j, k);

                        const auto faceIndex = curCell.FaceSeq(1);
                        auto &curFace = face(faceIndex);

                        curFace.atBdry = false;

                        curFace.type = FACE::QUADRILATERAL;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(1);
                        curFace.includedNode(2) = curCell.NodeSeq(5);
                        curFace.includedNode(3) = curCell.NodeSeq(8);
                        curFace.includedNode(4) = curCell.NodeSeq(4);

                        curFace.leftCell = adjCell.CellSeq();
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }

            /// Internal J direction
            for (size_t k = 1; k < nK; ++k)
                for (size_t i = 1; i < nI; ++i)
                    for (size_t j = 2; j < nJ; ++j)
                    {
                        const auto &curCell = b.cell(i, j, k);
                        const auto &adjCell = b.cell(i, j - 1, k);

                        const auto faceIndex = curCell.FaceSeq(3);
                        auto &curFace = face(faceIndex);

                        curFace.atBdry = false;

                        curFace.type = FACE::QUADRILATERAL;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(6);
                        curFace.includedNode(2) = curCell.NodeSeq(5);
                        curFace.includedNode(3) = curCell.NodeSeq(1);
                        curFace.includedNode(4) = curCell.NodeSeq(2);

                        curFace.leftCell = adjCell.CellSeq();
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }

            /// Internal k direction
            for (size_t i = 1; i < nI; ++i)
                for (size_t j = 1; j < nJ; ++j)
                    for (size_t k = 2; k < nK; ++k)
                    {
                        const auto &curCell = b.cell(i, j, k);
                        const auto &adjCell = b.cell(i, j, k - 1);

                        const auto faceIndex = curCell.FaceSeq(5);
                        auto &curFace = face(faceIndex);

                        curFace.atBdry = false;

                        curFace.type = FACE::QUADRILATERAL;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(4);
                        curFace.includedNode(2) = curCell.NodeSeq(3);
                        curFace.includedNode(3) = curCell.NodeSeq(2);
                        curFace.includedNode(4) = curCell.NodeSeq(1);

                        curFace.leftCell = adjCell.CellSeq();
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }

            /// I-MIN
            for (size_t k = 1; k < nK; ++k)
                for (size_t j = 1; j < nJ; ++j)
                {
                    const auto &curCell = b.cell(1, j, k);
                    const auto faceIndex = curCell.FaceSeq(1);
                    auto &curFace = face(faceIndex);
                    const auto &curSurf = b.surf(1);

                    if (visited[faceIndex - 1])
                    {
                        if (!curFace.atBdry)
                        {
                            /// Second-Round of Double-Sided Case
                            /// Update undetermined cell index.
                            if (curFace.leftCell == 0)
                                curFace.leftCell = curCell.CellSeq();
                            else if (curFace.rightCell == 0)
                                curFace.rightCell = curCell.CellSeq();
                            else
                                throw std::runtime_error("Double-Sided face should not appear more than twice!");
                        }
                        else
                            throw std::runtime_error("Boundary face shouldn't appear twice!");
                    }
                    else
                    {
                        curFace.type = FACE::QUADRILATERAL;

                        curFace.atBdry = !curSurf.neighbourSurf;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(1);
                        curFace.includedNode(2) = curCell.NodeSeq(5);
                        curFace.includedNode(3) = curCell.NodeSeq(8);
                        curFace.includedNode(4) = curCell.NodeSeq(4);

                        /// On I-MIN Surface, if current face is Single-Sided,
                        /// then index of left cell is set to 0 according to
                        /// right-hand convention; If current face is Double-Sided,
                        /// it is also set to 0 at this stage, and will be
                        /// updated further by loops of other blocks.
                        curFace.leftCell = 0;
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }
                }

            /// I-MAX
            for (size_t k = 1; k < nK; ++k)
                for (size_t j = 1; j < nJ; ++j)
                {
                    const auto &curCell = b.cell(nI - 1, j, k);
                    const auto faceIndex = curCell.FaceSeq(2);
                    auto &curFace = face(faceIndex);
                    const auto &curSurf = b.surf(2);

                    if (visited[faceIndex - 1])
                    {
                        if (!curFace.atBdry)
                        {
                            /// Second-Round of Double-Sided Case
                            /// Update undetermined cell index.
                            if (curFace.leftCell == 0)
                                curFace.leftCell = curCell.CellSeq();
                            else if (curFace.rightCell == 0)
                                curFace.rightCell = curCell.CellSeq();
                            else
                                throw std::runtime_error("Double-Sided face should not appear more than twice!");
                        }
                        else
                            throw std::runtime_error("Boundary face shouldn't appear twice!");
                    }
                    else
                    {
                        curFace.type = FACE::QUADRILATERAL;

                        curFace.atBdry = !curSurf.neighbourSurf;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(2);
                        curFace.includedNode(2) = curCell.NodeSeq(3);
                        curFace.includedNode(3) = curCell.NodeSeq(7);
                        curFace.includedNode(4) = curCell.NodeSeq(6);

                        /// On I-MAX Surface, if current face is Single-Sided,
                        /// then index of left cell is set to 0 according to
                        /// right-hand convention; If current face is Double-Sided,
                        /// it is also set to 0 at this stage, and will be
                        /// updated further by loops of other blocks.
                        curFace.leftCell = 0;
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }
                }

            /// J-MIN
            for (size_t k = 1; k < nK; ++k)
                for (size_t i = 1; i < nI; ++i)
                {
                    const auto &curCell = b.cell(i, 1, k);
                    const auto faceIndex = curCell.FaceSeq(3);
                    auto &curFace = face(faceIndex);
                    const auto &curSurf = b.surf(3);

                    if (visited[faceIndex - 1])
                    {
                        if (!curFace.atBdry)
                        {
                            /// Second-Round of Double-Sided Case
                            /// Update undetermined cell index.
                            if (curFace.leftCell == 0)
                                curFace.leftCell = curCell.CellSeq();
                            else if (curFace.rightCell == 0)
                                curFace.rightCell = curCell.CellSeq();
                            else
                                throw std::runtime_error("Double-Sided face should not appear more than twice!");
                        }
                        else
                            throw std::runtime_error("Boundary face shouldn't appear twice!");
                    }
                    else
                    {
                        curFace.type = FACE::QUADRILATERAL;

                        curFace.atBdry = !curSurf.neighbourSurf;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(6);
                        curFace.includedNode(2) = curCell.NodeSeq(5);
                        curFace.includedNode(3) = curCell.NodeSeq(1);
                        curFace.includedNode(4) = curCell.NodeSeq(2);

                        /// On J-MIN Surface, if current face is Single-Sided,
                        /// then index of left cell is set to 0 according to
                        /// right-hand convention; If current face is Double-Sided,
                        /// it is also set to 0 at this stage, and will be
                        /// updated further by loops of other blocks.
                        curFace.leftCell = 0;
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }
                }

            /// J-MAX
            for (size_t k = 1; k < nK; ++k)
                for (size_t i = 1; i < nI; ++i)
                {
                    const auto &curCell = b.cell(i, nJ - 1, k);
                    const auto faceIndex = curCell.FaceSeq(4);
                    auto &curFace = face(faceIndex);
                    const auto &curSurf = b.surf(4);

                    if (visited[faceIndex - 1])
                    {
                        if (!curFace.atBdry)
                        {
                            /// Second-Round of Double-Sided Case
                            /// Update undetermined cell index.
                            if (curFace.leftCell == 0)
                                curFace.leftCell = curCell.CellSeq();
                            else if (curFace.rightCell == 0)
                                curFace.rightCell = curCell.CellSeq();
                            else
                                throw std::runtime_error("Double-Sided face should not appear more than twice!");
                        }
                        else
                            throw std::runtime_error("Boundary face shouldn't appear twice!");
                    }
                    else
                    {
                        curFace.type = FACE::QUADRILATERAL;

                        curFace.atBdry = !curSurf.neighbourSurf;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(3);
                        curFace.includedNode(2) = curCell.NodeSeq(4);
                        curFace.includedNode(3) = curCell.NodeSeq(8);
                        curFace.includedNode(4) = curCell.NodeSeq(7);

                        /// On J-MAX Surface, if current face is Single-Sided,
                        /// then index of left cell is set to 0 according to
                        /// right-hand convention; If current face is Double-Sided,
                        /// it is also set to 0 at this stage, and will be
                        /// updated further by loops of other blocks.
                        curFace.leftCell = 0;
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }
                }

            /// K-MIN
            for (size_t i = 1; i < nI; ++i)
                for (size_t j = 1; j < nJ; ++j)
                {
                    const auto &curCell = b.cell(i, j, 1);
                    const auto faceIndex = curCell.FaceSeq(5);
                    auto &curFace = face(faceIndex);
                    const auto &curSurf = b.surf(5);

                    if (visited[faceIndex - 1])
                    {
                        if (!curFace.atBdry)
                        {
                            /// Second-Round of Double-Sided Case
                            /// Update undetermined cell index.
                            if (curFace.leftCell == 0)
                                curFace.leftCell = curCell.CellSeq();
                            else if (curFace.rightCell == 0)
                                curFace.rightCell = curCell.CellSeq();
                            else
                                throw std::runtime_error("Double-Sided face should not appear more than twice!");
                        }
                        else
                            throw std::runtime_error("Boundary face shouldn't appear twice!");
                    }
                    else
                    {
                        curFace.type = FACE::QUADRILATERAL;

                        curFace.atBdry = !curSurf.neighbourSurf;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(4);
                        curFace.includedNode(2) = curCell.NodeSeq(3);
                        curFace.includedNode(3) = curCell.NodeSeq(2);
                        curFace.includedNode(4) = curCell.NodeSeq(1);

                        /// On K-MIN Surface, if current face is Single-Sided,
                        /// then index of left cell is set to 0 according to
                        /// right-hand convention; If current face is Double-Sided,
                        /// it is also set to 0 at this stage, and will be
                        /// updated further by loops of other blocks.
                        curFace.leftCell = 0;
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }
                }

            /// K-MAX
            for (size_t i = 1; i < nI; ++i)
                for (size_t j = 1; j < nJ; ++j)
                {
                    const auto &curCell = b.cell(i, j, nK - 1);
                    const auto faceIndex = curCell.FaceSeq(6);
                    auto &curFace = face(faceIndex);
                    const auto &curSurf = b.surf(6);

                    if (visited[faceIndex - 1])
                    {
                        if (!curFace.atBdry)
                        {
                            /// Second-Round of Double-Sided Case
                            /// Update undetermined cell index.
                            if (curFace.leftCell == 0)
                                curFace.leftCell = curCell.CellSeq();
                            else if (curFace.rightCell == 0)
                                curFace.rightCell = curCell.CellSeq();
                            else
                                throw std::runtime_error("Double-Sided face should not appear more than twice!");
                        }
                        else
                            throw std::runtime_error("Boundary face shouldn't appear twice!");
                    }
                    else
                    {
                        curFace.type = FACE::QUADRILATERAL;

                        curFace.atBdry = !curSurf.neighbourSurf;

                        curFace.includedNode.resize(4);
                        curFace.includedNode(1) = curCell.NodeSeq(8);
                        curFace.includedNode(2) = curCell.NodeSeq(5);
                        curFace.includedNode(3) = curCell.NodeSeq(6);
                        curFace.includedNode(4) = curCell.NodeSeq(7);

                        /// On K-MAX Surface, if current face is Single-Sided,
                        /// then index of left cell is set to 0 according to
                        /// right-hand convention; If current face is Double-Sided,
                        /// it is also set to 0 at this stage, and will be
                        /// updated further by loops of other blocks.
                        curFace.leftCell = 0;
                        curFace.rightCell = curCell.CellSeq();

                        visited[faceIndex - 1] = true;
                    }
                }
        }

        /// Assign ZONE info.
        /// Here, only possible choice for cell is hex;
        /// only possible choice for face is quad.
        /// Index assignment convention:
        ///   1 : Nodal coordinates.
        ///   2 : Cell specification.
        ///   3 : Internal faces.
        ///   4 - N : Boundary faces. 

        /// Count total zones.
        m_totalZoneNum = 0;
        for (size_t i = 1; i <= NBLK; ++i)
        {
            const auto &b = nmf->block(i);
            for (short j = 1; j <= NMF::Block3D::NumOfSurf; ++j)
            {
                const auto &s = b.surf(j);
                if (s.neighbourSurf == nullptr)
                    ++m_totalZoneNum;
            }
        }
        m_totalZoneNum += 3;

        /// Allocate storage.
        m_zone.resize(m_totalZoneNum);

        /// Assign default index.
        for (size_t i = 1; i <= m_totalZoneNum; ++i)
        {
            auto &z = zone(i);
            z.ID = i;
            z.obj = nullptr;
        }

        /// ZONE-1
        auto &zone1 = zone(1);
        zone1.type = "";
        zone1.name = "NODE";

        /// ZONE-2
        auto &zone2 = zone(2);
        zone2.type = "fluid";
        zone2.name = "FLUID";

        /// ZONE-3
        auto &zone3 = zone(3);
        zone3.type = "interior";
        zone3.name = "int_FLUID";

        /// ZONE-4 and later
        size_t patch_idx = 4;
        for (size_t i = 1; i <= NBLK; ++i)
        {
            const auto &b = nmf->block(i);
            for (short j = 1; j <= NMF::Block3D::NumOfSurf; ++j)
            {
                const auto &s = b.surf(j);
                if (s.neighbourSurf == nullptr)
                {
                    auto &z = zone(patch_idx);
                    ++patch_idx;

                    /// Set to "wall" by default.
                    /// Can be mapped according to NMF specification or assigned mannually in FLUENT.
                    z.type = "wall";

                    z.name = "B" + std::to_string(i) + "F" + std::to_string(j);
                }
            }
        }

        /// Convert to primary form.
        clear_entry();

        add_entry(new HEADER("Block-Glue " + version_str()));
        add_entry(new DIMENSION(3));

        /// Nodal coordinates
        auto part1 = new NODE(1, 1, numOfNode(), NODE::ANY, 3);
        for (size_t i = 0; i < numOfNode(); ++i)
            part1->at(i) = node(i + 1).coordinate;
        zone(1).obj = part1;
        add_entry(part1);

        /// Cell specifications
        auto part2 = new CELL(2, 1, numOfCell(), CELL::FLUID, CELL::HEXAHEDRAL);
        zone(2).obj = part2;
        add_entry(part2);

        /// Internal faces connectivity
        size_t face_pos_L = 1;
        size_t face_pos_R = face_pos_L + innerFaceNum - 1;
        auto part3 = new FACE(3, face_pos_L, face_pos_R, BC::INTERIOR, FACE::QUADRILATERAL);
        for (size_t i = 0; i < innerFaceNum; ++i)
        {
            auto &raw_f = part3->at(i);
            const auto &derived_f = face(i + 1);

            raw_f.x = 4;

            raw_f.n[0] = derived_f.includedNode.at(0);
            raw_f.n[1] = derived_f.includedNode.at(1);
            raw_f.n[2] = derived_f.includedNode.at(2);
            raw_f.n[3] = derived_f.includedNode.at(3);

            raw_f.c[0] = derived_f.rightCell;
            raw_f.c[1] = derived_f.leftCell;
        }
        zone(3).obj = part3;
        add_entry(part3);

        /// Boundary faces connectivity
        patch_idx = 4;
        for (size_t i = 1; i <= NBLK; ++i)
        {
            const auto &b = nmf->block(i);
            for (short j = 1; j <= NMF::Block3D::NumOfSurf; ++j)
            {
                const auto &s = b.surf(j);
                if (s.neighbourSurf == nullptr)
                {
                    const size_t cfn = b.surface_face_num(j);

                    face_pos_L = face_pos_R + 1;
                    face_pos_R = face_pos_L + cfn - 1;

                    auto part_bfi = new FACE(patch_idx, face_pos_L, face_pos_R, BC::WALL, FACE::QUADRILATERAL);
                    for (size_t k = 0; k < cfn; ++k)
                    {
                        auto &raw_f = part_bfi->at(k);
                        const auto &derived_f = face(k + face_pos_L);

                        raw_f.x = 4;

                        raw_f.n[0] = derived_f.includedNode.at(0);
                        raw_f.n[1] = derived_f.includedNode.at(1);
                        raw_f.n[2] = derived_f.includedNode.at(2);
                        raw_f.n[3] = derived_f.includedNode.at(3);

                        raw_f.c[0] = derived_f.rightCell;
                        raw_f.c[1] = derived_f.leftCell;
                    }

                    zone(patch_idx).obj = part_bfi;
                    ++patch_idx;
                    add_entry(part_bfi);
                }
            }
        }

        /// Zone mapping relations.
        /// Here storage index is consistent with real index.
        m_zoneMapping.clear();
        for (size_t i = 1; i <= numOfZone(); ++i)
            m_zoneMapping.emplace(i, i);

        /// Zone specifications.
        /// No need to show NODE zone.
        add_entry(new COMMENT("Zone Sections"));
        for (size_t i = 2; i <= numOfZone(); ++i)
        {
            const auto &derived_z = zone(i);
            auto raw_z = new ZONE(i, derived_z.type, derived_z.name);
            add_entry(raw_z);
        }

        /// Finalize.
        delete nmf;
        delete p3d;
    }
}
