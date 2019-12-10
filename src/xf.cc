#include "../inc/xf.h"

namespace GridTool
{
	namespace XF
	{
		void MESH::writeToFile(const std::string &dst) const
		{
			if (numOfCell() == 0)
				throw std::runtime_error("Invalid num of cells.");
			if (numOfFace() == 0)
				throw std::runtime_error("Invalid num of faces.");
			if (numOfNode() == 0)
				throw std::runtime_error("Invalid num of nodes.");
			if (m_content.empty())
				throw std::runtime_error("Invalid num of contents.");

			// Open grid file
			std::ofstream fout(dst);
			if (fout.fail())
				throw std::runtime_error("Failed to open output grid file: " + dst);

			// Write until dimension declaration
			size_t i = 0;
			while (true)
			{
				m_content[i]->repr(fout);
				bool flag = dynamic_cast<DIMENSION*>(m_content[i]) != nullptr;
				++i;
				if (flag)
					break;
			}

			// Declaration of NODE, FACE, CELL
			fout << "(" << std::dec << SECTION::NODE << " (";
			fout << std::hex << 0 << " " << 1 << " " << m_totalNodeNum << " ";
			fout << std::dec << 0 << " " << (m_is3D ? 3 : 2) << "))" << std::endl;
			fout << "(" << std::dec << SECTION::CELL << " (";
			fout << std::hex << 0 << " " << 1 << " " << m_totalCellNum << " ";
			fout << std::dec << 0 << " " << 0 << "))" << std::endl;
			fout << "(" << std::dec << SECTION::FACE << " (";
			fout << std::hex << 0 << " " << 1 << " " << m_totalFaceNum << " ";
			fout << std::dec << 0 << " " << 0 << "))" << std::endl;

			// Contents
			for (; i < m_content.size(); ++i)
				m_content[i]->repr(fout);

			// Close grid file
			fout.close();
		}

		void MESH::cell_standardization(CELL_ELEM &c)
		{
			switch (c.type)
			{
			case CELL::TETRAHEDRAL:
				tet_standardization(c);
				break;
			case CELL::HEXAHEDRAL:
				hex_standardization(c);
				break;
			case CELL::PYRAMID:
				pyramid_standardization(c);
				break;
			case CELL::WEDGE:
				prism_standardization(c);
				break;
			case CELL::TRIANGULAR:
				triangle_standardization(c);
				break;
			case CELL::QUADRILATERAL:
				quad_standardization(c);
				break;
			default:
				throw std::runtime_error("Invalid cell element type: " + std::to_string(c.type) + ", maybe caused by internal error.");
			}
		}

		void MESH::tet_standardization(CELL_ELEM &tet)
		{
			// Check num of total faces
			if (tet.includedFace.size() != 4)
				throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::elem_type_idx2str(tet.type) + R"(" and num of faces: )" + std::to_string(tet.includedFace.size()));

			// Ensure all faces are triangular
			const auto &f0 = face(tet.includedFace.at(0));
			if (f0.type != FACE::TRIANGULAR)
				throw std::runtime_error("Internal error.");

			const auto &f1 = face(tet.includedFace.at(1));
			if (f1.type != FACE::TRIANGULAR)
				throw std::runtime_error("Internal error.");

			const auto &f2 = face(tet.includedFace.at(2));
			if (f2.type != FACE::TRIANGULAR)
				throw std::runtime_error("Internal error.");

			const auto &f3 = face(tet.includedFace.at(3));
			if (f3.type != FACE::TRIANGULAR)
				throw std::runtime_error("Internal error.");

			// Find all 4 nodes
			size_t n1 = f0.includedNode.at(0);
			if (n1 == 0)
				throw std::runtime_error("Internal error.");

			size_t n2 = f0.includedNode.at(1);
			if (n2 == 0)
				throw std::runtime_error("Internal error.");

			size_t n3 = f0.includedNode.at(2);
			if (n3 == 0)
				throw std::runtime_error("Internal error.");

			size_t n0 = 0;
			for (auto e : f2.includedNode)
				if (e != n1 && e != n2 && e != n3)
				{
					n0 = e;
					break;
				}
			if (n0 == 0)
				throw std::runtime_error("Internal error.");

			// Assign node index
			tet.includedNode.resize(4);
			tet.includedNode.at(0) = n0;
			tet.includedNode.at(1) = n1;
			tet.includedNode.at(2) = n2;
			tet.includedNode.at(3) = n3;
		}

		void MESH::pyramid_standardization(CELL_ELEM &pyramid)
		{
			// Check num of total faces
			if (pyramid.includedFace.size() != 5)
				throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::elem_type_idx2str(pyramid.type) + R"(" and num of faces: )" + std::to_string(pyramid.includedFace.size()));

			// Find the bottom quad and ensure other faces are triangular.
			size_t f0_idx = 0;
			for (auto e : pyramid.includedFace)
			{
				const auto &f = face(e);
				if (f.type == FACE::QUADRILATERAL)
				{
					if (f0_idx == 0)
						f0_idx = e;
					else
						throw std::runtime_error("Internal error.");
				}
				else if (f.type == FACE::TRIANGULAR)
					continue;
				else
					throw std::runtime_error("Internal error.");
			}
			if (f0_idx == 0)
				throw std::runtime_error("Internal error.");

			// Nodes at bottom
			const auto &f0 = face(f0_idx);

			const size_t n0 = f0.includedNode.at(0);
			if (n0 == 0)
				throw std::runtime_error("Internal error.");

			const size_t n1 = f0.includedNode.at(1);
			if (n1 == 0)
				throw std::runtime_error("Internal error.");

			const size_t n2 = f0.includedNode.at(2);
			if (n2 == 0)
				throw std::runtime_error("Internal error.");

			const size_t n3 = f0.includedNode.at(3);
			if (n3 == 0)
				throw std::runtime_error("Internal error.");

			// Find other 4 triangles
			size_t f1_idx = 0;
			for (auto e : pyramid.includedFace)
			{
				if (e == f0_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n0, n3))
				{
					f1_idx = e;
					break;
				}
			}
			if (f1_idx == 0)
				throw std::runtime_error("Internal error.");

			size_t f2_idx = 0;
			for (auto e : pyramid.includedFace)
			{
				if (e == f0_idx || e == f1_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n3, n2))
				{
					f2_idx = e;
					break;
				}
			}
			if (f2_idx == 0)
				throw std::runtime_error("Internal error.");

			size_t f3_idx = 0;
			for (auto e : pyramid.includedFace)
			{
				if (e == f0_idx || e == f1_idx || e == f2_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n2, n1))
				{
					f3_idx = e;
					break;
				}
			}
			if (f3_idx == 0)
				throw std::runtime_error("Internal error.");

			size_t f4_idx = 0;
			for (auto e : pyramid.includedFace)
			{
				if (e != f0_idx && e != f1_idx && e != f2_idx && e != f3_idx)
				{
					f4_idx = e;
					const auto &f = face(e);
					if (!f.includedNode.contains(n1, n0))
						throw std::runtime_error("Internal error.");

					break;
				}
			}
			if (f4_idx == 0)
				throw std::runtime_error("Internal error.");

			// Assign face index
			pyramid.includedFace.at(0) = f0_idx;
			pyramid.includedFace.at(1) = f1_idx;
			pyramid.includedFace.at(2) = f2_idx;
			pyramid.includedFace.at(3) = f3_idx;
			pyramid.includedFace.at(4) = f4_idx;

			// The last node
			size_t n4 = 0;
			const auto &f1 = face(f1_idx);
			for (auto e : f1.includedNode)
				if (e != n0 && e != n3)
				{
					n4 = e;
					break;
				}
			if (n4 == 0)
				throw std::runtime_error("Internal error.");

			// Assign node index
			pyramid.includedNode.resize(5);
			pyramid.includedNode.at(0) = n0;
			pyramid.includedNode.at(1) = n1;
			pyramid.includedNode.at(2) = n2;
			pyramid.includedNode.at(3) = n3;
			pyramid.includedNode.at(4) = n4;
		}

		void MESH::prism_standardization(CELL_ELEM &prism)
		{
			// Check num of total faces
			if (prism.includedFace.size() != 5)
				throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::elem_type_idx2str(prism.type) + R"(" and num of faces: )" + std::to_string(prism.includedFace.size()));

			// Ensure there're only 2 triangle and 3 quad
			size_t f0_idx = 0, f1_idx = 0;
			for (auto e : prism.includedFace)
			{
				const auto &f = face(e);
				if (f.type == FACE::TRIANGULAR)
				{
					if (f0_idx == 0)
						f0_idx = e;
					else if (f1_idx == 0)
						f1_idx = e;
					else
						throw std::runtime_error("There're more than 2 triangular faces in a prism cell.");
				}
				else if (f.type == FACE::QUADRILATERAL)
					continue;
				else
					throw std::runtime_error("Internal error.");
			}
			if (f0_idx == 0 || f1_idx == 0)
				throw std::runtime_error("Missing triangular faces in a prism cell.");

			// 2 triangular faces
			const auto &f0 = face(f0_idx);
			const auto &f1 = face(f1_idx);

			// 3 nodes on the bottom triangular face
			const size_t n0 = f0.includedNode.at(0);
			const size_t n1 = f0.includedNode.at(1);
			const size_t n2 = f0.includedNode.at(2);

			// Find face 4
			size_t f4_idx = 0;
			for (auto e : prism.includedFace)
			{
				const auto &f = face(e);
				if (f.type == FACE::QUADRILATERAL && f.includedNode.contains(n0, n1))
				{
					if (f4_idx == 0)
						f4_idx = e;
					else
						throw std::runtime_error("Inconsistent face composition.");
				}
			}
			if (f4_idx == 0)
				throw std::runtime_error("Missing face 4.");

			const auto &f4 = face(f4_idx);

			// Find face 3
			size_t f3_idx = 0;
			for (auto e : prism.includedFace)
			{
				const auto &f = face(e);
				if (f.type == FACE::QUADRILATERAL && f.includedNode.contains(n1, n2))
				{
					if (f3_idx == 0)
						f3_idx = e;
					else
						throw std::runtime_error("Inconsistent face composition.");
				}
			}
			if (f3_idx == 0 || f3_idx == f4_idx)
				throw std::runtime_error("Missing face 3.");

			const auto &f3 = face(f3_idx);

			// Find face 2
			size_t f2_idx = 0;
			for (auto e : prism.includedFace)
			{
				const auto &f = face(e);
				if (f.type == FACE::QUADRILATERAL)
				{
					if (e != f4_idx && e != f3_idx)
					{
						f2_idx = e;
						break;
					}
				}
			}
			if (f2_idx == 0)
				throw std::runtime_error("Missing face 2.");

			const auto &f2 = face(f2_idx);
			if (!f2.includedNode.contains(n2, n0))
				throw std::runtime_error("Inconsistent face composition.");

			// Assign face index
			prism.includedFace.at(0) = f0_idx;
			prism.includedFace.at(1) = f1_idx;
			prism.includedFace.at(2) = f2_idx;
			prism.includedFace.at(3) = f3_idx;
			prism.includedFace.at(4) = f4_idx;

			// 3 nodes on the top triangular face to be defined
			size_t n3 = 0, n5 = 0, n4 = 0;
			for (auto e : f2.includedNode)
			{
				if (e != n0 && e != n2)
				{
					if (n3 == 0)
						n3 = e;
					else if (n5 == 0)
						n5 = e;
					else
						throw std::runtime_error("Internal error.");
				}
			}
			if (n3 == 0 || n5 == 0 || n3 == n5)
				throw std::runtime_error("Missing nodes on the top");

			// Check if n3 and n5 needs to be swaped
			if (!(f1.includedNode.contains(n3) && f4.includedNode.contains(n3)))
				std::swap(n3, n5);

			// Ensure n5 located on f1, f2 and f3
			if (!(f1.includedNode.contains(n5) && f3.includedNode.contains(n5)))
				throw std::runtime_error("Internal error.");

			// Find n4
			for (auto e : f1.includedNode)
			{
				if (e != n3 && e != n5)
				{
					if (n4 == 0)
						n4 = e;
					else
						throw std::runtime_error("Internal error.");
				}
			}
			if (n4 == 0)
				throw std::runtime_error("Internal error.");

			// Assign node index
			prism.includedNode.resize(6);
			prism.includedNode.at(0) = n0;
			prism.includedNode.at(1) = n1;
			prism.includedNode.at(2) = n2;
			prism.includedNode.at(3) = n3;
			prism.includedNode.at(4) = n4;
			prism.includedNode.at(5) = n5;
		}

		void MESH::hex_standardization(CELL_ELEM &hex)
		{
			// Check num of total faces
			if (hex.includedFace.size() != 6)
				throw std::runtime_error(R"(Mismatch between cell type ")" + CELL::elem_type_idx2str(hex.type) + R"(" and num of faces: )" + std::to_string(hex.includedFace.size()));

			// Ensure all faces are quad
			for (auto e : hex.includedFace)
			{
				if (e == 0)
					throw std::runtime_error("Internal error.");
				const auto &f = face(e);
				if (f.type != FACE::QUADRILATERAL)
					throw std::runtime_error(R"(Inconsistent face type ")" + FACE::idx2str(f.type) + R"(" in a hex cell.)");
			}

			// Face 4 at bottom
			const size_t f4_idx = hex.includedFace.at(0);
			const auto &f4 = face(f4_idx);
			const size_t n0 = f4.includedNode.at(0);
			const size_t n1 = f4.includedNode.at(1);
			const size_t n2 = f4.includedNode.at(2);
			const size_t n3 = f4.includedNode.at(3);
			if (n0 == 0 || n1 == 0 || n2 == 0 || n3 == 0)
				throw std::runtime_error("Internal error.");

			// Face 0
			size_t f0_idx = 0;
			for (auto e : hex.includedFace)
			{
				if (e == f4_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n3, n0))
				{
					if (f0_idx == 0)
						f0_idx = e;
					else
						throw std::runtime_error("Duplicated face detected.");
				}
			}
			if (f0_idx == 0)
				throw std::runtime_error("Missing face 0");

			const auto &f0 = face(f0_idx);
			if (f0.includedNode.size() != 4)
				throw std::runtime_error("Internal error.");

			// Face 2
			size_t f2_idx = 0;
			for (auto e : hex.includedFace)
			{
				if (e == f4_idx || e == f0_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n0, n1))
				{
					if (f2_idx == 0)
						f2_idx = e;
					else
						throw std::runtime_error("Duplicated face detected.");
				}
			}
			if (f2_idx == 0)
				throw std::runtime_error("Missing face 2");

			const auto &f2 = face(f2_idx);
			if (f2.includedNode.size() != 4)
				throw std::runtime_error("Internal error.");

			// Face 1
			size_t f1_idx = 0;
			for (auto e : hex.includedFace)
			{
				if (e == f4_idx || e == f0_idx || e == f2_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n1, n2))
				{
					if (f1_idx == 0)
						f1_idx = e;
					else
						throw std::runtime_error("Duplicated face detected.");
				}
			}
			if (f1_idx == 0)
				throw std::runtime_error("Missing face 1");

			const auto &f1 = face(f1_idx);
			if (f1.includedNode.size() != 4)
				throw std::runtime_error("Internal error.");

			// Face 3
			size_t f3_idx = 0;
			for (auto e : hex.includedFace)
			{
				if (e == f4_idx || e == f0_idx || e == f2_idx || e == f1_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n2, n3))
				{
					if (f3_idx == 0)
						f3_idx = e;
					else
						throw std::runtime_error("Duplicated face detected.");
				}
			}
			if (f3_idx == 0)
				throw std::runtime_error("Missing face 3");

			const auto &f3 = face(f3_idx);
			if (f3.includedNode.size() != 4)
				throw std::runtime_error("Internal error.");

			// Face 5
			size_t f5_idx = 0;
			for (auto e : hex.includedFace)
			{
				if (e == f4_idx || e == f0_idx || e == f2_idx || e == f1_idx || e == f3_idx)
					continue;

				f5_idx = e;
				break;
			}
			if (f5_idx == 0)
				throw std::runtime_error("Missing face 5");

			const auto &f5 = face(f5_idx);
			if (f5.includedNode.size() != 4)
				throw std::runtime_error("Internal error.");

			// Assign face index
			hex.includedFace.at(0) = f0_idx;
			hex.includedFace.at(1) = f1_idx;
			hex.includedFace.at(2) = f2_idx;
			hex.includedFace.at(3) = f3_idx;
			hex.includedFace.at(4) = f4_idx;
			hex.includedFace.at(5) = f5_idx;

			// 4 nodes at top to be defined
			size_t n4 = 0, n5 = 0, n6 = 0, n7 = 0;

			// Node 4
			for (auto e : f5.includedNode)
			{
				if (f0.includedNode.contains(e) && f2.includedNode.contains(e))
				{
					n4 = e;
					break;
				}
			}
			if (n4 == 0 || n4 == n0)
				throw std::runtime_error("Missing node 4");

			// Node 5
			for (auto e : f2.includedNode)
			{
				if (e == n0 || e == n1 || e == n4)
					continue;

				n5 = e;
				break;
			}
			if (n5 == 0)
				throw std::runtime_error("Missing node 5");

			// Node 6
			for (auto e : f1.includedNode)
			{
				if (e == n1 || e == n2 || e == n5)
					continue;

				n6 = e;
				break;
			}
			if (n6 == 0)
				throw std::runtime_error("Missing node 6");

			// Node 7
			for (auto e : f0.includedNode)
			{
				if (e == n0 || e == n3 || e == n4)
					continue;

				n7 = e;
				break;
			}
			if (n7 == 0)
				throw std::runtime_error("Missing node 7");

			if (!f3.includedNode.contains(n6, n7))
				throw std::runtime_error("Inconsistent node composition.");

			// Assign node index
			hex.includedNode.resize(8);
			hex.includedNode.at(0) = n0;
			hex.includedNode.at(1) = n1;
			hex.includedNode.at(2) = n2;
			hex.includedNode.at(3) = n3;
			hex.includedNode.at(4) = n4;
			hex.includedNode.at(5) = n5;
			hex.includedNode.at(6) = n6;
			hex.includedNode.at(7) = n7;
		}

		void MESH::triangle_standardization(CELL_ELEM &tri)
		{
			// Check num of total faces
			if (tri.includedFace.size() != 3)
				throw std::runtime_error("Inconsitent face composition.");

			// Ensure all faces are lines
			for (auto e : tri.includedFace)
			{
				const auto &f = face(e);
				if (e == 0 || f.type != FACE::LINEAR || f.includedNode.size() != 2)
					throw std::runtime_error("Invalid face detected.");
			}

			// Faces
			// Keep the order of faces as it is.
			// Doesn't matter as any order of the 3 faces is valid.
			const auto f0_idx = tri.includedFace.at(0);
			const auto f1_idx = tri.includedFace.at(1);
			const auto f2_idx = tri.includedFace.at(2);

			const auto &f0 = face(f0_idx);
			const auto &f1 = face(f1_idx);
			const auto &f2 = face(f2_idx);

			// Nodes
			const size_t n0 = f0.includedNode.at(0);
			const size_t n1 = f0.includedNode.at(1);
			if (n0 == 0 || n1 == 0 || n0 == n1)
				throw std::runtime_error("Missing node detected.");

			size_t n2 = 0;
			for (auto e : f1.includedNode)
			{
				if (f0.includedNode.contains(e))
					continue;

				n2 = e;
				break;
			}
			if (n2 == 0 || n2 == n0 || n2 == n1 || !f2.includedNode.contains(n2))
				throw std::runtime_error("Node 2 is missing.");

			// Assign node index
			tri.includedNode.resize(3);
			tri.includedNode.at(0) = n0;
			tri.includedNode.at(1) = n1;
			tri.includedNode.at(2) = n2;
		}

		void MESH::quad_standardization(CELL_ELEM &quad)
		{
			// Check num of total faces
			if (quad.includedFace.size() != 4)
				throw std::runtime_error("Inconsitent face composition.");

			// Ensure all faces are lines
			for (auto e : quad.includedFace)
			{
				const auto &f = face(e);
				if (e == 0 || f.type != FACE::LINEAR || f.includedNode.size() != 2)
					throw std::runtime_error("Invalid face detected.");
			}

			// Face 0
			const auto f0_idx = quad.includedFace.at(0);
			const auto &f0 = face(f0_idx);

			// Node 0 and 1
			const auto n0 = f0.includedNode.at(0);
			const auto n1 = f0.includedNode.at(1);

			// Face 1
			size_t f1_idx = 0;
			for (auto e : quad.includedFace)
			{
				if (e == f0_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n1))
				{
					f1_idx = e;
					break;
				}
			}
			if (f1_idx == 0)
				throw std::runtime_error("Missing face 1");

			const auto &f1 = face(f1_idx);

			// Node 2
			const size_t n2 = f1.includedNode.at(0) == n1 ? f1.includedNode.at(1) : f1.includedNode.at(0);

			// Face 3
			size_t f3_idx = 0;
			for (auto e : quad.includedFace)
			{
				if (e == f0_idx || e == f1_idx)
					continue;

				const auto &f = face(e);
				if (f.includedNode.contains(n0))
				{
					f3_idx = e;
					break;
				}
			}
			if (f3_idx == 0)
				throw std::runtime_error("Missing face 3");

			const auto &f3 = face(f3_idx);

			// Node 3
			const size_t n3 = f3.includedNode.at(0) == n0 ? f3.includedNode.at(1) : f3.includedNode.at(0);

			// Check face 2
			size_t f2_idx = 0;
			for (auto e : quad.includedFace)
			{
				if (e == f0_idx || e == f1_idx || e == f3_idx)
					continue;

				f2_idx = e;
				break;
			}
			if (f2_idx == 0)
				throw std::runtime_error("Missing face 2");

			const auto &f2 = face(f2_idx);
			if (!f2.includedNode.contains(n2, n3))
				throw std::runtime_error("Inconsistent node and face includance on face 2");

			// Assign face index
			quad.includedFace.at(1) = f1_idx;
			quad.includedFace.at(2) = f2_idx;
			quad.includedFace.at(3) = f3_idx;

			// Assign node index
			quad.includedNode.resize(4);
			quad.includedNode.at(0) = n0;
			quad.includedNode.at(1) = n1;
			quad.includedNode.at(2) = n2;
			quad.includedNode.at(3) = n3;
		}
	}
}
