#include "xf_msh.h"

XF_NODE::XF_NODE(int zone, int first, int last, int type, int ND) :
	XF_MAIN_RECORD(XF_SECTION::NODE, zone, first, last)
{
	if (type == 0)
		m_type = XF_NODE::VIRTUAL;
	else if (type == 1)
		m_type = XF_NODE::ANY;
	else if (type == 2)
		m_type = XF_NODE::BOUNDARY;
	else
		throw("Invalid description of node type!");

	if (ND == 2)
		m_is3D = false;
	else if (ND == 3)
		m_is3D = true;
	else
		throw("Invalid specification of node dimension!");
	m_dim = ND;

	m_node.resize(ND * num());
	std::fill(m_node.begin(), m_node.end(), 0.0);
}

void XF_NODE::set_node_coordinate(size_t loc_idx, double x0, double x1, double x2)
{
	const size_t stx = loc_idx * 3;

	m_node[stx] = x0;
	m_node[stx + 1] = x1;
	m_node[stx + 2] = x2;
}

void XF_NODE::set_node_coordinate(size_t loc_idx, double x0, double x1)
{
	const size_t stx = loc_idx * 2;

	m_node[stx] = x0;
	m_node[stx + 1] = x1;
}

void XF_NODE::repr(std::ostream & out)
{
	const int n_dim = ND();
	const int N = num();

	out << "(" << std::dec << identity();
	out << " (" << std::hex << zone() << " " << first_index() << " " << last_index() << " ";
	out << std::dec << type() << " " << n_dim << ")(" << std::endl;

	size_t loc_idx = 0;
	out.precision(12);
	for (int i = 0; i < N; ++i)
	{
		for (int k = 0; k < n_dim; ++k)
			out << " " << m_node[loc_idx + k];
		out << std::endl;
		loc_idx += n_dim;
	}

	out << "))" << std::endl;
}

XF_CELL::XF_CELL(int zone, int first, int last, int type, int elem_type) :
	XF_MAIN_RECORD(XF_SECTION::CELL, zone, first, last)
{
	if (type == XF_CELL::FLUID)
		m_type = XF_CELL::FLUID;
	else if (type == XF_CELL::SOLID)
		m_type = XF_CELL::SOLID;
	else if (type == XF_CELL::DEAD)
		m_type = XF_CELL::DEAD;
	else
		throw("Invalid specification of cell type!");

	if (elem_type == XF_CELL::MIXED)
	{
		m_elem = XF_CELL::MIXED;
		m_mixedElemDesc.resize(num());
		std::fill(m_mixedElemDesc.begin(), m_mixedElemDesc.end(), XF_CELL::MIXED);
	}
	else if (elem_type == XF_CELL::TRIANGULAR)
	{
		m_elem = XF_CELL::TRIANGULAR;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == XF_CELL::TETRAHEDRAL)
	{
		m_elem = XF_CELL::TETRAHEDRAL;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == XF_CELL::QUADRILATERAL)
	{
		m_elem = XF_CELL::QUADRILATERAL;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == XF_CELL::HEXAHEDRAL)
	{
		m_elem = XF_CELL::HEXAHEDRAL;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == XF_CELL::PYRAMID)
	{
		m_elem = XF_CELL::PYRAMID;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == XF_CELL::WEDGE)
	{
		m_elem = XF_CELL::WEDGE;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == XF_CELL::POLYHEDRAL)
	{
		m_elem = XF_CELL::POLYHEDRAL;
		m_mixedElemDesc.resize(0);
	}
	else
		throw("Invalid specification of cell element type!");
}

void XF_CELL::repr(std::ostream & out)
{
	out << "(" << std::dec << identity() << " (";
	out << std::hex;
	out << zone() << " " << first_index() << " " << last_index() << " ";
	out << m_type << " " << m_elem << ")";

	if (m_elem != XF_CELL::MIXED)
		out << ")" << std::endl;
	else
	{
		out << "(";

		const size_t N = num();
		for (size_t i = 0; i < N; ++i)
		{
			if (i % 40 == 0)
				out << std::endl;
			out << " " << elem(i);
		}
		out << std::endl << "))" << std::endl;
	}
}

XF_FACE::XF_FACE(int zone, int first, int last, int bc, int face) :
	XF_MAIN_RECORD(XF_SECTION::FACE, zone, first, last)
{
	switch (bc)
	{
	case 2:
		m_bc = XF_BC::INTERIOR;
		break;
	case 3:
		m_bc = XF_BC::WALL;
		break;
	case 4:
		m_bc = XF_BC::PRESSURE_INLET;
		break;
	case 5:
		m_bc = XF_BC::PRESSURE_OUTLET;
		break;
	case 7:
		m_bc = XF_BC::SYMMETRY;
		break;
	case 8:
		m_bc = XF_BC::PERIODIC_SHADOW;
		break;
	case 9:
		m_bc = XF_BC::PRESSURE_FAR_FIELD;
		break;
	case 10:
		m_bc = XF_BC::VELOCITY_INLET;
		break;
	case 12:
		m_bc = XF_BC::PERIODIC;
		break;
	case 14:
		m_bc = XF_BC::FAN;
		break;
	case 20:
		m_bc = XF_BC::MASS_FLOW_INLET;
		break;
	case 24:
		m_bc = XF_BC::INTERFACE;
		break;
	case 31:
		m_bc = XF_BC::PARENT;
		break;
	case 36:
		m_bc = XF_BC::OUTFLOW;
		break;
	case 37:
		m_bc = XF_BC::AXIS;
		break;
	default:
		throw("Invalid specification of B.C. type!");
	}

	switch (face)
	{
	case 0:
		m_face = XF_FACE::MIXED;
		break;
	case 2:
		m_face = XF_FACE::LINEAR;
		break;
	case 3:
		m_face = XF_FACE::TRIANGULAR;
		break;
	case 4:
		m_face = XF_FACE::QUADRILATERAL;
		break;
	case 5:
		m_face = XF_FACE::POLYGONAL;
		throw("Currently not supported!");
	default:
		throw("Invalid face-type!");
	}

	m_connectivity.resize(num());
}

void XF_FACE::repr(std::ostream & out)
{
	out << "(" << std::dec << identity() << " (";
	out << std::hex;
	out << zone() << " " << first_index() << " " << last_index() << " ";
	out << bc_type() << " " << face_type() << ")(" << std::endl;

	const int N = num();
	if (m_face == XF_FACE::MIXED)
	{
		for (int i = 0; i < N; ++i)
		{
			const auto &loc_cnect = m_connectivity[i];
			out << " " << loc_cnect.x;
			for (int i = 0; i < loc_cnect.x; ++i)
				out << " " << loc_cnect.n[i];
			out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
		}
	}
	else
	{
		for (int i = 0; i < N; ++i)
		{
			const auto &loc_cnect = m_connectivity[i];
			for (int i = 0; i < loc_cnect.x; ++i)
				out << " " << loc_cnect.n[i];
			out << " " << loc_cnect.c[0] << " " << loc_cnect.c[1] << std::endl;
		}
	}

	out << "))" << std::endl;
}

void XF_ZONE::repr(std::ostream & out)
{
	out << std::dec;
	out << "(" << identity() << " (" << m_zoneID << " " << m_zoneType << " " << m_zoneName << ")())" << std::endl;
}

int XF_MSH::readFromFile(const std::string & src)
{
	std::ifstream fin(src);
	if (!fin)
		return -1;

	clear_entry();

	while (!fin.eof())
	{
		skipWhite(fin);
		eat(fin, '(');
		int ti;
		fin >> std::dec >> ti;
		if (ti == XF_SECTION::COMMENT)
		{
			eat(fin, '\"');
			std::string ts;
			char tc;
			while ((tc = fin.get()) != '\"')
				ts.push_back(tc);
			eat(fin, ')');
			add_entry(new XF_COMMENT(ts));
			skipWhite(fin);
		}
		else if (ti == XF_SECTION::HEADER)
		{
			eat(fin, '\"');
			std::string ts;
			char tc;
			while ((tc = fin.get()) != '\"')
				ts.push_back(tc);
			eat(fin, ')');
			add_entry(new XF_HEADER(ts));
			skipWhite(fin);
		}
		else if (ti == XF_SECTION::DIMENSION)
		{
			int nd = 0;
			fin >> std::dec >> nd;
			eat(fin, ')');
			add_entry(new XF_DIMENSION(nd));
			skipWhite(fin);
			m_dim = nd;
			m_is3D = (nd == 3);
		}
		else if (ti == XF_SECTION::NODE)
		{
			eat(fin, '(');
			int zone;
			fin >> std::hex >> zone;
			if (zone == 0)
			{
				// If zone-id is 0, indicating total number of nodes in the mesh.
				int tmp;
				fin >> std::hex;
				fin >> tmp;
				if (tmp != 1)
					throw("Invalid \"first-index\" in NODE declaration!");
				fin >> m_totalNodeNum;
				std::cout << "Total number of nodes: " << m_totalNodeNum << std::endl;
				fin >> tmp;
				if (tmp != 0)
					throw("Invalid \"type\" in NODE declaration!");
				char ndc = fin.get();
				if (ndc != ')')
				{
					fin >> tmp;
					eat(fin, ')');
				}
				eat(fin, ')');
			}
			else
			{
				// If zone-id is positive, it indicates the zone to which the nodes belong.
				int first, last, tp, nd;
				fin >> std::hex;
				fin >> first >> last;
				fin >> tp >> nd;
				auto e = new XF_NODE(zone, first, last, tp, nd);
				eat(fin, ')');
				eat(fin, '(');
				std::cout << "Reading " << e->num() << " nodes in zone " << zone << " (from " << first << " to " << last << ") ..." << std::endl;

				if (nd != dimension())
					throw("Inconsistent with previous DIMENSION declaration!");

				if (nd == 3)
				{
					double x, y, z;
					for (int i = first; i <= last; ++i)
					{
						size_t i_loc = i - first;
						fin >> x >> y >> z;
						e->set_node_coordinate(i_loc, x, y, z);
					}
				}
				else
				{
					double x, y;
					for (int i = first; i <= last; ++i)
					{
						size_t i_loc = i - first;
						fin >> x >> y;
						e->set_node_coordinate(i_loc, x, y);
					}
				}
				eat(fin, ')');
				eat(fin, ')');
				std::cout << "Done!" << std::endl;
				add_entry(e);
			}
			skipWhite(fin);
		}
		else if (ti == XF_SECTION::CELL)
		{
			eat(fin, '(');
			int zone;
			fin >> std::hex >> zone;
			if (zone == 0)
			{
				// If zone-id is 0, indicating total number of cells in the mesh.
				int tmp;
				fin >> std::hex;
				fin >> tmp;
				if (tmp != 1)
					throw("Invalid \"first-index\" in CELL declaration!");
				fin >> m_totalCellNum;
				std::cout << "Total number of cells: " << m_totalCellNum << std::endl;
				fin >> tmp;
				if (tmp != 0)
					throw("Invalid \"type\" in CELL declaration!");
				char ndc = fin.get();
				if (ndc != ')')
				{
					fin >> tmp;
					eat(fin, ')');
				}
				eat(fin, ')');
			}
			else
			{
				// If zone-id is positive, it indicates the zone to which the cells belong.
				int first, last, tp, elem;
				fin >> std::hex;
				fin >> first >> last;
				fin >> tp >> elem;
				auto e = new XF_CELL(zone, first, last, tp, elem);
				eat(fin, ')');

				if (elem == 0)
				{
					std::cout << "Reading " << e->num() << " mixed cells in zone " << zone << " (from " << first << " to " << last << ") ..." << std::endl;
					eat(fin, '(');
					for (int i = first; i <= last; ++i)
					{
						fin >> elem;
						size_t i_loc = i - first;
						switch (elem)
						{
						case 1:
							e->elem(i_loc) = XF_CELL::TRIANGULAR;
							break;
						case 2:
							e->elem(i_loc) = XF_CELL::TETRAHEDRAL;
							break;
						case 3:
							e->elem(i_loc) = XF_CELL::QUADRILATERAL;
							break;
						case 4:
							e->elem(i_loc) = XF_CELL::HEXAHEDRAL;
							break;
						case 5:
							e->elem(i_loc) = XF_CELL::PYRAMID;
							break;
						case 6:
							e->elem(i_loc) = XF_CELL::WEDGE;
							break;
						case 7:
							e->elem(i_loc) = XF_CELL::POLYHEDRAL;
							break;
						default:
							throw("Invalid cell type!");
						}
					}
					eat(fin, ')');
					std::cout << "Done!" << std::endl;
				}
				else
					std::cout << e->num() << " " << XF_CELL::cell_name(elem) << " in zone " << zone << " (from " << first << " to " << last << ")" << std::endl;

				eat(fin, ')');
				add_entry(e);
			}
			skipWhite(fin);
		}
		else if (ti == XF_SECTION::FACE)
		{
			eat(fin, '(');
			int zone;
			fin >> std::hex >> zone;
			if (zone == 0)
			{
				// If zone-id is 0, indicating total number of faces in the mesh.
				int tmp;
				fin >> tmp;
				if (tmp != 1)
					throw("Invalid \"first-index\" in FACE declaration!");
				fin >> m_totalFaceNum;
				std::cout << "Total number of faces: " << m_totalFaceNum << std::endl;
				fin >> tmp;
				char ndc = fin.get();
				if (ndc != ')')
				{
					fin >> tmp;
					eat(fin, ')');
				}
				eat(fin, ')');
			}
			else
			{
				// If zone-id is positive, it indicates a regular face section and will be
				// followed by a body containing information about the grid connectivity.
				size_t first, last;
				int bc, face;
				fin >> first >> last;
				fin >> bc >> face;
				auto e = new XF_FACE(zone, first, last, bc, face);
				eat(fin, ')');
				eat(fin, '(');

				std::cout << "Reading " << e->num() << " faces in zone " << zone << " (from " << first << " to " << last << ") ..." << std::endl;

				size_t tmp_n[4];
				size_t tmp_c[2];
				for (size_t i = first; i <= last; ++i)
				{
					// Local index
					size_t i_loc = i - first;

					// Read connectivity record
					for (int j = 0; j < face; ++j)
						fin >> tmp_n[j];
					fin >> tmp_c[0] >> tmp_c[1];

					// Store current connectivity info
					e->connectivity(i_loc).set(face, tmp_n, tmp_c);
				}
				eat(fin, ')');
				eat(fin, ')');
				std::cout << "Done!" << std::endl;
				add_entry(e);
			}
			skipWhite(fin);
		}
		else if (ti == XF_SECTION::ZONE)
		{
			std::cout << "Reading ZONE section ..." << std::endl;
			eat(fin, '(');
			int zone;
			fin >> std::dec >> zone;
			std::string ztp;
			fin >> ztp;
			skipWhite(fin);
			std::string zname;
			char t0;
			while ((t0 = fin.get()) != ')')
				zname.push_back(t0);
			eat(fin, '(');
			eat(fin, ')');
			eat(fin, ')');
			auto e = new XF_ZONE(zone, ztp, zname);
			add_entry(e);
			skipWhite(fin);
		}
		else
			throw("Unsupported section index: " + std::to_string(ti) + "!");
	}

	fin.close();
	return 0;
}

int XF_MSH::writeToFile(const std::string & dst) const
{
	std::ofstream fout(dst);
	if (!fout)
		return -1;

	const size_t N = m_content.size();
	if (!N)
		return -2;

	size_t i = 0;

	// Write until dimension declaration
	while (true)
	{
		m_content[i]->repr(fout);
		bool flag = dynamic_cast<XF_DIMENSION*>(m_content[i]) != nullptr;
		++i;
		if (flag)
			break;
	}

	// Declaration of NODE, FACE, CELL
	fout << "(" << std::dec << XF_SECTION::NODE << " (";
	fout << std::hex << 0 << " " << 1 << " " << m_totalNodeNum << " ";
	fout << std::dec << 0 << " " << (m_is3D ? 3 : 2) << "))" << std::endl;
	fout << "(" << std::dec << XF_SECTION::CELL << " (";
	fout << std::hex << 0 << " " << 1 << " " << m_totalCellNum << " ";
	fout << std::dec << 0 << " " << 0 << "))" << std::endl;
	fout << "(" << std::dec << XF_SECTION::FACE << " (";
	fout << std::hex << 0 << " " << 1 << " " << m_totalFaceNum << " ";
	fout << std::dec << 0 << " " << 0 << "))" << std::endl;

	// Contents
	for (; i < N; ++i)
		m_content[i]->repr(fout);

	// Finalize
	fout.close();
	return 0;
}

int XF_MSH::computeTopology_nodeCoordinates(std::vector<std::vector<double>>& dst) const
{
	// Check output array shape.
	if (dst.size() != numOfNode())
		return -1;
	if (dst[0].size() != dimension())
		return -2;

	// Parse related entries.
	for (size_t r = 0; r < m_content.size(); ++r)
	{
		auto curPtr = m_content[r];
		if (curPtr->identity() == XF_SECTION::NODE)
		{
			auto curObj = dynamic_cast<XF_NODE*>(curPtr);

			const int cur_first = curObj->first_index(); // 1-based index, logical
			const int cur_last = curObj->last_index(); // 1-based index, logical

			for (int i = cur_first; i <= cur_last; ++i)
				curObj->get_node_coordinate(i, dst[i - 1]); // 0-based index in storage
		}
	}

	return 0;
}

int XF_MSH::computeTopology_nodeBoundaryFlag(std::vector<bool>& dst) const
{
	// Check output array shape.
	if (dst.size() != numOfNode())
		return -1;

	// Parse related entries.
	for (size_t r = 0; r < m_content.size(); ++r)
	{
		auto curPtr = m_content[r];
		if (curPtr->identity() == XF_SECTION::NODE)
		{
			auto curObj = dynamic_cast<XF_NODE*>(curPtr);

			const int cur_first = curObj->first_index(); // 1-based index, logical
			const int cur_last = curObj->last_index(); // 1-based index, logical

			const bool flag = curObj->is_boundary_node(); // node type within this zone

			for (int i = cur_first; i <= cur_last; ++i)
				dst[i - 1] = flag;
		}
	}

	return 0;
}

int XF_MSH::computeTopology_nodeAdjacentNode(std::vector<std::vector<size_t>>& dst) const
{
	// Check output array shape.
	if (dst.size() != numOfNode())
		return -1;

	// Parse related entries.
	for (size_t r = 0; r < m_content.size(); ++r)
	{
		auto curPtr = m_content[r];
		if (curPtr->identity() == XF_SECTION::FACE)
		{
			auto curObj = dynamic_cast<XF_FACE*>(curPtr);

			const int cur_first = curObj->first_index(); // 1-based index, logical
			const int cur_last = curObj->last_index(); // 1-based index, logical

			for (int i = cur_first; i <= cur_last; ++i)
			{
				auto cnct = curObj->connectivity(i - 1);
				for (int j = 0; j < cnct.x; ++j)
				{
					int node_idx = cnct.n[j];

					int l = cnct.leftAdj(j); // 1-based index, logical
					int r = cnct.rightAdj(j); // 1-based index, logical

					dst[node_idx - 1].push_back(l);
					if (cnct.x > 2)
						dst[node_idx - 1].push_back(r);
				}
			}
		}
	}

	// Remove duplication
	for (size_t i = 0; i < dst.size(); ++i)
	{
		//Method 1
		//std::sort(dst[i].begin(), dst[i].end());
		//auto it = std::unique(dst[i].begin(), dst[i].end());
		//dst[i].resize(std::distance(dst[i].begin(), it));

		// Method 2
		std::set<size_t> st(dst[i].begin(), dst[i].end());
		dst[i].assign(st.begin(), st.end());
	}

	return 0;
}

int XF_MSH::computeTopology_nodeDependentFace(std::vector<std::vector<size_t>>& dst) const
{
	// Check output array shape.
	if (dst.size() != numOfNode())
		return -1;

	// Parse related entries.
	for (size_t r = 0; r < m_content.size(); ++r)
	{
		auto curPtr = m_content[r];
		if (curPtr->identity() == XF_SECTION::FACE)
		{
			auto curObj = dynamic_cast<XF_FACE*>(curPtr);

			// Face index, 1-based
			const int cur_first = curObj->first_index();
			const int cur_last = curObj->last_index();

			for (int i = cur_first; i <= cur_last; ++i)
			{
				auto cnct = curObj->connectivity(i - 1);
				for (int j = 0; j < cnct.x; ++j)
				{
					int node_idx = cnct.n[j]; // 1-based
					dst[node_idx - 1].push_back(i);
				}
			}
		}
	}

	return 0;
}

int XF_MSH::computeTopology_nodeDependentCell(std::vector<std::vector<size_t>>& dst) const
{
	// Check output array shape.
	if (dst.size() != numOfNode())
		return -1;

	// Parse related entries.
	for (size_t r = 0; r < m_content.size(); ++r)
	{
		auto curPtr = m_content[r];
		if (curPtr->identity() == XF_SECTION::FACE)
		{
			auto curObj = dynamic_cast<XF_FACE*>(curPtr);

			// Face index, 1-based
			const int cur_first = curObj->first_index();
			const int cur_last = curObj->last_index();

			for (int i = cur_first; i <= cur_last; ++i)
			{
				auto cnct = curObj->connectivity(i - 1);
				for (int j = 0; j < cnct.x; ++j)
				{
					int node_idx = cnct.n[j]; // 1-based
					dst[node_idx - 1].push_back(cnct.c[0]);
					dst[node_idx - 1].push_back(cnct.c[1]);
				}
			}
		}
	}

	// Remove duplication
	for (size_t i = 0; i < dst.size(); ++i)
	{
		std::set<size_t> st(dst[i].begin(), dst[i].end());
		dst[i].assign(st.begin(), st.end());
	}

	return 0;
}
