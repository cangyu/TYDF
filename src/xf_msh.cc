#include "../inc/xf_msh.h"
#include <typeinfo>

void XF_COMMENT::repr(std::ostream &out)
{
	out << "(" << std::dec << m_identity << " \"" << m_info << "\")" << std::endl;
}

void XF_HEADER::repr(std::ostream & out)
{
	out << "(" << std::dec << m_identity << " \"" << m_msg << "\")" << std::endl;
}

XF_DIMENSION::XF_DIMENSION(int dim) :
	XF_ENTRY(XF_SECTION::DIMENSION)
{
	if (dim == 2)
		m_is3D = false;
	else if (dim == 3)
		m_is3D = true;
	else
		throw("Invalid dimension!");
}

void XF_DIMENSION::repr(std::ostream & out)
{
	out << "(" << std::dec << m_identity << " " << ND() << ")" << std::endl;
}

XF_NODE::XF_NODE(int zone, int first, int last, int type, int ND) :
	XF_MAIN_RECORD(XF_SECTION::NODE, zone, first, last)
{
	if (type == 0)
		m_type = XF_NODE_TYPE::VIRTUAL;
	else if (type == 1)
		m_type = XF_NODE_TYPE::ANY;
	else if (type == 2)
		m_type = XF_NODE_TYPE::BOUNDARY;
	else
		throw("Invalid description of node type!");

	if (ND == 2)
		m_is3D = false;
	else if (ND == 3)
		m_is3D = true;
	else
		throw("Invalid specification of node dimension!");

	m_pnt.resize(ND * num());
	std::fill(m_pnt.begin(), m_pnt.end(), 0.0);
}

void XF_NODE::record_pnt_coordinate(size_t loc_idx, double x0, double x1, double x2)
{
	const size_t stx = loc_idx * 3;

	m_pnt[stx] = x0;
	m_pnt[stx + 1] = x1;
	m_pnt[stx + 2] = x2;
}

void XF_NODE::record_pnt_coordinate(size_t loc_idx, double x0, double x1)
{
	const size_t stx = loc_idx * 2;

	m_pnt[stx] = x0;
	m_pnt[stx + 1] = x1;
}

void XF_NODE::repr(std::ostream & out)
{
	const int n_dim = ND();
	const int N = num();

	out << "(" << std::dec << m_identity;
	out << " (" << std::hex << zone() << " " << first_index() << " " << last_index() << " ";
	out << std::dec << type() << " " << n_dim << ")(" << std::endl;

	size_t loc_idx = 0;
	out.precision(12);
	for (int i = 0; i < N; ++i)
	{
		for (int k = 0; k < n_dim; ++k)
			out << " " << m_pnt[loc_idx + k];
		out << std::endl;
		loc_idx += n_dim;
	}

	out << "))" << std::endl;
}

XF_CELL::XF_CELL(int zone, int first, int last, int type, int elem_type) :
	XF_MAIN_RECORD(XF_SECTION::CELL, zone, first, last)
{
	if (type == 1)
		m_type = XF_CELL_TYPE::FLUID;
	else if (type == 17)
		m_type = XF_CELL_TYPE::SOLID;
	else if (type == 0)
		m_type = XF_CELL_TYPE::DEAD;
	else
		throw("Invalid specification of cell type!");

	if (elem_type == 0)
	{
		m_elem = XF_CELL_ELEM_TYPE::MIXED;
		m_mixedElemDesc.resize(num());
		std::fill(m_mixedElemDesc.begin(), m_mixedElemDesc.end(), XF_CELL_ELEM_TYPE::MIXED);
	}
	else if (elem_type == 1)
	{
		m_elem = XF_CELL_ELEM_TYPE::TRIANGULAR;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == 2)
	{
		m_elem = XF_CELL_ELEM_TYPE::TETRAHEDRAL;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == 3)
	{
		m_elem = XF_CELL_ELEM_TYPE::QUADRILATERAL;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == 4)
	{
		m_elem = XF_CELL_ELEM_TYPE::HEXAHEDRAL;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == 5)
	{
		m_elem = XF_CELL_ELEM_TYPE::PYRAMID;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == 6)
	{
		m_elem = XF_CELL_ELEM_TYPE::WEDGE;
		m_mixedElemDesc.resize(0);
	}
	else if (elem_type == 7)
	{
		m_elem = XF_CELL_ELEM_TYPE::POLYHEDRAL;
		m_mixedElemDesc.resize(0);
	}
	else
		throw("Invalid specification of cell element type!");
}

void XF_CELL::repr(std::ostream & out)
{
	out << "(" << std::dec << m_identity << " (";
	out << std::hex;
	out << m_zone << " " << m_first << " " << m_last << " ";
	out << m_type << " " << m_elem << ")";

	if (m_elem != XF_CELL_ELEM_TYPE::MIXED)
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
		m_bc = XF_BC_TYPE::INTERIOR;
		break;
	case 3:
		m_bc = XF_BC_TYPE::WALL;
		break;
	case 4:
		m_bc = XF_BC_TYPE::PRESSURE_INLET;
		break;
	case 5:
		m_bc = XF_BC_TYPE::PRESSURE_OUTLET;
		break;
	case 7:
		m_bc = XF_BC_TYPE::SYMMETRY;
		break;
	case 8:
		m_bc = XF_BC_TYPE::PERIODIC_SHADOW;
		break;
	case 9:
		m_bc = XF_BC_TYPE::PRESSURE_FAR_FIELD;
		break;
	case 10:
		m_bc = XF_BC_TYPE::VELOCITY_INLET;
		break;
	case 12:
		m_bc = XF_BC_TYPE::PERIODIC;
		break;
	case 14:
		m_bc = XF_BC_TYPE::FAN;
		break;
	case 20:
		m_bc = XF_BC_TYPE::MASS_FLOW_INLET;
		break;
	case 24:
		m_bc = XF_BC_TYPE::INTERFACE;
		break;
	case 31:
		m_bc = XF_BC_TYPE::PARENT;
		break;
	case 36:
		m_bc = XF_BC_TYPE::OUTFLOW;
		break;
	case 37:
		m_bc = XF_BC_TYPE::AXIS;
		break;
	default:
		throw("Invalid specification of B.C. type!");
	}

	switch (face)
	{
	case 0:
		m_face = XF_FACE_TYPE::F_MIXED;
		break;
	case 2:
		m_face = XF_FACE_TYPE::F_LINEAR;
		break;
	case 3:
		m_face = XF_FACE_TYPE::F_TRIANGULAR;
		break;
	case 4:
		m_face = XF_FACE_TYPE::F_QUADRILATERAL;
		break;
	case 5:
		m_face = XF_FACE_TYPE::F_POLYGONAL;
		throw("Currently not supported!");
	default:
		throw("Invalid face-type!");
	}

	m_connectivity.resize(num());
}

void XF_FACE::repr(std::ostream & out)
{
	out << "(" << std::dec << m_identity << " (";
	out << std::hex;
	out << zone() << " " << first_index() << " " << last_index() << " ";
	out << bc_type() << " " << face_type() << ")(" << std::endl;

	const int N = num();
	if (m_face == XF_FACE_TYPE::F_MIXED)
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
	out << "(" << m_identity << " (" << m_zoneID << " " << m_zoneType << " " << m_zoneName << ")())" << std::endl;
}

int XF_MSH::readFromFile(const std::string & src)
{
	std::ifstream fin(src);
	if (!fin)
		return -1;

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
				if(nd == 3)
				{
					double x, y, z;
					for (int i = first; i <= last; ++i)
					{
						size_t i_loc = i - first;
						fin >> x >> y >> z;
						e->record_pnt_coordinate(i_loc, x, y, z);
					}
					m_is3D = true;
				}
				else
				{
					double x, y;
					for (int i = first; i <= last; ++i)
					{
						size_t i_loc = i - first;
						fin >> x >> y;
						e->record_pnt_coordinate(i_loc, x, y);
					}
					m_is3D = false;
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
							e->elem(i_loc) = XF_CELL_ELEM_TYPE::TRIANGULAR;
							break;
						case 2:
							e->elem(i_loc) = XF_CELL_ELEM_TYPE::TETRAHEDRAL;
							break;
						case 3:
							e->elem(i_loc) = XF_CELL_ELEM_TYPE::QUADRILATERAL;
							break;
						case 4:
							e->elem(i_loc) = XF_CELL_ELEM_TYPE::HEXAHEDRAL;
							break;
						case 5:
							e->elem(i_loc) = XF_CELL_ELEM_TYPE::PYRAMID;
							break;
						case 6:
							e->elem(i_loc) = XF_CELL_ELEM_TYPE::WEDGE;
							break;
						case 7:
							e->elem(i_loc) = XF_CELL_ELEM_TYPE::POLYHEDRAL;
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

int XF_MSH::writeToFile(const std::string & dst)
{
	std::ofstream fout(dst);
	const size_t N = m_content.size();

	// Dimension specification
	size_t i = 0;
	for(; i < N; ++i)
	{
		m_content[i]->repr(fout);

		if(dynamic_cast<XF_DIMENSION*>(m_content[i]))
		{
			// Node declaration
			fout << "(" << std::dec << XF_SECTION::NODE << " (";
			fout << std::hex << 0 << " " << 1 << " " << m_totalNodeNum << " ";
			fout << std::dec << 0 << " " <<  (m_is3D ? 3 : 2) << "))" << std::endl;

			// Cell declaration 
			fout << "(" << std::dec << XF_SECTION::CELL << " (";
			fout << std::hex << 0 << " " << 1 << " " << m_totalCellNum << " ";
			fout << std::dec << 0 << " " << 0 << "))" << std::endl;

			// Face declaration
			fout << "(" << std::dec << XF_SECTION::FACE << " (";
			fout << std::hex << 0 << " " << 1 << " " << m_totalFaceNum << " ";
			fout << std::dec << 0 << " " << 0 << "))" << std::endl;
			
			break;
		}
	}

	for (i+=1; i < N; ++i)
	{
		if(typeid(*m_content[i]) == typeid(XF_DIMENSION))
			continue;

		m_content[i]->repr(fout);
	}

	fout.close();
	return 0;
}
