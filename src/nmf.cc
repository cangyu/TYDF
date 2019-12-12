#include <iostream>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include "../inc/nmf.h"

static bool isWhite(char c)
{
	return c == '\n' || c == ' ' || c == '\t';
}

static bool isBlankLine(const std::string &s)
{
	for (const auto &e : s)
		if (!isWhite(e))
			return false;
	return true;
}

static bool checkStarting(const std::string &s, char c)
{
	for (const auto &e : s)
	{
		if (isWhite(e))
			continue;
		else
			return e == c;
	}
	return false;
}

static void distribute_index(size_t s, size_t e, std::vector<size_t> &dst)
{
	if (s > e)
	{
		const size_t n = s - e + 1;
		dst.resize(n);
		size_t val = s;

		for (size_t i = 0; i < n; ++i)
			dst[i] = val--;
	}
	else
	{
		const size_t n = e - s + 1;
		dst.resize(n);
		size_t val = s;
		for (size_t i = 0; i < n; ++i)
			dst[i] = val++;
	}
}

static void str_formalize(std::string &s)
{
	std::transform(s.begin(), s.end(), s.begin(), ::toupper);
	for (auto &e : s)
		if (e == '-')
			e = '_';
}

namespace GridTool
{
	namespace NMF
	{
		bool BC::isValidBCIdx(int x)
		{
			static const std::set<int> candidate_idx_set = BC_ENUM;
			return candidate_idx_set.find(x) != candidate_idx_set.end();
		}

		bool BC::isValidBCStr(const std::string &x)
		{
			static const std::set<std::string> candidate_str_set = BC_STR;
			std::string x_(x);
			str_formalize(x_);
			return candidate_str_set.find(x_) != candidate_str_set.end();
		}

		const std::string &BC::idx2str(int x)
		{
			static const std::vector<std::string> candidate_str_set = BC_STR;
			static const std::vector<int> candidate_idx_set = BC_ENUM;

			if (!isValidBCIdx(x))
				throw std::runtime_error("Not found!");
			else
			{
				auto it1 = candidate_idx_set.begin();
				auto it2 = candidate_str_set.begin();
				while (*it1 != x)
				{
					++it1;
					++it2;
				}
				return *it2;
			}
		}

		int BC::str2idx(const std::string &x)
		{
			static const std::vector<int> candidate_idx_set = BC_ENUM;
			static const std::vector<std::string> candidate_str_set = BC_STR;

			std::string x_(x);
			str_formalize(x_);

			if (!isValidBCStr(x_))
				throw std::runtime_error("Not found!");
			else
			{
				auto it1 = candidate_idx_set.begin();
				auto it2 = candidate_str_set.begin();
				while (x_ != *it2)
				{
					++it1;
					++it2;
				}
				return *it1;
			}
		}

		BLOCK::BLOCK(size_t nI, size_t nJ) :
			m_idx(0),
			m_name(""),
			m_dim{ nI, nJ, 1 }
		{
			if (nI < 2 || nJ < 2)
				throw std::invalid_argument("Invalid dimension.");
		}
		BLOCK::BLOCK(size_t nI, size_t nJ, size_t nK) :
			m_idx(0),
			m_name(""),
			m_dim{ nI, nJ, nK }
		{
			if (nI < 2 || nJ < 2 || nK < 2)
				throw std::invalid_argument("Invalid dimension.");
		}

		void Block3D::setup_dependence()
		{
			for (short i = 0; i < NumOfVertex; ++i)
			{
				auto &v = m_vertex[i];
				v.local_index = i + 1;
				v.dependentBlock = this;
			}
			for (short i = 0; i < NumOfFrame; ++i)
			{
				auto &e = m_frame[i];
				e.local_index = i + 1;
				e.dependentBlock = this;
			}
			for (short i = 0; i < NumOfSurf; ++i)
			{
				auto &s = m_surf[i];
				s.local_index = i + 1;
				s.dependentBlock = this;
			}
		}

		void Block3D::establish_connections()
		{
			// Connection between frame edges and surrounding surfaces.
			// The order matters a lot!
			m_surf[0].includedFrame = { &m_frame[4], &m_frame[11], &m_frame[7], &m_frame[8] };
			m_surf[1].includedFrame = { &m_frame[5], &m_frame[10], &m_frame[6], &m_frame[9] };
			m_surf[2].includedFrame = { &m_frame[8], &m_frame[3], &m_frame[9], &m_frame[0] };
			m_surf[3].includedFrame = { &m_frame[11], &m_frame[2], &m_frame[10], &m_frame[1] };
			m_surf[4].includedFrame = { &m_frame[0], &m_frame[5], &m_frame[1], &m_frame[4] };
			m_surf[5].includedFrame = { &m_frame[3], &m_frame[6], &m_frame[2], &m_frame[7] };

			m_frame[0].dependentSurf = { &m_surf[2], &m_surf[4] };
			m_frame[1].dependentSurf = { &m_surf[4], &m_surf[3] };
			m_frame[2].dependentSurf = { &m_surf[3], &m_surf[5] };
			m_frame[3].dependentSurf = { &m_surf[5], &m_surf[2] };
			m_frame[4].dependentSurf = { &m_surf[0], &m_surf[4] };
			m_frame[5].dependentSurf = { &m_surf[4], &m_surf[1] };
			m_frame[6].dependentSurf = { &m_surf[1], &m_surf[5] };
			m_frame[7].dependentSurf = { &m_surf[5], &m_surf[0] };
			m_frame[8].dependentSurf = { &m_surf[0], &m_surf[2] };
			m_frame[9].dependentSurf = { &m_surf[2], &m_surf[1] };
			m_frame[10].dependentSurf = { &m_surf[1], &m_surf[3] };
			m_frame[11].dependentSurf = { &m_surf[3], &m_surf[0] };

			// Connection between surrounding surfaces and block vertexes.
			// The order matters a lot!
			m_surf[0].includedVertex = { &m_vertex[0], &m_vertex[3], &m_vertex[7], &m_vertex[4] };
			m_surf[1].includedVertex = { &m_vertex[1], &m_vertex[2], &m_vertex[6], &m_vertex[5] };
			m_surf[2].includedVertex = { &m_vertex[0], &m_vertex[4], &m_vertex[5], &m_vertex[1] };
			m_surf[3].includedVertex = { &m_vertex[3], &m_vertex[7], &m_vertex[6], &m_vertex[2] };
			m_surf[4].includedVertex = { &m_vertex[0], &m_vertex[1], &m_vertex[2], &m_vertex[3] };
			m_surf[5].includedVertex = { &m_vertex[4], &m_vertex[5], &m_vertex[6], &m_vertex[7] };

			m_vertex[0].dependentSurf = { &m_surf[4], &m_surf[0], &m_surf[2] };
			m_vertex[1].dependentSurf = { &m_surf[4], &m_surf[2], &m_surf[1] };
			m_vertex[2].dependentSurf = { &m_surf[4], &m_surf[1], &m_surf[3] };
			m_vertex[3].dependentSurf = { &m_surf[4], &m_surf[3], &m_surf[0] };
			m_vertex[4].dependentSurf = { &m_surf[5], &m_surf[0], &m_surf[2] };
			m_vertex[5].dependentSurf = { &m_surf[5], &m_surf[2], &m_surf[1] };
			m_vertex[6].dependentSurf = { &m_surf[5], &m_surf[1], &m_surf[3] };
			m_vertex[7].dependentSurf = { &m_surf[5], &m_surf[3], &m_surf[0] };

			// Connection between frame edges and block vertexes.
			m_frame[0].includedVertex = { &m_vertex[0], &m_vertex[1] };
			m_frame[1].includedVertex = { &m_vertex[3], &m_vertex[2] };
			m_frame[2].includedVertex = { &m_vertex[7], &m_vertex[6] };
			m_frame[3].includedVertex = { &m_vertex[4], &m_vertex[5] };
			m_frame[4].includedVertex = { &m_vertex[0], &m_vertex[3] };
			m_frame[5].includedVertex = { &m_vertex[1], &m_vertex[2] };
			m_frame[6].includedVertex = { &m_vertex[5], &m_vertex[6] };
			m_frame[7].includedVertex = { &m_vertex[4], &m_vertex[7] };
			m_frame[8].includedVertex = { &m_vertex[0], &m_vertex[4] };
			m_frame[9].includedVertex = { &m_vertex[1], &m_vertex[5] };
			m_frame[10].includedVertex = { &m_vertex[2], &m_vertex[6] };
			m_frame[11].includedVertex = { &m_vertex[3], &m_vertex[7] };

			m_vertex[0].dependentFrame = { &m_frame[8], &m_frame[4], &m_frame[0] };
			m_vertex[1].dependentFrame = { &m_frame[9], &m_frame[0], &m_frame[5] };
			m_vertex[2].dependentFrame = { &m_frame[10], &m_frame[5], &m_frame[1] };
			m_vertex[3].dependentFrame = { &m_frame[11], &m_frame[1], &m_frame[4] };
			m_vertex[4].dependentFrame = { &m_frame[8], &m_frame[7], &m_frame[3] };
			m_vertex[5].dependentFrame = { &m_frame[9], &m_frame[3], &m_frame[6] };
			m_vertex[6].dependentFrame = { &m_frame[10], &m_frame[6], &m_frame[2] };
			m_vertex[7].dependentFrame = { &m_frame[11], &m_frame[2], &m_frame[7] };
		}

		size_t Block3D::surface_internal_node_num(short s) const
		{
			switch (s)
			{
			case 1:
			case 2:
				return (IDIM() - 2) * (JDIM() - 2);
			case 3:
			case 4:
				return (JDIM() - 2) * (KDIM() - 2);
			case 5:
			case 6:
				return (KDIM() - 2) * (IDIM() - 2);
			default:
				throw not_a_surface(s);
			}
		}

		size_t Block3D::frame_node_num(short idx) const
		{
			switch (idx - 1)
			{
			case 4:
			case 5:
			case 6:
			case 7:
				return IDIM();
			case 8:
			case 9:
			case 10:
			case 11:
				return JDIM();
			case 0:
			case 1:
			case 2:
			case 3:
				return KDIM();
			default:
				throw not_a_frame(idx);
			}
		}

		size_t Block3D::surface_node_num(short idx) const
		{
			switch (idx - 1)
			{
			case 0:
			case 1:
				return IDIM() * JDIM();
			case 2:
			case 3:
				return JDIM() * KDIM();
			case 4:
			case 5:
				return KDIM() * IDIM();
			default:
				throw not_a_surface(idx);
			}
		}

		size_t Block3D::surface_face_num(short idx) const
		{
			switch (idx - 1)
			{
			case 0:
			case 1:
				return (IDIM() - 1) * (JDIM() - 1);
			case 2:
			case 3:
				return (JDIM() - 1) * (KDIM() - 1);
			case 4:
			case 5:
				return (KDIM() - 1) * (IDIM() - 1);
			default:
				throw not_a_surface(idx);
			}
		}

		size_t &Block3D::surface_face_index(short f, size_t pri, size_t sec)
		{
			HEX_CELL *p = nullptr;

			switch (f)
			{
			case 1:
				p = &cell(pri, sec, 1);
				break;
			case 2:
				p = &cell(pri, sec, KDIM() - 1);
				break;
			case 3:
				p = &cell(1, pri, sec);
				break;
			case 4:
				p = &cell(IDIM() - 1, pri, sec);
				break;
			case 5:
				p = &cell(sec, 1, pri);
				break;
			case 6:
				p = &cell(sec, JDIM() - 1, pri);
				break;
			default:
				throw not_a_surface(f);
			}
			return p->FaceSeq(f);
		}

		size_t &Block3D::vertex_node_index(short v)
		{
			HEX_CELL *p = nullptr;
			switch (v)

			{
			case 1:
				p = &cell(1, 1, 1);
				break;
			case 2:
				p = &cell(1, 1, KDIM() - 1);
				break;
			case 3:
				p = &cell(IDIM() - 1, 1, KDIM() - 1);
				break;
			case 4:
				p = &cell(IDIM() - 1, 1, 1);
				break;
			case 5:
				p = &cell(1, JDIM() - 1, 1);
				break;
			case 6:
				p = &cell(1, JDIM() - 1, KDIM() - 1);
				break;
			case 7:
				p = &cell(IDIM() - 1, JDIM() - 1, KDIM() - 1);
				break;
			case 8:
				p = &cell(IDIM() - 1, JDIM() - 1, 1);
				break;
			default:
				throw not_a_vertex(v);
			}
			return p->NodeSeq(v);
		}

		void Block3D::surface_node_coordinate(short f, size_t pri_seq, size_t sec_seq, size_t &i, size_t &j, size_t &k)
		{
			switch (f)
			{
			case 1:
				k = 1;
				i = pri_seq;
				j = sec_seq;
				break;
			case 2:
				k = KDIM();
				i = pri_seq;
				j = sec_seq;
				break;
			case 3:
				i = 1;
				j = pri_seq;
				k = sec_seq;
				break;
			case 4:
				i = IDIM();
				j = pri_seq;
				k = sec_seq;
				break;
			case 5:
				j = 1;
				k = pri_seq;
				i = sec_seq;
				break;
			case 6:
				j = JDIM();
				k = pri_seq;
				i = sec_seq;
				break;
			default:
				throw not_a_surface(f);
			}
		}

		void Block3D::surface_internal_node_occurance(short f, size_t pri, size_t sec, std::vector<size_t*> &oc)
		{
			size_t i = 0, j = 0, k = 0;

			surface_node_coordinate(f, pri, sec, i, j, k);
			switch (f)
			{
			case 1:
				oc[0] = &cell(i, j, k).NodeSeq(1);
				oc[1] = &cell(i - 1, j, k).NodeSeq(4);
				oc[2] = &cell(i, j - 1, k).NodeSeq(5);
				oc[3] = &cell(i - 1, j - 1, k).NodeSeq(8);
				break;
			case 2:
				oc[0] = &cell(i, j, k - 1).NodeSeq(2);
				oc[1] = &cell(i - 1, j, k - 1).NodeSeq(3);
				oc[2] = &cell(i, j - 1, k - 1).NodeSeq(6);
				oc[3] = &cell(i - 1, j - 1, k - 1).NodeSeq(7);
				break;
			case 3:
				oc[0] = &cell(i, j, k).NodeSeq(1);
				oc[1] = &cell(i, j - 1, k).NodeSeq(5);
				oc[2] = &cell(i, j, k - 1).NodeSeq(2);
				oc[3] = &cell(i, j - 1, k - 1).NodeSeq(6);
				break;
			case 4:
				oc[0] = &cell(i - 1, j, k).NodeSeq(4);
				oc[1] = &cell(i - 1, j - 1, k).NodeSeq(8);
				oc[2] = &cell(i - 1, j, k - 1).NodeSeq(3);
				oc[3] = &cell(i - 1, j - 1, k - 1).NodeSeq(7);
				break;
			case 5:
				oc[0] = &cell(i, j, k).NodeSeq(1);
				oc[1] = &cell(i - 1, j, k).NodeSeq(4);
				oc[2] = &cell(i, j, k - 1).NodeSeq(2);
				oc[3] = &cell(i - 1, j, k - 1).NodeSeq(3);
				break;
			case 6:
				oc[0] = &cell(i, j - 1, k).NodeSeq(5);
				oc[1] = &cell(i - 1, j - 1, k).NodeSeq(8);
				oc[2] = &cell(i, j - 1, k - 1).NodeSeq(6);
				oc[3] = &cell(i - 1, j - 1, k - 1).NodeSeq(7);
				break;
			default:
				throw not_a_surface(f);
			}
		}

		void Block3D::frame_internal_node_occurace(short f, size_t idx, std::vector<size_t*> &oc)
		{
			switch (f - 1)
			{
			case 0:
				oc[0] = &cell(1, 1, idx).NodeSeq(1);
				oc[1] = &cell(1, 1, idx - 1).NodeSeq(2);
				break;
			case 1:
				oc[0] = &cell(IDIM() - 1, 1, idx).NodeSeq(4);
				oc[1] = &cell(IDIM() - 1, 1, idx - 1).NodeSeq(3);
				break;
			case 2:
				oc[0] = &cell(IDIM() - 1, JDIM() - 1, idx).NodeSeq(8);
				oc[1] = &cell(IDIM() - 1, JDIM() - 1, idx - 1).NodeSeq(7);
				break;
			case 3:
				oc[0] = &cell(1, JDIM() - 1, idx).NodeSeq(5);
				oc[1] = &cell(1, JDIM() - 1, idx - 1).NodeSeq(6);
				break;
			case 4:
				oc[0] = &cell(idx, 1, 1).NodeSeq(1);
				oc[1] = &cell(idx - 1, 1, 1).NodeSeq(4);
				break;
			case 5:
				oc[0] = &cell(idx, 1, KDIM() - 1).NodeSeq(2);
				oc[1] = &cell(idx - 1, 1, KDIM() - 1).NodeSeq(3);
				break;
			case 6:
				oc[0] = &cell(idx, JDIM() - 1, KDIM() - 1).NodeSeq(6);
				oc[1] = &cell(idx - 1, JDIM() - 1, KDIM() - 1).NodeSeq(7);
				break;
			case 7:
				oc[0] = &cell(idx, JDIM() - 1, 1).NodeSeq(5);
				oc[1] = &cell(idx - 1, JDIM() - 1, 1).NodeSeq(8);
				break;
			case 8:
				oc[0] = &cell(1, idx, 1).NodeSeq(1);
				oc[1] = &cell(1, idx - 1, 1).NodeSeq(5);
				break;
			case 9:
				oc[0] = &cell(1, idx, KDIM() - 1).NodeSeq(2);
				oc[1] = &cell(1, idx - 1, KDIM() - 1).NodeSeq(6);
				break;
			case 10:
				oc[0] = &cell(IDIM() - 1, idx, KDIM() - 1).NodeSeq(3);
				oc[1] = &cell(IDIM() - 1, idx - 1, KDIM() - 1).NodeSeq(7);
				break;
			case 11:
				oc[0] = &cell(IDIM() - 1, idx, 1).NodeSeq(4);
				oc[1] = &cell(IDIM() - 1, idx - 1, 1).NodeSeq(8);
				break;
			default:
				throw not_a_frame(f);
			}
		}

		void Mapping3D::readFromFile(const std::string &path)
		{
			std::string s;
			std::stringstream ss;

			// Open file
			std::ifstream mfp(path);
			if (mfp.fail())
				throw std::runtime_error("Can not open target input file: \"" + path + "\".");

			// Skip header
			do {
				std::getline(mfp, s);
			} while (isBlankLine(s) || checkStarting(s, '#'));

			// Read block nums
			static const std::regex pattern1(R"(\s*(\d+)\s*)");
			std::smatch res1;
			if (std::regex_match(s, res1, pattern1))
			{
				int NumOfBlk = std::stoi(res1[1].str());
				if (NumOfBlk > 0)
				{
					release_all(); // NOT release all existing resources until it is ensured that this input file is valid.
					m_blk.resize(NumOfBlk, nullptr); // Re-Allocate storage for new recordings.
				}
				else
					throw std::runtime_error("Invalid num of blocks: \"" + res1[1].str() + "\".");
			}
			else
				throw std::runtime_error("Failed to match the single line, where only the num of blocks is specified.");

			// Read dimension info of each block
			static const std::regex pattern2(R"(\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*)");
			const auto NumOfBlk = nBlock();
			for (size_t i = 0; i < NumOfBlk; i++)
			{
				std::getline(mfp, s);
				std::smatch res2;
				if (std::regex_match(s, res2, pattern2))
				{
					const size_t idx = std::stoi(res2[1].str());
					if (idx < 1 || idx > NumOfBlk)
						throw std::runtime_error("Invalid order of block: " + std::to_string(idx));

					const int i_max = std::stoi(res2[2].str());
					if (i_max < 1)
						throw std::runtime_error("Invalid I dimension: " + std::to_string(i_max));

					const int j_max = std::stoi(res2[3].str());
					if (j_max < 1)
						throw std::runtime_error("Invalid J dimension: " + std::to_string(j_max));

					const int k_max = std::stoi(res2[4].str());
					if (k_max < 1)
						throw std::runtime_error("Invalid K dimension: " + std::to_string(k_max));

					auto &e = m_blk(idx);
					e = new Block3D(i_max, j_max, k_max);
					e->index() = idx;
				}
				else
					throw std::runtime_error("Failed to match 4 integers.");
			}

			// Skip separators
			while (std::getline(mfp, s))
			{
				if (!isBlankLine(s) && !checkStarting(s, '#'))
					break;
			}

			// Read connections
			if (!mfp.eof())
			{
				do {
					str_formalize(s);
					ss.clear();
					ss << s;
					std::string bc_str;
					ss >> bc_str;
					if (BC::str2idx(bc_str) == BC::ONE_TO_ONE)
					{
						size_t cB[2];
						short cF[2];
						size_t cS1[2], cE1[2], cS2[2], cE2[2];
						std::string swp;
						for (int i = 0; i < 2; i++)
							ss >> cB[i] >> cF[i] >> cS1[i] >> cE1[i] >> cS2[i] >> cE2[i];
						ss >> swp;
						add_entry(bc_str, cB[0], cF[0], cS1[0], cE1[0], cS2[0], cE2[0], cB[1], cF[1], cS1[1], cE1[1], cS2[1], cE2[1], swp == "TRUE");
					}
					else
					{
						size_t cB;
						short cF;
						size_t cS1, cE1, cS2, cE2;
						ss >> cB >> cF >> cS1 >> cE1 >> cS2 >> cE2;
						add_entry(bc_str, cB, cF, cS1, cE1, cS2, cE2);
					}
				} while (std::getline(mfp, s));
			}

			// Close input file
			mfp.close();
		}

		void Mapping3D::compute_topology()
		{
			connecting();

			const int nsf = coloring_surface();
			if (nsf < Block3D::NumOfSurf)
				throw std::runtime_error("Internal error occured when counting surfaces.");
			m_surf.resize(nsf);
			for (auto &e : m_surf)
				e.clear();
			for (auto b : m_blk)
				for (short i = 1; i <= Block3D::NumOfSurf; ++i)
				{
					auto &s = b->surf(i);
					m_surf(s.global_index).push_back(&s);
				}

			const int nfm = coloring_frame();
			if (nfm < Block3D::NumOfFrame)
				throw std::runtime_error("Internal error occured when counting frames.");
			m_frame.resize(nfm);
			for (auto &e : m_frame)
				e.clear();
			for (auto b : m_blk)
				for (short i = 1; i <= Block3D::NumOfFrame; ++i)
				{
					auto &e = b->frame(i);
					m_frame(e.global_index).push_back(&e);
				}

			const int nvt = coloring_vertex();
			if (nvt < Block3D::NumOfVertex)
				throw std::runtime_error("Internal error occured when counting vertexes.");
			m_vertex.resize(nvt);
			for (auto &e : m_vertex)
				e.clear();
			for (auto b : m_blk)
				for (short i = 1; i <= Block3D::NumOfVertex; ++i)
				{
					auto &v = b->vertex(i);
					m_vertex(v.global_index).push_back(&v);
				}
		}

		void Mapping3D::summary(std::ostream &out = std::cout)
		{
			out << "======================================== SUMMARY =======================================" << std::endl;
			out << "Num of blocks: " << nBlock() << std::endl;
			size_t nSa = 0, nSi = 0, nSb = 0;
			nSurface(nSa, nSi, nSb);
			out << "Num of surfaces: " << nSa << ", among which " << nSi << " are internal, " << nSb << " at boundary" << std::endl;
			out << "Num of frames: " << nFrame() << std::endl;
			out << "Num of vertexs: " << nVertex() << std::endl;
			out << "Num of HEX cells: " << nCell() << std::endl;
			size_t nFa = 0, nFi = 0, nFb = 0;
			nFace(nFa, nFi, nFb);
			out << "Num of QUAD faces: " << nFa << ", among which " << nFi << " are internal, " << nFb << " at boundary" << std::endl;
			out << "Num of nodes: " << nNode() << " (duplication removed)" << std::endl;
			out << "------------------------------------------Blocks----------------------------------------";
			for (auto b : m_blk)
			{
				out << "\nIndex: " << b->index() << std::endl;
				out << "I=" << b->IDIM() << ", J=" << b->JDIM() << ", K=" << b->KDIM() << std::endl;
				out << "Num of HEX cells: " << b->cell_num() << std::endl;
				out << "Num of QUAD faces: " << b->face_num() << std::endl;
				out << "Num of nodes: " << b->node_num() << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Local Vertex Index:" << std::resetiosflags(std::ios::left);
				for (short i = 1; i <= Block3D::NumOfVertex; ++i)
					out << std::setw(5) << b->vertex(i).local_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Global Vertex Index:" << std::resetiosflags(std::ios::left);
				for (short i = 1; i <= Block3D::NumOfVertex; ++i)
					out << std::setw(5) << b->vertex(i).global_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Local Frame Index:" << std::resetiosflags(std::ios::left);
				for (short i = 1; i <= Block3D::NumOfFrame; ++i)
					out << std::setw(5) << b->frame(i).local_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Global Frame Index:" << std::resetiosflags(std::ios::left);
				for (short i = 1; i <= Block3D::NumOfFrame; ++i)
					out << std::setw(5) << b->frame(i).global_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Local Surface Index:" << std::resetiosflags(std::ios::left);
				for (short i = 1; i <= Block3D::NumOfSurf; ++i)
					out << std::setw(5) << b->surf(i).local_index;
				out << std::endl;

				out << std::setiosflags(std::ios::left) << std::setw(24) << "Global Surface Index:" << std::resetiosflags(std::ios::left);
				for (short i = 1; i <= Block3D::NumOfSurf; ++i)
					out << std::setw(5) << b->surf(i).global_index;
				out << std::endl;
			}
			out << "-----------------------------------------Surfaces---------------------------------------";
			for (size_t i = 1; i <= nSa; ++i)
			{
				out << "\nGlobal Index: " << i << std::endl;
				auto sf_rep = m_surf(i)[0];
				out << "Num of faces: " << sf_rep->dependentBlock->surface_face_num(sf_rep->local_index) << std::endl;
				out << "Num of nodes: " << sf_rep->dependentBlock->surface_node_num(sf_rep->local_index) << std::endl;
				out << "Occurance:";
				for (auto e : m_surf(i))
					out << " (" << e->dependentBlock->index() << ", " << e->local_index << ")";
				out << std::endl;
			}
			out << "------------------------------------------Frames----------------------------------------";
			for (size_t i = 1; i <= nFrame(); ++i)
			{
				out << "\nGlobal Index: " << i << std::endl;
				auto f_rep = m_frame(i)[0];
				out << "Num of nodes: " << f_rep->dependentBlock->frame_node_num(f_rep->local_index) << std::endl;
				out << "Occurance:";
				for (auto e : m_frame(i))
					out << " (" << e->dependentBlock->index() << ", " << e->local_index << ")";
				out << std::endl;
			}
			out << "-----------------------------------------Vertexes---------------------------------------";
			for (size_t i = 1; i <= nVertex(); ++i)
			{
				out << "\nGlobal Index: " << i << std::endl;
				out << "Occurance:";
				for (auto e : m_vertex(i))
					out << " (" << e->dependentBlock->index() << ", " << e->local_index << ")";
				out << std::endl;
			}
			out << "========================================== END =========================================" << std::endl;
		}

		void Mapping3D::numbering()
		{
			for (auto b : m_blk)
				b->allocate_cell_storage();

			numbering_cell();
			numbering_face();
			numbering_node();
		}

		void Mapping3D::writeToFile(const std::string &path)
		{
			// Open target file
			std::ofstream f_out(path);
			if (f_out.fail())
				throw std::runtime_error("Can not open target output file: \"" + path + "\".");

			// Header
			f_out << "# ============================= Neutral Map File generated by the Grid-Glue software =========================" << std::endl;
			f_out << "# ============================================================================================================" << std::endl;
			f_out << "# Block#    IDIM    JDIM    KDIM" << std::endl;
			f_out << "# ------------------------------------------------------------------------------------------------------------" << std::endl;

			// Block info
			const size_t NumOfBlk = nBlock();
			f_out << std::setw(8) << std::right << NumOfBlk << std::endl;
			for (size_t i = 0; i < NumOfBlk; i++)
			{
				f_out << std::setw(8) << std::right << i + 1;
				f_out << std::setw(8) << std::right << m_blk[i]->IDIM();
				f_out << std::setw(8) << std::right << m_blk[i]->JDIM();
				f_out << std::setw(8) << std::right << m_blk[i]->KDIM();
				f_out << std::endl;
			}

			// Interface and boundary info
			f_out << "# ============================================================================================================" << std::endl;
			f_out << "# Type           B1    F1       S1    E1       S2    E2       B2    F2       S1    E1       S2    E2      Swap" << std::endl;
			f_out << "# ------------------------------------------------------------------------------------------------------------" << std::endl;
			for (auto & e : m_entry)
			{
				f_out << std::setw(13) << std::left << BC::idx2str(e->Type());
				f_out << std::setw(6) << std::right << e->Range1().B();
				f_out << std::setw(6) << std::right << e->Range1().F();
				f_out << std::setw(9) << std::right << e->Range1().S1();
				f_out << std::setw(6) << std::right << e->Range1().E1();
				f_out << std::setw(9) << std::right << e->Range1().S2();
				f_out << std::setw(6) << std::right << e->Range1().E2();
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = static_cast<DoubleSideEntry*>(e);
					f_out << std::setw(9) << std::right << p->Range2().B();
					f_out << std::setw(6) << std::right << p->Range2().F();
					f_out << std::setw(9) << std::right << p->Range2().S1();
					f_out << std::setw(6) << std::right << p->Range2().E1();
					f_out << std::setw(9) << std::right << p->Range2().S2();
					f_out << std::setw(6) << std::right << p->Range2().E2();
					f_out << std::setw(10) << std::right << (p->Swap() ? "TRUE" : "FALSE");
				}
				f_out << std::endl;
			}

			// Close output file
			f_out.close();
		}

		void Mapping3D::connecting()
		{
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = dynamic_cast<DoubleSideEntry*>(e);
					auto B1 = m_blk(p->Range1().B());
					auto B2 = m_blk(p->Range2().B());
					auto F1 = &B1->surf(p->Range1().F());
					auto F2 = &B2->surf(p->Range2().F());

					// Surface connectivity
					F1->neighbourSurf = F2;
					F2->neighbourSurf = F1;

					// Counterpart concerning frames on each surface.
					// There're 4 possible mapping cases.
					if (p->Swap())
					{
						// When the primary directions of F1 and F2 are not aligned,
						// the primary direction of F1 goes parallel with the secondary
						// direction of F2, and the secondary direction of F1 goes
						// parallel with the primary direction of F2. However, under
						// the right-hand convention, there're 2 further possibilities:
						// one is the primary direction of F1 and the secondary direction
						// of F2 not only go parallel, but also run in the same direction,
						// in this case, the remaining pair MUST runs in different direction;
						// the other is the primary direction of F1 and the secondary
						// direction of F2 go parallel, but run in different direction,
						// in this case, the remaining pair MUST run in same direction.
						// The exact case is determined by checking the trend of
						// corresponding ranges.

						if (p->Range1().pri_trend() == p->Range2().sec_trend()) // case 1
						{
							F1->counterpartFrame[0] = F2->includedFrame[1];
							F1->counterpartFrame[1] = F2->includedFrame[2];
							F1->counterpartFrame[2] = F2->includedFrame[3];
							F1->counterpartFrame[3] = F2->includedFrame[0];

							F1->counterpartFrameIsOpposite[0] = false;
							F1->counterpartFrameIsOpposite[1] = true;
							F1->counterpartFrameIsOpposite[2] = false;
							F1->counterpartFrameIsOpposite[3] = true;

							F2->counterpartFrame[0] = F1->includedFrame[3];
							F2->counterpartFrame[1] = F1->includedFrame[0];
							F2->counterpartFrame[2] = F1->includedFrame[1];
							F2->counterpartFrame[3] = F1->includedFrame[2];

							F2->counterpartFrameIsOpposite[0] = true;
							F2->counterpartFrameIsOpposite[1] = false;
							F2->counterpartFrameIsOpposite[2] = true;
							F2->counterpartFrameIsOpposite[3] = false;

							F1->counterpartVertex[0] = F2->includedVertex[1];
							F1->counterpartVertex[1] = F2->includedVertex[2];
							F1->counterpartVertex[2] = F2->includedVertex[3];
							F1->counterpartVertex[3] = F2->includedVertex[0];

							F2->counterpartVertex[0] = F1->includedVertex[3];
							F2->counterpartVertex[1] = F1->includedVertex[0];
							F2->counterpartVertex[2] = F1->includedVertex[1];
							F2->counterpartVertex[3] = F1->includedVertex[2];
						}
						else // case 2
						{
							F1->counterpartFrame[0] = F2->includedFrame[3];
							F1->counterpartFrame[1] = F2->includedFrame[0];
							F1->counterpartFrame[2] = F2->includedFrame[1];
							F1->counterpartFrame[3] = F2->includedFrame[2];

							F1->counterpartFrameIsOpposite[0] = true;
							F1->counterpartFrameIsOpposite[1] = false;
							F1->counterpartFrameIsOpposite[2] = true;
							F1->counterpartFrameIsOpposite[3] = false;

							F2->counterpartFrame[0] = F1->includedFrame[1];
							F2->counterpartFrame[1] = F1->includedFrame[2];
							F2->counterpartFrame[2] = F1->includedFrame[3];
							F2->counterpartFrame[3] = F1->includedFrame[0];

							F2->counterpartFrameIsOpposite[0] = false;
							F2->counterpartFrameIsOpposite[1] = true;
							F2->counterpartFrameIsOpposite[2] = false;
							F2->counterpartFrameIsOpposite[3] = true;

							F1->counterpartVertex[0] = F2->includedVertex[3];
							F1->counterpartVertex[1] = F2->includedVertex[0];
							F1->counterpartVertex[2] = F2->includedVertex[1];
							F1->counterpartVertex[3] = F2->includedVertex[2];

							F2->counterpartVertex[0] = F1->includedVertex[1];
							F2->counterpartVertex[1] = F1->includedVertex[2];
							F2->counterpartVertex[2] = F1->includedVertex[3];
							F2->counterpartVertex[3] = F1->includedVertex[0];
						}
					}
					else
					{
						// Even the primary directions of F1 and F2 are aligned, 
						// they may be in opposite directions. This is further deteced 
						// by the compare the trend from S1 to E1 in F1 with that in F2. 
						// If these 2 trends are the same, the 2 primary directions are 
						// not only parallel but also in the same directions, otherwise, 
						// they are only parallel, but runs in different directions.

						if (p->Range1().pri_trend() != p->Range2().pri_trend()) // Parallel, but goes in opposite direction.
						{
							F1->counterpartFrame[0] = F2->includedFrame[2];
							F1->counterpartFrame[1] = F2->includedFrame[3];
							F1->counterpartFrame[2] = F2->includedFrame[0];
							F1->counterpartFrame[3] = F2->includedFrame[1];

							F1->counterpartFrameIsOpposite[0] = true;
							F1->counterpartFrameIsOpposite[1] = true;
							F1->counterpartFrameIsOpposite[2] = true;
							F1->counterpartFrameIsOpposite[3] = true;

							F2->counterpartFrame[0] = F1->includedFrame[2];
							F2->counterpartFrame[1] = F1->includedFrame[3];
							F2->counterpartFrame[2] = F1->includedFrame[0];
							F2->counterpartFrame[3] = F1->includedFrame[1];

							F2->counterpartFrameIsOpposite[0] = true;
							F2->counterpartFrameIsOpposite[1] = true;
							F2->counterpartFrameIsOpposite[2] = true;
							F2->counterpartFrameIsOpposite[3] = true;

							F1->counterpartVertex[0] = F2->includedVertex[2];
							F1->counterpartVertex[1] = F2->includedVertex[3];
							F1->counterpartVertex[2] = F2->includedVertex[0];
							F1->counterpartVertex[3] = F2->includedVertex[1];

							F2->counterpartVertex[0] = F1->includedVertex[2];
							F2->counterpartVertex[1] = F1->includedVertex[3];
							F2->counterpartVertex[2] = F1->includedVertex[0];
							F2->counterpartVertex[3] = F1->includedVertex[1];
						}
						else // Parallel, and goes in the same direction.
						{
							F1->counterpartFrame[0] = F2->includedFrame[0];
							F1->counterpartFrame[1] = F2->includedFrame[1];
							F1->counterpartFrame[2] = F2->includedFrame[2];
							F1->counterpartFrame[3] = F2->includedFrame[3];

							F1->counterpartFrameIsOpposite[0] = false;
							F1->counterpartFrameIsOpposite[1] = false;
							F1->counterpartFrameIsOpposite[2] = false;
							F1->counterpartFrameIsOpposite[3] = false;

							F2->counterpartFrame[0] = F1->includedFrame[0];
							F2->counterpartFrame[1] = F1->includedFrame[1];
							F2->counterpartFrame[2] = F1->includedFrame[2];
							F2->counterpartFrame[3] = F1->includedFrame[3];

							F2->counterpartFrameIsOpposite[0] = false;
							F2->counterpartFrameIsOpposite[1] = false;
							F2->counterpartFrameIsOpposite[2] = false;
							F2->counterpartFrameIsOpposite[3] = false;

							F1->counterpartVertex[0] = F2->includedVertex[0];
							F1->counterpartVertex[1] = F2->includedVertex[1];
							F1->counterpartVertex[2] = F2->includedVertex[2];
							F1->counterpartVertex[3] = F2->includedVertex[3];

							F2->counterpartVertex[0] = F1->includedVertex[0];
							F2->counterpartVertex[1] = F1->includedVertex[1];
							F2->counterpartVertex[2] = F1->includedVertex[2];
							F2->counterpartVertex[3] = F1->includedVertex[3];
						}
					}
				}
			}
		}

		int Mapping3D::coloring_surface()
		{
			int global_cnt = 0;
			for (auto b : m_blk)
			{
				for (short j = 1; j <= Block3D::NumOfSurf; ++j)
				{
					auto s = &b->surf(j);
					if (s->global_index != 0)
						continue;

					s->global_index = ++global_cnt;
					if (s->neighbourSurf)
						s->neighbourSurf->global_index = s->global_index;
				}
			}
			return global_cnt;
		}

		int Mapping3D::coloring_frame()
		{
			int global_cnt = 0;
			for (auto b : m_blk)
			{
				for (short j = 1; j <= Block3D::NumOfFrame; ++j)
				{
					auto e = &b->frame(j);
					if (e->global_index != 0)
						continue;

					++global_cnt;
					std::queue<Block3D::FRAME*> q;
					q.push(e);

					// BFS
					while (!q.empty())
					{
						auto ce = q.front();
						q.pop();
						ce->global_index = global_cnt;

						for (auto sf : ce->dependentSurf)
						{
							if (!sf)
								throw std::runtime_error("Dependent surface of a frame should NOT be empty.");

							if (sf->neighbourSurf)
							{
								Block3D::FRAME *t = nullptr;
								for (int ii = 0; ii < 4; ++ii)
									if (sf->includedFrame[ii] == ce)
									{
										t = sf->counterpartFrame[ii];
										break;
									}
								if (!t)
									throw std::runtime_error("Internal error.");

								if (t->global_index == 0)
									q.push(t);
							}
						}
					}
				}
			}
			return global_cnt; // The total num of block frames.
		}

		int Mapping3D::coloring_vertex()
		{
			int global_cnt = 0;
			for (auto b : m_blk)
			{
				for (short j = 1; j <= Block3D::NumOfVertex; ++j)
				{
					auto v = &b->vertex(j);
					if (v->global_index != 0)
						continue;

					++global_cnt;
					std::stack<Block3D::VERTEX*> s;
					s.push(v);

					// DFS
					while (!s.empty())
					{
						auto cv = s.top();
						s.pop();
						cv->global_index = global_cnt;

						for (auto sf : cv->dependentSurf)
						{
							if (!sf)
								throw std::runtime_error("Dependent surface of a vertex should NOT be empty.");

							if (sf->neighbourSurf)
							{
								Block3D::VERTEX *t = nullptr;
								for (int ii = 0; ii < 4; ++ii)
									if (sf->includedVertex[ii] == cv)
									{
										t = sf->counterpartVertex[ii];
										break;
									}
								if (!t)
									throw std::runtime_error("Counterpart vertex should exist.");

								if (t->global_index == 0)
									s.push(t);
							}
						}
					}
				}
			}
			return global_cnt;
		}

		void Mapping3D::numbering_cell()
		{
			const auto totalCellNum = nCell();

			size_t cnt = 0;
			for (auto b : m_blk)
				for (size_t k = 1; k < b->KDIM(); ++k)
					for (size_t j = 1; j < b->JDIM(); ++j)
						for (size_t i = 1; i < b->IDIM(); ++i)
							b->cell(i, j, k).CellSeq() = ++cnt;

			if (cnt != totalCellNum)
				throw std::length_error("Inconsistent num of cells.");
		}

		void Mapping3D::numbering_face()
		{
			size_t totalFaceNum = 0, innerFaceNum = 0, bdryFaceNum = 0;
			nFace(totalFaceNum, innerFaceNum, bdryFaceNum);

			size_t cnt = 0;
			for (auto b : m_blk)
			{
				/* Internal faces */
				// K - direction
				for (size_t k = 1; k <= b->KDIM() - 2; ++k)
					for (size_t j = 1; j <= b->JDIM() - 1; ++j)
						for (size_t i = 1; i <= b->IDIM() - 1; ++i)
							b->cell(i, j, k + 1).FaceSeq(1) = b->cell(i, j, k).FaceSeq(2) = ++cnt;

				// I - direction
				for (size_t i = 1; i <= b->IDIM() - 2; ++i)
					for (size_t k = 1; k <= b->KDIM() - 1; ++k)
						for (size_t j = 1; j <= b->JDIM() - 1; ++j)
							b->cell(i + 1, j, k).FaceSeq(3) = b->cell(i, j, k).FaceSeq(4) = ++cnt;

				// J - direction
				for (size_t j = 1; j <= b->JDIM() - 2; ++j)
					for (size_t i = 1; i <= b->IDIM() - 1; ++i)
						for (size_t k = 1; k <= b->KDIM() - 1; ++k)
							b->cell(i, j + 1, k).FaceSeq(5) = b->cell(i, j, k).FaceSeq(6) = ++cnt;

				/* External faces */
				// Single-Sided
				for (short f = 1; f <= Block3D::NumOfSurf; ++f)
				{
					auto &sf = b->surf(f);
					if (!sf.neighbourSurf)
					{
						if (sf.local_index == 1)
						{
							for (size_t j = 1; j <= b->JDIM() - 1; ++j)
								for (size_t i = 1; i <= b->IDIM() - 1; ++i)
									b->cell(i, j, 1).FaceSeq(1) = ++cnt;
						}
						else if (sf.local_index == 2)
						{
							for (size_t j = 1; j <= b->JDIM() - 1; ++j)
								for (size_t i = 1; i <= b->IDIM() - 1; ++i)
									b->cell(i, j, b->KDIM() - 1).FaceSeq(2) = ++cnt;
						}
						else if (sf.local_index == 3)
						{
							for (size_t k = 1; k <= b->KDIM() - 1; ++k)
								for (size_t j = 1; j <= b->JDIM() - 1; ++j)
									b->cell(1, j, k).FaceSeq(3) = ++cnt;
						}
						else if (sf.local_index == 4)
						{
							for (size_t k = 1; k <= b->KDIM() - 1; ++k)
								for (size_t j = 1; j <= b->JDIM() - 1; ++j)
									b->cell(b->IDIM() - 1, j, k).FaceSeq(4) = ++cnt;
						}
						else if (sf.local_index == 5)
						{
							for (size_t i = 1; i <= b->IDIM() - 1; ++i)
								for (size_t k = 1; k <= b->KDIM() - 1; ++k)
									b->cell(i, 1, k).FaceSeq(5) = ++cnt;
						}
						else if (sf.local_index == 6)
						{
							for (size_t i = 1; i <= b->IDIM() - 1; ++i)
								for (size_t k = 1; k <= b->KDIM() - 1; ++k)
									b->cell(i, b->JDIM() - 1, k).FaceSeq(6) = ++cnt;
						}
						else
							throw std::invalid_argument("Internal error: Wrong local index of block surface.");
					}
				}
			}

			// Double-Sided
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = static_cast<DoubleSideEntry*>(e);
					const auto &rg1 = p->Range1();
					const auto &rg2 = p->Range2();
					auto b1 = &block(rg1.B());
					const auto f1 = rg1.F();
					auto b2 = &block(rg2.B());
					const auto f2 = rg2.F();

					std::vector<size_t> b1_dim_pri, b1_dim_sec, b2_dim_pri, b2_dim_sec;
					distribute_index(rg1.S1(), rg1.E1(), b1_dim_pri);
					distribute_index(rg1.S2(), rg1.E2(), b1_dim_sec);
					distribute_index(rg2.S1(), rg2.E1(), b2_dim_pri);
					distribute_index(rg2.S2(), rg2.E2(), b2_dim_sec);
					if (p->Swap())
						std::swap(b2_dim_pri, b2_dim_sec);

					if (b1_dim_pri.size() != b2_dim_pri.size() || b1_dim_sec.size() != b2_dim_sec.size())
						throw std::runtime_error("Inconsistent num of nodes.");

					const auto n1 = rg1.pri_node_num();
					const auto n2 = rg1.sec_node_num();
					if (p->Swap())
					{
						for (size_t l1 = 1; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 1; l2 <= n2 - 1; ++l2)
							{
								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];
								b1->surface_face_index(f1, b1i1, b1i2) = b2->surface_face_index(f2, b2i2, b2i1) = ++cnt;
							}
					}
					else
					{
						for (size_t l1 = 1; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 1; l2 <= n2 - 1; ++l2)
							{
								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];
								b1->surface_face_index(f1, b1i1, b1i2) = b2->surface_face_index(f2, b2i1, b2i2) = ++cnt;
							}
					}
				}
			}

			if (cnt != totalFaceNum)
				throw std::length_error("Inconsistent num of faces detected.");
		}

		void Mapping3D::numbering_node()
		{
			const auto totalNodeNum = nNode();

			size_t cnt = 0;

			// Block interior
			for (auto b : m_blk)
			{
				std::vector<size_t*> boc(8, nullptr);
				for (size_t k = 2; k <= b->KDIM() - 1; ++k)
					for (size_t j = 2; j <= b->JDIM() - 1; ++j)
						for (size_t i = 2; i <= b->IDIM() - 1; ++i)
						{
							++cnt;
							b->interior_node_occurance(i, j, k, boc);
							for (auto r : boc)
								*r = cnt;
						}
			}

			// Vertex
			for (const auto &e : m_vertex)
			{
				++cnt;
				for (auto r : e)
					r->dependentBlock->vertex_node_index(r->local_index) = cnt;
			}

			// Interior of double-sided surface
			for (auto e : m_entry)
			{
				if (e->Type() == BC::ONE_TO_ONE)
				{
					auto p = static_cast<DoubleSideEntry*>(e);
					const auto &rg1 = p->Range1();
					const auto &rg2 = p->Range2();
					auto b1 = &block(rg1.B());
					auto b2 = &block(rg2.B());
					const auto f1 = rg1.F();
					const auto f2 = rg2.F();
					const auto n1 = rg1.pri_node_num();
					const auto n2 = rg1.sec_node_num();

					std::vector<size_t> b1_dim_pri, b1_dim_sec, b2_dim_pri, b2_dim_sec;
					distribute_index(rg1.S1(), rg1.E1(), b1_dim_pri);
					distribute_index(rg1.S2(), rg1.E2(), b1_dim_sec);
					distribute_index(rg2.S1(), rg2.E1(), b2_dim_pri);
					distribute_index(rg2.S2(), rg2.E2(), b2_dim_sec);
					if (p->Swap())
						std::swap(b2_dim_pri, b2_dim_sec);

					if (b1_dim_pri.size() != b2_dim_pri.size() || b1_dim_sec.size() != b2_dim_sec.size())
						throw std::runtime_error("Inconsistent num of nodes.");

					std::vector<size_t*> sioc(4, nullptr);
					if (p->Swap())
					{
						for (size_t l1 = 2; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 2; l2 <= n2 - 1; ++l2)
							{
								++cnt;

								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];

								b1->surface_internal_node_occurance(f1, b1i1, b1i2, sioc);
								for (auto r : sioc)
									*r = cnt;

								b2->surface_internal_node_occurance(f2, b2i2, b2i1, sioc);
								for (auto r : sioc)
									*r = cnt;
							}
					}
					else
					{
						for (size_t l1 = 2; l1 <= n1 - 1; ++l1)
							for (size_t l2 = 2; l2 <= n2 - 1; ++l2)
							{
								++cnt;

								const auto b1i1 = b1_dim_pri[l1 - 1];
								const auto b1i2 = b1_dim_sec[l2 - 1];
								const auto b2i1 = b2_dim_pri[l1 - 1];
								const auto b2i2 = b2_dim_sec[l2 - 1];

								b1->surface_internal_node_occurance(f1, b1i1, b1i2, sioc);
								for (auto r : sioc)
									*r = cnt;

								b2->surface_internal_node_occurance(f2, b2i1, b2i2, sioc);
								for (auto r : sioc)
									*r = cnt;
							}
					}
				}
			}

			// Interior of single-sided surface
			for (auto b : m_blk)
			{
				std::vector<size_t*> sioc(4, nullptr);
				for (short i = 1; i <= Block3D::NumOfSurf; ++i)
				{
					auto &f = b->surf(i);
					if (!f.neighbourSurf)
					{
						size_t pri_end = 0, sec_end = 0;
						switch (i)
						{
						case 1:
						case 2:
							pri_end = b->IDIM() - 1;
							sec_end = b->JDIM() - 1;
							break;
						case 3:
						case 4:
							pri_end = b->JDIM() - 1;
							sec_end = b->KDIM() - 1;
							break;
						case 5:
						case 6:
							pri_end = b->KDIM() - 1;
							sec_end = b->IDIM() - 1;
							break;
						default:
							break;
						}
						for (size_t sec = 2; sec <= sec_end; ++sec)
							for (size_t pri = 2; pri <= pri_end; ++pri)
							{
								++cnt;
								b->surface_internal_node_occurance(i, pri, sec, sioc);
								for (auto r : sioc)
									*r = cnt;
							}
					}
				}
			}

			// Frame
			for (const auto &e : m_frame)
			{
				// Identify directions
				std::map<Block3D::FRAME*, size_t> ptr2idx;
				for (size_t i = 0; i < e.size(); ++i)
					ptr2idx[e[i]] = i;
				if (ptr2idx.size() != e.size())
					throw std::runtime_error("Something wired happens...");
				std::vector<bool> swap_flag(e.size(), false);
				std::vector<bool> visited(e.size(), false);
				visited[0] = true;
				std::queue<Block3D::FRAME*> q;
				q.push(e[0]);

				while (!q.empty())
				{
					auto r = q.front();
					q.pop();
					auto par_idx = ptr2idx[r];

					for (auto s : r->dependentSurf)
					{
						if (s->neighbourSurf)
						{
							size_t i = 0;
							while (i < s->includedFrame.size() && s->includedFrame[i] != r)
								++i;
							if (i >= s->includedFrame.size())
								throw std::runtime_error("Internal error detected.");

							auto rr = s->counterpartFrame[i];
							auto loc_idx = ptr2idx[rr];
							if (!visited[loc_idx])
							{
								visited[loc_idx] = true;
								swap_flag[loc_idx] = s->counterpartFrameIsOpposite[i] ? !swap_flag[par_idx] : swap_flag[par_idx];
								q.push(rr);
							}
						}
					}
				}

				// Assign index
				auto r0 = e[0];
				auto b0 = r0->dependentBlock;
				const auto itn = b0->frame_internal_node_num(r0->local_index);
				for (size_t lidx = 0; lidx < itn; ++lidx)
				{
					auto cur_cnt = ++cnt;
					for (size_t i = 0; i < e.size(); ++i)
					{
						std::vector<size_t*> fnoc(2, nullptr);
						auto r = e[i];
						auto b = r->dependentBlock;
						const size_t loc_pos = swap_flag[i] ? (itn + 1 - lidx) : (lidx + 2);
						b->frame_internal_node_occurace(r->local_index, loc_pos, fnoc);
						*fnoc[0] = *fnoc[1] = cur_cnt;
					}
				}
			}

			// Ensure total amount is right
			if (cnt != totalNodeNum)
				throw std::length_error("Inconsistent num of nodes detected.");
		}
	}
}
