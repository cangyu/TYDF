#include "../inc/nmf.h"

static void str_formalize(std::string &s)
{
	std::transform(s.begin(), s.end(), s.begin(), ::toupper);
	for (auto &e : s)
		if (e == '-')
			e = '_';
}

namespace GridTool::NMF
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
					// NOT release all existing resources until 
					// it is ensured that this input file is valid.
					release_all();

					// Re-Allocate storage for new recordings
					m_blk.resize(NumOfBlk, nullptr);
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
					const int idx = std::stoi(res2[1].str());
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
	}