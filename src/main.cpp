#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <cctype>
#include <algorithm>
#include <cstring>

using namespace std;

const string swp_t("TRUE");
const string swp_f("FALSE");

const string bc_o2o("ONE_TO_ONE");
const string bc_sym("SYM");
const string bc_wal("WALL");
const string bc_far("FAR");

class NMF_BCRange
{
public:
	uint32_t blk_seq;
	uint32_t face_seq;
	uint32_t pri_start_seq, pri_end_seq;
	uint32_t sec_start_seq, sec_end_seq;

	NMF_BCRange(uint32_t *src) :
		blk_seq(src[0]), face_seq(src[1]),
		pri_start_seq(src[2]), pri_end_seq(src[3]),
		sec_start_seq(src[4]), sec_end_seq(src[5])
	{
		//Empty Initialization Body
	}

	NMF_BCRange() :
		blk_seq(0), face_seq(0),
		pri_start_seq(0), pri_end_seq(0),
		sec_start_seq(0), sec_end_seq(0)
	{
		//Empty Initialization Body
	}

	NMF_BCRange(uint32_t b, uint32_t f, uint32_t s1, uint32_t e1, uint32_t s2, uint32_t e2) :
		blk_seq(b), face_seq(f),
		pri_start_seq(s1), pri_end_seq(e1),
		sec_start_seq(s2), sec_end_seq(e2)
	{
		//Empty Initialization Body
	}

	uint32_t pri_node_num()
	{
		return pri_end_seq - pri_start_seq + 1;
	}

	uint32_t sec_node_num()
	{
		return sec_end_seq - sec_start_seq + 1;
	}

	uint32_t node_num()
	{
		return pri_node_num() * sec_node_num();
	}

	uint32_t face_num()
	{
		return (pri_node_num() - 1) * (sec_node_num() - 1);
	}
};

class NMF_Entry
{
public:
	string bc;
	NMF_BCRange *rg1, *rg2;
	bool swap;

	NMF_Entry(const string &t, uint32_t *s) :
		bc(t), rg1(new NMF_BCRange(s)), rg2(nullptr), swap(false)
	{
		//Empty Initialization Body
	}

	NMF_Entry(const string &t, uint32_t *s1, uint32_t *s2, bool f) :
		bc(t), rg1(new NMF_BCRange(s1)), rg2(new NMF_BCRange(s2)), swap(f)
	{
		//Empty Initialization Body
	}

	~NMF_Entry()
	{
		delete rg1;
		delete rg2;
	}
};

class NMF_BLKDim
{
public:
	uint32_t i_dim, j_dim, k_dim;
	uint32_t caste[6];

	void calc_caste()
	{
		if (i_dim == 0 || j_dim == 0 || k_dim == 0)
			for (uint8_t i = 0; i < 6; i++)
				caste[i] = 0;
		else
		{
			uint32_t u = i_dim, v = j_dim, w = k_dim;
			caste[0] = caste[1] = w * v;
			caste[2] = caste[3] = (u - 2) * w;
			caste[4] = caste[5] = (u - 2) * (v - 2);
			for (uint8_t i = 1; i < 6; i++)
				caste[i] += caste[i - 1];
		}
	}

	NMF_BLKDim() : i_dim(0), j_dim(0), k_dim(0)
	{
		calc_caste();
	}

	NMF_BLKDim(uint32_t in, uint32_t jn, uint32_t kn) :
		i_dim(in), j_dim(jn), k_dim(kn)
	{
		calc_caste();
	}

	uint32_t cell_num()
	{
		return (i_dim - 1) * (j_dim - 1) * (k_dim - 1);
	}

	uint32_t face_num()
	{
		uint32_t ret = 0;
		ret += i_dim * (j_dim - 1) * (k_dim - 1);
		ret += (i_dim - 1) * j_dim * (k_dim - 1);
		ret += (i_dim - 1) * (j_dim - 1) * k_dim;
		return ret;
	}

	uint32_t node_num()
	{
		return i_dim * j_dim * k_dim;
	}

	uint32_t internal_node_num()
	{
		return (i_dim - 2)*(j_dim - 2)*(k_dim - 2);
	}

	uint32_t shell_node_num()
	{
		return node_num() - internal_node_num();
	}

	uint32_t shell_node_idx(uint32_t i, uint32_t j, uint32_t k)
	{
		uint32_t u = i_dim, v = j_dim, w = k_dim;

		if (i == 0)
			return j + k * v;
		else if (i == u - 1)
			return caste[0] + j + k * v;
		else if (j == 0)
			return caste[1] + k + (i - 1) * w;
		else if (j == v - 1)
			return caste[2] + k + (i - 1) * w;
		else if (k == 0)
			return caste[3] + (i - 1) + (j - 1) * (u - 2);
		else if (k == w - 1)
			return caste[4] + (i - 1) + (j - 1) * (u - 2);
		else
			throw "Not a boundary node.";
	}

	void shell_node_coordinate(uint32_t idx, uint32_t &i, uint32_t &j, uint32_t &k)
	{
		uint32_t u = i_dim, v = j_dim, w = k_dim;

		if (idx < caste[0])
		{
			i = 0;
			j = idx % v;
			k = idx / v;
		}
		else if (idx < caste[1])
		{
			uint32_t cp = idx - caste[0];
			i = u - 1;
			j = cp % v;
			k = cp / v;
		}
		else if (idx < caste[2])
		{
			uint32_t cp = idx - caste[1];
			i = 1 + cp / w;
			j = 0;
			k = cp % w;
		}
		else if (idx < caste[3])
		{
			uint32_t cp = idx - caste[2];
			i = 1 + cp / w;
			j = v - 1;
			k = cp % w;
		}
		else if (idx < caste[4])
		{
			uint32_t cp = idx - caste[3];
			i = cp % (u - 2) + 1;
			j = cp / (u - 2) + 1;
			k = 0;
		}
		else if (idx < caste[5])
		{
			uint32_t cp = idx - caste[4];
			i = cp % (u - 2) + 1;
			j = cp / (u - 2) + 1;
			k = w - 1;
		}
		else
			throw "Invalid node index.";
	}
};

class NMF
{
public:
	uint32_t blk_num;
	vector<NMF_BLKDim> blk_dim;
	vector<NMF_Entry *> mapping_entry;

	NMF(uint32_t n = 0) :
		blk_num(n),
		blk_dim(vector<NMF_BLKDim>(n)),
		mapping_entry(vector<NMF_Entry *>())
	{
		//Empty Initialization Body
	}

	NMF(const string &path) : mapping_entry(vector<NMF_Entry *>())
	{
		//Load file
		ifstream mfp;
		mfp.open(path);
		if (mfp.fail())
		{
			cerr << "Can not open the neutral mapping file!" << endl;
			return;
		}

		//Skip header
		for (int i = 0; i < 4; i++)
			mfp.ignore(numeric_limits<streamsize>::max(), mfp.widen('\n'));

		//Read block dimension info
		mfp >> blk_num;
		blk_dim = vector<NMF_BLKDim>(blk_num);

		for (uint32_t i = 0; i < blk_num; i++)
		{
			uint32_t idx;
			mfp >> idx;
			mfp >> blk_dim[idx - 1].i_dim;
			mfp >> blk_dim[idx - 1].j_dim;
			mfp >> blk_dim[idx - 1].k_dim;
			blk_dim[idx - 1].calc_caste();
		}

		//Skip blank lines
		string tmp;
		while (getline(mfp, tmp) && tmp == "");

		//Skip seperators
		for (int i = 0; i < 3; i++)
			mfp.ignore(numeric_limits<streamsize>::max(), mfp.widen('\n'));

		//Read bc mappings
		string s;
		while (mfp >> s)
		{
			transform(s.begin(), s.end(), s.begin(), ::toupper);
			if (s == bc_o2o) {
				uint32_t tmp[2][6] = { 0 };
				for (uint8_t i = 0; i < 2; i++)
					for (uint8_t j = 0; j < 6; j++)
						mfp >> tmp[i][j];
				mfp >> s;
				transform(s.begin(), s.end(), s.begin(), ::toupper);
				bool flag = s == swp_t;

				mapping_entry.push_back(new NMF_Entry(bc_o2o, tmp[0], tmp[1], flag));
			}
			else {
				uint32_t tmp[6] = { 0 };
				for (uint8_t i = 0; i < 6; i++)
					mfp >> tmp[i];

				mapping_entry.push_back(new NMF_Entry(s, tmp));
			}
		}

		//Done
		mfp.close();
	}

	~NMF()
	{
		for (uint32_t i = 0; i < mapping_entry.size(); i++)
			delete mapping_entry[i];
	}

	int save(const string &path)
	{
		ofstream f_out;
		f_out.open(path);
		if (f_out.fail())
		{
			cerr << "Can not open target output file!" << endl;
			return -1;
		}

		f_out << "# ======================== Neutral Map File generated by the GridGlue software ===============================" << endl;
		f_out << "# ============================================================================================================" << endl;
		f_out << "# Block#    IDIM    JDIM    KDIM" << endl;
		f_out << "# ------------------------------------------------------------------------------------------------------------" << endl;
		f_out << setw(8) << right << blk_num << endl;
		for (uint32_t i = 0; i < blk_num; i++)
		{
			f_out << setw(8) << right << i + 1;
			f_out << setw(8) << right << blk_dim[i].i_dim;
			f_out << setw(8) << right << blk_dim[i].j_dim;
			f_out << setw(8) << right << blk_dim[i].k_dim << endl;
		}

		f_out << "# ============================================================================================================" << endl;
		f_out << "# Type           B1    F1       S1    E1       S2    E2       B2    F2       S1    E1       S2    E2      Swap" << endl;
		f_out << "# ------------------------------------------------------------------------------------------------------------" << endl;
		for (uint32_t i = 0; i < mapping_entry.size(); i++)
		{
			f_out << setw(13) << left << mapping_entry[i]->bc;
			f_out << setw(6) << right << mapping_entry[i]->rg1->blk_seq;
			f_out << setw(6) << right << mapping_entry[i]->rg1->face_seq;
			f_out << setw(9) << right << mapping_entry[i]->rg1->pri_start_seq;
			f_out << setw(6) << right << mapping_entry[i]->rg1->pri_end_seq;
			f_out << setw(9) << right << mapping_entry[i]->rg1->sec_start_seq;
			f_out << setw(6) << right << mapping_entry[i]->rg1->sec_end_seq;

			if (mapping_entry[i]->rg2)
			{
				f_out << setw(9) << right << mapping_entry[i]->rg2->blk_seq;
				f_out << setw(6) << right << mapping_entry[i]->rg2->face_seq;
				f_out << setw(9) << right << mapping_entry[i]->rg2->pri_start_seq;
				f_out << setw(6) << right << mapping_entry[i]->rg2->pri_end_seq;
				f_out << setw(9) << right << mapping_entry[i]->rg2->sec_start_seq;
				f_out << setw(6) << right << mapping_entry[i]->rg2->sec_end_seq;
				f_out << setw(10) << right << (mapping_entry[i]->swap ? swp_t : swp_f);
			}

			f_out << endl;
		}

		f_out.close();
		return 0;
	}
};

int main(int argc, char **argv)
{
	/*
	if (argc != 3)
	{
		cout << "Usage: gridglue map.nmf grid.fmt" << endl;
		return -1;
	}
	*/

	NMF bc_mp("../../test/mapping2.nmf");

	uint32_t total_cell_num = 0, total_face_num = 0, total_node_num = 0;
	for (uint32_t i = 0; i < bc_mp.blk_num; i++)
	{
		total_cell_num += bc_mp.blk_dim[i].cell_num();
		total_face_num += bc_mp.blk_dim[i].face_num();
		total_node_num += bc_mp.blk_dim[i].internal_node_num();
	}

	for (uint32_t i = 0; i < bc_mp.mapping_entry.size(); i++)
		if (bc_mp.mapping_entry[i]->bc == bc_o2o)
			total_face_num -= bc_mp.mapping_entry[i]->rg1->face_num();

	vector<vector<uint32_t>> shell_node_seq(bc_mp.blk_num);
	for (uint32_t i = 0; i < bc_mp.blk_num; i++)
	{
		shell_node_seq[i].resize(bc_mp.blk_dim[i].shell_node_num());
		memset(shell_node_seq[i].data(), 0, shell_node_seq[i].size() * sizeof(uint32_t));
	}
		


	cout << total_cell_num << endl;
	cout << total_face_num << endl;
	cout << total_node_num << endl;

	bc_mp.save("../../test/test2.nmf");

	return 0;
}
