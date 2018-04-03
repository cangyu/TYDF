#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <cctype>
#include <algorithm>

using namespace std;

const string bc_o2o("ONE_TO_ONE");
const string bc_sym("SYM");
const string bc_wal("WALL");
const string bc_far("FAR");

//unordered_map<string, int> nmf_bc_tbl({{bc_o2o, 0}, {bc_sym, 1}, {bc_wal, 2}, {bc_far, 3}});

class NMF_BCRange {
public:
    uint32_t blk_seq;
    uint32_t face_seq;
    uint32_t pri_start_seq, pri_end_seq;
    uint32_t sec_start_seq, sec_end_seq;

    NMF_BCRange(uint32_t *src) :
            blk_seq(src[0]), face_seq(src[1]),
            pri_start_seq(src[2]), pri_end_seq(src[3]),
            sec_start_seq(src[4]), sec_end_seq(src[5]) {
        //Empty Initialization Body
    }

    NMF_BCRange(uint32_t b = 0, uint32_t f = 0, uint32_t s1 = 0, uint32_t e1 = 0, uint32_t s2 = 0, uint32_t e2 = 0) :
            blk_seq(b), face_seq(f),
            pri_start_seq(s1), pri_end_seq(e1),
            sec_start_seq(s2), sec_end_seq(e2) {
        //Empty Initialization Body
    }
};

class NMF_Entry {
public:
    string bc;
    NMF_BCRange *rg1, *rg2;
    bool swap;

    NMF_Entry(string t, uint32_t *s) :
            bc(t),
            rg1(new NMF_BCRange(s)),
            rg2(nullptr),
            swap(false) {
        //Empty Initialization Body
    }

    NMF_Entry(string t, uint32_t *s1, uint32_t *s2, bool f) :
            bc(t),
            rg1(new NMF_BCRange(s1)),
            rg2(new NMF_BCRange(s2)),
            swap(f) {
        //Empty Initialization Body
    }

    ~NMF_Entry() {
        if (rg1)
            delete rg1;
        if (rg2)
            delete rg2;
    }
};

class NMF_BLKDim {
public:
    uint32_t i_dim, j_dim, k_dim;

    NMF_BLKDim(uint32_t in = 0, uint32_t jn = 0, uint32_t kn = 0) : i_dim(in), j_dim(jn), k_dim(kn) {
        //Empty Initialization Body
    }
};

class NMF {
public:
    uint32_t blk_num;
    vector<NMF_BLKDim> blk_dim;
    vector<NMF_Entry> mapping_entry;

    NMF() : blk_num(0), blk_dim(vector<NMF_BLKDim>()), mapping_entry(vector<NMF_Entry>()) {
        //Empty Initialization Body
    }

    NMF(uint32_t n) : blk_num(n), blk_dim(vector<NMF_BLKDim>(n)), mapping_entry(vector<NMF_Entry>()) {
        //Empty Initialization Body
    }
};

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage: gridglue map.nmf grid.fmt\n");
        return -1;
    }

    // Parse mapping file
    ifstream mfp;
    mfp.open(argv[1]);
    for (int i = 0; i < 4; i++)
        mfp.ignore(numeric_limits<streamsize>::max(), mfp.widen('\n'));


    uint32_t n = 0;
    mfp >> n;
    NMF bc_mp(n);
    for (uint32_t i = 0; i < n; i++) {
        uint32_t idx;
        mfp >> idx;
        mfp >> bc_mp.blk_dim[idx - 1].i_dim;
        mfp >> bc_mp.blk_dim[idx - 1].j_dim;
        mfp >> bc_mp.blk_dim[idx - 1].k_dim;
    }

    for (int i = 0; i < 4; i++)
        mfp.ignore(numeric_limits<streamsize>::max(), mfp.widen('\n'));

    string s;
    while (mfp >> s) {
        transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (!s.compare(bc_o2o)) {
            uint32_t tmp[2][6] = {0};
            for (uint8_t i = 0; i < 2; i++)
                for (uint8_t j = 0; j < 6; j++)
                    mfp >> tmp[i][j];
            bool flag = false;
            mfp >> flag;

            bc_mp.mapping_entry.push_back(NMF_Entry(s, tmp[0], tmp[1], flag));

        } else {
            uint32_t tmp[6] = {0};
            for (uint8_t i = 0; i < 6; i++)
                mfp >> tmp[i];

            bc_mp.mapping_entry.push_back(NMF_Entry(s, tmp));
        }
    }


    mfp.close();

    return 0;
}
